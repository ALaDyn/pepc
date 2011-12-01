
module module_pepcfields
    implicit none
    private

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    public pepc_fields
    public pepc_grid_fields

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  subroutine-implementation  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Calculate fields and potential for supplied particles, work is taken from t_particle array
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine pepc_fields(np_local, npart_total, particles, &
        np_mult_, cf_par, itime, weighted, curve_type, num_neighbours, neighbours, no_dealloc, no_restore)

        use treevars
        use module_interaction_specific
        use timings
        use module_tree_domains
        use module_calc_force
        use tree_walk_pthreads
        use tree_walk_communicator
        use module_allocation
        use module_tree
        use module_htable
        use module_branching
        use module_mirror_boxes
        implicit none
        include 'mpif.h'

        integer, intent(inout) :: np_local    !< # particles on this CPU
        integer, intent(in) :: npart_total !< total # simulation particles
        type(t_particle), allocatable, intent(inout) :: particles(:)
        real, intent(in) :: np_mult_
        type(t_calc_force_params), intent(in) :: cf_par
        integer, intent(in) :: itime  ! timestep
        integer, intent(in) :: weighted   ! TODO: put into cf_par
        integer, intent(in) :: curve_type ! type of space-filling curve TODO: put into cf_par
        integer, intent(in) :: num_neighbours !< number of shift vectors in neighbours list (must be at least 1 since [0, 0, 0] has to be inside the list)
        integer, intent(in) :: neighbours(3, num_neighbours) ! list with shift vectors to neighbour boxes that shall be included in interaction calculation, at least [0, 0, 0] should be inside this list
        logical, intent(in) :: no_dealloc, no_restore

        integer :: ierr, i
        integer :: npnew, npold
        integer, allocatable :: indxl(:),irnkl(:)
        integer :: islen(num_pe),irlen(num_pe)
        integer :: fposts(num_pe+1),gposts(num_pe+1)
        real*8 :: ttrav, ttrav_loc, tcomm(3) ! timing integrals
        integer :: ibox
        real*8 :: vbox(3)
        character(30) :: cfile
        integer*8, allocatable :: leaf_keys(:)
        integer :: nppmax
        integer :: neighbour_pe_particles !< number of particles that have been fetched from neighbouring PEs durin tree_domains

        call status('FIELDS')

        ! copy call parameters to treevars module
        npart      = npart_total
        np_mult    = np_mult_
        npp        = np_local

        call timer_start(t_all)

        ! workload per particle must be nonzero
        do i=1,np_local
          particles(i)%work = max(particles(i)%work, 1.)
        end do

        call timer_start(t_fields_tree)

        !TODO: make adjustable by user or find a good estimation. Additional Question: Does this value have to be globally constant?
        nppmax = int(1.25 * max(npart/num_pe,1000)) ! allow 25% fluctuation around average particle number per PE in sorting library for load balancing

        allocate(indxl(nppmax),irnkl(nppmax))
        ! Domain decomposition: allocate particle keys to PEs
        call tree_domains(particles, nppmax,indxl,irnkl,islen,irlen,fposts,gposts,npnew,npold, weighted, curve_type, neighbour_pe_particles)
        call allocate_tree(cf_par%theta)

        ! build local part of tree
        call timer_start(t_local)
        call htable_clear_and_insert_root()
        allocate(leaf_keys(npp+neighbour_pe_particles))
        call tree_build_from_particles(particles, npp+neighbour_pe_particles, leaf_keys)
        ! remove the boundary particles from the htable - we are not interested in them any more
        call htable_remove_keys(leaf_keys(npp+1:npp+neighbour_pe_particles), neighbour_pe_particles)
        neighbour_pe_particles = 0
        ! build tree from local particle keys up to root (the boundary particles are not included in the tree construction)
        call tree_build_upwards(leaf_keys(1:npp), npp)
        deallocate(leaf_keys)

        if (htable(1)%leaves .ne. npp) then
            write(*,*) 'PE', me, ' did not find all its particles inside the htable after local tree buildup: htable(1)%leaves =', htable(1)%leaves, ' but npp =', npp
            call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        endif
        call timer_stop(t_local)

        ! Should now have multipole information up to root list level(s) (only up to branch level, the information is correct)
        ! By definition, this is complete: each branch node is self-contained.
        ! This information has to be broadcast to the other PEs so that the top levels can be filled in.

        ! identification of branch nodes
        call timer_start(t_branches_find)
        call find_branches()
        call timer_stop(t_branches_find)

        ! exchange branch nodes
        nleaf_me = nleaf       !  Retain leaves and twigs belonging to local PE
        ntwig_me = ntwig
        call timer_start(t_exchange_branches)
        call tree_exchange(pebranch, nbranch, branch_key, nbranch_sum)
        call timer_stop(t_exchange_branches)
        ! build global part of tree
        call timer_start(t_global)
        call tree_build_upwards(branch_key, nbranch_sum)
        call timer_stop(t_global)

        if (htable(1)%leaves .ne. npart_total) then
            write(*,*) 'PE', me, ' did not find all particles inside the htable after global tree buildup: htable(1)%leaves =', htable(1)%leaves, ' but npp =', npp
            call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
        endif

        call timer_stop(t_fields_tree)
        call timer_reset(t_walk)
        call timer_reset(t_walk_local)
        call timer_reset(t_comm_total)
        call timer_reset(t_comm_recv)
        call timer_reset(t_comm_sendreqs)

        max_req_list_length = 0
        cum_req_list_length = 0
        comm_loop_iterations = 0
        sum_fetches=0      ! total # multipole fetches/iteration
        sum_ships=0      ! total # multipole shipments/iteration

        call timer_start(t_fields_passes)

        call particleresults_clear(particles(:))

        do ibox = 1,num_neighbours ! sum over all boxes within ws=1

            vbox = lattice_vect(neighbours(:,ibox))

            ! tree walk finds interaction partners and calls interaction routine for particles on short list
            call tree_walk(npp,particles,cf_par,ttrav,ttrav_loc, vbox, tcomm)

            call timer_add(t_walk, ttrav)           ! traversal time (until all walks are finished)
            call timer_add(t_walk_local, ttrav_loc) ! traversal time (local)
            call timer_add(t_comm_total,    tcomm(TIMING_COMMLOOP))
            call timer_add(t_comm_recv,     tcomm(TIMING_RECEIVE))
            call timer_add(t_comm_sendreqs, tcomm(TIMING_SENDREQS))

        end do ! ibox = 1,num_neighbours

        ! add lattice contribution
        call timer_start(t_lattice)
        ! add lattice contribution and other per-particle-forces
        ! TODO: do not call calc_force_per_particle here!
        call calc_force_per_particle(particles, npp, cf_par)
        call timer_stop(t_lattice)

        ! restore initial particle order specified by calling routine to reassign computed forces
        ! notice the swapped order of the index-fields -> less changes in restore.f90 compared to tree_domains.f90

        call timer_stop(t_fields_passes)
        call timer_start(t_restore)

        if (no_restore) then
          ! we have to inform the calling routine, that the particle number has changed and all fields have been reallocated
          np_local = npnew ! is equal to npp
        else
          call restore(npnew,npold,nppmax,irnkl,indxl,irlen,islen,gposts,fposts, particles)
        endif

        call timer_stop(t_restore)

        deallocate(indxl,irnkl)

        call timer_start(t_fields_stats)

        nkeys_total = nleaf+ntwig

        if (tree_debug) call tree_stats(itime)

        call timer_stop(t_fields_stats)

        if( load_file_debug ) then
            call system("mkdir -p " // "load")
            write(cfile,'("load/load_",i6.6,".dat")') me
            open(60, file=trim(cfile),STATUS='UNKNOWN', POSITION = 'APPEND')
            write(60,'(i5,2f20.10, i12)') itime,interactions_local, mac_evaluations_local,npp
            close(60)
        end if

        ! deallocate particle and result arrays
        call timer_start(t_deallocate)

        if (.not. no_dealloc) then
            call deallocate_tree()
        endif

        call timer_stop(t_deallocate)
        call timer_stop(t_all)

        call status('FIELDS DONE')


    end subroutine pepc_fields


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Calculate fields and potential for coordinates of supplied particles
    !> without inserting them into the tree or including them in the actual computation
    !> everything in the supplied t_particle structure except the coordinates is simply ignored
    !> i.e. these particles behave like neutral testparticles
    !>
    !> before invoking this function, pepc_fields() must have been call with no_dealloc=.true.
    !>
    !> the coordinates should overlap with the local tree(!) domain
    !> TODO: provide function to prepare such a grid, take care whether no_backsort=.true./.false.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine pepc_grid_fields(npoints_local, particles, cf_par, num_neighbours, neighbours)

        use treevars
        use module_interaction_specific
        use timings
        use module_tree_domains
        use module_calc_force
        use tree_walk_pthreads
        use tree_walk_communicator
        use module_allocation
        use module_tree
        use module_htable
        use module_branching
        use module_mirror_boxes
        implicit none
        include 'mpif.h'

        integer, intent(in) :: npoints_local    !< # points on this CPU
        type(t_particle), intent(inout) :: particles(:)
        type(t_calc_force_params), intent(in) :: cf_par
        integer, intent(in) :: num_neighbours !< number of shift vectors in neighbours list (must be at least 1 since [0, 0, 0] has to be inside the list)
        integer, intent(in) :: neighbours(3, num_neighbours) ! list with shift vectors to neighbour boxes that shall be included in interaction calculation, at least [0, 0, 0] should be inside this list

        real*8 :: ttrav, ttrav_loc, tcomm(3) ! timing integrals
        integer :: ibox
        real*8 :: vbox(3)


        call status('GRID FIELDS')

        call timer_start(t_walk_grid)

        call particleresults_clear(particles(:))

        do ibox = 1,num_neighbours ! sum over all boxes within ws=1

            vbox = lattice_vect(neighbours(:,ibox))

            ! tree walk finds interaction partners and calls interaction routine for particles on short list
            call tree_walk(npoints_local,particles,cf_par,ttrav,ttrav_loc, vbox, tcomm)

            call timer_add(t_walk, ttrav)           ! traversal time (until all walks are finished)
            call timer_add(t_walk_local, ttrav_loc) ! traversal time (local)
            call timer_add(t_comm_total,    tcomm(TIMING_COMMLOOP))
            call timer_add(t_comm_recv,     tcomm(TIMING_RECEIVE))
            call timer_add(t_comm_sendreqs, tcomm(TIMING_SENDREQS))

        end do ! ibox = 1,num_neighbours

        call timer_stop(t_walk_grid)

        ! add lattice contribution
        call timer_start(t_lattice_grid)
        ! add lattice contribution and other per-particle-forces
        ! TODO: do not call calc_force_per_particle here!
        call calc_force_per_particle(particles, npp, cf_par)
        call timer_stop(t_lattice_grid)

        call status('GRID FIELDS DONE')

    end subroutine pepc_grid_fields


end module module_pepcfields







