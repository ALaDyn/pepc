
module module_libpepc_main
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

    public libpepc_restore_particles
    public libpepc_traverse_tree
    public libpepc_grow_tree

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    integer, allocatable ::  indxl(:),  irnkl(:)
    integer, allocatable :: fposts(:), gposts(:)
    integer, allocatable ::  islen(:),  irlen(:)
    integer :: nppmax
    integer :: npnew !> particle number after domain decomposition
    integer :: npold !< original particle number (before domain decomposition)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  subroutine-implementation  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


contains


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Builds the tree from the given particles, redistributes particles
    !> to other MPI ranks if necessary (i.e. reallocates particles and changes np_local)
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine libpepc_grow_tree(np_local, npart_total, particles)
        use module_htable
        use module_branching
        use treevars
        use module_pepc_types
        use module_timings
        use module_tree
        use module_domains
        use module_debug, only : pepc_status
        use module_allocation, only : allocate_tree
        implicit none
        include 'mpif.h'

        integer, intent(inout) :: np_local    !< number of particles on this CPU, i.e. number of particles in particles-array
        integer, intent(in) :: npart_total !< total number of simulation particles (sum over np_local over all MPI ranks)
        type(t_particle), allocatable, intent(inout) :: particles(:) !< input particle data, initializes %x, %data, %work appropriately (and optionally set %label) before calling this function

        integer*8, allocatable :: leaf_keys(:)
        integer :: neighbour_pe_particles !< number of particles that have been fetched from neighbouring PEs durin tree_domains
        integer :: i, ierr

        call pepc_status('GROW TREE')

        call MPI_BARRIER( MPI_COMM_WORLD, ierr)  ! Wait for everyone to catch up

        ! copy call parameters to treevars module
        npart      = npart_total
        npp        = np_local

        call timer_start(t_all)

        ! workload per particle must be nonzero
        do i=1,np_local
          particles(i)%work = max(particles(i)%work, 1.)
        end do

        call timer_start(t_fields_tree)

        !TODO: make adjustable by user or find a good estimation. Additional Question: Does this value have to be globally constant?
        nppmax = int(1.25 * max(npart/num_pe,1000)) ! allow 25% fluctuation around average particle number per PE in sorting library for load balancing

        ! fields for sorting library results, we have to deallocate them in case someone did not call restore()
        if (allocated(indxl)) deallocate(indxl, irnkl, fposts, gposts, islen, irlen)
        allocate(indxl(nppmax), irnkl(nppmax), fposts(num_pe+1), gposts(num_pe+1), islen(num_pe), irlen(num_pe))

        ! Domain decomposition: allocate particle keys to PEs
        call tree_domains(particles, nppmax,indxl,irnkl,islen,irlen,fposts,gposts,npnew,npold, neighbour_pe_particles)
        np_local = npnew ! == npp, just to inform the calling routine about the current size of the particles-field
        call allocate_tree(0.09_8) !(theta2) TODO: the argument of this function should be theta^2 for BH-MAC, something else (what) for other cases, maybe, make this a call parameter for tree_Walk or make it interaction-specific

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

        nkeys_total = nleaf+ntwig

        call timer_stop(t_fields_tree)

        call pepc_status('TREE GROWN')

    end subroutine



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Traverses the complete tree for the given particles, i.e. computes
    !> the field values at their positions. Although missing information
    !> is automatically requested from remote MPI ranks, it is important
    !> that the particle coordinates fit to the local MPI ranks domain
    !> to avoid excessive communication
    !> If field values on some regular grid are needed, they can be
    !> generated using pepc_prepare_local_grid()
    !> Otherwise, it makes sense to provide the same particles as given/returned
    !> from to pepc_grow_tree()
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine libpepc_traverse_tree(nparticles, particles, clearresults_before_traversal)
        use module_pepc_types
        use treevars
        use module_walk
        use module_mirror_boxes
        use module_walk_communicator
        use module_timings
        use module_calc_force
        use module_debug, only : pepc_status
        implicit none
        integer, intent(in) :: nparticles    !< number of particles on this CPU, i.e. number of particles in particles-array
        type(t_particle), allocatable, intent(inout) :: particles(:) !< input particle data, initializes %x, %data, %work appropriately (and optionally set %label) before calling this function
        logical, intent(in), optional :: clearresults_before_traversal

        logical :: clearresults
        integer :: ibox
        real*8 :: ttrav, ttrav_loc, tcomm(3) ! timing integrals

        call pepc_status('TRAVERSE TREE')

        clearresults = .true.

        if (present(clearresults_before_traversal)) clearresults = clearresults_before_traversal

        call timer_reset(t_walk)
        call timer_reset(t_walk_local)
        call timer_reset(t_comm_total)
        call timer_reset(t_comm_recv)
        call timer_reset(t_comm_sendreqs)

        max_req_list_length  = 0
        cum_req_list_length  = 0
        comm_loop_iterations = 0
        sum_fetches          = 0 ! total # multipole fetches/iteration
        sum_ships            = 0 ! total # multipole shipments/iteration

        call timer_start(t_fields_passes)

        if (clearresults) call particleresults_clear(particles, npp)


        do ibox = 1,num_neighbour_boxes ! sum over all boxes within ws=1

            ! tree walk finds interaction partners and calls interaction routine for particles on short list
            call tree_walk(nparticles,particles,ttrav,ttrav_loc, lattice_vect(neighbour_boxes(:,ibox)), tcomm)

            call timer_add(t_walk, ttrav)           ! traversal time (until all walks are finished)
            call timer_add(t_walk_local, ttrav_loc) ! traversal time (local)
            call timer_add(t_comm_total,    tcomm(TIMING_COMMLOOP))
            call timer_add(t_comm_recv,     tcomm(TIMING_RECEIVE))
            call timer_add(t_comm_sendreqs, tcomm(TIMING_SENDREQS))

        end do ! ibox = 1,num_neighbour_boxes

        ! add lattice contribution
        call timer_start(t_lattice)
        ! add lattice contribution and other per-particle-forces
        ! TODO: do not call calc_force_per_particle here!
        call calc_force_per_particle(particles, npp)
        call timer_stop(t_lattice)

        call timer_stop(t_fields_passes)

        nkeys_total = nleaf+ntwig

        call pepc_status('TRAVERSAL DONE')

    end subroutine


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Restores the initial particle distribution (before calling pepc_grow_tree() ).
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine libpepc_restore_particles(np_local, particles)
        use module_pepc_types
        use module_timings
        use module_debug, only : pepc_status
        use module_domains, only : restore
        implicit none
        integer, intent(inout) :: np_local    !< number of particles on this CPU, i.e. number of particles in particles-array
        type(t_particle), allocatable, intent(inout) :: particles(:) !< input particle data on local MPI rank - is replaced by original particle data that was given before calling pepc_grow_tree()

        call pepc_status('RESTORE PARTICLES')

        ! restore initial particle order specified by calling routine to reassign computed forces
        ! notice the swapped order of the index-fields -> less changes in restore() compared to tree_domains()
        call timer_start(t_restore)
        call restore(npnew,npold,nppmax,irnkl,indxl,irlen,islen,gposts,fposts, particles)
        call timer_stop(t_restore)
      
        np_local = npold ! we have got the original particle order again

        deallocate(indxl, irnkl, fposts, gposts, islen, irlen)

        call pepc_status('RESTORATION DONE')

    end subroutine


end module module_libpepc_main







