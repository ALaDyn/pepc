
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
        !>   Calculate fields and potential for supplied particle coordinates p_x, p_y, p_z and charges p_q
        !>
        !>   Returns fields Ex, Ey, Ez and potential pot excluding external terms
        !>   @param[in] np_local local number of particles
        !>   @param[in] npart_total total particle number
        !>   @param[in] p_x dimension(1:np_local) - x-component of particle coordinates
        !>   @param[in] p_y dimension(1:np_local) - y-component of particle coordinates
        !>   @param[in] p_z dimension(1:np_local) - z-component of particle coordinates
        !>   @param[in] p_q dimension(1:np_local) - particle charge
        !>   @param[in] p_w dimension(1:np_local) - particle workload from previous iteration (should be set to 1.0 for the first timestep)
        !>   @param[in] p_label dimension(1:np_local) - particle label (may any number except zero)
        !>   @param[out] p_Ex dimension(1:np_local) - x-component of electric field
        !>   @param[out] p_Ey dimension(1:np_local) - y-component of electric field
        !>   @param[out] p_Ez dimension(1:np_local) - z-component of electric field
        !>   @param[out] p_pot dimension(1:np_local) - electric potential
        !>   @param[in] np_mult_ memory allocation parameter
        !>   @param[in] cf_par parameters for force summation
        !>   @param[in] itime current simulation timestep number
        !>   @param[in] weighted selector for load balancing
        !>   @param[in] curve_type selector for type of space filling curve
        !>   @param[in] num_neighbours number of neighbour boxes to be considered during tree walk
        !>   @param[in] neighbours shift vectors to neighbour boxes
        !>   @param[in] no_dealloc if set to .true., deallocation of tree-structures is prevented to allow for front-end triggered diagnostics
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine pepc_fields(np_local,npart_total,p_x, p_y, p_z, p_q, p_w, p_label, &
	     p_Ex, p_Ey, p_Ez, p_pot, np_mult_, cf_par, itime, weighted, curve_type, &
	     num_neighbours, neighbours, no_dealloc)

	  use treevars
      use module_multipole_helpers
	  use timings
	  use module_tree_domains
      use module_fmm_framework
      use module_calc_force
	  use tree_walk_pthreads
	  use tree_walk_communicator
	  use module_allocation
	  implicit none
	  include 'mpif.h'

	  integer, intent(in) :: np_local  ! # particles on this CPU
	  integer, intent(in) :: npart_total ! total # simulation particles
	  real, intent(in) :: np_mult_       ! multipole opening angle
	  type(t_calc_force_params), intent(in) :: cf_par
	  integer, intent(in) :: itime  ! timestep
	  integer, intent(in) :: weighted
	  real*8, intent(in), dimension(np_local) :: p_x, p_y, p_z  ! coords and velocities: x1,x2,x3, y1,y2,y3, etc
	  real*8, intent(in), dimension(np_local) :: p_q ! charges, masses
	  integer, intent(in), dimension(np_local) :: p_label  ! particle label
	  real*8, intent(out), dimension(np_local) :: p_ex, p_ey, p_ez, p_pot  ! fields and potential to return
	  integer, intent(in) :: num_neighbours !< number of shift vectors in neighbours list (must be at least 1 since [0, 0, 0] has to be inside the list)
	  integer, intent(in) :: neighbours(3, num_neighbours) ! list with shift vectors to neighbour boxes that shall be included in interaction calculation, at least [0, 0, 0] should be inside this list
	  real*8, dimension(np_local) :: p_w ! work loads
	  integer, intent(in) :: curve_type ! type of space-filling curve
	  logical, intent(in) :: no_dealloc
	  integer :: nppm_ori, ierr
	  integer :: npnew, npold
	  integer, allocatable :: indxl(:),irnkl(:)
	  integer :: islen(num_pe),irlen(num_pe)
	  integer :: fposts(num_pe+1),gposts(num_pe+1)
	  integer :: i
	  real*8 :: ttrav, ttrav_loc, tcomm(3) ! timing integrals
	  integer :: ibox
	  real*8 :: vbox(3)
	  character(30) :: cfile

      ! fields, potential and load weights returned by force-sum: allocated in pepc_fields:
      type(t_particle_results), allocatable :: particle_results(:)

      if (me==0 .and. tree_debug) write(*,'(a)') 'LPEPC | FIELDS..'

	  ! copy call parameters to treevars module
	  npart      = npart_total
	  np_mult    = np_mult_
	  npp        = np_local

	  call allocate_particles(nppm_ori)
	  allocate(particle_results(nppm_ori))

	  call timer_start(t_all)
	  call timer_start(t_fields_begin)

	  if (force_debug) then
	     write (*,'(a7,a50/2i5,4f15.2)') 'PEPC | ','Params: itime, mac, theta, eps, force_const:', &
				itime, cf_par%mac, cf_par%theta, cf_par%eps, cf_par%force_const
	     write (*,'(a7,a20/(i16,4f15.3,i8))') 'PEPC | ','Initial buffers: ',(p_label(i), p_x(i), p_y(i), p_z(i), p_q(i), &
				p_label(i),i=1,npp)
	  endif

          ! Copy particle buffers to tree arrays
          do i=1,npp
            particles(i) = t_particle([p_x(i), p_y(i), p_z(i)], &  ! position
                                      [    0.,     0.,     0.], &  ! velocity - not relevant for tree code
                                       p_q(i),                  &  ! charge
                                       max(p_w(i), 1._8),       &  ! workload from last step
                                       -1_8,                    &  ! key - will be assigned later
                                       p_label(i),              &  ! particle label for tracking purposes
                                       me )                        ! particle owner
          end do

	  ! Trap bad particle labels
	  if (any(p_label(1:npp) == 0)) then
	    write (*,*) '*** Error: particle labels must be nonzero:', i, p_label(i)
	    call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
	  endif

	  call timer_stop(t_fields_begin)
	  call timer_start(t_fields_tree)

	  allocate(indxl(nppm_ori),irnkl(nppm_ori))
	  ! Domain decomposition: allocate particle keys to PEs
	  call tree_domains(indxl,irnkl,islen,irlen,fposts,gposts,npnew,npold, weighted, curve_type)
	  call allocate_tree(cf_par%theta)

	  ! calculate spherical multipole expansion of central box
	  if (cf_par%include_far_field_if_periodic) call fmm_framework_timestep()

	  ! build local part of tree
	  call timer_stamp(t_stamp_before_local)
	  call tree_local
	  ! exchange branch nodes
	  call timer_stamp(t_stamp_before_exchange)
	  call tree_exchange
	  ! build global part of tree
	  call timer_stamp(t_stamp_before_global)
	  call tree_global

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

      particle_results(:) = t_particle_results([0., 0., 0.], 0., 1.)

	  call timer_stamp(t_stamp_before_walkloop)

	  do ibox = 1,num_neighbours ! sum over all boxes within ws=1

	    vbox = lattice_vect(neighbours(:,ibox))

	    ! tree walk finds interaction partners and calls interaction routine for particles on short list
	    call tree_walk(npp,particles,particle_results,cf_par,itime,ttrav,ttrav_loc, vbox, tcomm)

        call timer_add(t_walk, ttrav)           ! traversal time (until all walks are finished)
        call timer_add(t_walk_local, ttrav_loc) ! traversal time (local)
	    call timer_add(t_comm_total,    tcomm(TIMING_COMMLOOP))
	    call timer_add(t_comm_recv,     tcomm(TIMING_RECEIVE))
	    call timer_add(t_comm_sendreqs, tcomm(TIMING_SENDREQS))

	  end do ! ibox = 1,num_neighbours

	  call timer_stamp(t_stamp_after_walkloop)

	  ! add lattice contribution
	  call timer_start(t_lattice)
	  ! add lattice contribution and other per-particle-forces
	  call calc_force_per_particle(npp, particles, particle_results, cf_par)
	  call timer_stop(t_lattice)

	  ! restore initial particle order specified by calling routine to reassign computed forces
	  ! notice the swapped order of the index-fields -> less changes in restore.f90 compared to tree_domains.f90

	  call timer_stop(t_fields_passes)
	  call timer_start(t_restore)

	  call restore(npnew,npold,nppm_ori,irnkl,indxl,irlen,islen,gposts,fposts, particle_results,p_pot,p_ex,p_ey,p_ez,p_w)

	  call timer_stop(t_restore)

	  deallocate(indxl,irnkl)

	  call timer_start(t_fields_stats)

	  nkeys_total = nleaf+ntwig

	  if (force_debug) then
	     write (ipefile,'("Tree forces:"/"   p    q   m   ux   pot  ",f8.2)')
	     write (*,'("Tree forces:"/"   p    q   m   ux   pot  ",f8.2)')
	     write (ipefile,'("Tree forces:"/"   p    q   m   ux   pot  ",f8.2)') cf_par%force_const

	     do i=1,np_local
	        write (ipefile,'(1x,i7,4(1pe14.5))') particles(i)%label, particles(i)%q, particles(i)%u(1), p_pot(i), p_ex(i)
	        write (*,'(1x,i7,4(1pe14.5))') particles(i)%label, particles(i)%x(1), particles(i)%q, particles(i)%u(1), p_pot(i)
	     end do

	  endif

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

          deallocate(particle_results)
          if (.not. no_dealloc) then
            call deallocate_tree(nppm_ori)
            call deallocate_particles()
          endif

          call timer_stop(t_deallocate)
          call timer_stop(t_all)

	end subroutine pepc_fields


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !>   Calculate fields and potential from gridded coordinates x,y,z **using existing tree**
        !>   pepc_fields(...,no_dealloc+.true.) must have been called before
        !>   this call does not modify the tree excpet that it fetches additionally needed particles
        !>   ideally, the supplied grid coordinates on each PE should coincide with
        !>   the spatial volume that is covered by the local branch nodes
        !>   TODO: add a routine, that produces such a grid automatically
        !>
        !>   Returns fields Ex, Ey, Ez and potential pot excluding external terms
        !>
        !>   @param[in] npgrid local number of grid points
        !>   @param[in] npart_total total particle number
        !>   @param[in] p_x dimension(1:np_local) - x-component of particle coordinates
        !>   @param[in] p_y dimension(1:np_local) - y-component of particle coordinates
        !>   @param[in] p_z dimension(1:np_local) - z-component of particle coordinates
        !>   @param[in] p_label dimension(1:np_local) - particle label (may any number except zero)
        !>   @param[out] p_Ex dimension(1:np_local) - x-component of electric field
        !>   @param[out] p_Ey dimension(1:np_local) - y-component of electric field
        !>   @param[out] p_Ez dimension(1:np_local) - z-component of electric field
        !>   @param[out] p_pot dimension(1:np_local) - electric potential
        !>   @param[in] cf_par parameters for force summation
        !>   @param[in] itime current simulation timestep number
        !>   @param[in] weighted selector for load balancing
        !>   @param[in] curve_type selector for type of space filling curve
        !>   @param[in] num_neighbours number of neighbour boxes to be considered during tree walk
        !>   @param[in] neighbours shift vectors to neighbour boxes
        !>   @param[in] no_dealloc if set to .true., deallocation of tree-structures is prevented to allow for front-end triggered diagnostics
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine pepc_grid_fields(npgrid,p_x, p_y, p_z, p_label, &
	     p_Ex, p_Ey, p_Ez, p_pot, cf_par, itime,  num_neighbours, neighbours)

	  use treevars
      use module_multipole_helpers
	  use module_htable
	  use timings
	  use module_calc_force
	  use module_fmm_framework
	  use tree_walk_pthreads
	  use tree_walk_communicator

	  implicit none
	  include 'mpif.h'

	  integer, intent(in) :: npgrid  ! # particles on this CPU
	  type(t_calc_force_params), intent(in) :: cf_par
	  integer, intent(in) :: itime  ! timestep
	  real*8, intent(in), dimension(npgrid) :: p_x, p_y, p_z  ! coords: x1,x2,x3, y1,y2,y3, etc
	  integer, intent(in), dimension(npgrid) :: p_label  ! particle label
	  real*8, intent(out), dimension(npgrid) :: p_ex, p_ey, p_ez, p_pot  ! fields and potential to return
	  integer, intent(in) :: num_neighbours !< number of shift vectors in neighbours list (must be at least 1 since [0, 0, 0] has to be inside the list)
	  integer, intent(in) :: neighbours(3, num_neighbours) ! list with shift vectors to neighbour boxes that shall be included in interaction calculation, at least [0, 0, 0] should be inside this list

      type(t_particle), allocatable :: grid_particles(:)
      type(t_particle_results), allocatable :: grid_particle_results(:)

	  integer :: i
	  real*8 :: ttrav, ttrav_loc, tcomm(3) ! timing integrals
	  integer :: ibox
	  real*8 :: vbox(3)

      if (me==0 .and. tree_debug) write(*,'(a)') 'LPEPC | GRID-FIELDS..'

      if (.not. (allocated(htable) .and. allocated(tree_nodes))) then
         write(*,*) 'pepc_grid_fields(): pepc_fields() must have been called with no_dealloc=.true. before'
        return
      endif

	  allocate(grid_particles(npgrid), grid_particle_results(npgrid))

      ! Copy particle buffers to tree arrays
      do i=1,npgrid
        grid_particles(i) = t_particle([p_x(i), p_y(i), p_z(i)], &  ! position
                                       [    0.,     0.,     0.], &  ! velocity - not relevant for tree code
                                        0.,                      &  ! charge - set to zero
                                        1._8,                    &  ! workload
                                       -1_8,                     &  ! key - will be assigned later
                                        p_label(i),              &  ! particle label
                                        me )                        ! particle owner
          end do

      grid_particle_results(:) = t_particle_results([0., 0., 0.], 0., 1.)

	  do ibox = 1,num_neighbours ! sum over all boxes within ws=1
	    vbox = lattice_vect(neighbours(:,ibox))
	    ! tree walk finds interaction partners and calls interaction routine for particles on short list
	    call tree_walk(npgrid,grid_particles,grid_particle_results,cf_par,itime,ttrav,ttrav_loc, vbox, tcomm)
	  end do ! ibox = 1,num_neighbours

	  ! add lattice contribution and other per-particle-forces
!	  call calc_force_per_particle(npgrid, grid_particles, grid_particle_results, cf_par)

	  nkeys_total = nleaf+ntwig

      ! Copy results back to local arrays
      do i=1,npgrid
        p_Ex(i)  = grid_particle_results(i)%e(1)
        p_Ey(i)  = grid_particle_results(i)%e(2)
        p_Ez(i)  = grid_particle_results(i)%e(3)
        p_pot(i) = grid_particle_results(i)%pot
      end do

      ! deallocate particle and result arrays
      deallocate (grid_particles, grid_particle_results)

	end subroutine pepc_grid_fields


end module module_pepcfields







