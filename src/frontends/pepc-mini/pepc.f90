
! ==============================================================
!
!
!                  PEPC-MINI
!
!    Parallel Efficient Parallel Coulomb-solver: Electrostatics 
!
!   $ Revision $
!
!   Driver code for Coulomb-solver library lpepc
!
!  Major changes:
!   June 2005: Development begun
!
!   See README.compile for summary of program units
!  ==============================================================

program pepcmini

  use treetypes
  use physvars
  use particle_pusher
  use timings
  use module_fmm_framework
  use module_mirror_boxes
  use module_pepcfields
  use files
  use module_setup
  use module_calc_force, only : theta2, eps2, mac_select, force_law
  implicit none
  include 'mpif.h'

  integer :: ierr, ifile

  ! Allocate array space for tree
  call libpepc_setup("pepc-mini", my_rank, n_cpu)

  ! Set up O/P files
  call openfiles

  ! Time stamp
  if (my_rank==0) call stamp(6,1)
  if (my_rank==0) call stamp(15,1)

  ! Each CPU gets copy of initial data
  call pepc_setup()

  ! initialize framework for lattice contributions (is automatically ignored if periodicity = [false, false, false]
  call fmm_framework_init(my_rank, wellsep = 1)

  ! Set up particles
  call special_start(ispecial)

  ! initialize calc force params
  theta2      = theta**2
  mac_select  = mac
  eps2        = eps**2
  force_law   = 3

  ! Loop over all timesteps
  do while (itime < nt)
     itime = itime + 1
     trun  = trun  + dt

     if (my_rank==0 ) then
        ifile=6
           write(ifile,'(//a,i8,(3x,a,f12.3))') &
                ' Timestep ',itime &
                ,' total run time = ',trun 
     endif
     
     call MPI_BARRIER( MPI_COMM_WORLD, ierr)  ! Wait for everyone to catch up
     call timer_start(t_tot)

    call pepc_fields(np_local, npart_total, particles, &
        itime, num_neighbour_boxes, neighbour_boxes, .false., .false.)

     ! Integrator
     call velocities(1,np_local,dt)
     call push(1,np_local,dt)  ! update positions

     ! periodic systems demand periodic boundary conditions
     if (do_periodic) call constrain_periodic(particles(1:np_local)%x(1),particles(1:np_local)%x(2),particles(1:np_local)%x(3),np_local)

     ! timings dump
     call timer_stop(t_tot) ! total loop time without diags

     call timings_LocalOutput(itime)
     call timings_GatherAndOutput(itime)

     flush(6)

  end do

  ! finalize framework for lattice contributions
  call fmm_framework_finalize()

  ! deallocate array space for particles
  call cleanup(my_rank,n_cpu)
  
  ! Time stamp
  if (my_rank==0) call stamp(6,2)
  if (my_rank==0) call stamp(15,2)

  ! Tidy up O/P files
  call closefiles

  ! cleanup of lpepc static data
  call libpepc_finalize()

end program pepcmini
