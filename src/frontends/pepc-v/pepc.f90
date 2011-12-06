
! ==============================================================
!
!
!                  PEPC-V
!
!    Parallel Efficient Parallel Coulomb-solver: Vortex particles
!
!  ==============================================================

program pepcv

  use physvars
  use manipulate_particles
  use module_pepcfields
  use timings
  use files
  use diagnostics
  use module_setup
  implicit none
  include 'mpif.h'

  integer :: ierr
  real :: trun                     ! total run time including restarts and offset
  integer :: itime, stage
  type(t_calc_force_params) ::cf_par

  ! Allocate array space for tree
  call libpepc_setup("pepc-v", my_rank, n_cpu)

  ! Set up O/P files
  call openfiles

  ! Time stamp
  if (my_rank==0) call stamp(6,1)
  if (my_rank==0) call stamp(15,1)

  ! Each CPU gets copy of initial data
  call pepc_setup(itime, trun)

  ! Set up particles
  call special_start()

  ! initialize calc force params
  cf_par%theta2      = theta**2
  cf_par%mac         = mac
  cf_par%eps2        = eps**2
  cf_par%force_const = force_const
  cf_par%force_law   = 2

  ! Loop over all timesteps
  do while (itime < nt)

     itime = itime + 1
     if (my_rank==0 ) write(*,'(a5,i8,a3,i8,a7,i8,a)') 'Step',itime,' of',nt,', using',n,' particles -------------'
     
     ! Runge-Kutta time integration
     do stage = 1,rk_stages

        if (my_rank == 0) write(*,*) 'Time:',trun,'/',te

        call MPI_BARRIER( MPI_COMM_WORLD, ierr)  ! Wait for everyone to catch up
        call timer_start(t_tot)

        call pepc_fields(np, n, vortex_particles, cf_par, itime, &
                         1, [0, 0, 0], .false., .true.)

        call push_rk2(stage)

        call timer_stop(t_tot)   ! total loop time without diags

        call timings_LocalOutput(itime,stage)
        call timings_GatherAndOutput(itime,stage)

        flush(6)

        trun  = trun  + dt/rk_stages

     end do

     ! Some linear diagnostics
     call linear_diagnostics(itime,trun)

     ! dump, if needed
     call dump(itime, trun)

  end do

  ! deallocate array space for particles
  call cleanup(my_rank,n_cpu)
  
  ! Time stamp
  if (my_rank==0) call stamp(6,2)
  if (my_rank==0) call stamp(15,2)

  ! Tidy up O/P files
  call closefiles

  ! cleanup of lpepc static data
  call libpepc_finalize()

end program pepcv
