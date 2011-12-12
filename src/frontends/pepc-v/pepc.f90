
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
  use module_pepc
  use module_timings
  use files
  use diagnostics
  implicit none

  integer :: ierr, i
  real :: trun                     ! total run time including restarts and offset
  integer :: itime, stage

  ! Allocate array space for tree
  call pepc_initialize("pepc-v", my_rank, n_cpu, .true.)

  ! Set up O/P files
  call openfiles

  ! Time stamp
  if (my_rank==0) call stamp(6,1)
  if (my_rank==0) call stamp(15,1)

  ! Each CPU gets copy of initial data
  call pepc_setup(itime, trun)

  ! Set up particles
  call special_start()

  call pepc_prepare()

  ! Loop over all timesteps
  do while (itime < nt)

     itime = itime + 1
     if (my_rank==0 ) write(*,'(a5,i8,a3,i8,a7,i8,a)') 'Step',itime,' of',nt,', using',n,' particles -------------'
     
     ! Runge-Kutta time integration
     do stage = 1,rk_stages

        if (my_rank == 0) write(*,*) 'Time:',trun,'/',te

        call timer_start(t_tot)

        call pepc_grow_and_traverse(np, n, vortex_particles, itime, .true., .false., .true.)

        ! TODO: check, if this is correct (moved force_const out of calc_force)

        do i=1,np
          vortex_particles(i)%results%u( 1:3) = vortex_particles(i)%results%u( 1:3) * force_const
          vortex_particles(i)%results%af(1:3) = vortex_particles(i)%results%af(1:3) * force_const
        end do

        call push_rk2(stage)

        call timer_stop(t_tot)   ! total loop time without diags

        call timings_LocalOutput(itime,stage)
        call timings_GatherAndOutput(itime,stage)

        flush(6)

        trun  = trun  + dt/rk_stages

     end do

     call remeshing(itime)

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
  call pepc_finalize()

end program pepcv
