
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
  use module_interaction_specific, only: theta2
  implicit none

  integer :: i
  real :: trun                     ! total run time including restarts and offset
  integer :: itime, stage, t_flag
  integer, parameter :: t_remesh = t_userdefined_first + 1

  ! Allocate array space for tree
  call pepc_initialize("pepc-v", my_rank, n_cpu, .true.)
  call pepc_read_parameters_from_first_argument()

  ! Set up O/P files
  call openfiles

  ! Time stamp
  if (my_rank==0) call stamp(6,1)
  if (my_rank==0) call stamp(15,1)

  ! Each CPU gets copy of initial data
  call pepc_setup(itime, trun)

  ! Set up particles
  call special_start()

  call pepc_prepare(3)

  ! Loop over all timesteps
  do while (itime < nt)

     itime = itime + 1
     if (my_rank==0 ) write(*,'(a5,i8,a3,i8,a7,i8,a)') 'Step',itime,' of',nt,', using',n,' particles -------------'
     
     ! Runge-Kutta time integration
     do stage = 1,rk_stages

        if (my_rank == 0) write(*,*) 'Time:',trun,'/',te

        call timer_start(t_tot)

        call pepc_particleresults_clear(vortex_particles, np)

        if (theta2 .gt. 0.0) then
            call pepc_grow_and_traverse(np, n, vortex_particles, itime, .false., .true.)
        else
            call direct_sum(np, vortex_particles, vortex_particles%results, my_rank, n_cpu)
        end if

        do i=1,np
          vortex_particles(i)%results%u( 1:3) = vortex_particles(i)%results%u( 1:3) * force_const
          vortex_particles(i)%results%af(1:3) = vortex_particles(i)%results%af(1:3) * force_const
          vortex_particles(i)%results%div     = vortex_particles(i)%results%div * force_const
        end do

        call verify_direct()

        !if (stage == rk_stages)  call dump(itime, trun)

        call push_rk2(stage)

        if (stage .lt. rk_stages) then
            call timer_stop(t_tot)   ! total loop time without diags
            call timings_LocalOutput(itime,stage)
            call timings_GatherAndOutput(itime,stage)
        end if

        flush(6)

        trun  = trun  + dt/rk_stages

     end do

     ! dump, if needed (need to do this before remeshing or we will loose velocity information)
     call dump(itime, trun)

     ! if remeshing is requested at all and if it is time right now, do it!
     if ((rem_freq .gt. 0) .and. (mod(itime,rem_freq)==0)) then

        if (my_rank==0) write(*,'("PEPC-V | ", a)') 'Starting remeshing...'
        call timer_start(t_remesh)

        call remeshing()

        call timer_stop(t_remesh)
        if (my_rank==0) write(*,'("PEPC-V | ", a,f12.8,a)') 'Finished remeshing after ',timer_read(t_remesh),' seconds'
        t_flag = -rk_stages

     else

        t_flag = rk_stages

     end if

     call timer_stop(t_tot)   ! total loop time incl. remeshing if requested
     call timings_LocalOutput(itime,t_flag)
     call timings_GatherAndOutput(itime,t_flag)

     ! Some linear diagnostics
     call linear_diagnostics(itime,trun)
     call divergence_diag(itime,trun)

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