
! ==============================================================
!
!
!                  PEPC-MW
!
!    Parallel Efficient Parallel Coulomb-solver
!
!   $ Revision $
!
!   Driver code for Coulomb-solver library lpepc
!
!  ==============================================================

program pepc

  use treetypes
  use physvars
  use benchmarking
  use module_timings
  use module_pepc
  use module_pepc_wrappers
  use module_mirror_boxes, only : do_periodic, constrain_periodic
  use module_laser
  use module_pusher
  use module_io
  use module_fields
  use module_acf
  use module_diagnostics
  use module_workflow
  use module_units
  use module_param_dump
  use module_treediags
  use module_vtk
  use module_directsum
  use module_energy
  use module_particle_setup
  use module_calc_force, only : theta2, eps2, mac_select, force_law
  implicit none
  include 'mpif.h'

  integer :: vtk_step

  integer :: ifile
  type(acf) :: momentum_acf

  ! Initialization of signal handler - deactivated atm since it outputs the call stack for every mpi rank which results in a very messy output
  !call InitSignalHandler()

  ! Allocate array space for tree
  call pepc_initialize("pepc-mw", my_rank, n_cpu, .true.)

  call benchmark_pre

  ! Set up O/P files
  call openfiles

  ! Time stamp
  if (my_rank==0) call stamp(6,1)
  if (my_rank==0) call stamp(15,1)

  ! Each CPU gets copy of initial data
  call pepc_setup()

  ! Set up particles
  call particle_setup(ispecial)

  ! parameter output
  if (my_rank == 0) then
    call PrintPhysicalParameters(6)
    call PrintPhysicalParameters(24)
  endif

  ! initial particle output
  ! no initial checkpoint since this would override the current checkpoint if in resume-mode
  call write_particles(.false.)
  if (( idump .gt. 0 ) .and. ((ispecial==9).or.(ispecial==10).or.(ispecial==11))) call sum_radial(itime)

  ! time-dependent setup stuff
  call workflow(my_rank, 0, 0._8, dt)

  call momentum_acf%initialize(nt-momentum_acf_from_timestep, dt*unit_t0_in_fs, my_rank, n_cpu, MPI_COMM_PEPC)

  call benchmark_inner

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
           write(ifile,'(//a)') "==================================================================="
           write(ifile,'(//a,i8,3x,a,f12.3,"  (",f12.3," fs)")') &
                ' Timestep ',itime &
                ,' total run time = ',trun, trun*unit_t0_in_fs
     endif
     
     ! time-dependent setup stuff
     call workflow(my_rank, itime, trun, dt)

     ! dump trajectory
     if (my_rank == 0 .and. itime == nt) call dump_trajectory()

     call timer_start(t_tot)

     ! laser propagation according to beam_config
     call laser_update()
     call PrintLaserParameters()

     call pepc_fields_coulomb_wrapper(np_local,npart_total,x(1:np_local),y(1:np_local),z(1:np_local), &
                 q(1:np_local),work(1:np_local),pelabel(1:np_local), &
                  ex(1:np_local),ey(1:np_local),ez(1:np_local),pot(1:np_local), &
                      itime, treediags, .false., force_const)

  !   call verifydirect(x, y, z, q, ex, ey, ez, pot, np_local, [1, 2, np_local-1, np_local], &
  !                     0, my_rank, n_cpu, MPI_COMM_PEPC)

     ! output of tree diagnostics
     if (treediags) then
       if (itime == 1) then
         vtk_step = VTK_STEP_FIRST
       else if (itime == nt) then
         vtk_step = VTK_STEP_LAST
       else
         vtk_step = VTK_STEP_NORMAL
       endif

       call write_branches_to_vtk(itime,   trun*unit_t0_in_fs, vtk_step)
       call write_spacecurve_to_vtk(itime, trun*unit_t0_in_fs, vtk_step)
     endif

     ! add any external forces (laser field etc)
     call force_laser(1, np_local)

     if (itime == nt) call gather_particle_diag()

     ! Velocity and position update - explicit schemes only
     call integrator(1, np_local, integrator_scheme)

     ! periodic systems demand periodic boundary conditions
     if (do_periodic) call constrain_periodic(x(1:np_local),y(1:np_local),z(1:np_local),np_local)

     call energies(Ukine,Ukini)

     ! periodic particle dump
     call write_particles(.true.)

     if ( idump .gt. 0 ) then
       if ( mod(itime, idump ) .eq. 0) then
         if ((ispecial==9).or.(ispecial==10).or.(ispecial==11)) call sum_radial(itime)

         call field_dump(itime)
       end if
     endif

     !> special diagnostics for [J. Phys. A: Math. Theor 42 (2009), 214048] Th. Raitza et al: Collision frequency of electrons in laser excited small clusters
     if (workflow_setup == 3) call cluster_diagnostics(itime, trun*unit_t0_in_fs, momentum_acf)

     ! timings dump
     call timer_stop(t_tot) ! total loop time without diags

     call timings_LocalOutput(itime)
     call timings_GatherAndOutput(itime)

     call flushfiles()

  end do

  call benchmark_post

  ! final particle dump
  call write_particles(.true.)

  ! deallocate array space for particles
  call cleanup(my_rank,n_cpu)
  
  ! Time stamp
  if (my_rank==0) call stamp(6,2)
  if (my_rank==0) call stamp(15,2)

  ! Tidy up O/P files
  call closefiles

  call benchmark_end

  ! cleanup of lpepc static data
  call pepc_finalize()

end program pepc
