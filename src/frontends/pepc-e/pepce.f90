
! ==============================================================
!
!
!                  PEPC-E
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

program pepce

  use treetypes
  use physvars
  use module_pepc
  use module_pepc_wrappers
  use module_mirror_boxes, only : do_periodic, constrain_periodic
  use particle_pusher
  use benchmarking
  use module_timings
  use files
  use energies
  use module_diagnostics
  use module_pepc_wrappers
  use module_calc_force, only : theta2, eps2, mac_select, force_law
  implicit none

  integer :: ierr, ifile

  ! Allocate array space for tree
  call pepc_initialize("pepc-e", my_rank, n_cpu, .true.)

  call benchmark_pre

  ! Set up O/P files
  call openfiles

  ! Time stamp
  if (my_rank==0) call stamp(6,1)
  if (my_rank==0) call stamp(15,1)

  ! Each CPU gets copy of initial data
  call pepc_setup()

  ! Set up particles
  call configure()

  ! initial particle output
  ! no initial checkpoint since this would override the current checkpoint if in resume-mode
  call write_particles(.false.)

  call benchmark_inner

  flush(6)

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
     
     call timer_start(t_tot)

     call pepc_fields_coulomb_wrapper(np_local,npart_total,x(1:np_local),y(1:np_local),z(1:np_local), &
                  q(1:np_local),work(1:np_local),pelabel(1:np_local), &
                  ex(1:np_local),ey(1:np_local),ez(1:np_local),pot(1:np_local), &
                      itime, .false., .false., force_const)

     if (itime == nt) then
        call gather_particle_diag()
        if (my_rank == 0) call benchmarking_dump_diagnostics()
     end if

     ! Integrator
     call velocities(1,np_local,dt)  ! pure ES, NVT ensembles
     call push(1,np_local,dt)  ! update positions

     ! periodic systems demand periodic boundary conditions
     if (do_periodic) call constrain_periodic(x(1:np_local),y(1:np_local),z(1:np_local),np_local)

     call energy_cons(Ukine,Ukini)

     ! periodic particle dump
     call write_particles(.true.)

     ! timings dump
     call timer_stop(t_tot) ! total loop time without diags

     call timings_LocalOutput(itime)
     call timings_GatherAndOutput(itime)

     flush(6)

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

end program pepce
