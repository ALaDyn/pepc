
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
!
!  ==============================================================

program pepce

  use physvars
  use utils
  implicit none
  include 'mpif.h'

  ! timing stuff
  real :: t0, t_key, t_domain=0., t_build=0., t_branch=0., t_fill=0., t_props=0.
  real :: t_walk=0., t_walkc=0., t_en, t_force=0.
  real :: t_push, t_diag, t_start_push, t_prefetch=0., Tpon, ttot, t_laser
  integer :: tremain ! remaining wall_clock seconds
  integer :: ierr, lvisit_active, ifile, debug

  ! Initialize the MPI system
  call MPI_INIT(ierr)

  ! Get the id number of the current task
  call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)

  ! Get the number of MPI tasks
  call MPI_COMM_size(MPI_COMM_WORLD, n_cpu, ierr)

  ! Time stamp
  if (my_rank==0) call stamp(6,1)
  if (my_rank==0) call stamp(15,1)


  call openfiles       ! Set up O/P files

  call setup           ! Each CPU gets copy of initial data

!  debug=2

  call pepc_setup(my_rank,n_cpu,npart_total,theta,db_level,np_mult,fetch_mult)  ! Allocate array space for tree

  call param_dump      ! Dump initial data
  call configure       ! Set up particles

  do itime = 1,nt
     trun = trun + dt


     if (my_rank==0 ) then
        do ifile = 6,15,9
           write(ifile,'(//a,i8,(3x,a,f8.2))') &
                ' Timestep ',itime+itime_start &
                ,' total run time = ',trun 

        end do

     endif
     call cputime(t0)
     call MPI_BARRIER( MPI_COMM_WORLD, ierr)  ! Wait for everyone to catch up

     call pepc_fields(npp,x(1:npp),y(1:npp),z(1:npp), &
! ux(1:npp),uy(1:npp),uz(1:npp), &   ! will eventually need velocities for B-fields -> 2nd interface
                 q(1:npp),m(1:npp),work(1:npp),pelabel(1:npp), &
                 ex(1:npp),ey(1:npp),ez(1:npp),pot(1:npp), &
                 mac, theta, eps, force_const, err_f, xl, yl, zl, itime, &
                 t_domain,t_build,t_prefetch,t_walk,t_walkc,t_force)   ! Compute Coulomb fields and pot using lpepc

! TODO: need proper mac selection instead of beam_config

     call error_test(npp)
     stop
     ! Integrator
     call cputime(t_start_push)

     call velocities(1,npp,dt)  ! pure ES, NVT ensembles

     call push_x(1,npp,dt)  ! update positions


     call cputime(t_push)



     call energy_cons(Ukine,Ukini)

     call cputime(t_diag)



     if (my_rank==0 .and. db_level .ge.1) then
        ttot = t_push-t0 ! total loop time without diags

        if (itime ==1 .or. mod(itime,iprot).eq.0) then
           ifile = 6
        else 
           ifile = 15
        endif
        write(ifile,'(//a/)') 'Timing:  Routine   time(s)  percentage'
        write(ifile,'(a20,2f12.3,a1)') 'Domains: ',t_domain,100*t_domain/ttot
        write(ifile,'(a20,2f12.3,a1)') 'Build: ',t_build,100*t_build/ttot

        write(ifile,'(a20,2f12.3,a1)') 'Prefetch: ',t_prefetch,100*t_prefetch/ttot
        write(ifile,'(a20,2f12.3,a1)') 'Walk serial: ',t_walk,100*t_walk/ttot
        write(ifile,'(a20,2f12.3,a1)') 'Walk comm: ',t_walkc,100*t_walkc/ttot
        write(ifile,'(a20,2f12.3,a1)') 'Forces: ',t_force,100*t_force/ttot
        write(ifile,'(a20,2f12.3,a1)') 'Pusher: ',t_push-t_start_push,100*(t_push-t_start_push)/ttot
        write(ifile,'(a20,2f12.3,a1)') 'Diagnostics: ',t_diag-t_push,100*(t_diag-t_push)/ttot

        write(ifile,'(a20,2f12.3,a1)') 'Total: ',ttot,100.
        write(ifile,'(a20,i4,6f12.3)') 'Timing format: ',n_cpu,t_domain,t_build,t_prefetch,t_walk+t_walkc,t_force,ttot

     endif
     call MPI_BARRIER( MPI_COMM_WORLD, ierr)  ! Wait for everyone to catch up

  end do


  call closefiles      ! Tidy up O/P files


  ! Time stamp
  if (my_rank==0) call stamp(6,2)
  if (my_rank==0) call stamp(15,2)
  ! End the MPI run
  call MPI_FINALIZE(ierr)


end program pepce




