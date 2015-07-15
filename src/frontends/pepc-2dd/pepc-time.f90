! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
!
! Copyright (C) 2002-2014 Juelich Supercomputing Centre,
!                         Forschungszentrum Juelich GmbH,
!                         Germany
!
! PEPC is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! PEPC is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with PEPC.  If not, see <http://www.gnu.org/licenses/>.
!

program pepc

  ! pepc modules
  use module_pepc
  use module_pepc_types
  use module_timings
  use module_debug
  use newton_krylov
  use tintegration
  use zufall
  use encap
  use field_helper, only: compute_field,setup_field_grid,pepc_setup,prepare_grid,field_grid_update_grid

  ! frontend helper routines
  use helper
  use module_mpi_io
  use module_shortcut
  implicit none

  ! timing variables
  real*8    :: timer(5)
  real*8    :: t1, t2,t,toll
  type(t_particle), allocatable              :: q(:),r(:)
  type(physics_pars_t)                       :: physics_pars
  type(field_grid_t)                         :: field_grid
  type(pepc_pars_t)                          :: pepc_pars
  real(kind_particle), allocatable           :: x_out(:),x_0(:)
  real(kind_particle)                        :: P(3),static_gmres(3),gamma,errnk,m,gmerrtol,&
                                                pmin,pmax,emax,xmin,new_extent(2),new_offset(2)
  integer(kind_particle)                     :: iter,rc,ip,jp,nn,maxit,errmsg,gmkmax

  ! control variable
  logical :: doDiag,dogrid,dorestart,dointerp, explicit =.false.

  !!! initialize pepc library and MPI
  gmerrtol = 1.0e-8
  gmkmax   = 100

  call pepc_initialize("pepc-2dd", my_rank, n_ranks, .true., comm=communicator)
!  call pepc_initialize("pepc-2dd", my_rank, n_ranks, .true.)


  root = my_rank.eq.0
  timer(1) = get_time()
  call set_parameter()

    select case (initial_setup)
          case (1)  ! Load State
              call read_restart_2d(particles,"restart_0_")
          case (2)  !  Default initial setup - 2D random particles with thermal velocity
              call init_particles(particles)
          case (3)  !  2D random particles in a square
              call init_particles_square(particles)
          case (4)  ! 2D random particles in a ring
              call init_particles_ring_2D(particles)
          case (5)  ! wire random distribution
              call init_particles_wire(particles)
          case (6)  ! stripe random distribution
              call init_particles_stripe(particles)
          case (7)  ! disc random distribution
              call init_particles_disc(particles)
          case (8)  ! 2D Two wires distribution
              call init_particles_two_wires(particles)
          case (9)  ! Default initial setup - 3D random particles with thermal velocity
              call init_particles3D(particles)
          case (10)  ! 3D Solenoid
              call init_particles_solenoid(particles)
          case (11)  ! two bodies problem - grav
              call two_bodies_grav(particles)
          case (12)  ! 2D Landau Damping
              call init_landau_damping2d(particles)
          case (13)  ! 1D Landau Damping
              call init_landau_damping1d(particles)
          case (14)  ! 1D Langmuir waves
              call init_langmuir_waves(particles)


          case default

              call init_particles(particles)

  end select

  if (.not. allocated(x_out) )        allocate(x_out(12*np)  , stat=rc )
  if (.not. allocated(x_0) )        allocate(x_0(12*np)  , stat=rc )

  x_out = 0.0_8
  x_0   = 0.0_8

  call pepc_setup(pepc_pars)
  call setup_field_grid(field_grid, pepc_pars%pepc_comm)

  timer(2) = get_time()
  t        = 0

  if(root) write(*,'(a,es12.4)') " === init time [s]: ", timer(2) - timer(1)

  if (ischeme=="leapfrog") explicit = .true.

  new_extent = extent
  new_offset = offset

  do step=0, nt
    if(root) then
      write(*,*) " "
      write(*,'(a,i12)')    " ====== computing step  :", step
      write(*,'(a,es12.4)') " ====== simulation time :", step * dt
    end if

    timer(3) = get_time()

    call pepc_particleresults_clear(particles)
    timer(1) = get_time()
    t1 = get_time()

    call pepc_grow_tree(particles)

    np = size(particles, kind=kind(np))
    t2 = get_time()
    if(root) write(*,'(a,es12.4)') " ====== tree grow time  :", t2-t1
    t1 = get_time()
    call pepc_traverse_tree(particles)
    t2 = get_time()
    if(root) write(*,'(a,es12.4)') " ====== tree walk time  :", t2-t1
    call pepc_timber_tree()
    !call pepc_restore_particles(np, particles)

    timer(5) = get_time()
    t        = t + timer(5) - timer(1)



    if (explicit) then
        call leapfrog(particles)
    else
!        call f(np,dt,particles,x_0,merda)
!        call gmres(np,dt,merda,inv_m,particles,x_0,gmerrtol,gmkmax,x_out,iter)

        pmax  = 0.0_8
        emax  = 0.0_8
        xmin  = 100.0_8
        pmin  = 100.0_8

        do ip = 1,np

            pmax  = max(pmax, sqrt( dot_product( particles(ip)%data%m*particles(ip)%data%v , particles(ip)%data%m*particles(ip)%data%v ) ) )
            emax  = max(emax, sqrt( dot_product( particles(ip)%results%E , particles(ip)%results%E ) ) )
            xmin  = min(xmin, sqrt( dot_product( particles(ip)%x , particles(ip)%x ) )  )
            pmin  = min(xmin, sqrt( dot_product(particles(ip)%data%m*particles(ip)%data%v , particles(ip)%data%m*particles(ip)%data%v ) )  )

        enddo

!        write(*,*) "time condition === ", dt, pmin/emax, xmin/pmax

        call nsolgm(np,dt,particles,midpoint,x_out,errmsg,iter,errnk,static_gmres)
        if(root) then
            write(*,'(a,i12)')    " ====== Newton - Krylov Error Message  :", errmsg
            write(*,'(a,i12)')    " ====== Newton - Krylov Iterations     :", iter
            write(*,'(a,es12.4)') " ====== Newton - Krylov Final Residual :", errnk
            write(*,'(a,es12.4)') " ====== GMRES Statistics Min Iteration :", static_gmres(1)
            write(*,'(a,es12.4)') " ====== GMRES Statistics Ave Iteration :", static_gmres(2)
            write(*,'(a,es12.4)') " ====== GMRES Statistics Max Iteration :", static_gmres(3)

        endif

        do ip = 1,np
            jp                      = (ip-1)*12
            particles(ip)%x         = x_out(jp+1:jp+3)
!            particles(ip)%x(2)      = zero
            P                       = x_out(jp+4:jp+6)
            m                       = particles(ip)%data%m
            gamma                   = 1.0_8!sqrt( 1.0_8 + dot_product(P/m,P/m) )
            particles(ip)%data%v    = P/m/gamma !+ particles(ip)%data%q*particles(ip)%results%A
!            particles(ip)%data%v(2) = zero

        enddo
    endif

    dointerp = .false.!mod( step , diag_interval) .eq. 0

    if (dointerp) then
        call prepare_grid(particles,new_extent,new_offset)
!        write(*,*) "offset/extent: ", field_grid%offset,field_grid%extent
        dogrid      = new_extent(1) .ne. extent(1) .or. new_extent(2) .ne. extent(2)
        dogrid      = dogrid.or.( new_offset(1) .ne. offset(1) .or. new_offset(2) .ne. offset(2) )
        if ( dogrid  ) call field_grid_update_grid( field_grid,pepc_pars%pepc_comm,new_extent,new_offset )
!        write(*,*) "new offset/extent: ", field_grid%offset,field_grid%extent

        call compute_field(pepc_pars, field_grid, particles)
        call write_field_on_grid_ascii(field_grid,int(step, kind = kind_particle))
        call write_particle_ascii(particles,int(step, kind = kind_particle) )
        call electric_energy(field_grid,real(step, kind = kind_particle)*dt)
        call energy(particles,real(step, kind = kind_particle)*dt)

    endif

!    call write_particle_ascii(particles,step)

!    call write_ascii(np,particles(:)%x(1),'prova.dat',pepc_pars%pepc_comm)
    dorestart = .false.!mod( step , restart_step) .eq. 0
    if (dorestart)   call write_restart_2d(particles,int(step, kind=kind_particle))

    timer(4) = get_time()
    if(root) write(*,'(a,es12.4)') " == time in step [s]                              : ", timer(4) - timer(2)

    call timings_GatherAndOutput(step, 0)


  end do


  deallocate(particles)
  !deallocate(particles0)

  timer(5) = get_time()

  if(root) then
    write(*,*)            " "
    write(*,'(a)')        " ===== finished pepc simulation"
    write(*,'(a,es12.4)') " ===== total run time [s]: ",timer(5) - timer(2)
  end if

  !!! cleanup pepc and MPI
  call pepc_finalize()

end program pepc

