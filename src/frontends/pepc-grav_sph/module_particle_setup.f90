! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2017 Juelich Supercomputing Centre, 
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>  Encapsulates anything necessary for particle setup,
!>  i.e. the old configure() and special_start() routines
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_particle_setup
  use physvars, only: &
       particles, &
       n_cpu, &
       my_rank, &
       np_local, &
       npart_total, &
       idim, &
       do_sph

!  use module_files, only: &
!       read_particles


  implicit none
  save
  private

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  public particle_setup

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains



  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>
  !>
  !>
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine particle_setup(iconf)

    use files, only: &
         read_particles

    implicit none
    include 'mpif.h'

    integer, intent(in) :: iconf
    integer :: fences(-1:n_cpu-1)
    integer :: ierr

    call MPI_SCAN(np_local, fences(my_rank), 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER(MPI_IN_PLACE, 1, MPI_INTEGER, fences(0), 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    fences(-1) = 0

    select case (iconf)
    case (-1)
       if (my_rank == 0) write(*,*) "Using special start... case -1 (reading data file)"
       call read_particles()

    case(1)
       if (my_rank == 0) write(*,*) "Using special start... case 1 (3D homogeneous distribution)"
       call init_generic(fences)
       call particle_setup_homogeneous(fences, 3)

    case(2)
       if (my_rank == 0) write(*,*) "Using special start... case 2 (one sphere benchmark)"
       call init_generic(fences)
       call particle_setup_sphere(fences)

    case(3)
       if (my_rank == 0) write(*,*) "Using special start... case 3 (two sphere benchmark)"
!       call particle_setup_doublesphere(fences)
!       call rescale_coordinates_spherical()
!       call init_generic(fences)
!       call init_fields_zero()

    case(342)
       if (my_rank == 0) write(*,*) "Using special start... case 342 (42 spheres benchmark)"
!       call particle_setup_manyspheres(fences)
!       call rescale_coordinates_spherical()
!       call init_generic(fences)
!       call init_fields_zero()

    case(4)
       if (my_rank == 0) write(*,*) "Using special start... case 4: Plummer distribution (core cut)"
!       call particle_setup_plummer(fences)
!       call rescale_coordinates_spherical()
!       call init_generic(fences)
!       call init_fields_zero()

    case(5)
       if (my_rank == 0) write(*,*) "Using special start... case 5: 2D disc"
!       call particle_setup_flatdisc(fences)
!       z_plasma = 0
!       call rescale_coordinates_spherical()
!       call init_generic(fences)
!       uz(1:np_local) = 0.
!       call init_fields_zero()

    case(6)
       if (my_rank == 0) write(*,*) "Using special start... case 6 (2D homogeneous distribution)"
       call init_generic(fences)
       call particle_setup_homogeneous(fences, 2)

    case(7)
       if (my_rank == 0) write(*,*) "Using special start... case 7 (3D homogeneous distribution xyz periodic)"
!       call init_generic(fences)
!       call particle_setup_homogeneous(fences, 3)

 
    case(8)
       if (my_rank == 0) write(*,*) "Using special start... case 8 (fast 3D homogeneous distribution)"
       call init_generic(fences)
       call particle_setup_homogeneous_fast(fences, 3)

    case(12)
       if (my_rank == 0) write(*,*) "Using special start... case 12 (1D wave)"
       call init_generic(fences)
       call particle_setup_1d_wave(fences)

    case(13)
       if (my_rank == 0) write(*,*) "Using special start... case 13 (1st 1D shock test)"
       call init_generic(fences)
       call particle_setup_1d_shock_1(fences)

    case(14)
       if (my_rank == 0) write(*,*) "Using special start... case 14 (2nd 1D shock test)"
       call init_generic(fences)
       call particle_setup_1d_shock_2(fences)

    case(15)
       if (my_rank == 0) write(*,*) "Using special start... case 15 (3rd 1D shock test)"
       call init_generic(fences)
       call particle_setup_1d_shock_3(fences)

    case(16)
       if (my_rank == 0) write(*,*) "Using special start... case 16 (3D sphere collapse (grav only) )"
       call init_generic(fences)
       call particle_setup_3d_sphere_grid(fences)
       do_sph = .false.

    case(17)
       if (my_rank == 0) write(*,*) "Using special start... case 16 (3D sphere collapse (grav + sph) )"
       call init_generic(fences)
       call particle_setup_3d_sphere_grid(fences)

    end select

  end subroutine particle_setup




  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>
  !>
  !>
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine particle_setup_homogeneous(fences, dim)

    implicit none
    
    integer, intent(in) :: fences(-1:n_cpu-1)
    integer, intent(in) :: dim

    integer :: mpi_cnt, p, i
    real*8, dimension(3) :: tmp
    
    tmp = [ 0._8, 0._8, 0._8 ]

    idim = dim
    
    do mpi_cnt = 0, n_cpu-1
       do p = 1, (fences(mpi_cnt) - fences(mpi_cnt-1))
          
          do i = 1, dim
             tmp(i) = par_rand()
          end do

          if ( my_rank .eq. mpi_cnt .and. p .le. np_local ) then
             particles(p)%x = tmp
!             particles(p)%data%v_minus_half = [ 0._8, 0._8, 0._8 ]
          end if

       end do
    end do
    
  end subroutine particle_setup_homogeneous


  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! >
  ! >
  ! >
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine particle_setup_sphere(fences)

    use physvars, only: &
         initialized_v_minus_half, &
         thermal_constant, &
         kappa, &
         dt

    use module_interaction_specific_types, only: &
         num_neighbour_particles

    use module_mirror_boxes, only: &
         periodicity

    implicit none
    integer, intent(in) :: fences(-1:n_cpu-1)
    integer :: mpi_cnt, p
    real*8 :: xt, yt, zt
    
    do mpi_cnt = 0, n_cpu-1
       do p = 1, (fences(mpi_cnt) - fences(mpi_cnt-1))
          
          xt = 1.0_8
          yt = 1.0_8
          zt = 1.0_8
          
          do while ( (xt*xt + yt*yt + zt*zt) > 0.5_8**2)
             xt = par_rand() - 0.5_8
             yt = par_rand() - 0.5_8
             zt = par_rand() - 0.5_8
          end do

          if ( my_rank .eq. mpi_cnt .and. p .le. np_local ) then
             particles(p)%x = [ xt, yt, zt ]
             particles(p)%data%q = 1._8/real(npart_total, 8)
             particles(p)%data%temperature = 100.
                      
             
          end if
       end do
    end do

    initialized_v_minus_half = .true.

    num_neighbour_particles = 40


    ! set number of dimension. important for factor for sph kernel
    idim = 3
    
    ! set periodicity
    periodicity = [.false., .false., .false.]
    
    dt = 0.001

    kappa = 1.66666666

    thermal_constant = 1. /(kappa * 100.) ! soundspeed^2 /( kappa * temperature)


  end subroutine particle_setup_sphere


  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! > A simple 1D sound wave in a periodic medium for testing sph.
  ! > Wavenumber is currently k=1, soundspeed=1.
  ! >
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine particle_setup_1d_wave(fences)

    use physvars, only: &
         dt, &
         n_cpu, &
         thermal_constant, &
         initialized_v_minus_half

    use module_mirror_boxes, only: &
         periodicity, &
         t_lattice_1

    use module_interaction_specific_types, only: &
         num_neighbour_particles

    
    implicit none
    include 'mpif.h'

    integer, intent(in) :: fences(-1:n_cpu-1)

    real*8 :: medium_temperature
    integer :: part_counter
    integer :: part_before_me
    integer :: part_including_mine
    real*8 :: offset
    real*8 :: left_y
    real*8 :: right_y
    real*8 :: area1, area2
    real*8 :: sound_speed
    real*8 :: setup_rho0, setup_rho1
    real*8 :: dx
    integer :: i
    integer :: actual_particle
    real*8 :: actual_x
    integer :: all_part
    real*8 :: omega_t_for_half_velocity
    real*8 :: const_pi = acos(-1.0)

    ! set number of dimension. important for factor for sph kernel
    idim = 1

    ! set periodicity
    periodicity = [.true., .false., .false.]

    ! if periodic set boxlength
    t_lattice_1 = [1, 0, 0]
    
    sound_speed = 1._8
      
    ! temperature of the medium the wave moves in
    medium_temperature = 100._8
    thermal_constant = sound_speed**2 /medium_temperature
    
    ! n_nn = 6 does not work! try shepard correction rho = (sum m_i *W )/(sum W)
    ! num_neighbour_particles = 8

    ! timestep length
    dt = 0.0001
    !* sqrt(thermal_constant * medium_temperature /10.)

    all_part = fences(n_cpu-1)
    part_before_me = fences(my_rank-1)
    part_including_mine = fences(my_rank)

 
    setup_rho0 = 1000._8
    setup_rho1 = 1._8
    
    dx = 0.00000001_8
    
    area1 = 0._8
    
    do i= 0, int(t_lattice_1(1)/dx) -1
       left_y  = setup_rho0 + setup_rho1 * sin(dx * real(i,   8) /t_lattice_1(1) * 2._8 * const_pi)
       right_y = setup_rho0 + setup_rho1 * sin(dx * real(i+1, 8) /t_lattice_1(1) * 2._8 * const_pi)
         
       area1 = area1 + (left_y + right_y ) /2._8 * dx
    end do

 
    offset = area1 /(all_part*2)
 
    write (*,*) 'area1:' , area1, 'offset:', offset
 
    actual_x = 0._8
    area2 = 0._8


    ! omega*t for the velocity for t = 0+1/2 (leap frog)
    omega_t_for_half_velocity = const_pi * sound_speed* dt /t_lattice_1(1)
    
    do part_counter = 1, all_part
    
       do while(area2 < area1/real(all_part,8)*real(part_counter-1,8)+offset )
          left_y  = setup_rho0 + setup_rho1 * sin( actual_x      /t_lattice_1(1) * 2._8 * const_pi)
          right_y = setup_rho0 + setup_rho1 * sin((actual_x +dx) /t_lattice_1(1) * 2._8 * const_pi)
          
          area2 = area2 +(left_y + right_y )/2._8*dx
          actual_x = actual_x + dx
       end do
    
       if( (part_counter > part_before_me) .and. (part_counter <= part_including_mine) ) then
          
          actual_particle = part_counter - part_before_me
          particles(actual_particle)%x = [ actual_x, 0._8, 0._8 ]
          
          ! v not necessary for integration. no need to initialize it (besides the initial output).

          ! v is initialized and v_minus_half computed before entering the loop over timesteps
          particles(actual_particle)%data%v = [ -sound_speed * setup_rho1/setup_rho0 * sin( particles(actual_particle)%x(1) * 2._8* const_pi / t_lattice_1(1) ) , 0._8, 0._8 ]

          ! v_minus_half is the velocity used for the leapfrog.
          ! TODO correct this v_minus_half
!          particles(actual_particle)%data%v_minus_half = [ -sound_speed * setup_rho1/setup_rho0 * sin( particles(actual_particle)%x(1) * 2._8* const_pi / t_lattice_1(1) + omega_t_for_half_velocity ) , 0._8, 0._8 ]

!          initialized_v_minus_half = .true.

          
          particles(actual_particle)%data%temperature =  medium_temperature
       
       end if
       
    end do
    
  end subroutine particle_setup_1d_wave
  

  
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! > 1D shock Problem 1, see Springel, V. 2010, ARA&A, 48, 391
  ! >
  ! > Initial conditions given by density (rho), velocity (v) and pressure (P)
  ! > on the left and right side of a discontinuity:
  ! > | rho_L | v_L | P_L | rho_R | v_R | P_R |
  ! > | 1.0   | 0.0 | 1.0 | 0.125 | 0.0 | 0.1 |
  ! >
  ! > Particles are distributed from -0.5 to 1.5 to provide a stable medium
  ! > outside of the test box from 0.0 to 1.0.
  ! >
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine particle_setup_1d_shock_1(fences)
    
    use physvars, only: &
         dt, &
         n_cpu, &
         my_rank, &
         thermal_constant, &
         nt, &
         initialized_v_minus_half, &
         do_gravity
    
    use module_mirror_boxes, only: &
         periodicity
    
    use module_interaction_specific_types, only: &
         num_neighbour_particles, &
         PARTICLE_TYPE_FIXED
    
    
    implicit none
    include 'mpif.h'
    
    integer, intent(in) :: fences(-1:n_cpu-1)

    integer :: part_counter
    integer :: part_before_me
    integer :: part_including_mine
    real*8 :: offset
    real*8 :: left_y
    real*8 :: right_y
    real*8 :: area1, area2
    real*8 :: setup_rho0, setup_rho1
    real*8 :: dx
    integer :: actual_particle
    real*8 :: actual_x
    integer :: all_part

    if( my_rank .eq. 0 ) then
       write(*,*) "Creating setup for the 1st 1D shock test."
       write(*,*) "Using some particles as boundary particles."
       write(*,*) "Use only 0.0 < x< 1.0 for comparison."
    end if

    ! set number of dimension. important for factor for sph kernel
    idim = 1

    ! switch gravity off
    do_gravity = .false.
    
    ! set periodicity
    periodicity = [.false., .false., .false.]
    
    all_part = fences(n_cpu-1)
    part_before_me = fences(my_rank-1)
    part_including_mine = fences(my_rank)
    
    
    setup_rho0 = 1._8
    setup_rho1 = 0.125_8
    
    dx = 0.00000001
    
    ! fill from -0.5 to 1.5 to avoid outflow out of box
    area1 = setup_rho0 * 1. + setup_rho1 * 1.
    
    
    offset = area1/(all_part*2)
    
    actual_x = -0.5_8
    area2 = 0._8
    
    
    do part_counter = 1, all_part

       do while(area2 < area1/real(all_part,8)*real(part_counter-1,8)+offset )
          if(actual_x <= 0.5) then
             left_y = setup_rho0
          else
             left_y = setup_rho1
          end if
          
          if( (actual_x + dx ) <= 0.5) then
             right_y = setup_rho0
          else
             right_y = setup_rho1
          end if
          
          area2 = area2 +(left_y + right_y )/2._8*dx
          actual_x = actual_x + dx
       end do
       
       
       if( (part_counter > part_before_me) .and. (part_counter <= part_including_mine) ) then
          
          actual_particle = part_counter - part_before_me
          particles(actual_particle)%x = [ actual_x, 0._8, 0._8 ]
          
          particles(actual_particle)%data%q = offset*2
          
          
          particles(actual_particle)%data%v = [ 0._8, 0._8, 0._8 ]
          ! this is the best possible initialisation for v_minus_half 
          particles(actual_particle)%data%v_minus_half = [ 0._8, 0._8, 0._8 ]

          if(actual_x <= 0.5) then
             particles(actual_particle)%data%temperature = 1.0 /(thermal_constant * 1.0)
          else
             particles(actual_particle)%data%temperature = 0.1 /(thermal_constant * 0.125)
          end if
          
          if(part_counter < 2* num_neighbour_particles .or. (all_part - part_counter) < 2 * num_neighbour_particles) then
             particles(actual_particle)%data%type = ibset(particles(actual_particle)%data%type, PARTICLE_TYPE_FIXED)
          end if
          
       end if
       
    end do

    initialized_v_minus_half = .true.
    
    ! timestep length
    dt = 0.0001

    ! num timesteps for comparison with analytical results
    nt = 2500

    ! number of neighbour particles known to produce reasonable results
    num_neighbour_particles = 8
    
  end subroutine particle_setup_1d_shock_1
  
  
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! > 1D shock Problem 2, see Springel, V. 2010, ARA&A, 48, 391
  ! >
  ! > Initial conditions given by density (rho), velocity (v) and pressure (P)
  ! > on the left and right side of a discontinuity:
  ! > | rho_L | v_L  | P_L | rho_R | v_R | P_R |
  ! > | 1.0   | −2.0 | 0.4 | 1.0   | 2.0 | 0.4 |
  ! >
  ! > Particles are distributed from -0.5 to 1.5 to provide a stable medium
  ! > outside of the test box from 0.0 to 1.0.
  ! >
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine particle_setup_1d_shock_2(fences)
    
    use physvars, only: &
         dt, &
         nt, &
         n_cpu, &
         my_rank, &
         thermal_constant, &
         initialized_v_minus_half, &
         do_gravity
    
    use module_mirror_boxes, only: &
         periodicity
    
    use module_interaction_specific_types, only: &
         num_neighbour_particles
    
    
    implicit none
    include 'mpif.h'
    
    integer, intent(in) :: fences(-1:n_cpu-1)

    integer :: part_before_me
    integer :: part_including_mine
    real*8 :: offset
    real*8 :: dx
    integer :: actual_particle
    integer :: all_part

    
    if( my_rank .eq. 0 ) then
       write(*,*) "Creating setup for the 2nd 1D shock test."
       write(*,*) "Using some particles as boundary particles."
       write(*,*) "Use only 0.0 < x< 1.0 for comparison."
    end if

    ! set number of dimension. important for factor for sph kernel
    idim = 1
    
    ! switch gravity off
    do_gravity = .false.

    ! set periodicity
    periodicity = [.false., .false., .false.]
    
    all_part = fences(n_cpu-1)
    part_before_me = fences(my_rank-1)
    part_including_mine = fences(my_rank)

    dx = 2./real(all_part)
    offset = dx/2. + dx * part_before_me - 0.5
    
    do actual_particle = 1, np_local
       
       particles(actual_particle)%x = [ offset + dx * (actual_particle-1), 0._8, 0._8 ]
       particles(actual_particle)%data%temperature = 0.4 /(thermal_constant * 1.0)
       particles(actual_particle)%data%q = 2./real(all_part)
       
       if(particles(actual_particle)%x(1) <= 0.5 ) then
          particles(actual_particle)%data%v = [ -2._8, 0._8, 0._8 ]
          ! v_minus_half will be computed from v and the force in pepc.f90
          particles(actual_particle)%data%v_minus_half = [ -2._8, 0._8, 0._8 ]
       else
          particles(actual_particle)%data%v = [ 2._8, 0._8, 0._8 ]
          particles(actual_particle)%data%v_minus_half = [ 2._8, 0._8, 0._8 ]
       end if
       
    end do

    initialized_v_minus_half = .true.
    
    ! timestep length
    dt = 0.0001

    ! num timesteps for comparison with analytical results
    nt = 2500

    ! number of neighbour particles known to produce reasonable results
    num_neighbour_particles = 8

  end subroutine particle_setup_1d_shock_2


  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! > 1D shock Problem 3, see Springel, V. 2010, ARA&A, 48, 391
  ! >
  ! > Initial conditions given by density (rho), velocity (v) and pressure (P)
  ! > on the left and right side of a discontinuity:
  ! > | rho_L | v_L | P_L  | rho_R | v_R | P_R  |
  ! > | 1.0   | 0.0 | 1000 | 1.0   | 0.0 | 0.01 |
  ! >
  ! > Particles are distributed from -0.5 to 1.5 to provide a stable medium
  ! > outside of the test box from 0.0 to 1.0.
  ! >
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine particle_setup_1d_shock_3(fences)
    
    use physvars, only: &
         dt, &
         nt, &
         n_cpu, &
         my_rank, &
         thermal_constant, &
         initialized_v_minus_half, &
         do_gravity
   
    use module_mirror_boxes, only: &
         periodicity

    use module_interaction_specific_types, only: &
         num_neighbour_particles, &
         PARTICLE_TYPE_FIXED
    
    
    implicit none
    include 'mpif.h'
    
    integer, intent(in) :: fences(-1:n_cpu-1)

    integer :: part_before_me
    integer :: part_including_mine
    real*8 :: offset
    real*8 :: dx
    integer :: actual_particle
    integer :: all_part

    if( my_rank .eq. 0 ) then
       write(*,*) "Creating setup for the 3rd 1D shock test."
       write(*,*) "Using some particles as boundary particles."
       write(*,*) "Use only 0.0 < x< 1.0 for comparison."
    end if

    ! set number of dimension. important for factor for sph kernel
    idim = 1
    
    ! switch gravity off
    do_gravity = .false.

    ! set periodicity
    periodicity = [.false., .false., .false.]
    
    all_part = fences(n_cpu-1)
    part_before_me = fences(my_rank-1)
    part_including_mine = fences(my_rank)
    
    dx = 2./real(all_part)
    offset = dx/2. + dx * part_before_me -0.5
    
    do actual_particle = 1, np_local
       
       particles(actual_particle)%x = [ offset + dx * (actual_particle-1), 0._8, 0._8 ]
       particles(actual_particle)%data%v = [ 0._8, 0._8, 0._8 ]
! v_minus_half will be computed from v and the force in pepc.f90
       particles(actual_particle)%data%v_minus_half = [ 0._8, 0._8, 0._8 ]

       particles(actual_particle)%data%q = dx ! to get a density of 1.0: 1/(particles between 0 and 1), this are half of the particles, so 1/real(all_part/2), which equals dx
       
       if(particles(actual_particle)%x(1) <= 0.5 ) then
          particles(actual_particle)%data%temperature = 1000.0 /(thermal_constant * 1.0)
       else
          particles(actual_particle)%data%temperature = 0.01 /(thermal_constant * 1.0)
       end if
          
       if( (actual_particle + part_before_me) < 2* num_neighbour_particles .or. (all_part - (part_before_me + actual_particle ) ) < 2 * num_neighbour_particles) then
          particles(actual_particle)%data%type = ibset(particles(actual_particle)%data%type, PARTICLE_TYPE_FIXED)
       end if
       
    end do
    
    initialized_v_minus_half = .true.

    ! timestep length
    dt = 0.00001

    ! num timesteps for comparison with analytical results
    nt = 2500

    ! number of neighbour particles known to produce reasonable results
    num_neighbour_particles = 8
    
  end subroutine particle_setup_1d_shock_3


  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! > 3D sphere collapse test for gravity and sph
  ! >
  ! >
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine particle_setup_3d_sphere_grid(fences)
    
    use physvars, only: &
         dt, &
         n_cpu, &
         my_rank, &
         thermal_constant, &
         initialized_v_minus_half
   
    use module_mirror_boxes, only: &
         periodicity

    use module_interaction_specific_types, only: &
         num_neighbour_particles, &
         PARTICLE_TYPE_FIXED
    
    
    implicit none
    include 'mpif.h'
    
    integer, intent(in) :: fences(-1:n_cpu-1)

    integer :: particles_per_dimension
    real*8, dimension(3) :: sphere_center
    real*8 :: sphere_radius
    real*8 :: sphere_radius2
    real*8, dimension(3) :: xyz_min_corner
    real*8 :: particle_separation
    integer :: grid_x, grid_y, grid_z
    real*8, dimension(3) :: coord
    real*8 :: dist_from_center2
    integer :: local_particle_counter = 0
    integer :: my_n_x_grid
    integer :: n_rest
    integer :: my_start_x
    integer :: ierr
    integer :: n_grid
    real :: pi = acos(-1.0)
    real, dimension(1000) :: area
    real :: total_area
    real :: actual_area_sum
    integer :: actual_cpu
    integer :: actual_level_counter
    integer :: already_distributed_levels
    


    if( my_rank .eq. 0 ) then
       write(*,*) "Creating setup for the 3D sphere collapse test."
       write(*,*) "Placing particles on a grid inside a sphere with given center (0,0,0) and radius (0.5 AU)."
    end if

    ! for a 3D grid only 6, 14, 20, ?? neighbours are valid.
    ! but since h is determined as half of the distance to the furthest neighbour and so the density- and force-contribution of the furthest particle is 0,
    ! the number of neighbour has to be one more: 7, 15, 21, ??
    num_neighbour_particles = 20


    ! set number of dimension. important for factor for sph kernel
    idim = 3
    
    ! set periodicity
    periodicity = [.false., .false., .false.]
    
    ! sphere parameter
    sphere_center = [0., 0., 0.]
    sphere_radius = 0.5

    ! Estimate the number of grid-points needed to get approximately npart_total particles:
    ! The volume of the sphere inside the cube is pi/6 * the volume of the cube.
    ! We want n_grid^3 * pi/6 > npart_total
    n_grid = int( (6.*npart_total/pi)**(1./3.)*0.6 )

    particles_per_dimension = n_grid

!    write(*,*) "n_grid:", n_grid

    ! derived values
    sphere_radius2= sphere_radius**2
    xyz_min_corner = sphere_center - [sphere_radius, sphere_radius, sphere_radius]
    particle_separation = sphere_radius * 2. / real(particles_per_dimension-1)

!    write(*,*) "particle separation:", particle_separation

    
    ! distribute the particles to the processes with approximately the same number of particles per process.
    ! we go through the sphere grid-level by grid-level and estimate the number of particles on this level inside the sphere.
    do actual_level_counter = 1, particles_per_dimension
       area(actual_level_counter) = pi * (sphere_radius**2 - (sphere_radius - (actual_level_counter-1)*particle_separation)**2)
    end do

!    write(*,*) "area:", area(1:particles_per_dimension)

    total_area = sum(area(1:particles_per_dimension))

!    write(*,*) "total area:", total_area


    
    
    ! Now the grid-levels are distributed to the processes that the area per process is approximately area/n_cpu

    already_distributed_levels = 0
    
    do actual_cpu = 0, n_cpu-2 ! the last process gets the rest
       actual_area_sum = 0
       actual_level_counter = 1
       
       do while((actual_area_sum + area(already_distributed_levels + actual_level_counter))  < total_area/n_cpu)
          actual_area_sum = actual_area_sum + area(already_distributed_levels + actual_level_counter)
          actual_level_counter = actual_level_counter + 1
       end do
       
       if(my_rank .eq. actual_cpu) then
          my_n_x_grid = actual_level_counter - 1
          my_start_x = already_distributed_levels + 1
       end if
       
       already_distributed_levels = already_distributed_levels + (actual_level_counter - 1)
    end do
       
    if(my_rank .eq. n_cpu -1) then
       my_n_x_grid = particles_per_dimension - already_distributed_levels
       my_start_x = already_distributed_levels + 1
    end if

    ! ! for positions on a 3d regular grid it is tested whether they are inside the sphere with sphere_radius around sphere_center.
    ! ! for n_cpu processes the grid is split in x direction.
    ! ! first, all cpus get the same number of gridpoints in x direction.
    ! my_n_x_grid = int(particles_per_dimension / n_cpu)

    ! ! the rest
    ! n_rest = mod(particles_per_dimension, n_cpu)

    ! ! is distributed to the first n_rest cpus
    ! if( my_rank < n_rest) then
    !    my_n_x_grid = my_n_x_grid +1
    ! end if
    
    ! ! now calculate the starting x coordinate
    ! if( my_rank < n_rest) then ! this process and all processes before got one more
    !    my_start_x = my_rank * my_n_x_grid
    ! else ! n_rest processes got one more and the others the same as this process
    !    my_start_x = n_rest * (my_n_x_grid + 1) + ( my_rank - n_rest ) * my_n_x_grid
    ! end if

    write(*,*) "3d grid distribution:", my_rank, my_n_x_grid, my_start_x
!    call flush()
    
!    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    
    ! if( my_n_x_grid * particles_per_dimension * particles_per_dimension > np_local) then
    !    write(*,*) "3d grid, local particle number could be to big. Exiting."
    !    ! I know that this is a very rough estimation but its enough for now.
    !    call flush()
       
    !    call MPI_ABORT(MPI_COMM_WORLD, ierr)
    ! end if


    do grid_x = my_start_x, my_start_x + my_n_x_grid
       do grid_y = 1, particles_per_dimension
          do grid_z = 1, particles_per_dimension

             coord = xyz_min_corner + particle_separation * [grid_x-1, grid_y-1, grid_z-1]
             dist_from_center2 = ( coord(1)**2 + coord(2)**2 + coord(3)**2)

             if( dist_from_center2 <= sphere_radius2 ) then ! use this particle
                particles(local_particle_counter + 1)%x = coord
                
                local_particle_counter = local_particle_counter + 1
             end if
          end do
       end do
    end do

    np_local = local_particle_counter

    call MPI_ALLREDUCE(np_local, npart_total, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

    initialized_v_minus_half = .true.

!    call MPI_ABORT(MPI_COMM_WORLD, ierr)


    particles(1:np_local)%data%q           = 1._8/real(npart_total, 8)
    particles(1:np_local)%data%temperature = 1.





  end subroutine particle_setup_3d_sphere_grid





  



  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>
  !>
  !>
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! subroutine particle_setup_doublesphere(fences)
  !   implicit none
  !   integer, intent(in) :: fences(-1:n_cpu-1)
  !   integer :: mpi_cnt, p
  !   real*8 :: xt, yt, zt

  !   do mpi_cnt = 0, n_cpu-1
  !      do p = 1, (fences(mpi_cnt) - fences(mpi_cnt-1))

  !         xt = 1.0_8
  !         yt = 1.0_8
  !         zt = 1.0_8

  !         do while ( (xt*xt + yt*yt + zt*zt) > 0.5_8)
  !            xt = par_rand() - 0.5_8
  !            yt = par_rand() - 0.5_8
  !            zt = par_rand() - 0.5_8
  !         end do

  !         if (par_rand() < 0.5) then
  !            xt = xt + 2.
  !         else
  !            xt = xt - 2.
  !         endif

  !         if ( my_rank == mpi_cnt .and. p <= np_local ) then
  !            z(p) = zt
  !            y(p) = yt
  !            x(p) = xt
  !         end if
  !      end do
  !   end do

  ! end subroutine particle_setup_doublesphere


  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>
  !>
  !>
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! subroutine particle_setup_manyspheres(fences)
  !   implicit none
  !   integer, intent(in) :: fences(-1:n_cpu-1)
  !   integer :: mpi_cnt, p
  !   real*8 :: xt, yt, zt
  !   real*8 :: delta(3)

  !   do mpi_cnt = 0, n_cpu-1
  !      do p = 1, (fences(mpi_cnt) - fences(mpi_cnt-1))

  !         xt = 1.0_8
  !         yt = 1.0_8
  !         zt = 1.0_8

  !         do while ( (xt*xt + yt*yt + zt*zt) > 0.5_8)
  !            xt = par_rand() - 0.5_8
  !            yt = par_rand() - 0.5_8
  !            zt = par_rand() - 0.5_8
  !         end do

  !         delta = 20._8*GetSphereCenter(nint(42.*par_rand())) - 10._8

  !         if ( my_rank == mpi_cnt .and. p <= np_local ) then
  !            z(p) = zt + delta(1)
  !            y(p) = yt + delta(2)
  !            x(p) = xt + delta(3)
  !         end if
  !      end do
  !   end do

  ! end subroutine particle_setup_manyspheres


  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>
  !>
  !>
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! subroutine particle_setup_plummer(fences)
  !   use module_units
  !   implicit none
  !   integer, intent(in) :: fences(-1:n_cpu-1)
  !   integer :: mpi_cnt, p
  !   real*8 :: xt, yt, zt
  !   real*8 :: r1, r2

  !   do mpi_cnt = 0, n_cpu-1
  !      do p = 1, (fences(mpi_cnt) - fences(mpi_cnt-1))

  !         r1 = 0.
  !         do while (.not.(r1 > 0.1 .and. r1 < 3))
  !            r1 = (par_rand()**(-0.2D01/0.3D01)-0.1D01)**(-0.5D00)
  !         end do

  !         zt = (0.1D01-0.2D01*par_rand())*r1
  !         r2 = par_rand()
  !         xt = (r1**0.2D01-zt**0.2D01)**(0.5D00)*cos(0.2D01*pi*r2)
  !         yt = (r1**0.2D01-zt**0.2D01)**(0.5D00)*sin(0.2D01*pi*r2)

  !         if ( my_rank == mpi_cnt .and. p <= np_local ) then
  !            z(p) = zt
  !            y(p) = yt
  !            x(p) = xt
  !         end if
  !      end do
  !   end do

  ! end subroutine particle_setup_plummer



  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>
  !>
  !>
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! subroutine particle_setup_flatdisc(fences)
  !   implicit none
  !   integer, intent(in) :: fences(-1:n_cpu-1)
  !   integer :: mpi_cnt, p
  !   real*8 :: xt, yt, zt

  !   do mpi_cnt = 0, n_cpu-1
  !      do p = 1, (fences(mpi_cnt) - fences(mpi_cnt-1))

  !         xt = 1.0_8
  !         yt = 1.0_8
  !         zt = 0.0_8

  !         do while ( (xt*xt + yt*yt) > 0.5_8)
  !            xt = par_rand() - 0.5_8
  !            yt = par_rand() - 0.5_8
  !         end do

  !         if ( my_rank == mpi_cnt .and. p <= np_local ) then
  !            z(p) = zt
  !            y(p) = yt
  !            x(p) = xt
  !         end if
  !      end do
  !   end do

  ! end subroutine particle_setup_flatdisc


  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>
  !>
  !>
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! subroutine particle_setup_madelung()
  !   use module_mirror_boxes
  !   implicit none
  !   real*8 :: delta(3)
  !   integer :: i,j,k,n(3), myidx, globalidx

  !   n(1) = NINT((npart_total/8)**(1./3.))*2
  !   n(2) = NINT((npart_total/n(1)/4)**(1./2.))*2
  !   n(3) = npart_total/n(1)/n(2)

  !   if (my_rank == 0) write(*,*) "Using ", n, "particles per edge"

  !   delta(1) = sqrt(dot_product(t_lattice_1,t_lattice_1))/real(n(1))
  !   delta(2) = sqrt(dot_product(t_lattice_2,t_lattice_2))/real(n(2))
  !   delta(3) = sqrt(dot_product(t_lattice_3,t_lattice_3))/real(n(3))

  !   myidx     = 0
  !   globalidx = 0

  !   do i = 0, n(1)-1
  !      do j = 0, n(2)-1
  !         do k = 0, n(3)-1

  !            globalidx = globalidx + 1

  !            if ( mod((globalidx-1)/2,n_cpu) == my_rank) then ! distribute pairs of electron and ion, since np_local is constructed a bit weird
  !               myidx = myidx + 1

  !               x(myidx)  = (i + 0.5)*delta(1)
  !               y(myidx)  = (j + 0.5)*delta(2)
  !               z(myidx)  = (k + 0.5)*delta(3)
  !               ux(myidx) = 0
  !               uy(myidx) = 0
  !               uz(myidx) = 0

  !               if (mod(i+j+k, 2) == 0) then
  !                  q(myidx)       = qe
  !                  m(myidx)       = mass_e
  !                  pelabel(myidx) = -globalidx
  !               else
  !                  q(myidx)       = qi
  !                  m(myidx)       = mass_i
  !                  pelabel(myidx) = globalidx
  !               end if

  !            end if

  !         end do
  !      end do
  !   end do

  !   work(1:np_local) = 1.

  !   if (myidx .ne. np_local) write(*,*) "ERROR in special_start(7): PE", my_rank, "set up", myidx, &
  !        "particles, but np_local=", np_local, "globalidx=", globalidx, "npart_total=",npart_total

  !   !write(*,*) "Particle positions: "
  !   !do i=1,np_local
  !   !  write(*,*) my_rank,i,x(i), y(i), z(i), q(i), pelabel(i)
  !   !end do

  ! end subroutine particle_setup_madelung



  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! >
  ! >
  ! >
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine particle_setup_homogeneous_fast(fences, dim)

    implicit none
    integer, intent(in) :: fences(-1:n_cpu-1)
    integer :: p, dim, i
    real*8, dimension(3) :: tmp
    
    ! initialize random number generator with some arbitrary seed
    tmp(1) = par_rand(my_rank + 13)

    tmp = [ 0._8, 0._8, 0._8 ]
    
    do p = 1, (fences(my_rank) - fences(my_rank-1))
       do i = 1, dim
          tmp(i) = par_rand()
       end do

       particles(p)%x = tmp
    end do

  end subroutine particle_setup_homogeneous_fast



  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>
  !>
  !>
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! subroutine particle_setup_icosahedron(fences)
  !   use module_units
  !   use module_icosahedron
  !   implicit none
  !   integer, intent(in) :: fences(-1:n_cpu-1)
  !   integer :: currlayer, particletype
  !   integer :: p
  !   real*8 :: r(3), xt, yt, zt

  !   ! "random" initialization of par_rand
  !   xt = par_rand(my_rank + rngseed)

  !   do p=1,np_local/2
  !      ! get particle position inside the cluster
  !      r = get_particle(p + fences(my_rank-1)/2-1, currlayer, particletype)

  !      ! put an ion there
  !      x(np_local-p+1)   = r(1)
  !      y(np_local-p+1)   = r(2)
  !      z(np_local-p+1)   = r(3)
  !      ux(np_local-p+1)  = 0.
  !      uy(np_local-p+1)  = 0.
  !      uz(np_local-p+1)  = 0.
  !      q(np_local-p+1)  = qi
  !      m(np_local-p+1)  = mass_i
  !      !pelabel(np_local-p+1)  = 1 + p + fences(my_rank-1)/2-1
  !      pelabel(np_local-p+1)  = 1 + particletype

  !      ! and put an electron into near proximity
  !      xt = 2*pi*par_rand()
  !      yt =   pi*par_rand()
  !      zt = 0.025*par_rand()

  !      x(p) = x(np_local-p+1) + zt * cos(xt) * sin(yt)
  !      y(p) = y(np_local-p+1) + zt * sin(xt) * sin(yt)
  !      z(p) = z(np_local-p+1) + zt * cos(yt)

  !      ! chose random velocity
  !      xt = par_rand()*2*pi
  !      yt = par_rand()  *pi
  !      zt = par_rand()  *vte

  !      ux(p)  = zt * cos(xt) * sin(yt)
  !      uy(p)  = zt * cos(xt) * sin(yt)
  !      uz(p)  = zt * cos(yt)
  !      q(p)  = qe
  !      m(p)  = mass_e
  !      pelabel(p)  = -(1 + particletype)
  !   end do

  ! end subroutine particle_setup_icosahedron



  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>
  !>
  !>
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! subroutine particle_setup_ionlattice()
  !   implicit none
  !   integer :: n, i, j, k
  !   real*8 :: delta(3)
  !   integer :: globalidx, myidx

  !   delta(1) = 1./real(n)
  !   delta(2) = 1./real(n)
  !   delta(3) = 1./real(n)

  !   myidx     = 0
  !   globalidx = 0

  !   do i = 0, n-1
  !      do j = 0, n-1
  !         do k = 0, n-1

  !            globalidx = globalidx + 1

  !            if ( mod((globalidx-1)/2,n_cpu) == my_rank) then ! distribute pairs of electron and ion, since np_local is constructed a bit weird
  !               myidx = myidx + 1

  !               x(np_local-myidx+1)  = (i + 0.5)*delta(1)
  !               y(np_local-myidx+1)  = (j + 0.5)*delta(2)
  !               z(np_local-myidx+1)  = (k + 0.5)*delta(3)
  !               ux(np_local-myidx+1) = 0
  !               uy(np_local-myidx+1) = 0
  !               uz(np_local-myidx+1) = 0

  !               q(np_local-myidx+1) = qi
  !               m(np_local-myidx+1) = mass_i
  !            end if

  !         end do
  !      end do
  !   end do

  !   do i = 1,ne
  !      globalidx = globalidx + 1

  !      if ( mod((globalidx-1)/2,n_cpu) == my_rank) then ! distribute pairs of electron and ion, since np_local is constructed a bit weird
  !         myidx = myidx + 1

  !         x(np_local-myidx+1)  = par_rand()
  !         y(np_local-myidx+1)  = par_rand()
  !         z(np_local-myidx+1)  = par_rand()

  !         ux(np_local-myidx+1) = 0
  !         uy(np_local-myidx+1) = 0
  !         uz(np_local-myidx+1) = 0

  !         q(np_local-myidx+1) = qe
  !         m(np_local-myidx+1) = mass_e
  !      end if
  !   end do

  !   if (myidx .ne. np_local) write(*,*) "ERROR in special_start(8): PE", my_rank, "set up", myidx, &
  !        "particles, but np_local=", np_local, "globalidx=", globalidx, "npart_total=",npart_total
  !   if (globalidx .ne. npart_total) write(*,*) "ERROR in special_start(8): PE", my_rank, "set up globalidx=", globalidx, ", but npart_total=",npart_total

  ! end subroutine particle_setup_ionlattice




  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>
  !> generic initialisation of charges(here masses), labels, work, velocities, temperatures
  !> all set to default values, so real particle setup has to be called afterwards
  !>
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_generic(fences)

    use module_interaction_specific_types, only: &
         PARTICLE_TYPE_DEFAULT

    implicit none
    integer, intent(in) :: fences(-1:n_cpu-1)
    integer :: p

    ! TODO: clean this quick'n'dirty setup stuff
    do p = 1, np_local
       particles(p)%data%q           = 1._8
       particles(p)%label            = fences(my_rank-1) + p
       particles(p)%work             = 1._8
       particles(p)%data%v_minus_half  = [ 0._8, 0._8, 0._8 ]
       particles(p)%data%v           = [ 0._8, 0._8, 0._8 ]
       particles(p)%data%temperature = 0._8
       particles(p)%data%type        = PARTICLE_TYPE_DEFAULT   ! moving, gas..., see module_interaction_specific_types
    end do

  end subroutine init_generic




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>
  !> initialisation of fields to zero
  !>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! subroutine init_fields_zero()
  !   implicit none

  !   ex(1:np_local) = 0.
  !   ey(1:np_local) = 0.
  !   ez(1:np_local) = 0.
  !   pot(1:np_local) = 0.

  ! end subroutine init_fields_zero




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>
  !> portable random number generator, see numerical recipes
  !> check for the random numbers:
  !> the first numbers should be 0.2853809, 0.2533582 and 0.0934685
  !> the parameter iseed is optional
  !>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function par_rand(iseed)
    implicit none
    real :: par_rand
    integer, intent(in), optional :: iseed

    integer, parameter :: IM1  = 2147483563
    integer, parameter :: IM2  = 2147483399
    real,    parameter :: AM   = 1.0/IM1
    integer, parameter :: IMM1 = IM1-1
    integer, parameter :: IA1  = 40014
    integer, parameter :: IA2  = 40692
    integer, parameter :: IQ1  = 53668
    integer, parameter :: IQ2  = 52774
    integer, parameter :: IR1  = 12211
    integer, parameter :: IR2  = 3791
    integer, parameter :: NTAB = 32
    integer, parameter :: NDIV = 1+IMM1/NTAB
    real,    parameter :: eps_ = 1.2e-7 ! epsilon(eps_)
    real,    parameter :: RNMX = 1.0 - eps_

    integer :: j, k
    integer, volatile, save :: idum  = -1
    integer, volatile, save :: idum2 =  123456789
    integer, volatile, save :: iy    =  0
    integer, volatile, save :: iv(NTAB)


    if (idum <=0 .or. present(iseed)) then
       if (present(iseed)) then
          idum = iseed
       else
          if (-idum < 1) then
             idum = 1
          else
             idum = -idum
          endif
       endif

       idum2 = idum

       do j = NTAB+7,0,-1
          k = idum/IQ1
          idum = IA1 * (idum-k*IQ1) - k*IR1
          if (idum < 0 ) idum = idum + IM1

          if (j<NTAB) iv(j+1) = idum

       end do
       iy = iv(1)
    end if

    k = idum/IQ1
    idum = IA1 * (idum-k*IQ1) - k*IR1
    if (idum < 0) idum = idum + IM1

    k = idum2/IQ2
    idum2 = IA2 * (idum2-k*IQ2) - k*IR2
    if (idum2 < 0) idum2 = idum2 + IM2

    j = iy/NDIV + 1
    iy = iv(j)-idum2
    iv(j) = idum

    if (iy < 1) iy = iy + IMM1
    par_rand = AM*iy
    if (par_rand > RNMX) par_rand = RNMX

  end function par_rand


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>
  !> returns one (of 50 fixed) random 3D vectors
  !>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! function GetSphereCenter(idx)
  !   implicit none
  !   integer, intent(in) :: idx
  !   real*8, dimension(3) :: GetSphereCenter

  !   real*8, dimension(3, 50) :: random_vectors = reshape([  0.7011 ,  0.0942 ,  0.0012, &
  !        0.6663 ,  0.5985 ,  0.4624,&
  !        0.5391 ,  0.4709 ,  0.4243,&
  !        0.6981 ,  0.6959 ,  0.4609,&
  !        0.6665 ,  0.6999 ,  0.7702,&
  !        0.1781 ,  0.6385 ,  0.3225,&
  !        0.1280 ,  0.0336 ,  0.7847,&
  !        0.9991 ,  0.0688 ,  0.4714,&
  !        0.1711 ,  0.3196 ,  0.0358,&
  !        0.0326 ,  0.5309 ,  0.1759,&
  !        0.5612 ,  0.6544 ,  0.7218,&
  !        0.8819 ,  0.4076 ,  0.4735,&
  !        0.6692 ,  0.8200 ,  0.1527,&
  !        0.1904 ,  0.7184 ,  0.3411,&
  !        0.3689 ,  0.9686 ,  0.6074,&
  !        0.4607 ,  0.5313 ,  0.1917,&
  !        0.9816 ,  0.3251 ,  0.7384,&
  !        0.1564 ,  0.1056 ,  0.2428,&
  !        0.8555 ,  0.6110 ,  0.9174,&
  !        0.6448 ,  0.7788 ,  0.2691,&
  !        0.3763 ,  0.4235 ,  0.7655,&
  !        0.1909 ,  0.0908 ,  0.1887,&
  !        0.4283 ,  0.2665 ,  0.2875,&
  !        0.4820 ,  0.1537 ,  0.0911,&
  !        0.1206 ,  0.2810 ,  0.5762,&
  !        0.5895 ,  0.4401 ,  0.6834,&
  !        0.2262 ,  0.5271 ,  0.5466,&
  !        0.3846 ,  0.4574 ,  0.4257,&
  !        0.5830 ,  0.8754 ,  0.6444,&
  !        0.2518 ,  0.5181 ,  0.6476,&
  !        0.2904 ,  0.9436 ,  0.6790,&
  !        0.6171 ,  0.6377 ,  0.6358,&
  !        0.2653 ,  0.9577 ,  0.9452,&
  !        0.8244 ,  0.2407 ,  0.2089,&
  !        0.9827 ,  0.6761 ,  0.7093,&
  !        0.7302 ,  0.2891 ,  0.2362,&
  !        0.3439 ,  0.6718 ,  0.1194,&
  !        0.5841 ,  0.6951 ,  0.6073,&
  !        0.1078 ,  0.0680 ,  0.4501,&
  !        0.9063 ,  0.2548 ,  0.4587,&
  !        0.8797 ,  0.2240 ,  0.6619,&
  !        0.8178 ,  0.6678 ,  0.7703,&
  !        0.2607 ,  0.8444 ,  0.3502,&
  !        0.5944 ,  0.3445 ,  0.6620,&
  !        0.0225 ,  0.7805 ,  0.4162,&
  !        0.4253 ,  0.6753 ,  0.8419,&
  !        0.3127 ,  0.0067 ,  0.8329,&
  !        0.1615 ,  0.6022 ,  0.2564,&
  !        0.1788 ,  0.3868 ,  0.6135,&
  !        0.4229 ,  0.9160 ,  0.5822], shape(random_vectors))

  !   GetSphereCenter = 1._8*random_vectors(:, idx + 1)

  ! end function GetSphereCenter

end module module_particle_setup
