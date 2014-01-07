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

!>
!> helper module
!>
module pepca_helper
  use module_pepc_types
  use module_timings
  use pepca_units
  implicit none
  private
  save

  public pepca_init
  public pepca_nml_t

  ! dimension of space
  integer(kind_dim), public, parameter :: dim = 2
  ! side lengths of particle box for particle_config = 0
  real*8 :: boxdims(3) = [10, 10, 10]

  integer, public, parameter :: OI_PARTICLES_VTK    = 1
  integer, public, parameter :: OI_PARTICLES_ASC    = 2
  integer, public, parameter :: OI_PARTICLES_MPI    = 3
  integer, public, parameter :: OI_DENSITIES_VTK    = 4
  integer, public, parameter :: OI_DOMAIN_VTK       = 5
  integer, public, parameter :: OI_VERIFY_PARTICLES = 6
  integer, public, parameter :: OI_MAXIDX = OI_VERIFY_PARTICLES

  integer, public, parameter :: PARAMS_VTE = 1
  integer, public, parameter :: PARAMS_VTI = 2
  integer, public, parameter :: PARAMS_KX  = 3
  integer, public, parameter :: PARAMS_A   = 4
  integer, public, parameter :: PARAMS_MAXIDX = PARAMS_A

  !> parameter collection for pepca
  type pepca_nml_t
    ! grid for density output
    integer :: Ngrid(1:3) = [500, 1000, 1]
    ! MPI variables
    integer(kind_pe) :: rank, nrank
    integer(kind_default) :: comm
    ! time variables
    real*8  :: dt  !< timestep (in simunits), set via pfasst parameters
    integer :: nt  !< number of timesteps, set via pfasst parameters
    ! configuration variables
    integer :: particle_config = 0
    real*8 :: setup_params(PARAMS_MAXIDX) = 0.
    ! total number of particles per species
    integer(kind_particle) :: numparts_total = 2500 !296568
    ! number of particles per species and rank, will be set automatically later
    integer(kind_particle) :: numparts
    ! output control intervals - set to 0 to deactivate output, see above for meaning of the fields
    integer :: output_interval(OI_MAXIDX) = [1, 1, 1, 1, 1, 0]
    ! use PFASST
    logical :: use_pfasst = .true.
  end type

  type(pepca_nml_t), public :: pepca_nml

  contains

  subroutine pepca_init(nml, particles, dt, nt)
    implicit none
    type(pepca_nml_t), intent(inout) :: nml
    type(t_particle), allocatable, target, intent(out) :: particles(:)
    real*8, intent(in)  :: dt
    integer, intent(in) :: nt

    call set_parameter(nml, dt, nt)
    call setup_particles(particles, nml)

  end subroutine


  !> set initial values for particle positions and velocities in y0
  subroutine setup_particles(particles, nml)
    use module_debug
    use module_mirror_boxes
    implicit none
    type(pepca_nml_t), intent(in) :: nml
    type(t_particle), allocatable, target, intent(out) :: particles(:)
    real*8 :: n0

    select case (nml%particle_config)
      case (0)
        ! homogeneously filled cube
        call generate_particles_homogeneous(particles, nml%numparts, nml%setup_params(PARAMS_VTE), nml%numparts, nml%setup_params(PARAMS_VTI), nml%comm)
      case (1)
        ! homogeneously filled cube with position displacement for inducing plasma waves
        call generate_particles_homogeneous(particles, nml%numparts, nml%setup_params(PARAMS_VTE), nml%numparts, nml%setup_params(PARAMS_VTI), nml%comm)
        n0 = nml%numparts / (product(boxdims(1:dim)) / unit_length_micron_per_simunit**dim)
        call displace_particles(particles, nml%setup_params(PARAMS_KX) / boxdims(1) * unit_length_micron_per_simunit, &
                                           nml%setup_params(PARAMS_A) / n0)
        periodicity = [.true., .true., .false.]
        t_lattice_1 = [boxdims(1) / unit_length_micron_per_simunit, 0._8, 0._8]
        t_lattice_2 = [0._8, boxdims(2) / unit_length_micron_per_simunit, 0._8]
        t_lattice_3 = [0._8, 0._8, boxdims(3) / unit_length_micron_per_simunit]
        call init_movement_constraint()
        call constrain_periodic(particles)
      case (-1)
        call read_particles(particles, 'E_phase_space.dat', nml%numparts, 'I_phase_space.dat', nml%numparts, nml%comm)
      case default
        DEBUG_ERROR(*, 'setup_particles() - invalid value for particle_config:', nml%particle_config)
    end select

  end subroutine



  subroutine set_parameter(nml, dt, nt)
    use module_pepc
    use module_debug
    use module_interaction_specific, only : theta2, eps2
    use treevars, only : num_threads, np_mult
    implicit none

    type(pepca_nml_t), intent(inout) :: nml
    real*8, intent(in)  :: dt
    integer, intent(in) :: nt

    integer, parameter :: fid = 12
    character(255)     :: para_file
    logical            :: read_para_file

    integer :: Ngrid(1:3)
    integer :: particle_config
    integer(kind_particle) :: numparts_total
    logical :: use_pfasst
    integer :: output_interval(OI_MAXIDX)
    real*8 :: setup_params(PARAMS_MAXIDX)

    real*8 :: eps = 1.e-5 ! interaction cutoff parameter

    namelist /pepcapfasst/ eps, Ngrid, particle_config, setup_params, numparts_total, use_pfasst, output_interval

    ! frontend parameters
    Ngrid           = nml%Ngrid
    particle_config = nml%particle_config
    setup_params    = nml%setup_params
    numparts_total  = nml%numparts_total
    use_pfasst      = nml%use_pfasst
    output_interval = nml%output_interval

    ! pepc parameters
    theta2      = 0.36
    num_threads = 8
    np_mult     = -500
    ! read in namelist file
    call pepc_read_parameters_from_first_argument(read_para_file, para_file)
    eps2 = (eps/unit_length_micron_per_simunit)**2

    if (read_para_file) then
      if(nml%rank==0) write(*,'(a)') ' == reading parameter file, section pepcapfasst: ', para_file
      open(fid,file=para_file)
      read(fid,NML=pepcapfasst)
      close(fid)
    else
      if(nml%rank==0) write(*,*) ' == no param file, using default parameters '
    end if

    ! frontend parameters
    nml%Ngrid           = Ngrid
    nml%particle_config = particle_config
    nml%setup_params    = setup_params
    nml%numparts_total  = numparts_total
    nml%numparts        = numparts_total/nml%nrank
    if (mod(numparts_total, nml%nrank) > nml%rank) nml%numparts = nml%numparts + 1
    nml%use_pfasst      = use_pfasst
    nml%output_interval = output_interval

    if (nml%output_interval(OI_PARTICLES_MPI)>0 .and. nml%output_interval(OI_VERIFY_PARTICLES)>0) then
      DEBUG_WARNING(*, 'You activated particle output and verification in the same run. This will overwrite reference data after comparison and is usually not what you want to do.')
    endif

    ! derived from pfasst parameters
    nml%dt      = dt
    nml%nt      = nt

    if(nml%rank==0) then
      write(*,'(a,i12)')       ' == total  of particles per spec          : ', nml%numparts_total
      write(*,'(a,i12)')       ' == number of particles per spec and rank : ', nml%numparts
      write(*,'(a,i12)')       ' == number of time steps                  : ', nml%nt
      write(*,'(a,es12.4)')    ' == time step (simunits)                  : ', nml%dt
      write(*,'(a,es12.4)')    ' == final time (simunits)                 : ', nml%dt*nml%nt
      write(*,'(a,es12.4)')    ' == time step (fs)                        : ', nml%dt*unit_time_fs_per_simunit
      write(*,'(a,es12.4)')    ' == final time (fs)                       : ', nml%dt*nml%nt*unit_time_fs_per_simunit
    end if

    call pepc_prepare(dim)
  end subroutine set_parameter


  subroutine read_particles(p, file_el, nel, file_ion, nion, comm)
    implicit none

    type(t_particle), allocatable, intent(inout) :: p(:)
    integer(kind_particle), intent(in) :: nel, nion
    character(*), intent(in) :: file_el, file_ion
    integer(kind_default), intent(in) :: comm

    integer(kind_pe) :: rank
    integer, parameter :: filehandle = 47
    integer(kind_particle) :: ip
    integer(kind_default) :: ierr

    allocate(p(nel+nion))
    call MPI_COMM_RANK(comm, rank, ierr)

    ! FIXME: currently, we read all particles on the root rank of MPI_COMM_SPACE; thus, particle numbers should be zero on all others
    if(rank==0) then
      write(*,'(a, 2(x,a))') ' == [load] reading particles from files', file_el, file_ion

      ! we read all particles the root rank

      ! Script for file preprocessing: sed -i 's/\([0-9]\)-/\1 -/g' filename
      !call system("sed -i 's/\([0-9]\)-/\1 -/g' " // trim(file_el))
      open(filehandle, file=trim(file_el), STATUS='OLD', POSITION = 'REWIND', ACTION='READ')
      do ip=1,nel
        read(filehandle, *)   p(ip)%x(1:dim), p(ip)%data%v(1:dim)
        p(ip)%x(dim+1:)      = 0.
        p(ip)%data%v(dim+1:) = 0.
        ! other stuff
        p(ip)%label       = -ip
        p(ip)%data%q      =  unit_qe
        p(ip)%data%m      =  unit_me
        p(ip)%work        =  1.0

        p(ip)%x = p(ip)%x / unit_length_micron_per_simunit
      end do
      close(filehandle)

      !call system("sed -i 's/\([0-9]\)-/\1 -/g' " // trim(file_ion))
      open(filehandle, file=trim(file_ion), STATUS='OLD', POSITION = 'REWIND', ACTION='READ')
      do ip=nel+1,nel+nion
        read(filehandle, *)   p(ip)%x(1:dim), p(ip)%data%v(1:dim)
        p(ip)%x(dim+1:)        = 0.
        p(ip)%data%v(dim+1:)   = 0.
        ! other stuff
        p(ip)%label       =  ip-nel
        p(ip)%data%q      =  unit_qp
        p(ip)%data%m      =  unit_mp
        p(ip)%work        =  1.0

        p(ip)%x = p(ip)%x / unit_length_micron_per_simunit
      end do
      close(filehandle)
    endif

  end subroutine read_particles


  subroutine generate_particles_homogeneous(p, nel, vte, nion, vti, comm)
    use pepca_units
    implicit none
    include 'mpif.h'

    type(t_particle), allocatable, intent(inout) :: p(:)
    integer(kind_particle), intent(in) :: nel, nion
    real*8, intent(in) :: vte, vti
    integer(kind_default), intent(in) :: comm

    integer(kind_pe) :: rank, nrank
    integer(kind_default) :: ierr
    real*8 :: pos(1:3), vel(1:3)
    integer(kind_particle) :: i,l
    integer(kind_dim) :: k

    integer(kind_particle) :: nel_tot, npart_left, npart_tot

    allocate(p(nel+nion))
    call MPI_COMM_RANK(comm,  rank, ierr)
    call MPI_COMM_SIZE(comm, nrank, ierr)

    call MPI_ALLREDUCE(nel,        nel_tot,  1, MPI_KIND_PARTICLE, MPI_SUM, comm, ierr)
    call MPI_ALLREDUCE(nel+nion, npart_tot,  1, MPI_KIND_PARTICLE, MPI_SUM, comm, ierr)
    call MPI_EXSCAN(   nel+nion, npart_left, 1, MPI_KIND_PARTICLE, MPI_SUM, comm, ierr)

    ! stupid parallel random number generation
    if(rank==0) write(*,'(a, 2(x,a))') ' == [generate] generating particles'

    l   = 0
    do i=1,npart_tot
      do k=1,dim
        pos(k) = par_rand() * boxdims(k)
      end do

      if (i<=nel_tot) then
        vel(1:3) = velocity_2d_maxwell(vte)
      else
        vel(1:3) = velocity_2d_maxwell(vti)
      endif


      if ((i>npart_left) .and. (l<nel+nion)) then
        l = l + 1
        p(l)%x(:)          = 0.
        p(l)%data%v(:)     = 0.
        p(l)%x(1:dim)      = pos(1:dim) / unit_length_micron_per_simunit
        p(l)%data%v(1:dim) = vel(1:dim)
        p(l)%label         = i
        p(l)%work          =  1.0

        if (i<=nel_tot) then ! first half of all particles are electrons, second half ions
          p(l)%data%q      =  unit_qe
          p(l)%data%m      =  unit_me
        else
          p(l)%data%q      =  unit_qp
          p(l)%data%m      =  unit_mp
        endif

      endif
    end do
  end subroutine generate_particles_homogeneous

  ! see Birdsall & Langdon pp 81-85
  subroutine displace_particles(particles, k, A_over_n0)
    use module_pepc_types
    implicit none
    type(t_particle), allocatable, intent(inout) :: particles(:)
    real*8, intent(in) :: k
    real*8, intent(in) :: A_over_n0

    integer(kind_particle) :: p

    do p=1, size(particles, kind=kind(p))
      particles(p)%x(1) = particles(p)%x(1) + A_over_n0 / k * cos(k*particles(p)%x(1))
    end do

  end subroutine


  function velocity_2d_maxwell(vt) result(v)
    implicit none

    real*8, intent(in) :: vt
    real*8, dimension(3) :: v

    real*8 :: xi, v0
    real(kind=8), parameter :: pi = 3.1415926535897932384626434D0

    xi = par_rand()
    v0 = vt * sqrt(-2.0D0 * log(max(xi, 10.0_8**(-15))))
    xi = par_rand()
    v(1) = v0 * cos(2 * pi * xi)
    v(2) = v0 * sin(2 * pi * xi)
    v(3) = 0.0D0

  end function velocity_2d_maxwell

  !>
  !> portable random number generator, see numerical recipes
  !> check for the random numbers:
  !> the first numbers should be 0.2853809, 0.2533582 and 0.0934685
  !> the parameter iseed is optional
  !>
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

end module
