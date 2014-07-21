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
module pepcboris_helper
  use module_pepc_kinds
  use module_pepc_types
  use module_timings
  implicit none
  private
  save

  public pepcboris_init
  public pepcboris_nml_t
  public get_magnetic_field
  public cross_prod
  public cross_prod_plus
  public get_mpi_rank

  ! dimension of space
  integer(kind_dim), public, parameter :: dim = 3

  integer, public, parameter :: PARAMS_X0  = 1
  integer, public, parameter :: PARAMS_Y0  = 2
  integer, public, parameter :: PARAMS_Z0  = 3
  integer, public, parameter :: PARAMS_VX0 = 4
  integer, public, parameter :: PARAMS_VY0 = 5
  integer, public, parameter :: PARAMS_VZ0 = 6
  integer, public, parameter :: PARAMS_OMEGAE = 7
  integer, public, parameter :: PARAMS_OMEGAB = 8
  integer, public, parameter :: PARAMS_RADIUS = 9
  integer, public, parameter :: PARAMS_VELOCITY_SPREAD = 10
  integer, public, parameter :: PARAMS_MAXIDX = PARAMS_VELOCITY_SPREAD

  integer, public, parameter :: WM_BORIS            = 1
  integer, public, parameter :: WM_BORIS_SDC_REF    = 2
  integer, public, parameter :: WM_ANALYTIC         = 3
  integer, public, parameter :: WM_CYCLOTRONIC      = 4
  integer, public, parameter :: WM_BORIS_PATACCHINI = 5
  integer, public, parameter :: WM_BORIS_LEAP_FROG  = 6
  integer, public, parameter :: WM_BORIS_TANALPHA   = 7
  integer, public, parameter :: WM_TAJIMA_LEAP_FROG_IMPLICIT = 8
  integer, public, parameter :: WM_TAJIMA_LEAP_FROG_EXPLICIT = 9
  integer, public, parameter :: WM_MATRIX_VERLET             = 10
  integer, public, parameter :: WM_CYCLOTRONIC_NOTAN         = 11
  integer, public, parameter :: WM_BORIS_PATACCHINI_NOTAN    = 12
  integer, public, parameter :: WM_BORIS_MLSDC               = 13
  integer, public, parameter :: WM_BORIS_SDC                 = 14

  integer, public, parameter :: IFILE_SUMMAND_PARTICLES     =  46
  integer, public, parameter :: IFILE_SUMMAND_ENERGY        = 146
  integer, public, parameter :: IFILE_SUMMAND_PARTICLES_AVG = 246
  integer, public, parameter :: IFILE_SUMMAND_NFEVAL        = 346
  integer, public, parameter :: IFILE_SUMMAND_NITER         = 446


  !> parameter collection for pepcboris
  type pepcboris_nml_t
    ! MPI variables
    integer(kind_pe) :: rank_space, nrank_space
    integer(kind_pe) :: rank_world, nrank_world !< rank/num_ranks in MPI_COMM_WORLD (used for globally unique output)
    integer(kind_pe) :: rank_time,  nrank_time
    integer(kind_default) :: comm_space, comm_time
    ! time variables
    real(kind_physics) :: dt  !< timestep (in simunits), set via pfasst parameters
    integer :: nt  !< number of timesteps, set via pfasst parameters
    ! dump step modulus
    integer :: dumpstep = 1
    ! configuration variables
    integer :: particle_config = 0
    integer :: dumptype = 0
    real(kind_physics) :: setup_params(PARAMS_MAXIDX) = 0.
    ! number of particles per species and rank, will be set automatically later
    integer(kind_particle) :: numparts = 100
    ! use PFASST
    integer :: workingmode = WM_BORIS_SDC
    !
    logical :: root_stdio
    logical :: root_file
  end type

  type(pepcboris_nml_t), public :: pepcboris_nml

  contains

  subroutine get_mpi_rank(comm, rank, size)
    use module_pepc_types
    implicit none
    include 'mpif.h'
    integer(kind_default), intent(in) :: comm
    integer(kind_pe), intent(out) :: rank, size
    integer(kind_default) :: mpi_err

    call MPI_COMM_RANK( comm, rank, mpi_err )
    call MPI_COMM_SIZE( comm, size, mpi_err )

  end subroutine


  subroutine pepcboris_init(particles, dt, nt)
    implicit none
    type(t_particle), allocatable, target, intent(out) :: particles(:)
    real(kind_physics), intent(in)  :: dt
    integer, intent(in) :: nt

    call set_parameter(dt, nt)
    call setup_particles(particles)

  end subroutine

  !> computes a x b
  pure function cross_prod(a, b)
    implicit none

    real(kind_physics), dimension(3), intent(in) :: a, b
    real(kind_physics), dimension(3) :: cross_prod

    cross_prod(1) = a(2) * b(3) - a(3) * b(2)
    cross_prod(2) = a(3) * b(1) - a(1) * b(3)
    cross_prod(3) = a(1) * b(2) - b(2) * a(1)
  end function cross_prod

  !> computes a x b + c
  pure function cross_prod_plus(a, b, c)
    implicit none

    real(kind_physics), dimension(3), intent(in) :: a, b, c
    real(kind_physics), dimension(3) :: cross_prod_plus

    cross_prod_plus(1) = a(2) * b(3) - a(3) * b(2) + c(1)
    cross_prod_plus(2) = a(3) * b(1) - a(1) * b(3) + c(2)
    cross_prod_plus(3) = a(1) * b(2) - b(2) * a(1) + c(3)
  end function cross_prod_plus

  !> set initial values for particle positions and velocities in y0
  subroutine setup_particles(particles)
    use module_debug
    use module_mirror_boxes
    implicit none
    type(t_particle), allocatable, target, intent(out) :: particles(:)
    integer(kind_particle) :: i
    real(kind_physics) :: myrand(3)

    ! FIXME: Random numbers used here are not guaranteed to be repeatable in subsequent runs (but at least this works for cub, juropa, zam862)
    select case (pepcboris_nml%particle_config)
      case (0)
        ! single particle at specified position
        pepcboris_nml%numparts = 1
        allocate(particles(pepcboris_nml%numparts))
        particles(1)%x      = pepcboris_nml%setup_params(PARAMS_X0:PARAMS_Z0)
        particles(1)%data%q = 1.
        particles(1)%data%m = 1.
        particles(1)%data%v = pepcboris_nml%setup_params(PARAMS_VX0:PARAMS_VZ0)
        particles(1)%label  = 1
        particles(1)%work   = 1.
      case (1)
        ! a bunch of particles around specified position
        allocate(particles(pepcboris_nml%numparts))

        do i=1,pepcboris_nml%numparts
          call random_number(myrand)
          particles(i)%x      = pepcboris_nml%setup_params(PARAMS_X0:PARAMS_Z0) &
                              + pepcboris_nml%setup_params(PARAMS_RADIUS)*(myrand - 0.5_8)
          particles(i)%data%q = 1.
          particles(i)%data%m = 1.
          call random_number(myrand)
          particles(i)%data%v = pepcboris_nml%setup_params(PARAMS_VX0:PARAMS_VZ0) &
                              + pepcboris_nml%setup_params(PARAMS_VELOCITY_SPREAD)*(myrand - 0.5_8)
          particles(i)%label  = i
          particles(i)%work   = 1.
        end do
      case (2)
        ! randomly placed particles
        allocate(particles(pepcboris_nml%numparts))

        do i=1,pepcboris_nml%numparts
          call random_number(myrand)
          particles(i)%x      = pepcboris_nml%setup_params(PARAMS_RADIUS)*(myrand - 0.5_8)
          particles(i)%data%q = 1.
          particles(i)%data%m = 1.
          call random_number(myrand)
          particles(i)%data%v = pepcboris_nml%setup_params(PARAMS_VELOCITY_SPREAD)*(myrand - 0.5_8)
          particles(i)%label  = i
          particles(i)%work   = 1.
        end do
      case default
        DEBUG_ERROR(*, 'setup_particles() - invalid value for particle_config:', pepcboris_nml%particle_config)
    end select

  end subroutine

  function get_magnetic_field()
    implicit none
    real(kind_physics), dimension(3) :: get_magnetic_field

    get_magnetic_field = [0.0_kind_physics, 0.0_kind_physics, pepcboris_nml%setup_params(PARAMS_OMEGAB)] ! FIXME * m/Q
  end function

  subroutine set_parameter(dt, nt)
    use module_pepc
    use module_debug
    use module_interaction_specific, only : theta2, eps2
    use treevars, only : num_threads, np_mult
    implicit none

    real(kind_physics), intent(in)  :: dt
    integer, intent(in) :: nt

    integer, parameter :: fid = 12
    character(255)     :: para_file
    logical            :: read_para_file

    integer :: particle_config, workingmode, dumptype, dumpstep
    integer(kind_particle) :: numparts
    real(kind_physics) :: setup_params(PARAMS_MAXIDX)

    namelist /pepcborispfasst/ particle_config, setup_params, workingmode, dumptype, numparts, dumpstep

    ! frontend parameters
    particle_config = pepcboris_nml%particle_config
    setup_params    = pepcboris_nml%setup_params
    workingmode     = pepcboris_nml%workingmode
    dumptype        = pepcboris_nml%dumptype
    numparts        = pepcboris_nml%numparts
    dumpstep        = pepcboris_nml%dumpstep

    ! pepc parameters
    theta2      = 0.36
    num_threads = 8
    np_mult     = -500
    ! read in namelist file
    call pepc_read_parameters_from_first_argument(read_para_file, para_file)

    if (read_para_file) then
      if(pepcboris_nml%rank_world==0) write(*,'(a)') ' == reading parameter file, section pepcborispfasst: ', para_file
      open(fid,file=para_file)
      read(fid,NML=pepcborispfasst)
      close(fid)
    else
      if(pepcboris_nml%rank_world==0) write(*,*) ' == no param file, using default parameters '
    end if

    ! frontend parameters
    pepcboris_nml%particle_config = particle_config
    pepcboris_nml%setup_params    = setup_params
    pepcboris_nml%workingmode     = workingmode
    pepcboris_nml%dumptype        = dumptype
    pepcboris_nml%numparts        = numparts
    pepcboris_nml%dumpstep        = dumpstep

    ! derived from pfasst parameters
    pepcboris_nml%dt = dt
    pepcboris_nml%nt = nt

    if(pepcboris_nml%rank_world==0) then
      write(*,'(a,i12)')       ' == particle config      : ', pepcboris_nml%particle_config
      write(*,'(a,i12)')       ' == number of time steps : ', pepcboris_nml%nt
      write(*,'(a,es12.4)')    ' == time step            : ', pepcboris_nml%dt
      write(*,'(a,es12.4)')    ' == final time           : ', pepcboris_nml%dt*pepcboris_nml%nt
      write(*,'(a,i12)')       ' == dumptype             : ', pepcboris_nml%dumptype
      write(*,'(a,i12)')       ' == dumpstep             : ', pepcboris_nml%dumpstep
    end if

    call pepc_prepare(dim)
  end subroutine set_parameter

end module
