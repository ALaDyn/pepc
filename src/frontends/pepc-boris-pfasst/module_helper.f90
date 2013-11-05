! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
!
! Copyright (C) 2002-2013 Juelich Supercomputing Centre,
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
  integer, public, parameter :: PARAMS_MAXIDX = PARAMS_OMEGAB

  integer, public, parameter :: WM_PFASST   = 1
  integer, public, parameter :: WM_VERLET   = 2
  integer, public, parameter :: WM_ANALYTIC = 3

  !> parameter collection for pepcboris
  type pepcboris_nml_t
    ! MPI variables
    integer(kind_pe) :: rank, nrank
    integer(kind_default) :: comm
    ! time variables
    real*8  :: dt  !< timestep (in simunits), set via pfasst parameters
    integer :: nt  !< number of timesteps, set via pfasst parameters
    ! configuration variables
    integer :: particle_config = 0
    real*8 :: setup_params(PARAMS_MAXIDX) = 0.
    ! number of particles per species and rank, will be set automatically later
    integer(kind_particle) :: numparts
    ! use PFASST
    integer :: workingmode = WM_PFASST
  end type

  type(pepcboris_nml_t), public :: pepcboris_nml

  contains

  subroutine pepcboris_init(nml, particles, dt, nt)
    implicit none
    type(pepcboris_nml_t), intent(inout) :: nml
    type(t_particle), allocatable, target, intent(out) :: particles(:)
    real*8, intent(in)  :: dt
    integer, intent(in) :: nt

    call set_parameter(nml, dt, nt)
    call setup_particles(particles, nml)

  end subroutine

  !> computes a x b
  pure function cross_prod(a, b)
    implicit none

    real*8, dimension(3), intent(in) :: a, b
    real*8, dimension(3) :: cross_prod

    cross_prod(1) = a(2) * b(3) - a(3) * b(2)
    cross_prod(2) = a(3) * b(1) - a(1) * b(3)
    cross_prod(3) = a(1) * b(2) - b(2) * a(1)
  end function cross_prod

  !> computes a x b + c
  pure function cross_prod_plus(a, b, c)
    implicit none

    real*8, dimension(3), intent(in) :: a, b, c
    real*8, dimension(3) :: cross_prod_plus

    cross_prod_plus(1) = a(2) * b(3) - a(3) * b(2) + c(1)
    cross_prod_plus(2) = a(3) * b(1) - a(1) * b(3) + c(2)
    cross_prod_plus(3) = a(1) * b(2) - b(2) * a(1) + c(3)
  end function cross_prod_plus

  !> set initial values for particle positions and velocities in y0
  subroutine setup_particles(particles, nml)
    use module_debug
    use module_mirror_boxes
    implicit none
    type(pepcboris_nml_t), intent(inout) :: nml
    type(t_particle), allocatable, target, intent(out) :: particles(:)

    select case (nml%particle_config)
      case (0)
        ! single particle at specified position
        nml%numparts = 1
        allocate(particles(1))
        particles(1)%x      = nml%setup_params(PARAMS_X0:PARAMS_Z0)
        particles(1)%data%q = 1.
        particles(1)%data%m = 1.
        particles(1)%data%v = nml%setup_params(PARAMS_VX0:PARAMS_VZ0)
        particles(1)%label  = 1
        particles(1)%work   = 1.
      case default
        DEBUG_ERROR(*, 'setup_particles() - invalid value for particle_config:', nml%particle_config)
    end select

  end subroutine

  function get_magnetic_field()
    implicit none
    real*8, dimension(3) :: get_magnetic_field

    get_magnetic_field = [0.0_8, 0.0_8, pepcboris_nml%setup_params(PARAMS_OMEGAB)] ! FIXME * m/Q
  end function

  subroutine set_parameter(nml, dt, nt)
    use module_pepc
    use module_debug
    use module_interaction_specific, only : theta2, eps2
    use treevars, only : num_threads, np_mult
    implicit none

    type(pepcboris_nml_t), intent(inout) :: nml
    real*8, intent(in)  :: dt
    integer, intent(in) :: nt

    integer, parameter :: fid = 12
    character(255)     :: para_file
    logical            :: read_para_file

    integer :: particle_config, workingmode
    real*8 :: setup_params(PARAMS_MAXIDX)

    namelist /pepcborispfasst/ particle_config, setup_params, workingmode

    ! frontend parameters
    particle_config = nml%particle_config
    setup_params    = nml%setup_params
    workingmode     = nml%workingmode

    ! pepc parameters
    theta2      = 0.36
    num_threads = 8
    np_mult     = -500
    ! read in namelist file
    call pepc_read_parameters_from_first_argument(read_para_file, para_file)

    if (read_para_file) then
      if(nml%rank==0) write(*,'(a)') ' == reading parameter file, section pepcborispfasst: ', para_file
      open(fid,file=para_file)
      read(fid,NML=pepcborispfasst)
      close(fid)
    else
      if(nml%rank==0) write(*,*) ' == no param file, using default parameters '
    end if

    ! frontend parameters
    nml%particle_config = particle_config
    nml%setup_params    = setup_params
    nml%workingmode     = workingmode

    ! derived from pfasst parameters
    nml%dt      = dt
    nml%nt      = nt

    if(nml%rank==0) then
      write(*,'(a,i12)')       ' == particle config      : ', nml%particle_config
      write(*,'(a,i12)')       ' == number of time steps : ', nml%nt
      write(*,'(a,es12.4)')    ' == time step            : ', nml%dt
      write(*,'(a,es12.4)')    ' == final time           : ', nml%dt*nml%nt
    end if

    call pepc_prepare(dim)
  end subroutine set_parameter

end module
