! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
!
! Copyright (C) 2002-2016 Juelich Supercomputing Centre,
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
module helper
  use module_pepc_kinds
  use module_pepc_types
  use module_timings
  use module_random
  implicit none

  ! timing variables
  integer, parameter :: t_user_total       = t_userdefined_first
  integer, parameter :: t_user_init        = t_userdefined_first + 1
  integer, parameter :: t_user_step        = t_userdefined_first + 2
  integer, parameter :: t_user_directsum   = t_userdefined_first + 3
  integer, parameter :: t_user_particleio  = t_userdefined_first + 4

  ! MPI variables
  integer(kind_pe) :: my_rank, n_ranks
  logical :: root

  ! time variables
  real*8 :: dt
  integer :: step

  ! control variables
  integer :: nt                   ! number of timesteps
  integer(kind_particle) :: tnp   ! total number of particles
  integer(kind_particle) :: np    ! local number of particles
  logical :: particle_output      ! turn vtk output on/off
  logical :: domain_output        ! turn vtk output on/off
  logical :: particle_test        ! check tree code results against direct summation
  logical :: reflecting_walls     ! reflect particles at walls
  integer :: diag_interval        ! number of timesteps between all diagnostics and IO
  real(kind_physics) :: plasma_dimensions(3) ! size of the simulation box

  integer, parameter :: particle_direct = -1 ! number of particle for direct summation

  ! particle data (position, velocity, mass, charge)
  type(t_particle), allocatable   :: particles(:)
  real(kind_physics), allocatable :: direct_L2(:)
  
  ! PRNG state, PER MPI RANK, NOT THREAD!
  type(random_state_t) :: rng_state

  contains

  subroutine set_parameter()

    use module_pepc
    use module_interaction_specific, only : theta2, eps2, force_law, include_far_field_if_periodic
    implicit none

    integer, parameter :: fid = 12
    character(255)     :: para_file
    logical            :: read_para_file

    namelist /pepcbenchmark/ tnp, nt, dt, particle_output, domain_output, reflecting_walls, particle_test, diag_interval, plasma_dimensions

    ! set default parameter values
    tnp               = 10000
    nt                = 25
    dt                = 1e-2
    particle_test     = .false.
    particle_output   = .false.
    domain_output     = .false.
    reflecting_walls  = .false.
    diag_interval     = 1
    plasma_dimensions = (/ 1.0_8, 1.0_8, 1.0_8 /)

    ! read in namelist file
    call pepc_read_parameters_from_first_argument(read_para_file, para_file)

    if (read_para_file) then
      if(root) write(*,'(a)') " == reading parameter file, section pepcbenchmark: ", para_file
      open(fid,file=para_file)
      read(fid,NML=pepcbenchmark)
      close(fid)
    else
      if(root) write(*,*) " == no param file, using default parameter "
    end if

    if(root) then
      write(*,'(a,i12)')       " == total number of particles : ", tnp
      write(*,'(a,i12)')       " == number of time steps      : ", nt
      write(*,'(a,es12.4)')    " == time step                 : ", dt
      write(*,'(a,i12)')       " == diag & IO interval        : ", diag_interval
      write(*,'(a,l12)')       " == particle test             : ", particle_test
      write(*,'(a,l12)')       " == particle output           : ", particle_output
      write(*,'(a,l12)')       " == domain output             : ", domain_output
      write(*,'(a,l12)')       " == reflecting walls          : ", reflecting_walls
      write(*,'(a,3(es12.4))') " == plasma dimensions         : ", plasma_dimensions
    end if

    call pepc_prepare(3_kind_dim)
  end subroutine set_parameter


  subroutine init_particles(p)
    implicit none

    type(t_particle), allocatable, intent(inout) :: p(:)
    integer(kind_particle) :: ip
    integer :: rc
    real*8 :: dummy
    real(kind_physics) :: pi, phi, theta, rnd(1:2), s

    if(root) write(*,'(a)') " == [init] init particles "

    ! set initially number of local particles
    np = tnp / n_ranks
    if (my_rank < MOD(tnp, 1_kind_particle*n_ranks)) np = np + 1

    allocate(particles(np), stat=rc)
    if (rc.ne.0) write(*,*) " === particle allocation error!"

    allocate(direct_L2(np), stat=rc)
    if (rc.ne.0) write(*,*) " === direct_L2 allocation error!"
    direct_L2 = -1.0_8

    ! set random seed to MPI rank
    rng_state%idum = -(my_rank + 1)
    call random(dummy, rng_state)

    pi = 2.0_kind_physics*acos(0.0_kind_physics)

    ! setup random qubic particle cloud
#define COUEX
#ifdef COUEX
    do ip=1, np
       p(ip)%label       = my_rank * (tnp / n_ranks) + ip - 1
       p(ip)%data%q      = -1.0_8 * 2.0_8 * &
                            plasma_dimensions(1) * plasma_dimensions(2) * &
                            plasma_dimensions(3) / tnp
       p(ip)%data%m      = 1.d-2
       if(p(ip)%data%q .gt. 0.0) p(ip)%data%m = p(ip)%data%m * 100.0_8

       ! obtain 2 random numbers
       call random(rnd, rng_state)

       ! get angles to point into sphere
       phi               = rnd(1) * 2*pi
       theta             = acos(rnd(2) * 2 - 1)

       ! compute a radial coordinate from distribution Eq.(9) in Kaplan's paper, \mu = 1
       s                 = ( 1./(1-real(p(ip)%label)/tnp) - 1 )**(1./3.)

       ! fill cartesian coordinates
       p(ip)%x           = s * [cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta)] * plasma_dimensions

       p(ip)%data%v      = 0.0_8

       p(ip)%work        = 1.0_8
    end do
#else
    do ip=1, np
       p(ip)%label       = my_rank * (tnp / n_ranks) + ip - 1
       p(ip)%data%q      = (-1.0_8 + 2.0_8*MOD(p(ip)%label,2_kind_particle)) * 2.0_8 * &
                            plasma_dimensions(1) * plasma_dimensions(2) * &
                            plasma_dimensions(3) / tnp
       p(ip)%data%m      = 1.0_8
       if(p(ip)%data%q .gt. 0.0) p(ip)%data%m = p(ip)%data%m * 100.0_8

       call random(p(ip)%x, rng_state)
       p(ip)%x           = p(ip)%x * plasma_dimensions

       call random_gauss(p(ip)%data%v, rng_state)
       p(ip)%data%v      = p(ip)%data%v / sqrt(p(ip)%data%m)

       p(ip)%work        = 1.0_8
    end do
#endif
  end subroutine init_particles


  subroutine push_particles(p)
    use module_mirror_boxes
    implicit none

    type(t_particle), allocatable, intent(inout) :: p(:)
    integer(kind_particle) :: ip
    real*8  :: fact

    if(root) write(*,'(a)') " == [pusher] push particles "

    fact = dt

    do ip=1, np
      p(ip)%data%v = p(ip)%data%v + fact * p(ip)%data%q / p(ip)%data%m * p(ip)%results%e
      p(ip)%x      = p(ip)%x      + dt   * p(ip)%data%v
    end do
  end subroutine push_particles


  subroutine filter_particles(p)
    implicit none
    include 'mpif.h'

    type(t_particle), allocatable, intent(inout) :: p(:)
    integer(kind_particle) :: ip
    integer :: id, ncoll, ncoll_total, ierr

    ncoll = 0

    do ip=1, np
      do id=1, 3
        if(p(ip)%x(id) < 0.0_8) then
          p(ip)%x(id) = - p(ip)%x(id)
          p(ip)%data%v(id) = - p(ip)%data%v(id)
          ncoll = ncoll + 1
        else if(p(ip)%x(id) > plasma_dimensions(id)) then
          p(ip)%x(id) = 2 * plasma_dimensions(id) - p(ip)%x(id)
          p(ip)%data%v(id) = - p(ip)%data%v(id)
          ncoll = ncoll + 1
        end if
      end do
    end do

    call mpi_reduce(ncoll, ncoll_total, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    if(root) write(*,'(a,i12)')    " == [filter] total number of wall collisions      : ", ncoll_total
  end subroutine filter_particles


  subroutine test_particles()
    use module_pepc_types
    use module_directsum
    implicit none
    include 'mpif.h'

    integer(kind_particle), allocatable   :: tindx(:)
    real*8, allocatable                   :: trnd(:)
    type(t_particle_results), allocatable :: trslt(:)
    integer(kind_particle)                :: tn, tn_global, ti
    integer                               :: rc
    real(kind_physics)                    :: L2sum_local, L2sum_global, L2

    call timer_start(t_user_directsum)

    if(allocated(direct_L2)) then
      deallocate(direct_L2)
    end if
    allocate(direct_L2(np))
    direct_L2 = -1.0_8

    if (particle_direct .eq. -1) then
       tn = np
    else
       tn = particle_direct / n_ranks
       if(my_rank.eq.(n_ranks-1)) tn = tn + MOD(particle_direct, n_ranks)
    endif

    allocate(tindx(tn), trnd(tn), trslt(tn))

    if (particle_direct .eq. -1) then
       do ti = 1, tn
          tindx(ti) = ti
       enddo
    else
       call random(trnd(1:tn), rng_state)

       tindx(1:tn) = int(trnd(1:tn) * (np-1)) + 1
    endif

    call directforce(particles, tindx, tn, trslt, MPI_COMM_WORLD)

    L2sum_local  = 0.0
    L2sum_global = 0.0
    do ti = 1, tn
      L2          = &
                    (particles(tindx(ti))%results%e(1) - trslt(ti)%e(1))**2+ &
                    (particles(tindx(ti))%results%e(2) - trslt(ti)%e(2))**2+ &
                    (particles(tindx(ti))%results%e(3) - trslt(ti)%e(3))**2
      L2sum_local = L2sum_local + L2
      direct_L2(tindx(ti)) = L2
    end do

    call MPI_ALLREDUCE(tn, tn_global, 1, MPI_KIND_PARTICLE, MPI_SUM, MPI_COMM_WORLD, rc)
    call MPI_ALLREDUCE(L2sum_local, L2sum_global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)

    L2sum_global = sqrt(L2sum_global) / tn_global

    call timer_stop(t_user_directsum)
    if(root) then
      write(*,'(a,i12)')    " == [direct test] number tested particles         : ", tn_global
      write(*,'(a,es12.4)') " == [direct test] L2 error in probed particles    : ", L2sum_global
      write(*,'(a,es12.4)') " == [direct test] time in test [s]                : ", timer_read(t_user_directsum)
    end if

    deallocate(tindx)
    deallocate(trnd)
    deallocate(trslt)
  end subroutine test_particles


  integer function vtk_step_of_step(step) result(vtk_step)
    use module_vtk
    implicit none

    integer, intent(in) :: step

    if (step .eq. 0) then
      vtk_step = VTK_STEP_FIRST
    else if (nt - 1 - step < diag_interval ) then
      vtk_step = VTK_STEP_LAST
    else
      vtk_step = VTK_STEP_NORMAL
    endif
  end function vtk_step_of_step


  subroutine write_particles(p)
    use module_vtk_helpers
    implicit none

    include 'mpif.h'

    type(t_particle), intent(in) :: p(:)

    integer :: vtk_step

    call timer_start(t_user_particleio)
    vtk_step = vtk_step_of_step(step)
    call vtk_write_particles("particles", MPI_COMM_WORLD, step, dt * step, vtk_step, p, coulomb_and_l2)
    call timer_stop(t_user_particleio)
    if(root) write(*,'(a,es12.4)') " == [write particles] time in vtk output [s]      : ", timer_read(t_user_particleio)

    contains

    subroutine coulomb_and_l2(d, r, vtkf)
      use module_vtk
      use module_interaction_specific_types
      implicit none

      type(t_particle_data), intent(in) :: d(:)
      type(t_particle_results), intent(in) :: r(:)
      type(vtkfile_unstructured_grid), intent(inout) :: vtkf

      call vtk_write_particle_data_results(d, r, vtkf)
      if(particle_test) call vtkf%write_data_array("L2 error", direct_L2(:))
    end subroutine
  end subroutine write_particles


  subroutine write_domain(p)
    use module_vtk
    use module_vtk_helpers
    use module_pepc, only : global_tree
    implicit none

    type(t_particle), allocatable, intent(in) :: p(:)

    integer :: vtk_step

    ! output of tree diagnostics
    vtk_step = vtk_step_of_step(step)
    call vtk_write_branches(step,  dt * step, vtk_step, global_tree)
    call vtk_write_leaves(step, dt * step, vtk_step, global_tree)
    call vtk_write_spacecurve(step, dt * step, vtk_step, p)
  end subroutine write_domain


  subroutine random_gauss(list, state)
    implicit none

    real(kind_physics), intent(inout) :: list(:)
    type(random_state_t), intent(inout) :: state

    real(kind_physics) :: v(2), pi, r, p
    integer :: n, i

    pi = 2.0_kind_physics*acos(0.0_kind_physics)
    n  = size(list)

    do i=1, n, 2

      call random(v, state)

      r = sqrt(-2.0_8 * log(v(1)))
      p = 2.0_8*pi*v(2)

      list(i)                = r * sin(p)
      if((i+1)<=n) list(i+1) = r * cos(p)

    end do
  end subroutine

end module
