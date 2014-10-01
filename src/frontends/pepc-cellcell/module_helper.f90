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
module helper
  use module_pepc_kinds
  use module_pepc_types
  use module_timings
  implicit none

  ! timing variables
  integer, parameter :: t_user_total              = t_userdefined_first
  integer, parameter :: t_user_init               = t_userdefined_first + 1
  integer, parameter :: t_user_grow_and_traverse  = t_userdefined_first + 2
  integer, parameter :: t_user_calculate_internal = t_userdefined_first + 3
  integer, parameter :: t_user_directsum          = t_userdefined_first + 4

  ! MPI variables
  integer(kind_pe) :: my_rank, n_ranks
  logical :: root

  integer :: step

  ! control variables
  integer(kind_particle) :: tnp   ! total number of particles
  integer(kind_particle) :: np    ! local number of particles
  real(kind_physics) :: dimensions(3) ! size of the simulation box
  integer :: method ! used methods. 0: pepc_grow_and_traverse, 1: pepc_calculate_internal, 2: both
  real(kind_physics) :: theta_stepsize
  integer :: theta_stepcount

  ! output variables
  logical :: vtk_output ! turn vtk output on/off
  character (len=255) :: output_dir
  character (len=255) :: output_filename
  integer :: times_unit

  interface random
    module procedure random8, random16
  end interface

  contains

  subroutine set_parameter()
    use module_pepc
    use module_utils, only: create_directory
    implicit none

    integer :: params_unit
    character(255)     :: para_file
    logical            :: read_para_file
    namelist /pepccellcell/ tnp, dimensions, method, theta_stepsize, theta_stepcount, vtk_output

    ! set default parameter values
    tnp = 10000
    dimensions = (/ 1.0_8, 1.0_8, 1.0_8 /)
    method = 1
    theta_stepsize = 0.1
    theta_stepcount = 10
    vtk_output = .false.

    ! read in namelist file
    call pepc_read_parameters_from_first_argument(read_para_file, para_file)

    if (read_para_file) then
      if(root) write(*,'(a)') " == reading parameter file, section pepcessential: ", para_file
      open(newunit=params_unit, file=para_file)
      read(params_unit,NML=pepccellcell)
      close(params_unit)
    else
      if(root) write(*,*) " == no param file, using default parameter "
    end if

    if(root) then
      write(*,'(a,i12)')       " == total number of particles : ", tnp
      write(*,'(a,3(es12.4))') " == dimensions                : ", dimensions
      write(*,'(a,i12)')       " == method                    : ", method
    end if

    call pepc_prepare(3_kind_dim)

    write (output_dir, '(a,i0)') "tnp", tnp
    call create_directory(trim(output_dir))
    call create_directory("vtk")
    call create_directory("vtk/"//trim(output_dir))
    write (output_filename, '( a , "/", "output_tnp", i0, ".dat" )' ) trim(output_dir), tnp
    open (newunit=times_unit,file=output_filename)


  end subroutine set_parameter


  subroutine init_particles(p)
    implicit none

    type(t_particle), allocatable, intent(inout) :: p(:)
    integer(kind_particle) :: ip
    integer :: rc
    real*8 :: dummy

    if(root) write(*,'(a)') " == [init] init particles "

    ! set initially number of local particles
    np = tnp / n_ranks
    if (my_rank < MOD(tnp, 1_kind_particle*n_ranks)) np = np + 1

    allocate(p(np), stat=rc)
    if (rc.ne.0) write(*,*) " === particle allocation error!"

    ! set random seed
    dummy = par_rand(1*my_rank)

    ! setup random cubic particle cloud
    do ip=1, np
      p(ip)%label       = my_rank * (tnp / n_ranks) + ip - 1
      p(ip)%data%q      = (-1.0_8 + 2.0_8*MOD(p(ip)%label,2_kind_particle)) * 2.0_8 * &
                            dimensions(1) * dimensions(2) * &
                            dimensions(3) / tnp
      call random(p(ip)%x)
      p(ip)%x           = p(ip)%x * dimensions
      p(ip)%work        = 1.0_8
    end do

  end subroutine init_particles


  subroutine calculate_errors(p, direct_results, tindx, m, mean_relerrs, relerrs)
    implicit none
    include 'mpif.h'

    type(t_particle), intent(in) :: p(:)
    type(t_particle_results), intent(in) :: direct_results(:)
    integer(kind_particle), intent(in) :: tindx(:)
    real(kind_physics), intent(out) :: mean_relerrs(2)
    real(kind_physics), intent(out) :: relerrs(:,:)
    integer, intent(in) :: m

    integer(kind_particle) :: ip
    integer :: rc
    real(kind_physics) :: mean_relerr_local

    mean_relerr_local = 0
    do ip = 1, size(tindx)
      relerrs(tindx(ip), m) = sqrt(dot_product( &
        p(tindx(ip))%results%e - direct_results(ip)%e, &
        p(tindx(ip))%results%e - direct_results(ip)%e)) / sqrt(dot_product(direct_results(ip)%e, direct_results(ip)%e))
      mean_relerr_local = mean_relerr_local + relerrs(tindx(ip), m)
    end do

    call MPI_ALLREDUCE(mean_relerr_local, mean_relerrs(m), 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
    mean_relerrs(m) = mean_relerrs(m) / tnp

  end subroutine calculate_errors


  subroutine calculate_internal(p)
    use module_pepc, only: pepc_calculate_internal
    use module_interaction_specific, only : theta2
    implicit none

    type(t_particle), allocatable, intent(inout) :: p(:)

    call pepc_calculate_internal(p, .true.)
    call timer_copy(t_fields_passes, t_user_calculate_internal)
    if(root) write(*,'(a,es12.4,a,es12.4)') " === theta: ", sqrt(theta2), ", pepc_calculate_internal time [s]: ", timer_read(t_user_calculate_internal)

  end subroutine calculate_internal


  subroutine grow_and_traverse(p)
    use module_pepc, only: pepc_grow_and_traverse
    use module_interaction_specific, only : theta2
    implicit none

    type(t_particle), allocatable, intent(inout) :: p(:)

    call pepc_grow_and_traverse(p, 0)
    call timer_copy(t_fields_passes, t_user_grow_and_traverse)
    if(root) write(*,'(a,es12.4,a,es12.4)') " === theta: ", sqrt(theta2), ", tree grow & traversal time [s]: ", timer_read(t_user_grow_and_traverse)

  end subroutine grow_and_traverse


  subroutine gather_results(relerrs_local, relerrs_global)
    implicit none
    include 'mpif.h'

    real(kind_physics), intent(in) :: relerrs_local(:,:)
    real(kind_physics), intent(out) :: relerrs_global(:,:)

    real(kind_physics), allocatable :: relerrs_recvbuf(:)
    integer(kind_default) :: sendcount
    integer(kind_default) :: recvcounts(n_ranks), displs(n_ranks)
    integer :: i, ierror
    logical, allocatable :: mask(:,:)

    sendcount = size(relerrs_local, 1) * size(relerrs_local, 2)
    call MPI_GATHER(sendcount, 1, MPI_INTEGER, recvcounts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)

    allocate(relerrs_recvbuf(2 * tnp))
    if (root) then
      displs(1) = 0
      do i=2, size(displs)
       displs(i) = displs(i-1) + recvcounts(i-1)
      end do
    end if
    call MPI_GATHERV(reshape(relerrs_local, (/ sendcount, 1 /)), sendcount, MPI_KIND_PHYSICS, relerrs_recvbuf,  recvcounts, displs, MPI_KIND_PHYSICS, 0, MPI_COMM_WORLD, ierror)

    if (root) then
      allocate(mask(2 * tnp, 2))
      mask(:,1) = .false.
      do i=1, n_ranks
        mask(displs(i)+1:displs(i)+recvcounts(i)/2, 1) = .true.
      end do
      mask(:,2) = .not. mask(:,1)

      relerrs_global(:,1) = pack(relerrs_recvbuf, mask(:,1))
      relerrs_global(:,2) = pack(relerrs_recvbuf, mask(:,2))
    end if

    deallocate(relerrs_recvbuf, mask)

  end subroutine


  subroutine write_results(mean_relerrs, relerrs)
    use module_interaction_specific, only : theta2
    implicit none

    real(kind_physics), intent(in) :: mean_relerrs(2)
    real(kind_physics), intent(in) :: relerrs(:,:)

    integer :: i
    integer :: unit1, unit2
    character (len=255) :: output_filename_calc_internal_per_theta, output_filename_grow_traverse_per_theta

    write (times_unit,*) sqrt(theta2), timer_read(t_user_directsum), timer_read(t_user_grow_and_traverse), timer_read(t_user_calculate_internal), mean_relerrs

    write (output_filename_grow_traverse_per_theta, '( a, "/", "output_tnp", i0, "_m", i0, "_theta0", f0.2, ".dat" )' ) trim(output_dir), tnp, 0, sqrt(theta2)
    write (output_filename_calc_internal_per_theta, '( a, "/", "output_tnp", i0, "_m", i0, "_theta0", f0.2, ".dat" )' ) trim(output_dir), tnp, 1, sqrt(theta2)

    !
    select case (method)
      case (0)
        open (newunit=unit1,file=output_filename_grow_traverse_per_theta)
          do i=1, size(relerrs, 1)
            write (unit1,*) relerrs(i, 2)
          end do
        close(unit1)
      case (1)
        open (newunit=unit2,file=output_filename_calc_internal_per_theta)
          do i=1, size(relerrs, 1)
            write (unit2,*) relerrs(i, 1)
          end do
        close(unit2)
      case (2)
        open (newunit=unit1,file=output_filename_grow_traverse_per_theta)
          do i=1, size(relerrs, 1)
            write (unit1,*) relerrs(i, 2)
          end do
        close(unit1)
        open (newunit=unit2,file=output_filename_calc_internal_per_theta)
          do i=1, size(relerrs, 1)
            write (unit2,*) relerrs(i, 1)
          end do
        close(unit2)
    end select

  end subroutine write_results


  subroutine finalize(p)
    use module_pepc, only: pepc_finalize
    implicit none

    type(t_particle), allocatable, intent(inout) :: p(:)

    call timer_stop(t_user_total)
    if(root) then
      write(*,*)            " "
      write(*,'(a)')        " ===== finished pepc simulation"
      write(*,'(a,es12.4)') " ===== total run time [s]: ", timer_read(t_user_total)
    end if

    close(times_unit)
    deallocate(p)

    ! cleanup pepc and MPI
    call pepc_finalize()

  end subroutine finalize


  integer function vtk_step_of_step(step) result(vtk_step)
    use module_vtk
    implicit none

    integer, intent(in) :: step

    if (step .eq. 0) then
      vtk_step = VTK_STEP_FIRST
    else if (theta_stepcount -1 -step == 0) then
      vtk_step = VTK_STEP_LAST
    else
      vtk_step = VTK_STEP_NORMAL
    endif
  end function vtk_step_of_step


  subroutine write_particles(p, relerrs)
    use module_vtk_helpers
    use module_interaction_specific, only : theta2
    implicit none

    include 'mpif.h'

    type(t_particle), intent(in) :: p(:)
    real(kind_physics), intent(in) :: relerrs(:,:)

    integer :: vtk_step

    vtk_step = vtk_step_of_step(step)

    ! TODO: A few of the vtk files (particles.timeseries.pvd, particles.timeseries.visit) are written into the wrong directory and relative paths are broken.
    call vtk_write_particles(trim(output_dir)//"/particles", MPI_COMM_WORLD, step, sqrt(theta2), vtk_step, p, coulomb_and_relerrs)
    if(root) write(*,*) "== written particles in vtk output"

    contains

    subroutine coulomb_and_relerrs(d, r, vtkf)
      use module_vtk
      use module_interaction_specific_types
      implicit none

      type(t_particle_data), intent(in) :: d(:)
      type(t_particle_results), intent(in) :: r(:)
      type(vtkfile_unstructured_grid), intent(inout) :: vtkf

      call vtk_write_particle_data_results(d, r, vtkf)
      call vtkf%write_data_array("rel. error Dual", relerrs(:,1))
      call vtkf%write_data_array("rel. error BH", relerrs(:,2))
    end subroutine coulomb_and_relerrs
  end subroutine write_particles


  subroutine random8(array)
    implicit none
    real*8 :: array(:)
    integer :: i

    do i = 1,size(array)
       array(i) = par_rand()
    end do
  end subroutine random8


  subroutine random16(array)
    implicit none
    real*16 :: array(:)
    integer :: i

    do i = 1,size(array)
       array(i) = par_rand()
    end do
  end subroutine random16


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
