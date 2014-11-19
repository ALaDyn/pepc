! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
!
! Copyright (C) 2002-2014 Juelich Supercomputing Centre,
!                         Forschungszentrum Juelich GmbH,
!                         Germany
!
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

!!!!!!!!!!!!!!!!!!!!
!! helper module
!!!!!!!!!!!!!!!!!!!!

module helper

  use module_pepc_types
  use variables
  use module_cmdline

  public :: QsortC
  private :: Partition

  contains

    recursive subroutine QsortC(A)
        real*8, intent(in out), dimension(:) :: A
        integer :: iq

        if(size(A) > 1) then
            call Partition(A, iq)
            call QsortC(A(:iq-1))
            call QsortC(A(iq:))
        endif
    end subroutine QsortC

    subroutine Partition(A, marker)
        real*8, intent(in out), dimension(:) :: A
        integer, intent(out) :: marker
        integer :: i, j
        real*8 :: temp
        real*8 :: x      ! pivot point
        x = A(1)
        i= 0
        j= size(A) + 1

        DO
            j = j-1
            DO
                IF (A(j) <= x) exit
                j = j-1
            END DO
            i = i+1
            DO
                IF (A(i) >= x) exit
                i = i+1
            END DO
            IF (i < j) THEN
            ! exchange A(i) and A(j)
                temp = A(i)
                A(i) = A(j)
                A(j) = temp
            ELSEIF (i == j) THEN
                marker = i+1
                return
            ELSE
                marker = i
                return
            ENDIF
        END DO

    end subroutine Partition

!======================================================================================
  subroutine linspace(first, last, arr)
    implicit none

    real(KIND=8), intent(in) :: first, last
    real(KIND=8), intent(out) :: arr(:)

    integer :: i,n
    real(KIND=8) :: delta

    n=size(arr)
    delta = (last-first)/(n-1)

    do i=1,n
        arr(i) = first + (i-1)*delta
    end do
  end subroutine

!======================================================================================

    real*8 function norm(a)
        implicit none
        real*8, intent(in),dimension(3) :: a

        norm = sqrt(dotproduct(a,a))

    end function norm

!======================================================================================

    real*8 function dotproduct(a,b)
        implicit none
        real*8, intent(in),dimension(3) :: a,b

        dotproduct=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)

    end function dotproduct

!======================================================================================
  subroutine init_files()
      implicit none

      !if(root) open(unit=recycling_out,file='recycling.out',status='UNKNOWN',position='APPEND')
      if(root) open(unit=out,file='pepcf.out',status='UNKNOWN',position='APPEND')

  end subroutine init_files

!======================================================================================
  subroutine close_files()
      implicit none

      !if(root) close(recycling_out)
      if(root) close(out)

  end subroutine close_files

!======================================================================================
  subroutine flush_files()
      implicit none

      call close_files()
      call init_files()

  end subroutine flush_files

!======================================================================================
  real*8 function get_time()
    implicit none
    include 'mpif.h'

    get_time = MPI_WTIME()

  end function get_time

!======================================================================================
!Doesn't work if a or b are exactly zero. Compare to zero with func real_equal_zero
  logical function real_equal(a,b,epsi)
    implicit none
    real(KIND=8),intent(in) :: a,b,epsi

    if (abs(a-b) <= epsi * max(abs(a),abs(b))) then
       real_equal = .true.
    else
       real_equal = .false.
    end if

  end function real_equal

!======================================================================================
  logical function real_equal_zero(a,epsi)
    implicit none
    real(KIND=8),intent(in) :: a,epsi

    if (abs(a) <= epsi ) then
       real_equal_zero = .true.
    else
       real_equal_zero = .false.
    end if

  end function real_equal_zero

!======================================================================================
!Doesn't work if a or b are exactly zero. Compare to zero with func real_equal_zero
  logical function real_unequal(a,b,epsi)
    implicit none
    real(KIND=8),intent(in) :: a,b,epsi

    real_unequal = (.not. real_equal(a,b,epsi))

  end function real_unequal

!======================================================================================
  logical function real_unequal_zero(a,epsi)
    implicit none
    real(KIND=8),intent(in) :: a,epsi

    real_unequal_zero = (.not. real_equal_zero(a,epsi))

  end function real_unequal_zero

!======================================================================================
  logical function hit_r(p)
    implicit none
    type(t_particle), intent(in) :: p

    if (p%x(1) > xmax .and. p%data%species /= 0) then
        hit_r = .true.
    else
        hit_r = .false.
    end if

  end function hit_r

!======================================================================================
  logical function hit_l(p)
    implicit none
    type(t_particle), intent(in) :: p

    if (p%x(1) < xmin .and. p%data%species /= 0) then
        hit_l = .true.
    else
        hit_l = .false.
    end if

  end function hit_l


!======================================================================================
  logical function hit_side(p)
    implicit none
    type(t_particle), intent(in) :: p

    if (((p%x(2) < ymin).or.(p%x(2) > ymax).or.(p%x(3) < zmin).or.(p%x(3) > zmax)) .and. p%label > 0) then
        hit_side = .true.
    else
        hit_side = .false.
    end if

  end function hit_side



!======================================================================================
  subroutine reallocate_particles(list, oldn, newn)
    implicit none

    type(t_particle), allocatable, intent(inout) :: list(:)
    integer, intent(in) :: oldn, newn

    type(t_particle) :: tmp_p(oldn)

    tmp_p(1:oldn) = list(1:oldn)
    deallocate(list)
    allocate(list(newn))
    list(1:oldn) = tmp_p

  end subroutine


!======================================================================================
  subroutine find_yinx(y,x,l,u,x_l,x_u)
    implicit none

    real(KIND=8), intent(in) :: y
    real(KIND=8), intent(in), dimension(:) :: x
    integer, intent(out) :: l,u
    real(KIND=8), intent(out) :: x_l,x_u

    integer :: n,uold

    n = size(x)
    if ((y<x(1)) .or. (y>x(n))) then
        l=-1
        u=-1
        x_l=0._8
        x_u=0._8
    else
        u = n
        uold = u
        l = 1
        do while (.True.)
            if ((u - l == 1) .and. (x(u) > y)) exit
            if (x(u) > y) then
                uold = u
                u = l +(u-l)/2
            else
                l = u
                u = uold
            end if
        end do
        x_l = x(l)
        x_u = x(u)
    end if
  end subroutine


!=======================================================================================
  subroutine pepc_tree_diagnostics()
      use module_vtk
      use module_vtk_helpers
      use module_pepc

      implicit none

      integer :: vtk_step

      if (step == 0) then
          vtk_step = VTK_STEP_FIRST
      else if (step == nt + startstep) then
          vtk_step = VTK_STEP_LAST
      else
          vtk_step = VTK_STEP_NORMAL
      endif

      call pepc_statistics(step)
      call vtk_write_branches(step,   step*dt, vtk_step, global_tree)
      call vtk_write_spacecurve(step, step*dt, vtk_step, particles)

  end subroutine pepc_tree_diagnostics

end module
