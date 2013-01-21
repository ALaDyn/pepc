! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2013 Juelich Supercomputing Centre, 
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
  use zufall
  use module_cmdline

  contains

!======================================================================================

    real*8 function dotproduct(a,b)
        implicit none
        real*8, intent(in),dimension(3) :: a,b

        dotproduct=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)

    end function dotproduct

!======================================================================================
  subroutine init_files()
      implicit none

      !if(root) open(unit=timing_out,file='timing.out',status='UNKNOWN',position='APPEND')
      !if(root) open(unit=recycling_out,file='recycling.out',status='UNKNOWN',position='APPEND')
      if(root) open(unit=out,file='pepcf.out',status='UNKNOWN',position='APPEND')

  end subroutine init_files

!======================================================================================
  subroutine close_files()
      implicit none

      !if(root) close(timing_out)
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
!sortiert wallparticles nach hinten
  subroutine sort_particles(p)
    implicit none

    type(t_particle), allocatable, intent(inout) :: p(:)
    type(t_particle)                               :: p_help_vorn


    integer :: i,iback

    iback=np

    DO i=1,np
        if (p(i)%data%species == 0) then
            p_help_vorn=p(i)
            do while (p(iback)%data%species==0)
                iback=iback-1
            end do
            if (iback<=i) exit
            p(i)=p(iback)
            p(iback)=p_help_vorn
        end if
    END DO
  end subroutine

!=======================================================================================
  subroutine pepc_tree_diagnostics()
      use module_vtk
      use module_treediags
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
      call write_branches_to_vtk(step,   step*dt, vtk_step)
      call write_spacecurve_to_vtk(step, step*dt, vtk_step, particles)

  end subroutine pepc_tree_diagnostics

end module
