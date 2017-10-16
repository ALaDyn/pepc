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

!>
!> trajectory integration routines
!>
module pepca_integrator
  use module_pepc_kinds
  implicit none
  private
  
    public update_velocities
    public push_particles
  
  contains

    subroutine update_velocities(p, dt)
      use module_pepc_types
      use pepca_units
      use pepca_helper, only: dim
      implicit none

      type(t_particle), allocatable, intent(inout) :: p(:)
      real*8, intent(in) :: dt
      integer(kind_particle) :: ip
      real*8  :: force(1:3)

      do ip=1, size(p, kind=kind_particle)
        force(1:dim)  = p(ip)%data%q/unit_4piepsilon0 * p(ip)%results%e(1:dim)
        force(dim+1:) = 0.
        p(ip)%data%v  = p(ip)%data%v + dt * force / p(ip)%data%m
      end do
    end subroutine update_velocities


    subroutine push_particles(p, dt)
      use module_pepc_types
      use pepca_units
      implicit none

      type(t_particle), allocatable, intent(inout) :: p(:)
      real*8, intent(in) :: dt
      integer(kind_particle) :: ip

      do ip=1, size(p, kind=kind_particle)
        p(ip)%x = p(ip)%x + p(ip)%data%v * dt
      end do
    end subroutine push_particles

end module
