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
!> trajectory integration routines
!>
module pepca_integrator
  implicit none
  private
  
    public update_velocities
    public push_particles
  
  contains

    subroutine update_velocities(p, dt)
      use module_pepc_types
      use pepca_units
      use pepca_globals, only: dim
      implicit none

      type(t_particle), allocatable, intent(inout) :: p(:)
      real*8, intent(in) :: dt
      integer(kind_particle) :: ip
      real*8  :: force(1:3)

      do ip=1, size(p, kind=kind_particle)
        ! v represents momentum p/mc = gamma*v/c
        force(1:dim)  = p(ip)%data%q/unit_4piepsilon0 * p(ip)%results%e(1:dim)
        force(dim+1:) = 0.
        p(ip)%data%v  = p(ip)%data%v + dt * force / ( p(ip)%data%m * unit_c )
      end do
    end subroutine update_velocities


    subroutine push_particles(p, dt)
      use module_pepc_types
      use pepca_units
      implicit none

      type(t_particle), allocatable, intent(inout) :: p(:)
      real*8, intent(in) :: dt
      integer(kind_particle) :: ip
      real*8  :: gam

      do ip=1, size(p, kind=kind_particle)
        ! v represents momentum p/mc = gamma*v/c
        ! gamma = sqrt(1+ (p/mc)^2)
        gam     = sqrt( 1.0 + dot_product(p(ip)%data%v, p(ip)%data%v) )
        p(ip)%x = p(ip)%x + p(ip)%data%v * unit_c / gam * dt
      end do
    end subroutine push_particles

end module
