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

    public push_particles_boris

  contains

   subroutine push_particles_boris(nml, p, dt)
      use module_pepc_types
      use pepca_helper
      implicit none

      type(pepca_nml_t), intent(in) :: nml
      real*8, intent(in) :: dt
      type(t_particle), intent(inout) :: p(:)

      integer(kind_particle) :: ip
      real*8 :: beta, gam
      real*8, dimension(3) :: uminus, uprime, uplus, t, s

      real*8, dimension(3) :: B0

      B0 = [0.0_8, 0.0_8, nml%setup_params(PARAMS_B0)]

      do ip = 1, size(p)
        ! charge/mass*time-constant
        beta   = p(ip)%data%q / (2._8 * p(ip)%data%m) * dt
        ! first half step with electric field
        uminus(:) = p(ip)%data%v(:) + beta * p(ip)%results%e(:)
        ! gamma factor
        !gam    = sqrt( 1._8 + ( dot_product(uminus, uminus) ) / unit_c2 )
        gam    = 1._8
        ! rotation with magnetic field
        t      = beta/gam * B0
        uprime = uminus + cross_product(uminus, t)
        s      = 2._8 * t / (1._8 + dot_product(t, t))
        uplus  = uminus + cross_product(uprime, s)
        ! second half step with electric field
        p(ip)%data%v(:) = uplus(:) + beta * p(ip)%results%e(:)
        ! gam = sqrt(1._8 + dot_product(p(ip)%data%v * p(ip)%data%v) / unit_c2)
        gam = 1._8
        p(ip)%x = p(ip)%x + p(ip)%data%v / gam * dt
      end do

      contains

      pure function cross_product(a, b)
        implicit none

        real*8, dimension(3), intent(in) :: a, b
        real*8, dimension(3) :: cross_product

        cross_product(1) = a(2) * b(3) - a(3) * b(2)
        cross_product(2) = a(3) * b(1) - a(1) * b(3)
        cross_product(3) = a(1) * b(2) - b(2) * a(1)
      end function cross_product

   end subroutine push_particles_boris


end module
