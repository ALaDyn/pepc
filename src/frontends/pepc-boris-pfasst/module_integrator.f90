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
module pepcboris_integrator
  implicit none
  private

    ! these routines assume that x and v are lagging behind each other by one half step
    public update_velocities_boris
    public push_particles
    ! for these routines, x and v are assumed to live on full and identical time steps
    public update_velocities_velocity_verlet_boris
    public push_particles_velocity_verlet_boris
    ! dito, cyclotronic integrator from J.Comp.Phys. 228 (2009), 2604-2615
    public drift_cyclotronic
    public kick_cyclotronic

    real*8, allocatable, public :: eold(:,:)

  contains


    subroutine push_particles_velocity_verlet_boris(p, dt)
      use module_pepc_types
      use pepcboris_helper
      implicit none
      type(t_particle), intent(inout) :: p(:)
      real*8, intent(in) :: dt
      integer(kind_particle) :: ip
      real*8, dimension(3) :: B0
      real*8 :: beta

      B0 = get_magnetic_field()

      if (.not. allocated(eold)) allocate(eold(1:3, size(p)))

      do ip = 1, size(p, kind=kind_particle)
        ! charge/mass*time-constant
        beta   = p(ip)%data%q / (2._8 * p(ip)%data%m) * dt

        p(ip)%x(:) = p(ip)%x(:) + dt * ( p(ip)%data%v(:) + beta * cross_prod_plus(p(ip)%data%v, B0, p(ip)%results%e) )

        eold(1:3, ip) = p(ip)%results%e(:)
      end do
    end subroutine

   subroutine update_velocities_velocity_verlet_boris(p, dt)
      use module_pepc_types
      use pepcboris_helper
      implicit none

      real*8, intent(in) :: dt
      type(t_particle), intent(inout) :: p(:)

      integer(kind_particle) :: ip
      real*8 :: beta, gam
      real*8, dimension(3) :: uminus, uprime, uplus, t, s

      real*8, dimension(3) :: B0, Eavg

      B0 = get_magnetic_field()

      do ip = 1, size(p, kind=kind_particle)
        Eavg(:) = (p(ip)%results%e(:) + eold(:,ip)) / 2._8

        ! charge/mass*time-constant
        beta   = p(ip)%data%q / (2._8 * p(ip)%data%m) * dt
        ! first half step with electric field
        uminus(:) = p(ip)%data%v(:) + beta * Eavg(:)
        ! gamma factor
        !gam    = sqrt( 1._8 + ( dot_product(uminus, uminus) ) / unit_c2 )
        gam    = 1._8
        ! rotation with magnetic field
        t      = beta/gam * B0
        uprime = cross_prod_plus(uminus, t, uminus)
        s      = 2._8 * t / (1._8 + dot_product(t, t))
        uplus  = cross_prod_plus(uprime, s, uminus)
        ! second half step with electric field
        p(ip)%data%v(:) = uplus(:) + beta * Eavg(:)
      end do

   end subroutine update_velocities_velocity_verlet_boris


   subroutine push_particles(p, dt)
      use module_pepc_types
      implicit none
      type(t_particle), intent(inout) :: p(:)
      real*8, intent(in) :: dt
      real*8 :: gam
      integer(kind_particle) :: ip

      do ip = 1, size(p, kind=kind_particle)
        ! gam = sqrt(1._8 + dot_product(p(ip)%data%v * p(ip)%data%v) / unit_c2)
        gam = 1._8
        p(ip)%x = p(ip)%x + p(ip)%data%v / gam * dt
      end do
   end subroutine

   subroutine update_velocities_boris(p, dt)
      use module_pepc_types
      use pepcboris_helper
      implicit none

      real*8, intent(in) :: dt
      type(t_particle), intent(inout) :: p(:)

      integer(kind_particle) :: ip
      real*8 :: beta, gam
      real*8, dimension(3) :: uminus, uprime, uplus, t, s

      real*8, dimension(3) :: B0

      B0 = get_magnetic_field()

      do ip = 1, size(p, kind=kind_particle)
        ! charge/mass*time-constant
        beta   = p(ip)%data%q / (2._8 * p(ip)%data%m) * dt
        ! first half step with electric field
        uminus(:) = p(ip)%data%v(:) + beta * p(ip)%results%e(:)
        ! gamma factor
        !gam    = sqrt( 1._8 + ( dot_product(uminus, uminus) ) / unit_c2 )
        gam    = 1._8
        ! rotation with magnetic field
        t      = beta/gam * B0
        uprime = cross_prod_plus(uminus, t, uminus)
        s      = 2._8 * t / (1._8 + dot_product(t, t))
        uplus  = cross_prod_plus(uprime, s, uminus)
        ! second half step with electric field
        p(ip)%data%v(:) = uplus(:) + beta * p(ip)%results%e(:)
      end do

   end subroutine update_velocities_boris


   subroutine kick_cyclotronic(p, dt)
      use module_pepc_types
      implicit none
      type(t_particle), intent(inout) :: p(:)
      real*8, intent(in) :: dt
      integer(kind_particle) :: ip

      do ip = 1, size(p, kind=kind_particle)
        p(ip)%data%v = p(ip)%data%v + p(ip)%data%q / p(ip)%data%m * p(ip)%results%e * dt
      end do
   end subroutine

   subroutine drift_cyclotronic(p, dt)
      use module_pepc_types
      use pepcboris_helper
      use module_debug
      implicit none
      type(t_particle), intent(inout) :: p(:)
      real*8, intent(in) :: dt
      integer(kind_particle) :: ip
      real*8 :: vprime(2), phi, Omega, B0(3), c, s

      B0 = get_magnetic_field()

      if (any(abs(B0(1:2))>0.)) then
        DEBUG_WARNING(*, 'Cyclotronic integrator does only work if B=B*e_z, but B=', B0)
      endif

      do ip = 1, size(p, kind=kind_particle)
        Omega = p(ip)%data%q / p(ip)%data%m * B0(3)
        phi = Omega*dt
        c = cos(phi)
        s = sin(phi)

        vprime(1) =  p(ip)%data%v(1)*c + p(ip)%data%v(2)*s
        vprime(2) = -p(ip)%data%v(1)*s + p(ip)%data%v(2)*c

        p(ip)%x(1) = p(ip)%x(1) + (p(ip)%data%v(2) - vprime(2))/Omega
        p(ip)%x(2) = p(ip)%x(2) + (vprime(1) - p(ip)%data%v(1))/Omega
        p(ip)%x(3) = p(ip)%x(3) + p(ip)%data%v(3)*dt

        p(ip)%data%v(1:2) = vprime(1:2)
      end do
   end subroutine


end module
