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

    public update_velocities_boris
    public push_particles
    public halfbackstep_velocities_boris
    public halfforwardstep_velocities_boris
    public update_velocities_velocity_verlet_boris
    public push_particles_velocity_verlet_boris

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

   !> velocity half step back: first rotate v by q/m B0*dt/2 and then accelerate one half step back, i.e. with -dt/2
   !> compare Birdsall/Langdon, chapter 2-4, pg. 14
   subroutine halfbackstep_velocities_boris(p, dt)
      use module_pepc_types
      use pepcboris_helper
      implicit none

      real*8, intent(in) :: dt
      type(t_particle), intent(inout) :: p(:)

      integer(kind_particle) :: ip
      real*8 :: beta, gam
      real*8, dimension(3) :: uminus, uprime, uplus, t, s

      real*8, dimension(3) :: B0
      real*8 :: absB0, ttilda

      B0 = get_magnetic_field()
      absB0 = sqrt(dot_product(B0, B0))

      do ip = 1, size(p, kind=kind_particle)
        ! charge/mass*time-constant
        beta   = p(ip)%data%q / (2._8 * p(ip)%data%m) * dt
        ! no initial half step with electric field
        uminus(:) = p(ip)%data%v(:)
        ! gamma factor
        !gam    = sqrt( 1._8 + ( dot_product(uminus, uminus) ) / unit_c2 )
        gam    = 1._8
        ! half rotation with magnetic field backwards
        ttilda = beta/gam * absB0
        t      = tan(atan(ttilda)/2._8) * B0 / absB0
        uprime = cross_prod_plus(uminus, t, uminus)
        s      = 2._8 * t / (1._8 + dot_product(t, t))
        uplus  = cross_prod_plus(uprime, s, uminus)
        ! second half step with electric field
        p(ip)%data%v(:) = uplus(:) - beta * p(ip)%results%e(:)
     end do

   end subroutine halfbackstep_velocities_boris

   subroutine halfforwardstep_velocities_boris(p, dt)
      use module_pepc_types
      use pepcboris_helper
      implicit none

      real*8, intent(in) :: dt
      type(t_particle), intent(inout) :: p(:)

      integer(kind_particle) :: ip
      real*8 :: beta, gam
      real*8, dimension(3) :: uminus, uprime, uplus, t, s

      real*8, dimension(3) :: B0
      real*8 :: absB0, ttilda

      B0 = get_magnetic_field()
      absB0 = sqrt(dot_product(B0, B0))

      do ip = 1, size(p, kind=kind_particle)
        ! charge/mass*time-constant
        beta   = p(ip)%data%q / (2._8 * p(ip)%data%m) * dt
        ! initial half step with electric field
        uminus(:) = p(ip)%data%v(:) + beta * p(ip)%results%e(:)
        ! gamma factor
        !gam    = sqrt( 1._8 + ( dot_product(uminus, uminus) ) / unit_c2 )
        gam    = 1._8
        ! half rotation with magnetic field backwards
        ttilda = beta/gam * absB0
        t      = -tan(atan(ttilda)/2._8) * B0 / absB0
        uprime = cross_prod_plus(uminus, t, uminus)
        s      = 2._8 * t / (1._8 + dot_product(t, t))
        uplus  = cross_prod_plus(uprime, s, uminus)
        ! no second half step with electric field
        p(ip)%data%v(:) = uplus(:)! + beta * p(ip)%results%e(:)
     end do

   end subroutine halfforwardstep_velocities_boris

end module
