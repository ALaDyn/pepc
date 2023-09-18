! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
!
! Copyright (C) 2002-2023 Juelich Supercomputing Centre,
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

module time_integration

   use module_pepc_kinds
   use manipulate_particles
   implicit none

contains
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !>
   !>  Call the chosen integration scheme for particle updtaing
   !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine update_rk(stage, trun)

      use physvars
      integer, intent(in)    :: stage   ! In which RK stage are we?
      real,    intent(inout) :: trun

      if (rk_stages == 2) then

         trun = trun + dt / rk_stages
         call push_rk2(stage)

      elseif (rk_stages == 4) then

         if (stage == 4) trun = trun + dt
         call push_rk4(stage)

      endif

   end subroutine update_rk
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !>
   !>   Update particle positions - using 2nd order RK
   !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine push_rk2(stage)

      use physvars
      integer, intent(in) :: stage   ! In which RK stage are we?
      integer(kind_particle) :: i

      do i = 1, np

         !&<
         if (stage .eq. 1) then
            ! Euler predictor
            vortex_particles(i)%data%x_rk(1:3)     = vortex_particles(i)%x(1:3)
            vortex_particles(i)%data%alpha_rk(1:3) = vortex_particles(i)%data%alpha(1:3)
            vortex_particles(i)%data%u_rk(1:3)     = vortex_particles(i)%results%u(1:3)
            vortex_particles(i)%data%af_rk(1:3)    = vortex_particles(i)%results%af(1:3)
            vortex_particles(i)%x(1:3)             = vortex_particles(i)%x(1:3) &
                                                   + dt * vortex_particles(i)%results%u(1:3)
            vortex_particles(i)%data%alpha(1:3)    = vortex_particles(i)%data%alpha(1:3) &
                                                   + dt * vortex_particles(i)%results%af(1:3)
         else
            ! Trapezoidal corrector
            vortex_particles(i)%x(1:3)             = vortex_particles(i)%data%x_rk(1:3) &
                                                   + 0.5 * dt * (vortex_particles(i)%data%u_rk(1:3) + vortex_particles(i)%results%u(1:3))
            vortex_particles(i)%data%alpha(1:3)    = vortex_particles(i)%data%alpha_rk(1:3) &
                                                   + 0.5 * dt * (vortex_particles(i)%data%af_rk(1:3) + vortex_particles(i)%results%af(1:3))
         end if
         !&>

      end do

      call kick_out_particles()

   end subroutine push_rk2

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !>
   !>   Update particle positions - using 4th order RK
   !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine push_rk4(stage)
      use physvars
      integer, intent(in) :: stage   ! In which RK stage are we?
      integer(kind_particle) :: i
      real(kind_physics), parameter :: rk_d(4) = [1./6.,1./3.,1./3.,1./6.]
      real(kind_physics), parameter :: rk_s(3) = [0.5,0.5,1.0]

      do i = 1, np

         !&<
         if (stage == 1) then
            ! store initial status and initialize derivatives
            vortex_particles(i)%data%x_rk     = vortex_particles(i)%x
            vortex_particles(i)%data%alpha_rk = vortex_particles(i)%data%alpha
            vortex_particles(i)%data%u_rk     = 0.d0
            vortex_particles(i)%data%af_rk    = 0.d0
         end if

         ! store derivatives
         vortex_particles(i)%data%u_rk  = vortex_particles(i)%data%u_rk  + rk_d(stage) * vortex_particles(i)%results%u
         vortex_particles(i)%data%af_rk = vortex_particles(i)%data%af_rk + rk_d(stage) * vortex_particles(i)%results%af

         ! update status
         if (stage /= 4) then
            vortex_particles(i)%x          = vortex_particles(i)%data%x_rk     + rk_s(stage) * dt * vortex_particles(i)%results%u
            vortex_particles(i)%data%alpha = vortex_particles(i)%data%alpha_rk + rk_s(stage) * dt * vortex_particles(i)%results%af
         else
            vortex_particles(i)%x          = vortex_particles(i)%data%x_rk     + dt * vortex_particles(i)%data%u_rk
            vortex_particles(i)%data%alpha = vortex_particles(i)%data%alpha_rk + dt * vortex_particles(i)%data%af_rk
            vortex_particles(i)%results%u  = vortex_particles(i)%data%u_rk
            ! TODO: Psi function should be updated with another call to tree
         end if
         !&>

      end do

      call kick_out_particles()
   end subroutine push_rk4

end module time_integration
