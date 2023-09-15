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

end module time_integration
