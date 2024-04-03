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

!!!!!!!!!!!!!!!!!!!!
!! helper module
!!!!!!!!!!!!!!!!!!!!

module module_diagnostic

!  use module_pepc_kinds
!  use module_pepc_types
!  use module_interaction_Specific_types
   use module_globals    !,only: root
   use module_shortcut, only: zero, one, half
   implicit none

contains

   subroutine hamiltonian(np, p, pold, t)
      use module_globals, only: lorentz_tilde, nsp, ns_tot, uth, vth, tnp
      use mpi
      implicit none
      type(t_particle), allocatable, intent(in) :: p(:), pold(:)
      real(kind_particle), intent(in) :: t
      integer(kind_particle), intent(in) :: np
      integer(kind_particle)                    :: ip, jp, species
      integer                                   :: rc, rd
      real(kind_particle)                       :: upot(1:nsp), ukin(1:nsp), udar(1:nsp), gam, v2, ploc(1:3), &
                                                   uk_mean(1:nsp), uk_std(1:nsp), count_species(1:nsp), rtnp, etot, pmax(1:3), divA!,vthermal(1:3),g

      upot = zero
      ukin = zero
      udar = zero
      ploc = zero
      pmax = zero
!      vthermal = zero
      uk_mean = zero
      uk_std = zero
      etot = zero
      divA = zero
      count_species = zero

      rtnp = real(tnp, kind=kind_particle)

      do ip = 1, np

         species = p(ip)%label
!          write(*,*) ip,species
         count_species(species) = count_species(species) + 1

                !!!!   ENERGY !!!!!!!!!!

         v2 = dot_product(p(ip)%data%v, p(ip)%data%v)
         gam = v2 / (p(ip)%data%g + one)

         upot(species) = upot(species) + half * p(ip)%data%q * p(ip)%results%pot
         udar(species) = udar(species) + half / lorentz_tilde * p(ip)%data%q * dot_product(p(ip)%results%A(1:3), p(ip)%data%v(1:3)) / p(ip)%data%g

         ukin(species) = ukin(species) + p(ip)%data%m * gam
         uk_std(species) = uk_std(species) + (p(ip)%data%m * gam)**2
!                !!!!!! CANONICAL MOMENTUM - Max Norm !!!!!!
!
         pmax(1) = max(abs(p(ip)%data%m * (p(ip)%data%v(1) - pold(ip)%data%v(1)) + &
                           p(ip)%data%q / lorentz_tilde * (p(ip)%results%A(1) - pold(ip)%results%A(1))), pmax(1))

         pmax(2) = max(abs(p(ip)%data%m * (p(ip)%data%v(2) - pold(ip)%data%v(2)) + &
                           p(ip)%data%q / lorentz_tilde * (p(ip)%results%A(2) - pold(ip)%results%A(2))), pmax(2))

         pmax(3) = max(abs(p(ip)%data%m * (p(ip)%data%v(3) - pold(ip)%data%v(3)) + &
                           p(ip)%data%q / lorentz_tilde * (p(ip)%results%A(3) - pold(ip)%results%A(3))), pmax(3))

         divA = max(divA, abs(p(ip)%results%dxA(1) + p(ip)%results%dyA(2) + p(ip)%results%dzA(3)))

!          g            = one/sqrt( one - ( uth(1)**2 + vth(1)**2 )**2 )
!          vthermal(1)  = vthermal(1)   + p(ip)%data%m*g*uth(1)
!          vthermal(2)  = vthermal(2)   + p(ip)%data%m*g*vth(1)
!          ploc(1:2)    = ploc(1:2)     + p(ip)%data%m*p(ip)%data%v(1:2)

         ploc(1:3) = ploc(1:3) + p(ip)%data%m * p(ip)%data%v(1:3) + p(ip)%data%q / lorentz_tilde * p(ip)%results%A(1:3)

      end do

      call MPI_ALLREDUCE(MPI_IN_PLACE, count_species, nsp, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
      call MPI_ALLREDUCE(MPI_IN_PLACE, upot(1:nsp), nsp, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
      call MPI_ALLREDUCE(MPI_IN_PLACE, ukin(1:nsp), nsp, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
      call MPI_ALLREDUCE(MPI_IN_PLACE, udar(1:nsp), nsp, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
      call MPI_ALLREDUCE(MPI_IN_PLACE, uk_std(1:nsp), nsp, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
      call MPI_ALLREDUCE(MPI_IN_PLACE, ploc(1:3), 3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
!      call MPI_ALLREDUCE(MPI_IN_PLACE , ploc(1:2)     , 2  , MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
!      call MPI_ALLREDUCE(MPI_IN_PLACE , vthermal(1:2) , 2  , MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
      call MPI_ALLREDUCE(MPI_IN_PLACE, pmax(1), 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, rc)
      call MPI_ALLREDUCE(MPI_IN_PLACE, pmax(2), 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, rc)
      call MPI_ALLREDUCE(MPI_IN_PLACE, pmax(3), 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, rc)
      call MPI_ALLREDUCE(MPI_IN_PLACE, divA, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, rc)

!      if ( vthermal(1) .gt. zero )  ploc(1)  = ploc(1)/vthermal(1)
!      if ( vthermal(2) .gt. zero )  ploc(2)  = ploc(2)/vthermal(2)

      do ip = 1, nsp

         if (count_species(ip) .gt. zero) then
            uk_mean(ip) = sqrt(uk_std(ip) / count_species(ip))
            uk_std(ip) = count_species(ip) / (count_species(ip) - one) * (uk_std(ip) / count_species(ip) - uk_mean(ip)**2)
         end if

      end do

      if (root) then
         etot = sum(upot(1:nsp))
         etot = etot + sum(udar(1:nsp))
         etot = etot + sum(ukin(1:nsp))
         open (newunit=rc, file=trim(folder)//trim("energy_")//trim(adjustl(ischeme))//".dat", form='formatted', status='unknown', position='append')
         open (newunit=rd, file=trim(folder)//trim("momentum_")//trim(adjustl(ischeme))//".dat", form='formatted', status='unknown', position='append')
         write (rc, *) t, ukin(1:nsp), uk_mean(1:nsp), uk_std(1:nsp), upot(1:nsp), udar(1:nsp), etot, divA
         write (rd, *) t, ploc, pmax
         close (rc)
         close (rd)
      end if

   end subroutine hamiltonian

   subroutine densities(np, p, t)
      use module_globals, only: nsp, ns_tot, percentages
      use mpi
      implicit none
      type(t_particle), allocatable, intent(in) :: p(:)
      real(kind_particle), intent(in) :: t
      integer(kind_particle), intent(in) :: np
      integer(kind_particle)                    :: ip, species, i, j
      integer                                   :: rc
      real(kind_particle)                       :: E2(1:nsp, 1:3), B2(1:nsp, 1:3), E2_rms(1:nsp, 1:3), B2_rms(1:nsp, 1:3), count_species(1:nsp)

      E2 = zero
      B2 = zero
      E2_rms = zero
      B2_rms = zero
      count_species = zero

      do ip = 1, np

         i = p(ip)%label
         count_species(i) = count_species(i) + 1

         do j = 1, 3
            E2(i,j)     = E2(i,j)     + p(ip)%results%E(j)**2                 !&
            B2(i,j)     = B2(i,j)     + p(ip)%results%B(j)**2                 !&
            E2_rms(i,j) = E2_rms(i,j) + p(ip)%results%E(j)**4                 !&
            B2_rms(i,j) = B2_rms(i,j) + p(ip)%results%B(j)**4                 !&
         end do

      end do

      call MPI_ALLREDUCE(MPI_IN_PLACE, count_species  , nsp, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc) !&
      call MPI_ALLREDUCE(MPI_IN_PLACE, E2(1:nsp,1)    , nsp, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc) !&
      call MPI_ALLREDUCE(MPI_IN_PLACE, E2(1:nsp,2)    , nsp, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc) !&
      call MPI_ALLREDUCE(MPI_IN_PLACE, E2(1:nsp,3)    , nsp, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc) !&
      call MPI_ALLREDUCE(MPI_IN_PLACE, B2(1:nsp,1)    , nsp, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc) !&
      call MPI_ALLREDUCE(MPI_IN_PLACE, B2(1:nsp,2)    , nsp, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc) !&
      call MPI_ALLREDUCE(MPI_IN_PLACE, B2(1:nsp,3)    , nsp, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc) !&
      call MPI_ALLREDUCE(MPI_IN_PLACE, E2_rms(1:nsp,1), nsp, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc) !&
      call MPI_ALLREDUCE(MPI_IN_PLACE, E2_rms(1:nsp,2), nsp, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc) !&
      call MPI_ALLREDUCE(MPI_IN_PLACE, E2_rms(1:nsp,3), nsp, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc) !&
      call MPI_ALLREDUCE(MPI_IN_PLACE, B2_rms(1:nsp,1), nsp, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc) !&
      call MPI_ALLREDUCE(MPI_IN_PLACE, B2_rms(1:nsp,2), nsp, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc) !&
      call MPI_ALLREDUCE(MPI_IN_PLACE, B2_rms(1:nsp,3), nsp, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc) !&

      do i = 1, nsp

         if (count_species(i) .gt. zero) then
            E2(i, :) = E2(i, :) / count_species(i)
            B2(i, :) = B2(i, :) / count_species(i)
            E2_rms(i, :) = sqrt(E2_rms(i, :) / count_species(i))
            B2_rms(i, :) = sqrt(B2_rms(i, :) / count_species(i))
         end if

      end do

      if (root) then
         open (newunit=rc, file=trim(folder)//trim("density_")//trim(adjustl(ischeme))//".dat", form='formatted', status='unknown', position='append')
         write (rc, *) t, E2(1:nsp, 1), E2(1:nsp, 2), E2(1:nsp, 3), B2(1:nsp, 1), B2(1:nsp, 2), B2(1:nsp, 3), &
            E2_rms(1:nsp, 1), E2_rms(1:nsp, 2), E2_rms(1:nsp, 3), B2_rms(1:nsp, 1), B2_rms(1:nsp, 2), B2_rms(1:nsp, 3)
         close (rc)

      end if

   end subroutine densities

   subroutine beam_rnv(tnp, p, t)
      use mpi
      implicit none
      type(t_particle), allocatable, intent(in) :: p(:)
      real(kind_particle), intent(in) :: t
      integer(kind_particle), intent(in) :: tnp
      integer(kind_particle)                    :: ip, jp, np
      integer                                   :: rc
      real(kind_particle)                       :: rloc, rglo, vloc(1:3), vglo(1:3), vrloc, vrglo, rtnp

      rloc = zero
      vloc(1:3) = zero
      vrloc = zero
      rglo = zero
      vglo(1:3) = zero
      vrglo = zero

      np = size(p, kind=kind_particle)
      rtnp = real(tnp, kind=kind_particle)

      do ip = 1, np
         rloc = rloc + dot_product(p(ip)%x(1:2), p(ip)%x(1:2))
         vrloc = vrloc + dot_product(p(ip)%data%v(1:2), p(ip)%x(1:2)) / p(ip)%data%g
         do jp = 1, 3
            vloc(jp) = vloc(jp) + p(ip)%data%v(jp)**2 / p(ip)%data%g**2
         end do
      end do

      call MPI_ALLREDUCE(rloc , rglo , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc) !&
      call MPI_ALLREDUCE(vloc , vglo , 3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc) !&
      call MPI_ALLREDUCE(vrloc, vrglo, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc) !&

      rglo = sqrt(rglo / rtnp)
      vglo(1:3) = sqrt(vglo(1:3) / rtnp)
      vrglo = vrglo / rtnp / rglo

      if (root) then
         open (newunit=rc, file=trim(folder)//trim("beam_rnv_")//trim(adjustl(ischeme))//".dat", form='formatted', status='unknown', position='append')
         write (rc, *) t, rglo, vglo(1:3), vrglo
         close (rc)
      end if

   end subroutine

end module
