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

module files
   implicit none
   private

   public openfiles
   public closefiles

contains

   subroutine openfiles
      use physvars

      if (my_rank .eq. 0) then
         !  master diagnostics output
         open (15, file='run.out')
         open (70, file='domains.dat')
      end if

   end subroutine openfiles

   subroutine closefiles
      use physvars

      if (my_rank .eq. 0) then
         close (15)
         close (70)
      end if

      close (80)  ! initial particle data

   end subroutine closefiles

end module files
