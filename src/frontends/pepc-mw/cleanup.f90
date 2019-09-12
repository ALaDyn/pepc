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

!  This subroutine de-allocates the array space for the physvars module.
!  @param my_rank_l rank of the cpu
!  @param n_cpu_l number of cpu`s
subroutine cleanup(my_rank_l,n_cpu_l)
  
  use physvars
  implicit none

  integer, intent(in) :: my_rank_l ! MPI cpu rank
  integer, intent(in) :: n_cpu_l  ! MPI # CPUs

! copy call parameters to physvars module

  my_rank = my_rank_l
  n_cpu = n_cpu_l
 
  ! particle array deallocation in physvars

  deallocate ( particles, energy)

end subroutine cleanup
