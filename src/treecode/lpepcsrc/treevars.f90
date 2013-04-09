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
!>  Encapsulates all global variables for lpepc
!>
module treevars
  implicit none

  ! I/O units
  integer, parameter :: stats_u = 60 !< statistics

  !  Associated parallelization stuff
  integer :: me       !< Rank of current task
  integer :: num_pe   !< # cpus used by program
  integer :: MPI_COMM_lpepc !< communicator that has been supplied to or created by pepc_initialize
  integer :: num_threads = 3 !< number of threads to be used for hybrid parallelization (Pthreads, OpenMP, etc.), for compatibility, we set it to num_walk_threads in tree_walk_read_parameters() for now

  integer :: nlev, & !< max refinement level
             idim !< dimension of the system

! Memory control
  real    :: np_mult = 1.5
  integer :: interaction_list_length_factor = 1 !< factor for increasing todo_list_length and defer_list_length in case of respective warning (e.g. for very inhomogeneous or 2D cases set to 2..8)

  contains

  subroutine treevars_prepare(dim)
    implicit none
    integer, intent(in) :: dim

    idim = dim
    nlev = 60 / idim
  end subroutine treevars_prepare


  subroutine treevars_finalize()
    implicit none
  end subroutine treevars_finalize
end module treevars



