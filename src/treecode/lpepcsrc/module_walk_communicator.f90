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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>  Encapsulates helper functions that simplify communication during tree traversal
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_walk_communicator
  use module_pepc_types, only: t_tree_node
  implicit none

  private

    integer, public :: walk_status

    ! internal status
    integer, public, parameter :: WALK_STILL_RUNNING = 0 !< value of walk_status: there are still particles left to be processed
    integer, public, parameter :: WALK_IAM_FINISHED  = 1 !< value of walk_status: all particles on this PE have finished their walk
    integer, public, parameter :: WALK_ALL_MSG_DONE  = 2 !< value of walk_status: there are no more pending requests or answers in the message queue
    integer, public, parameter :: WALK_I_NOTIFIED_0  = 3 !< value of walk_status: rank 0 has been informed that this PE is finished
    integer, public, parameter :: WALK_ALL_FINISHED  = 4 !< value of walk_status: rank 0 told this pe, that every PE is finished with walking, i.e. no more communication is necessary

    public init_comm_data
    public uninit_comm_data
    public notify_walk_finished
    public send_request
    public send_data
    public unpack_data
    public send_walk_finished
    public broadcast_walk_finished

  contains

end module module_walk_communicator
