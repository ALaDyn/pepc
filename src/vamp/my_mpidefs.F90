 ! Machine-specific MPI stuff
 ! Cray uses 8-byte integers throughout
 ! MPI standard specifies 4-byte, eg: IBM

module my_mpidefs

  implicit none
  include 'mpif.h'
  integer :: me, &      ! Rank of current task (PE)
       num_pe, &  ! # PEs used by program
       lastpe, &  ! Rank of last PE
       destpe, &  ! Rank of destination
       me_minus_one, & ! me-1
       me_plus_one, &  ! me+1
       ierr, &    ! MPI error return code
       nbuf, &   ! Buffer length
       status(MPI_STATUS_SIZE), &
       index, &  ! Index of completed comm
       tag1 = 1, &  ! Tags
       tag2 = 2

  integer, parameter :: one = 1, &   ! MPI call constants
       root = 0

  integer, allocatable :: recv_counts(:), & ! arrays for gather operations
       recv_strides(:), &
       send_counts(:), &  ! arrays for scatter
       send_strides(:), &
       stat_pe(:,:), & ! status
       pe_handle(:), &  ! Handles for non-blocking comm  2*num_pe
       send_key_handle(:), &  !(num_pe)
       recv_key_handle(:), &
       send_child_handle(:), &  !(6*nppm)
       recv_child_handle(:)


end module my_mpidefs

