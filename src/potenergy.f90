


!  ===================================================================
!
!                              POTENERGY
!
!   Calculate potential energy:
!     nterm(p) = terms in particle interaction list
!   Pseudoparticles are given by: inlis(j,p), j=1,nterm(p)
!
!  ===================================================================

subroutine potenergy(epot_total)
  use treevars
  use utils
  implicit none
  integer :: key2addr        ! Mapping function to get hash table address from key

  integer :: p, i, j 

  real :: upartial, epot_total
  logical :: pot_debug=.false.

  epot_total = 0.  ! Global potential energy

  !  Sum Potential Energy

  upartial = 0. ! Single PE partial potential energy sum
  do p=1, npp
     upartial = upartial + 0.5*q(p)*pot(p)

     if (pot_debug) then
        write (ipefile,'(a,i5,a,i5,3f10.3,a,f12.4)') 'local particle ',p,' label ',pelabel(p),x(p),y(p),q(p),' pot ',pot(p)
     endif
  end do

  if (pot_debug) write (ipefile,'(a,1pe11.4)') 'partial PE sum',upartial


  call MPI_ALLREDUCE(upartial, epot_total, one, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)


end subroutine potenergy
