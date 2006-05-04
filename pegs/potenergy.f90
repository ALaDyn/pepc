


!  ===================================================================
!
!                              POTENERGY
!
!   Calculate potential energy:
!     nterm(p) = terms in particle interaction list
!   Pseudoparticles are given by: inlis(j,p), j=1,nterm(p)
!
!  ===================================================================

subroutine potenergy(epot_dust,epot_star)
  use physvars
  use treevars
  use utils
  implicit none
  include 'mpif.h'


  integer :: key2addr        ! Mapping function to get hash table address from key

  integer :: p, i, j, ierr 

  real :: upartial, epot_dust, epot_total, epot_star
  logical :: pot_debug=.false.

  epot_total = 0.  ! Global potential energy

  !  Sum Potential Energy

  upartial = 0. ! Single PE partial potential energy sum
                ! - pot(p) includes contrib from stars
  do p=1, npp
     upartial = upartial + 0.5*m(p)*pot(p)

     if (pot_debug) then
        write (ipefile,'(a,i5,a,i5,3f10.3,a,f12.4)') 'local particle ',p,' label ',pelabel(p),x(p),y(p),q(p),' pot ',pot(p)
     endif
  end do

  if (pot_debug) write (ipefile,'(a,1pe11.4)') 'partial PE sum',upartial


  call MPI_ALLREDUCE(upartial, epot_dust, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

  epot_star = 0. 
  do i=1, ni
     epot_star = epot_star + 0.5*m_star(i)*pot_star(i)
  end do


  epot_total = epot_dust + epot_star



end subroutine potenergy
