


!  ===================================================================
!
!                              POTENERGY
!
!   Calculate electrostatic and magnetic potential energies:
!
!  ===================================================================

subroutine potenergy(epot_total,emag_total)
  use physvars
  use treevars
  implicit none
  include 'mpif.h'


  integer :: p, ierr 

  real*8 :: upartial, umag, gamma
  real*8, intent(out) :: epot_total, emag_total
  logical :: pot_debug=.false.

  epot_total = 0.  ! Global potential energy
  emag_total = 0.  ! Global magnetic energy

  !  Sum Potential Energy

  upartial = 0. ! Single PE partial potential energy sum
  umag = 0. ! Single PE partial potential energy sum

  do p=1, npp
     upartial = upartial + 0.5*q(p)*pot(p)
     if (bfields) then
        gamma = sqrt(1.0+ux(p)**2+uy(p)**2+uz(p)**2)
        umag = umag + 0.5*q(p)/gamma*( ux(p)*ax(p) + uy(p)*ay(p) + uz(p)*az(p) )
     endif
     if (pot_debug) then
        write (ipefile,'(a,i5,a,i5,3f10.3,a,f12.4,a,3f12.5)') & 
	'local particle ',p,' label ',pelabel(p),x(p),y(p),q(p),' pot ',pot(p),' mag ',ax(p), ay(p), az(p)
     endif
  end do

  if (pot_debug) write (ipefile,'(a,1pe11.4)') 'partial PE sum',upartial


  call MPI_ALLREDUCE(upartial, epot_total,1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(umag, emag_total, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)


end subroutine potenergy
