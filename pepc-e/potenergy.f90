


!  ===================================================================
!
!                              POTENERGY
!
!   Calculate E.S. and magnetic potential energies:
!
!  ===================================================================

subroutine potenergy(epot_total)
  use physvars
  implicit none
  include 'mpif.h'

  integer :: key2addr        ! Mapping function to get hash table address from key

  integer :: p, i, j, ierr 

  real*8 :: upartial, umag, gamma
  real*8, intent(out) :: epot_total
  logical :: pot_debug=.false.
!  logical :: pot_debug=.true.

  epot_total = 0.  ! Global potential energy

  !  Sum Potential Energy

  upartial = 0. ! Single PE partial potential energy sum

  do p=1, np_local
     upartial = upartial + 0.5*q(p)*pot(p)

     if (my_rank == 0 .and. pot_debug) then
        write (ifile_cpu,'(a,i5,a,i5,3f10.3,a,f12.4)') & 
	'local particle ',p,' label ',pelabel(p),x(p),y(p),q(p),' pot ',pot(p)
     endif
  end do

  if (pot_debug) write (ifile_cpu,'(a,1pe11.4,i2)') 'partial PE sum',upartial,my_rank


  call MPI_ALLREDUCE(upartial, epot_total,1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)


end subroutine potenergy
