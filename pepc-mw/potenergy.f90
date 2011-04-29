


!  ===================================================================
!
!                              POTENERGY
!
!   Calculate E.S. and magnetic potential energies:
!
!  ===================================================================

subroutine potenergy(epot_total)
  use physvars
  use module_fmm_framework
  implicit none
  include 'mpif.h'

  integer :: p, ierr

  real*8 :: upartial
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


  call MPI_ALLREDUCE(upartial, epot_total,1, MPI_REAL8, MPI_SUM, MPI_COMM_PEPC, ierr)

  ! this can also be done in fields.f90, but thematically it fits better here :-)
  potfarfield  = potfarfield/2.
  potnearfield = potnearfield/2.

  call MPI_ALLREDUCE(MPI_IN_PLACE, potfarfield,  1, MPI_REAL8, MPI_SUM, MPI_COMM_PEPC, ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE, potnearfield, 1, MPI_REAL8, MPI_SUM, MPI_COMM_PEPC, ierr)



end subroutine potenergy
