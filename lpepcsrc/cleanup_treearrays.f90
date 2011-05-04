subroutine pepc_cleanup(my_rank,n_cpu)
  use treevars
  implicit none
  include 'mpif.h'

  integer, intent(in) :: my_rank  ! MPI cpu rank
  integer, intent(in) :: n_cpu  ! MPI # CPUs

! copy call parameters to treevars module

  me = my_rank
  num_pe = n_cpu
 
  if (me==0) then
     write(*,'(a)') 'LPEPC | De-allocating particle and tree arrays ...'
  endif

  ! particle array deallocation

  deallocate ( x, y, z, ux, uy, uz, & 
       q, m, work, &
       Ex, Ey, Ez, pot, &
       Ax, Ay, Az, &
       Bx, By, Bz,  &
       Axo, Ayo, Azo, &
       pepid, pelabel, pekey )    

  deallocate ( nbranches, igap )

end subroutine pepc_cleanup






