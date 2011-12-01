!  This subroutine de-allocates the array space for the physvars module.
!  @param my_rank_l rank of the cpu
!  @param n_cpu_l number of cpu`s
subroutine cleanup(my_rank_l,n_cpu_l)
  
  use physvars
  implicit none
  include 'mpif.h'

  integer, intent(in) :: my_rank_l ! MPI cpu rank
  integer, intent(in) :: n_cpu_l  ! MPI # CPUs

! copy call parameters to physvars module

  my_rank = my_rank_l
  n_cpu = n_cpu_l
 
  ! particle array deallocation in physvars

  deallocate ( particles )

end subroutine cleanup
