subroutine openfiles

  use physvars
  implicit none
  character(30) :: cfile
  character(1) :: csnap

  if (my_rank == 0) then
     !  master diagnostics output
     open(15,file='run.out')
     open(81,file='parts_all.dat')

     open(70,file='domains.dat')
     open(75,file='energy.dat')      ! energies
  endif

  !  stdout for PE my_rank
  !  must first create subdirectory 'data/peXXXX' in run directory

  csubme =   achar(my_rank/1000+48) &
       // achar(mod(my_rank/100,10)+48) &
       // achar(mod(my_rank/10,10)+48) &
       // achar(mod(my_rank,10)+48)  ! Convert 4-digit PE number into character string
  cfile="data/pe"//csubme//"/dump."//csubme
!  write (*,'(a3,i6,a15,a30)') 'PE ',my_rank,' opening ',cfile
  if (debug_level>1) open(20,file=cfile) ! suppress I/O for debug=0
  ifile_cpu = 20

end subroutine openfiles
