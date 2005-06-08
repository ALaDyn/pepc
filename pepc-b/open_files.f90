subroutine openfiles

  use physvars
  character(30) :: cfile
  character(1) :: csnap
  character(3) :: cme

  if (my_rank == 0) then
     !  master diagnostics output
     open(15,file='run.out')
     open(81,file='parts_all.dat')

     open(70,file='domains.dat')
     open(75,file='energy.dat')      ! energies
  endif

  !  stdout for PE my_rank
  !  must first create subdirectory 'peXXX' in run directory

  cme = achar(my_rank/100+48) // achar(mod(my_rank/10,10)+48) // achar(mod(my_rank,10)+48)  ! Convert 3-digit PE number into character string
  cfile="pe"//cme//"/dump."//cme
  open(20,file=cfile)
  ifile_cpu = 20

end subroutine openfiles
