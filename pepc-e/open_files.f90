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
  endif


end subroutine openfiles
