subroutine openfiles

  use treevars
  character(30) :: cfile
  character(1) :: csnap
  character(3) :: cme

  if (me == 0) then
     !  master diagnostics output
     open(15,file='run.out')
     open(81,file='parts_all.dat')

     open(70,file='domains.dat')
     open(75,file='energy.dat')      ! energies
  endif

  !  stdout for PE me
  !  must first create subdirectory 'peXXX' in run directory

  cme = achar(me/100+48) // achar(mod(me/10,10)+48) // achar(mod(me,10)+48)  ! Convert 3-digit PE number into character string
  cfile="pe"//cme//"/dump.out"
  open(20,file=cfile)
  ipefile = 20

end subroutine openfiles
