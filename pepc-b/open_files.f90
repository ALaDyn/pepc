subroutine openfiles

  use physvars
  implicit none
  character(30) :: cfile
  character(1) :: csnap

  if (my_rank == 0) then
     !  master diagnostics output
     open(15,file='tree.out')  ! Tree stats
     open(24,file='pepcb.out') ! Physics log

     open(70,file='domains.dat')
     open(71,file='laser.dat')   ! laser parameters
     open(75,file='energy.dat')      ! energies
     write(*,*) 'debug level: ',debug_level,' idump',idump
  endif

  !  stdout for PE my_rank
  !  must first create subdirectory 'data' in run directory

  csubme =   achar(my_rank/1000+48) &
       // achar(mod(my_rank/100,10)+48) &
       // achar(mod(my_rank/10,10)+48) &
       // achar(mod(my_rank,10)+48)  ! Convert 4-digit PE number into character string
  cfile="data/out."//csubme
  if (debug_level>2 .and. idump>0) then
      open(90,file='openfiles.out')
      write (90,'(a3,i6,a15,a30)') 'PE ',my_rank,' opening ',cfile
      open(20,file=cfile) ! suppress I/O for debug=0
  endif
  ifile_cpu = 20

end subroutine openfiles
