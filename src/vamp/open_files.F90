subroutine openfiles

  use treevars
!VAMPINST include
      INCLUDE 'VTcommon.h'
      INTEGER VTIERR
!
  character(30) :: cfile
  character(1) :: csnap
  character(3) :: cme

!VAMPINST subroutine_start
       CALL VTENTER(IF_openfiles,VTNOSCL,VTIERR)
!      write(*,*) 'VT: openfiles S>',VTIERR,
!     *    IF_openfiles,ICLASSH
!
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

!VAMPINST subroutine_end
       CALL VTLEAVE(ICLASSH,VTIERR)
!      write(*,*) 'VT: openfiles S<',VTIERR,ICLASSH
!
end subroutine openfiles
