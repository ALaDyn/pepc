subroutine closefiles
  use treevars
!VAMPINST include
      INCLUDE 'VTcommon.h'
      INTEGER VTIERR
!

!VAMPINST subroutine_start
       CALL VTENTER(IF_closefiles,VTNOSCL,VTIERR)
!      write(*,*) 'VT: closefiles S>',VTIERR,
!     *    IF_closefiles,ICLASSH
!
  if (me == 0) then
     close(15)
     close(81)  ! particle dump 
     close(70)
     close(75)
  endif
  close(20)
  close(80)  ! initial particle data


!VAMPINST subroutine_end
       CALL VTLEAVE(ICLASSH,VTIERR)
!      write(*,*) 'VT: closefiles S<',VTIERR,ICLASSH
!
end subroutine closefiles
