
! ==============================================================
!
!                    KEYTEST
!
!   Builds binary keys from random particle data
!
!  ==============================================================

program keytest

  use treevars
!VAMPINST include
      INCLUDE 'VTcommon.h'
      INTEGER VTIERR
!


!VAMPINST program_start
      CALL VTdefs()
      CALL VTENTER(IF_keytest,VTNOSCL,VTIERR)
      write(*,*) 'VT: keytest P>',VTIERR, &
         IF_keytest,ICLASSH
!
  call randion
  call makekeys

  call draw2d
  call maketable

!VAMPINST program_end
      CALL VTLEAVE(ICLASSH,VTIERR)
      write(*,*) 'VT: keytest P<',VTIERR,ICLASSH
!
end program
