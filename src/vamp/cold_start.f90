! ==============================================
!
!                COLD_START
!
!  Initialise set of particles with zero velocity
!
! ==============================================

subroutine cold_start(i1,n)

  use treevars
  use utils
!VAMPINST include
      INCLUDE 'VTcommon.h'
      INTEGER VTIERR
!
  implicit none


  integer :: i1,n

!VAMPINST subroutine_start
       CALL VTENTER(IF_cold_start,VTNOSCL,VTIERR)
!      write(*,*) 'VT: cold_start S>',VTIERR,
!     *    IF_cold_start,ICLASSH
!
  ux(i1:i1+n) = 0.
  uy(i1:i1+n) = 0.
  uz(i1:i1+n) = 0.


!VAMPINST subroutine_end
       CALL VTLEAVE(ICLASSH,VTIERR)
!      write(*,*) 'VT: cold_start S<',VTIERR,ICLASSH
!
end subroutine cold_start
