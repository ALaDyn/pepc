! ==============================================
!
!                SCRAMBLE_V
!
!  Sets up physical system: particle positions, velocities
!
! ==============================================

subroutine scramble_v(i1,n)

  use treevars
  use utils
!VAMPINST include
      INCLUDE 'VTcommon.h'
      INTEGER VTIERR
!
  implicit none

 integer :: dum1, dum2, dum3  
real :: uxt, uyt, uzt
  integer :: i, j, k, kk, p, i1, n, n1

!VAMPINST subroutine_start
       CALL VTENTER(IF_scramble_v,VTNOSCL,VTIERR)
!      write(*,*) 'VT: scramble_v S>',VTIERR,
!     *    IF_scramble_v,ICLASSH
!
  dum1 = -71 - 10*me
  dum2 = -113301 - 10*me
  dum3 = -8651 - 10*me
  !  exclude odd one out
  if (mod(n,2).ne.0) then
     n1=n-1
  else
    n1 = n
  endif
 write(15,*) 'scrambling ',n1,' particles'
  !  scramble indices to remove correlation between ux,uy,uz
  do i=1,n1
     p=i+i1-1
     j=n1*rano(dum1)+i1
     k=n1*rano(dum2)+i1
     kk=n1*rano(dum3)+i1
     uxt=ux(p)
     uyt=uy(p)
     uzt=uz(p)
     ux(p)=ux(kk)
     uy(p)=uy(j)
     uz(p)=uz(k)
     ux(kk)=uxt
     uy(j)=uyt
     uz(k)=uzt
  end do

!VAMPINST subroutine_end
       CALL VTLEAVE(ICLASSH,VTIERR)
!      write(*,*) 'VT: scramble_v S<',VTIERR,ICLASSH
!
end subroutine scramble_v
