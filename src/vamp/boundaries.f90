 
!  ===============================================================
!
!                           PMOVE
!
!   Update particle positions - used with leap-frog scheme
!
!  ===============================================================

subroutine push(ips,ipf,delt)

  use treevars
!VAMPINST include
      INCLUDE 'VTcommon.h'
      INTEGER VTIERR
!
  integer, intent(in) :: ips, ipf  ! 1st and last particle numbers
  real, intent(in) :: delt
  integer :: i,p
  real :: gamma

  !  relativistic particle push in space

!VAMPINST subroutine_start
       CALL VTENTER(IF_push,VTNOSCL,VTIERR)
!      write(*,*) 'VT: push S>',VTIERR,
!     *    IF_push,ICLASSH
!
  do p=ips,ipf
     gamma = sqrt(1.0 + ux(p)**2 + uy(p)**2 + uz(p)**2)
!     gamma = 1.0
     x(p)=x(p)+ux(p)/gamma*delt
     y(p)=y(p)+uy(p)/gamma*delt
     z(p)=z(p)+uz(p)/gamma*delt

  end do

!VAMPINST subroutine_end
       CALL VTLEAVE(ICLASSH,VTIERR)
!      write(*,*) 'VT: push S<',VTIERR,ICLASSH
!
end subroutine push
