 
!  ===============================================================
!
!                           RESET_IONS
!
!   Reset ion mass and velocity after const-temp eqm phase
!
!  ===============================================================

subroutine reset_ions

  use treevars
!VAMPINST include
      INCLUDE 'VTcommon.h'
      INTEGER VTIERR
!

  integer :: i,p
  real :: ratio_clamp  ! mass ratio used for NVT dynamics

!VAMPINST subroutine_start
       CALL VTENTER(IF_reset_ions,VTNOSCL,VTIERR)
!      write(*,*) 'VT: reset_ions S>',VTIERR,
!     *    IF_reset_ions,ICLASSH
!
  ratio_clamp = mass_i/mass_e  ! (typically 10)
  do p=1,npp

     if (q(p) > 0) then   ! only constrain ions to target boundaries
        m(p) = mass_ratio*mass_e  ! mass_ratio preserved from initial inputs
        ux(p) = ux(p)/sqrt(mass_ratio/ratio_clamp)  ! scale velocities back so that T_i = mv^2 conserved
        uy(p) = uy(p)/sqrt(mass_ratio/ratio_clamp) 
        uz(p) = uz(p)/sqrt(mass_ratio/ratio_clamp) 
     endif
  end do
!VAMPINST subroutine_end
       CALL VTLEAVE(ICLASSH,VTIERR)
!      write(*,*) 'VT: reset_ions S<',VTIERR,ICLASSH
!
end subroutine reset_ions
