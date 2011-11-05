!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>  Encapsulates anything concerning the particle velocity update
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_integration_scheme
 
      implicit none
      save
      private

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer, public :: integrator_scheme = 1
      real*8, public :: Te0 = 0., Te_uncor = 0., chie = 0., delta_Te = 0.
      real*8, public :: Ti0 = 0., Ti_uncor = 0., chii = 0., delta_Ti = 0.

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      public push_TE	!< 2v Boris rotation for TE fields (Ex, Ey, Bz)
      public push_TM	!< 2v Boris roation for TM fields (Ey, Bx, By)
      public push_x	!< Position update with relativistic velocities
      public push_restrict	!< Position update with restricted delta-x
      public push_nonrel !< Non-rel. position update
      public push_full3v !< 3V integrator with all field components
      public velocities  !< MD velocity update with constraints

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      contains




!  ==============================================
!
!     2v particle pusher
!
!  TE fields only (p-pol):  (Ex,Ey,0)  (0,0,Bz)
!  nonrelativistic for now
!  ==============================================


subroutine push_TE(p_start,p_finish,dts)
  use module_physvars
  use module_particle_props
!  use module_laser

  implicit none
  integer, intent(in) :: p_start, p_finish
  real, intent(in) :: dts
  integer :: p
  real*8 :: beta, gam1,tt, ss
  real*8 :: uxd,uyp,uxp,uxm,uym, uzm
  real*8 :: exi, eyi, ezi, bxi, byi, bzi



  do p = p_start, p_finish
     beta=q(p)/m(p)*dts*0.5  ! charge/mass constant

     !  Sum internal and external fields
     exi = ex(p)
     eyi = ey(p)
     ezi = 0.
     bxi = 0.
     byi = 0.
     bzi = bz(p)

     !   half-accn

     uxm = ux(p) + beta*exi
     uym = uy(p) + beta*eyi
     uzm = uz(p)

     !   rotation

!     gam1=sqrt(1.0+uxm**2+uym**2+uzm**2)
     gam1=1.d0
     tt=beta*bzi/gam1
     ss=2.d0*tt/(1.d0+tt**2)

     uxd = uxm + uym*tt
     uyp = uym - uxd*ss
     uxp = uxd + uyp*tt

     !   half-accn

     ux(p)=uxp+beta*exi
     uy(p)=uyp+beta*eyi

  end do

end subroutine push_TE


!  ==============================================
!
!     2v particle pusher
!
!  TM fields only (s-pol):  (0,0,Ez)  (Bx,By,0)
!
!  ==============================================


subroutine push_TM(p_start,p_finish,dts)
  use module_physvars
  use module_particle_props
  use module_laser

  implicit none
  integer, intent(in) :: p_start, p_finish
  real, intent(in) :: dts
  integer :: p
  real :: beta, gam1,tt, sy, sz, tz, ty
  real :: uxd,uyp,uzp,uxp,uxm,uym,uzm
  real :: exi, eyi, ezi, bxi, byi, bzi
  real :: phipon, ez_em, by_em, bx_em, az_em, xd, yd, zd


  if (itime>0 .and. beam_config==4) focus(1) = x_crit  ! laser tracks n_c

  do p = p_start, p_finish
     beta=q(p)/m(p)*dts*0.5  ! charge/mass constant

     xd = x(p)-focus(1)
     yd = y(p)-focus(2)
     zd = z(p)-focus(3)

     ! evaluate external fields at particle positions

     if (beam_config.eq.4) then
! pond. standing wave on step-profile
        call empond(tlaser,tpulse,sigma,vosc,omega,xd,yd,zd,ez_em,by_em,bx_em,az_em,phipon)

     else if (beam_config.eq.6) then
! plane wave with Gaussian spot
        call emplane(tlaser,tpulse,sigma,vosc,omega,xd,yd,zd,ez_em,by_em,bx_em,az_em,phipon)
     endif

     !  Sum internal and external fields
     exi = ex(p)
     eyi = ey(p)
     ezi = ez(p)+ez_em
     bxi = 0.
     byi = by_em
     bzi = 0.

     ! transverse momentum from pz=az including thermal motion 
     !	uzi = -q(p)/m(p)*az_em + uz(p)


     !   first half-accn
     uxm = ux(p) + beta*exi
     uym = uy(p) + beta*eyi
     uzm = uz(p) + beta*ezi
     !     uzm = uz(p) - q(p)/m(p)*az_em
     !   rotation
     gam1=dsqrt(1.d0 + uxm**2 + uym**2 + uzm**2)
     ty = beta*byi/gam1
     tz = beta*bzi/gam1
     tt = 1.0 + ty**2+tz**2
     sy = 2.0*ty/tt
     sz = 2.0*tz/tt

     uxd = uxm + uym*tz - uzm*ty
     uyp = uym - uxd*sz
     uzp = uzm + uxd*sy
     uxp = uxd + uyp*tz - uzp*ty

     !   second half-accn
     ux(p) = uxp + beta*exi
     uy(p) = uyp + beta*eyi
     uz(p) = uzp + beta*ezi

  end do


end subroutine push_TM

 
!  ==============================================
!
!     Full 3v particle pusher (Boris rotation scheme)
!
!  TM fields only (s-pol):  (Ex,0,Ez)  (0,By,0)
!  TE fields  (0,Ey,0) (Bx,0,Bz)
!  ==============================================


subroutine push_full3v(p_start,p_finish,dts)
  use module_physvars
  use module_particle_props
  implicit none
  integer, intent(in) :: p_start, p_finish
  real, intent(in) :: dts
  integer :: p
  real :: beta, gam1,tt, sx, sy, sz, tz, ty, tx
  real :: uxd, uyd, uzd, uyp, uzp, uxp, uxm, uym, uzm
  real :: exi, eyi, ezi, bxi, byi, bzi
  real :: xd, yd, zd


  if (itime>0 .and. beam_config==4) focus(1) = x_crit  ! laser tracks n_c

  do p = p_start, p_finish
     beta=q(p)/m(p)*dts*0.5  ! charge/mass constant

     xd = x(p)-focus(1)
     yd = y(p)-focus(2)
     zd = z(p)-focus(3)



     exi = ex(p)
     eyi = ey(p)
     ezi = ez(p)
     bxi = bx(p)
     byi = by(p)
     bzi = bz(p)



     !   first half-accn
     uxm = ux(p) + beta*exi
     uym = uy(p) + beta*eyi
     uzm = uz(p) + beta*ezi

     !   rotation
!     gam1=dsqrt(1.d0 + uxm**2 + uym**2 + uzm**2)
	gam1=1.
     tx = beta*bxi/gam1
     ty = beta*byi/gam1
     tz = beta*bzi/gam1
     tt = 1.0 + tx**2 + ty**2 + tz**2

     sx = 2.0*tx/tt
     sy = 2.0*ty/tt
     sz = 2.0*tz/tt

     uxd = uxm + uym*tz - uzm*ty
     uyd = uym + uzm*tx - uxm*tz
     uzd = uzm + uxm*ty - uym*tx

     uxp = uxm + uyd*sz - uzd*sy
     uyp = uym + uzd*sx - uxd*sz
     uzp = uzm + uxd*sy - uyd*sx

     !   second half-accn
     ux(p) = uxp + beta*exi
     uy(p) = uyp + beta*eyi
     uz(p) = uzp + beta*ezi

  end do


end subroutine push_full3v

!  ===============================================================
!
!                           PUSH_NONREL
!
!   Nonrelativistic particle position update - used with leap-frog scheme
!
!  ===============================================================

subroutine push_nonrel(ips,ipf,delt)

  use module_physvars
  use module_particle_props
  integer, intent(in) :: ips, ipf  ! 1st and last particle numbers
  real, intent(in) :: delt
  integer :: p

  do p=ips,ipf

     x(p)=x(p)+ux(p)*delt
     y(p)=y(p)+uy(p)*delt
     if (idim.eq.3) z(p)=z(p)+uz(p)*delt
     if (p==1) write(70,*) itime*delt,x(p),y(p),ux(p),uy(p)
  end do


end subroutine push_nonrel

!  ===============================================================
!
!                           PMOVE
!
!   Update particle positions - used with relativistic leap-frog scheme
!
!  ===============================================================

subroutine push_x(ips,ipf,delt)

  use module_physvars
  use module_particle_props
  integer, intent(in) :: ips, ipf  ! 1st and last particle numbers
  real, intent(in) :: delt
  integer :: p
  real :: gamma, r_glue2, dx, dy, dz, dr2

  r_glue2 = (max(xl,yl,zl)*glue_radius)**2

  !  relativistic particle push in space

  do p=ips,ipf
! find radius from plasma_centre
     dx = x(p)-plasma_centre(1)
     dy = y(p)-plasma_centre(2)
     dz = z(p)-plasma_centre(3)
     dr2 = dx**2+dy**2+dz**2
     if (dr2 < r_glue2 ) then 
       gamma = sqrt(1.0 + ux(p)**2 + uy(p)**2 + uz(p)**2)
       x(p)=x(p)+ux(p)/gamma*delt
       if (idim > 1) y(p)=y(p)+uy(p)/gamma*delt
       if (idim == 3) z(p)=z(p)+uz(p)/gamma*delt
     else
!	if (mod(me,200).eq.0) write(*,*) "particle ",p," glued at ",sqrt(dr2)
	! leave particle where it is (should flag it to remove from force lists)
     endif
  end do

end subroutine push_x

!  ===============================================================
!
!>      Position update for ion quiet start
!> 	Restrict movement to fraction of aii
!  ===============================================================

subroutine push_restrict(ips,ipf,delt,drmax)

  use module_physvars
  use module_particle_props
  integer, intent(in) :: ips, ipf  ! 1st and last particle numbers
  real, intent(in) :: delt, drmax
  integer :: p

  do p=ips,ipf
       if (abs(ux(p)*delt).lt.drmax) then
	 x(p)=x(p)+ux(p)*delt
       else
	 x(p)=x(p)+drmax*sign(1.d0,ux(p))
       endif
       if (abs(uy(p)*delt).lt.drmax) then
	 y(p)=y(p)+uy(p)*delt
       else
	 y(p)=y(p)+drmax*sign(1.d0,uy(p))
       endif

       if (idim == 3) z(p)=z(p)+uz(p)*delt
  end do

end subroutine push_restrict


!  ===================================================================
!
!                              VELOCITIES
!
!   Calculate velocities from accelerations
! 
!   apply thermodynamic constraint according to ensemble
!
!
!
!  ===================================================================


subroutine velocities(p_start,p_finish,delta_t)


  use module_physvars
  use module_particle_props
  implicit none
  include 'mpif.h'

  real, intent(in) :: delta_t
  integer, intent(in) :: p_start,p_finish  ! min, max particle nos.

  integer p, i, ne_loc, ierr
  real*8, dimension(np_local) :: uhx, uhy, uhz, accx, accy, accz
  real*8 :: sum_vxe, sum_vye, sum_vze, sum_v2e, sum_2ve, acmax
  real :: Te0, Te_uncor, Ti0, Ti_uncor, chie, chii, delta_u
  real*8 :: sum_vxi, sum_vyi, sum_vzi, sum_v2i, sum_2vi, mass_eqm
  real*8 :: global_v2e, global_v2i, gammah, delta_Te, delta_Ti, Te_local, ac_norm

!  Available ensemble modes
!      1 = NVE - total energy conserved
!      2 = NVT - global Te, Ti conserved
!      3 = global NVT electron Te conserved; ions frozen
!      4 = local NVT: each PE keeps Te clamped; ions frozen
!      5 = local NVT, ions only; electrons added at end of run

! Accelerations
  acmax=0.
  do i=1,np_local
     accx(i) = q(i)*ex(i)/m(i)
     accy(i) = q(i)*ey(i)/m(i)
     accz(i) = q(i)*ez(i)/m(i)
     acmax = max(abs(accx(i)),abs(accy(i)),abs(accz(i)),acmax)
  end do

  delta_u = acmax*delta_t

 pusher: select case(scheme)

 case(2)
     ! Conserve kinetic energies of electrons and ions (initial Te const)
     ! adapted from
     !  Allen and Tildesley p230, Brown & Clark, Mol. Phys. 51, 1243 (1984)

     !  Definitions:
     !
     !   unconstrained velocities    uh(x,y,z) = v*(x,y,z)

     !  1)  Unconstrained half-step for electrons

     sum_vxe=0.0  ! partial sums
     sum_vye=0.0
     sum_vze=0.0
     sum_v2e=0.0
     sum_vxi=0.0  ! partial sums (ions)
     sum_vyi=0.0
     sum_vzi=0.0
     sum_v2i=0.0

     do p=1,np_local
           uhx(p) = ux(p) + 0.5*delta_t*accx(p)
           uhy(p) = uy(p) + 0.5*delta_t*accy(p)
           uhz(p) = uz(p) + 0.5*delta_t*accz(p)
           gammah = sqrt(1.0 +uhx(p)**2 + uhy(p)**2 + uhz(p)**2) 
        if (pelabel(p)<=ne) then
           ! electrons
           sum_vxe  = sum_vxe  + uhx(p)/gammah
           sum_vye  = sum_vye  + uhy(p)/gammah
           sum_vze  = sum_vze  + uhz(p)/gammah
           sum_v2e = sum_v2e + gammah-1.

        else if (pelabel(p)<=ne+ni) then
           ! ions
           sum_vxi  = sum_vxi  + uhx(p)/gammah
           sum_vyi  = sum_vyi  + uhy(p)/gammah
           sum_vzi  = sum_vzi  + uhz(p)/gammah
           sum_v2i = sum_v2i + gammah-1.
        endif
     end do
     sum_2ve = sum_vxe**2 + sum_vye**2 + sum_vze**2
     sum_2vi = sum_vxi**2 + sum_vyi**2 + sum_vzi**2

     ! Find global KE sums
     call MPI_ALLREDUCE(sum_v2e, global_v2e, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
     call MPI_ALLREDUCE(sum_v2i, global_v2i, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)


    ! Te_uncor = 0.5*(global_v2/ne - global_2v/ne**2)      !  uncorrected temperature
     Te_uncor = 511*2./3.*global_v2e/ne  ! This should equal 3/2 kT for 3v Maxwellian
     Te0 = Te_keV  ! normalised electron temp

     Ti_uncor = 511*mass_i/mass_e*2./3.*global_v2i/ni  ! This should equal 3/2 kT for 3v Maxwellian
     Ti0 = Ti_keV  ! normalised electron temp


     chie = sqrt(abs(Te0/Te_uncor))     ! multipliers from Temperature ratio - reset once every cycle
     chii = sqrt(abs(Ti0/Ti_uncor)) 

     if (my_rank==0)	write (*,*) 'Te_unc ',Te_uncor,' Te0 ', Te0, ' chie ',chie
     if (my_rank==0)	write (*,*) 'Ti_unc ',Ti_uncor,' Ti0 ', Ti0, ' chii ',chii

     !  3)  Complete full step

     do p=1,np_local
        if (pelabel(p)<=ne) then
           ux(p) = (2*chie-1.)*ux(p) + chie*delta_t*accx(p)
           uy(p) = (2*chie-1.)*uy(p) + chie*delta_t*accy(p)
           if (idim==3) uz(p) = (2*chie-1.)*uz(p) + chie*delta_t*accz(p)

        elseif (pelabel(p)<=ne+ ni) then
           ux(p) = (2*chii-1.)*ux(p) + chii*delta_t*accx(p)
           uy(p) = (2*chii-1.)*uy(p) + chii*delta_t*accy(p)
           if (idim==3) uz(p) = (2*chii-1.)*uz(p) + chii*delta_t*accz(p)
        endif
     end do




     delta_Ti=0.


     delta_Te = 2*Te0*(1.0/chie**2-1.0)       !  heating
     !     heate = heate + delta_Te


  case(3)

! electrons clamped, ions frozen

     sum_vxe=0.0  ! partial sums
     sum_vye=0.0
     sum_vze=0.0
     sum_v2e=0.0
     ne_loc = 0  ! # local electrons

     do p=1,np_local

        if (pelabel(p)<=ne) then
           ! electrons
           ne_loc = ne_loc + 1
           uhx(p) = ux(p) + 0.5*delta_t*accx(p)
           uhy(p) = uy(p) + 0.5*delta_t*accy(p)
           uhz(p) = uz(p) + 0.5*delta_t*accz(p)
           gammah = sqrt(1.0 +uhx(p)**2 + uhy(p)**2 + uhz(p)**2) 
           sum_vxe  = sum_vxe  + uhx(p)/gammah
           sum_vye  = sum_vye  + uhy(p)/gammah
           sum_vze  = sum_vze  + uhz(p)/gammah
           sum_v2e = sum_v2e + gammah-1.

        endif
     end do
     sum_2ve = sum_vxe**2 + sum_vye**2 + sum_vze**2

     ! Find global KE sums
     call MPI_ALLREDUCE(sum_v2e, global_v2e, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)


    ! Te_uncor = 0.5*(global_v2/ne - global_2v/ne**2)      !  uncorrected temperature
     Te_uncor = 511*2./3.*global_v2e/ne  ! This should equal 3/2 kT for 3v Maxwellian
     Te_local = 511*2./3.*sum_v2e/ne_loc
     Te0 = Te_keV  ! normalised electron temp
     chie = sqrt(abs(Te0/Te_uncor))     ! multipliers from Temperature ratio - reset once every cycle
     chie = min(1.25,max(chie,0.75))  ! Set bounds of +- 50%

     if (my_rank==0) write (*,*) 'Te_unc ',Te_uncor,' Te0 ', Te0, ' chie ',chie

     !  3)  Complete full step

     do p=1,np_local
        if (pelabel(p)<=ne) then
           ux(p) = (2*chie-1.)*ux(p) + chie*delta_t*accx(p)
           uy(p) = (2*chie-1.)*uy(p) + chie*delta_t*accy(p)
           if (idim==3) uz(p) = (2*chie-1.)*uz(p) + chie*delta_t*accz(p)
        endif
     end do

     delta_Ti=0.
     delta_Te = 2*Te0*(1.0/chie**2-1.0)       !  heating


  case(4)

! electrons clamped locally, ions frozen
! - require T=T_e on each PE to avoid local drifts

     sum_vxe=0.0  ! partial sums
     sum_vye=0.0
     sum_vze=0.0
     sum_v2e=0.0
     ne_loc = 0  ! # local electrons

     do p=1,np_local

        if (pelabel(p)<=ne) then
           ! electrons
           ne_loc = ne_loc+1  
           uhx(p) = ux(p) + 0.5*delta_t*accx(p)
           uhy(p) = uy(p) + 0.5*delta_t*accy(p)
           uhz(p) = uz(p) + 0.5*delta_t*accz(p)
           gammah = sqrt(1.0 +uhx(p)**2 + uhy(p)**2 + uhz(p)**2) 
           sum_vxe  = sum_vxe  + uhx(p)/gammah
           sum_vye  = sum_vye  + uhy(p)/gammah
           sum_vze  = sum_vze  + uhz(p)/gammah
           sum_v2e = sum_v2e + gammah-1.

        endif
     end do
     sum_2ve = sum_vxe**2 + sum_vye**2 + sum_vze**2

     ! Find local KE temp

    ! Te_uncor = 0.5*(global_v2/ne - global_2v/ne**2)      !  uncorrected temperature
     Te_uncor = 511*2./3.*sum_v2e/ne_loc  ! This should equal 3/2 kT for 3v Maxwellian
     Te0 = Te_keV  ! normalised electron temp
     ! exponent of chie should be 1/2 - take 1/3 to soften oscillations
     chie = (abs(Te0/Te_uncor))**0.5     ! multipliers from Temperature ratio - reset once every cycle

     chie = min(1.25,max(chie,0.75))  ! Set bounds of +- 50%

     if (my_rank.eq.0) write (*,*) 'Te_unc ',Te_uncor,' Te0 ', Te0, ' chie ',chie

     !  3)  Complete full step

     do p=1,np_local
        if (pelabel(p)<=ne) then
           ux(p) = (2*chie-1.)*ux(p) + chie*delta_t*accx(p)
           uy(p) = (2*chie-1.)*uy(p) + chie*delta_t*accy(p)
           if (idim==3) uz(p) = (2*chie-1.)*uz(p) + chie*delta_t*accz(p)
        endif
     end do

     delta_Ti=0.
     delta_Te = 2*Te0*(1.0/chie**2-1.0)       !  heating


  case(5)

     ! Conserve kinetic energy of ions only (initial Ti const)
     mass_eqm = 20.  ! artificial ion mass for eqm stage
     sum_vxi=0.0  ! partial sums (ions)
     sum_vyi=0.0
     sum_vzi=0.0
     sum_v2i=0.0

!  Normalize accelerations according to initial velocities
     ac_norm = 15*vti/dt

!     do p=1,np_local
!	accx(p)=min(accx(p),accx(p)/max(1.d0,acmax)*ac_norm)
!	accy(p)=min(accy(p),accy(p)/max(1.d0,acmax)*ac_norm)
!	accz(p)=min(accz(p),accz(p)/max(1.d0,acmax)*ac_norm)
 !    end do

     do p=1,np_local
           uhx(p) = ux(p) + 0.5*delta_t*accx(p)
           uhy(p) = uy(p) + 0.5*delta_t*accy(p)
           uhz(p) = uz(p) + 0.5*delta_t*accz(p)

        if (pelabel(p)>=ne) then
           ! ions
           sum_vxi  = sum_vxi  + uhx(p)
           sum_vyi  = sum_vyi  + uhy(p)
           sum_vzi  = sum_vzi  + uhz(p)
           sum_v2i = sum_v2i + 0.5*m(p)*(uhx(p)**2+uhy(p)**2+uhz(p)**2)
        endif
     end do
     sum_2vi = sum_vxi**2 + sum_vyi**2 + sum_vzi**2

     ! Find global KE sums
!     call MPI_ALLREDUCE(sum_v2i, global_v2i, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
!     Ti_uncor = 511*2./3.*sum_v2i/np_local  ! This should equal 3/2 kT for 3v Maxwellian
     Ti_uncor = sum_v2i/np_local  ! This should equal 3/2 kT for 3v Maxwellian
     Ti0 = Ti_keV  ! normalised target ion temp


     chii = sqrt(abs(Ti0/Ti_uncor)) 
     chii = min(1.2,max(chii,0.3))  ! Set bounds of +- 50%


     !  3)  Complete full step

     do p=1,np_local

        if (pelabel(p)>=ne) then
       ! make ions lighter for eqm phase
           ux(p) = (2*chii-1.)*ux(p) + chii*delta_t*accx(p)
           uy(p) = (2*chii-1.)*uy(p) + chii*delta_t*accy(p)
           if (idim==3) uz(p) = (2*chii-1.)*uz(p) + chii*delta_t*accz(p)
!           ux(p) = ux(p)*sqrt(Ti0/Ti_uncor)
!           uy(p) = uy(p)*sqrt(Ti0/Ti_uncor)
!           uz(p) = uz(p)*sqrt(Ti0/Ti_uncor)

        endif
     end do
     delta_Ti = 2*Ti0*(1.0/chii**2-1.0)       !  heating
     if (my_rank==0)	then
	write (*,*) 'Ti_unc ',Ti_uncor,' Ti0 ', Ti0, ' chii ',chii,' heating:',delta_Ti,' bond',bond_const
        write (*,*) 'Max delta-ux: ',delta_u
     endif


  case default
  if (my_rank==0) then
	write(6,*) 'PEPC | velocities scheme',scheme
	write(24,*) 'PEPC | velocities scheme',scheme
  endif
     ! unconstrained motion by default (scheme=1,7)
   if (idim==3) then
     do p = p_start, p_finish
	ux(p) = ux(p) + delta_t * accx(p)
	uy(p) = uy(p) + delta_t * accy(p)
	uz(p) = uz(p) + delta_t * accz(p)
     end do
   else if (idim==2) then
     do p = p_start, p_finish
	ux(p) = ux(p) + delta_t * accx(p)
	uy(p) = uy(p) + delta_t * accy(p)
     end do
   else
     do p = p_start, p_finish
	ux(p) = ux(p) + delta_t * accx(p)
     end do
   endif
 
  end select pusher

end subroutine velocities



end module module_integration_scheme
