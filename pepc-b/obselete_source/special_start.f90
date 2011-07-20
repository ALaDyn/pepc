! ==============================================
!
!                SPECIAL_START
!
!  Initialise set of particles for special test configs
!
! ==============================================

subroutine special_start(iconf)

  use physvars
  use treevars
  use utils
  implicit none

  integer, intent(in) :: iconf
  integer :: i,iseed
  real :: rs,v0,gamma0,vt,xt,yt,zt,thetvel,phivel
  real :: dpx, vx_beam, vy_beam, vz_beam, st, ct, sp, cp, volb
  integer :: npb
  real :: gamma, Ndebye

  iseed = -17-me

  ! Define default container parameters in absence of plasma target
  Vplas = x_plasma * y_plasma * z_plasma
  Aplas = x_plasma * y_plasma
  focus = (/ xl/4., yl/2., zl/2. /) ! Centre of laser focal spot
  plasma_centre =  (/ xl/2., yl/2., zl/2. /) ! Centre of plasma
  Qplas = ne

  config: select case(iconf)

  case(1)        ! electron disc at x=0 with 2pi phase spread
     qe=-1.0
     mass_e=1.0
     qi=1.0
     mass_i=mass_ratio

     r_beam=sigma/2.
     gamma0=sqrt(1+vosc**2/2.)

     do while (i < nep)
        yt = r_beam*(2*rano(iseed)-1.)

        zt = r_beam*(2*rano(iseed)-1.)
        if (yt**2 + zt**2 <= r_beam**2 ) then
           i = i+1
           x(i) = 2*pi*rano(iseed)
           y(i) = yt + focus(2)
           z(i) = zt + focus(3)
           ux(i) = 0.
           uy(i) = 0.
           uz(i) = 0.

        endif
     end do
     q(1:nep) = qe           ! plasma electrons
     m(1:nep) = mass_e            ! electron mass

  case(2)        ! electron disc at laser focus with 2pi phase spread

     qe=-1.0
     mass_e=1.0
     qi=1.0
     mass_i=mass_ratio

     r_beam=sigma/2.
     gamma0=sqrt(1+vosc**2/2.)


     x_crit = focus(1)
     do while (i < nep)
        yt = r_beam*(2*rano(iseed)-1.)

        zt = r_beam*(2*rano(iseed)-1.)
        if (yt**2 + zt**2 <= r_beam**2 ) then
           i = i+1
           x(i) = focus(1)
           y(i) = yt + focus(2)
           z(i) = zt + focus(3)
           ux(i) = 0.
           uy(i) = 0.
           uz(i) = 0.

        endif
     end do
     q(1:nep) = qe           ! plasma electrons
     m(1:nep) = mass_e            ! electron mass


  case(3)       ! random placement in box with vte in x,y plane for B test 

     qe=-1.0
     mass_e=1.0
     qi=1.0
     mass_i=mass_ratio

     do i=1,npp
        if (i<=nep) then
           vt=vte
           q(i) = qe        ! plasma electrons
           m(i) = mass_e 
        else
           vt = vti
           q(i) = qi        ! plasma ions (need Z* here)
           m(i) = mass_i    ! ion mass

        endif
        xt =  xl * rano(iseed) 
        yt =  yl * rano(iseed)
!        zt =  zl * rano(iseed)
        zt = zl/2.-z_plasma/2.+z_plasma*i/1./npp
        x(i) = xt
        y(i) = yt
        z(i) = zt
        rs = rano(iseed)
        v0 = vt*sqrt(-2.0*log(rs))
        thetvel = 2*pi*rano(iseed)
        phivel = pi*rano(iseed)
        ux(i) = v0*cos(thetvel)
        uy(i) = v0*sin(thetvel)
!        uz(i) = v0/10.
        uz(i) = v0*cos(phivel)
     end do

  case(4)      !  Special config for FMM comparison
     open(30,file='fmm_config.data')
     do i=1,npp
        read(30,*) x(i), y(i), z(i), q(i)
        m(i) = 1.
        ux(i)=0.
        uy(i)=0.
        uz(i)=0.
     end do
     force_const=1.
     dt=0.
     eps=0.
     close(30)

  case(5)     ! beam initialised along x-axis and rotated by theta, phi


     write (*,*) 'Preparing particle beam'
     npb = max(ni,ne)  ! Take whichever species switched on
     np_beam=0 ! Discard normal beam mode (plasma+beam)
     npp = max(nep,nip)

     ct = cos(theta_beam)
     st = sin(theta_beam)
     cp = cos(phi_beam)
     sp = sin(phi_beam)
     vz_beam = u_beam*st
     vx_beam = u_beam*ct*cp
     vy_beam = u_beam*ct*sp
     Volb = pi*r_beam**2*x_beam
     qe = Volb*rho_beam/npb    ! charge
     mass_e = mass_beam
     mass_i = mass_beam
     qi=-qe
     dpx=x_beam/npb           ! x-axis spacing
     i=0
     do while (i < npp)
        yt = r_beam*(2*rano(iseed)-1.)
        zt = r_beam*(2*rano(iseed)-1.)

        if (yt**2 + zt**2 <= r_beam**2 ) then
           i = i+1
           x(i)= start_beam + x_beam*rano(iseed)
           y(i) = yt + focus(2)
           z(i) = zt + focus(3)
           ux(i) = vx_beam
           uy(i) = vy_beam
           uz(i) = vz_beam

        endif
     end do
     q(1:npp) = qe           ! plasma electrons
     m(1:npp) = abs(qe*mass_beam)  ! electron mass

     !     write (*,'(a30/(5f13.5))') 'Beam config:', &
     !          (x(i),y(i),z(i),q(i),m(i),i=1,npp)


  case(6)       ! random placement in box xl*xl*xl for periodic system 

     ! Use Debye units 
     ! - lengths normalized to lambda_De
     ! - velocities to v_te
     ! - time to wpe^-1
     ! Effective density adjusted via box length
     ! TODO:   gamma, Ndebye need to be made global (physvars.f90)

     yl=xl
     zl=xl
     qe=-1.0
     mass_e=1.0
     qi=1.0
     mass_i=mass_ratio
     !  Inter-electron spacing
     a_ee = xl/(4*pi/3.*ne)**(1./3.)
     !  Inter-ion spacing
     a_ii = a_ee*Zion**(1./3.) 

     !  Renormalise softening parameter
     eps = eps * a_ee

     !  Equivalent Gamma
     gamma= a_ee**2/3.
     !  # electrons in Debye sphere
     Ndebye = (3*gamma)**(-1.5)
     !  Adjust force constant
     force_const = 1./3./Ndebye

     write (*,*) "Ewald Setup:"
     write (*,'(5(a30,f15.5/))') "Particle spacing ",a_ee, &
          "Softening parameter:",eps, &
          "Gamma:",gamma, &
          "# electrons in Debye sphere:",Ndebye, &
          "force const:",force_const

     do i=1,npp
        if (i<=nep) then
           q(i) = qe        ! plasma electrons
           m(i) = mass_e 
        else
           q(i) = qi        ! plasma ions (need Z* here)
           m(i) = mass_i    ! ion mass

        endif
        xt =  xl * rano(iseed) 
        yt =  xl * rano(iseed)
        zt =  xl * rano(iseed)
        x(i) = xt
        y(i) = yt
        z(i) = zt
     end do


     ! Setup 3v Maxwellian electrons
     ! TODO: max random seed processor-dependent

     call maxwell1(ux,nppm,1,nep,vte)
     call maxwell1(uy,nppm,1,nep,vte)
     call maxwell1(uz,nppm,1,nep,vte)
     call scramble_v(1,nep)   ! remove x,y,z correlations
     ! Ions cold    
     call cold_start(nep+1,nip)


  end select config

  ex(1:npp) = 0.
  ey(1:npp) = 0.
  ez(1:npp) = 0.
  bx(1:npp) = 0.
  by(1:npp) = 0.
  bz(1:npp) = 0.
  ax(1:npp) = 0.
  ay(1:npp) = 0.
  az(1:npp) = 0.
  axo(1:npp) = 0.
  ayo(1:npp) = 0.
  azo(1:npp) = 0.
  pot(1:npp) = 0.


  pelabel(1:nep)       = me * nep + (/(i, i = 1, nep)/)      ! Electron labels
  pelabel(nep + 1:npp) = ne + me * nip + (/(i, i = 1, nip)/) ! Ion labels
  pepid(1:npp) = me                ! processor ID
  work(1:npp) = 1.
  call push_full3v(1,npp,dt/2.)

end subroutine special_start



