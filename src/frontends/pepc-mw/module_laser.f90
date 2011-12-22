!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>  Encapsulates anything that is concerned with laser setup and laser-particle interaction
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_laser
    implicit none
    save
    private

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  public parameter declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! possible choices for 1st digit of beam_config_in
    integer, public, parameter :: LASER_FOCUS_STATIC           = 0
    integer, public, parameter :: LASER_FOCUS_FOLLOW_XCRIT     = 1
    ! possible choices for 2nd digit of beam_config_in
    integer, public, parameter :: LASER_ENVELOPE_OFF           = 9
    integer, public, parameter :: LASER_ENVELOPE_CONSTANT      = 1
    integer, public, parameter :: LASER_ENVELOPE_SINSQUARED    = 2
    integer, public, parameter :: LASER_ENVELOPE_LINEAR_RISE   = 3
    integer, public, parameter :: LASER_ENVELOPE_LINEAR_FALL   = 4
    ! possible choices for 3rd digit of beam_config_in
    integer, public, parameter :: LASER_MODEL_UNIFORM_SINUSOID = 1
    integer, public, parameter :: LASER_MODEL_UNIFORM_SINUSOID_ONLYELECTRONS = 2
    ! possible choices for 4th digit of beam_config_in
    integer, public, parameter :: LASER_POLARIZATION_X         = 1
    integer, public, parameter :: LASER_POLARIZATION_Y         = 2
    integer, public, parameter :: LASER_POLARIZATION_Z         = 3

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    integer, public :: beam_config_in = 0 !< Particle or laser beam switch including variations

    integer, public :: beam_focus        = 0 !< 1st digit of beam_config_in
    integer, public :: beam_envelope     = 0 !< 2nd digit of beam_config_in
    integer, public :: beam_model        = 0 !< 3rd digit of beam_config_in
    integer, public :: beam_polarization = 0 !< 4th digit of beam_config_in
    character(40), public :: config_names(4)

    real*8, public :: omega     =  0.5    !< frequency
    real*8, public :: omega_wpl =  0.0    !< frequency omega in wpl_e
    real, public :: omega_hz  =  0.0    !< frequency omega in Hz
    real, public :: lambda    =  0.0    !< laser wavelength
    real, public :: lambda_nm =  0.0    !< laser wavelength in nm
    real*8, public :: rhocrit_nm3 = 0.      !< critical electron density in electrons per nm^3
    real*8, public :: I0_Wpercm2 = 0.   !< initial intensity in W/cm^2
    real*8, public :: E0 = 0.   !< laser field strength amplitude
    real*8, public :: vosc      =  0.1    !< pump strength

    real*8, public :: I_laser           !< Laser intensity (= amplitude**2), updated by each call to laser()
    real*8, public :: E_laser           !< Laser amplitude (= intensity**1/2), updated by each call to laser()
    real, public :: x_offset  =  0.     !< coordinate offset
    real, public :: z_offset  =  0.     !< coordinate offset
    real*8, public :: focus(3)  = [0., 0., 0.] !< centre of focal spot
    real*8, public :: t_pulse      = 0.     !< pulse duration (in simulation time units)
    real*8, public :: t_pulse_fs   = 10.     !< pulse duration (in fs)
    real*8, public :: sigma   =  1.0    !< 1/e pulse width (c/omega_p)
    real, public :: theta_inc =  0.     !< angle of incidence
    real, public :: rho_track =  1.5    !< tracking density for x_crit (/nc)
    real*8, public :: t_laser   = 0.      !< run time after laser switched on (in simulation time units)
    real*8, public :: elaser        !< deposited laser energy
    real, public :: propag_laser  !< distance travelled by laser after rezoning
    real, public :: intensity     !< normalised intensity = 0.5*vosc^2*omega^2
    real, public :: window_min    !< start of wakefield plasma
    real, public :: rezone_frac=0.75     !< Fraction of box to cross before rezoning switched on
    real*8, public :: fpon_max, ampl_max !< max amplitudes
    real*8, public :: navcycle     ! # timesteps in a laser cycle

    real*4, public, allocatable :: rho_helm(:)  !< Helmholtz density
    complex*16, public, allocatable :: Az_helm(:)   !< Helmholtz vector potential

    character*7, public :: beam_configs(0:9)=(/ &
    'off    ','beam   ','i-beam ','laser-u','ES pond','LWFA   ', &
    'EMplane','EM pond',' dust  ','       ' /)

    ! declarations from pepcb physvars module
    integer, public :: nxh=100 !< 1D Helmholtz grid dimension
    real, public :: dxh !< HH grid spacing
    real, public :: xh_start=0. !< Start point of Helmholtz grid
    real, public :: xh_end=10.  !< End point of Helmholtz grid
    real, public :: rho_upper      !< shelf/profile density above x_crit (/nc)
    real*8, public :: x_crit         !< critical surface
    real*8, public :: theta_beam = 0. !< beam elevation angle


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    public laser_setup
    public laser_update
    public force_laser
    public force_laser_at
    public laser_hist
    public emplane
    public empond


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    real, parameter :: pi=3.141592654


contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Setup all local module variables and data that depends on them
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine laser_setup()
        use physvars
        use module_units
        implicit none

        beam_focus        = modulo(beam_config_in / 1000, 10)
        beam_envelope     = modulo(beam_config_in /  100, 10)
        beam_model        = modulo(beam_config_in /   10, 10)
        beam_polarization = modulo(beam_config_in /    1, 10)
        navcycle          = 2*pi/dt/omega  ! # timesteps in a laser cycle
        t_pulse           = t_pulse_fs / unit_t0_in_fs

    end subroutine laser_setup


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Laser propagation according to beam_config
    !>
    !> Calculates time-dependend laser intensity
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function laser_envelope(ienvelope, tpulse, tlaser, vosc, envelope_name)
        implicit none
        real*8 :: laser_envelope
        integer, intent(in) :: ienvelope
        real*8, intent(in) :: tpulse, tlaser, vosc
        character(30), intent(out) :: envelope_name

        if (tlaser < 0.) then
            laser_envelope = 0.
            return
        endif

        select case(ienvelope)
            case(LASER_ENVELOPE_OFF)
                envelope_name = 'LASER_ENVELOPE_OFF'
                laser_envelope = 0.

            case(LASER_ENVELOPE_CONSTANT)
                envelope_name = 'LASER_ENVELOPE_CONSTANT'
                laser_envelope = 1.

            case(LASER_ENVELOPE_SINSQUARED)
                envelope_name = 'LASER_ENVELOPE_SINSQUARED'
                if (tlaser < tpulse) then
                    laser_envelope = max(0._8,sin(pi*tlaser/tpulse)**2)
                else
                    laser_envelope = 0.
                endif

            case(LASER_ENVELOPE_LINEAR_RISE)
                envelope_name = 'LASER_ENVELOPE_LINEAR_RISE'
                laser_envelope = min(1._8,tlaser/tpulse)

            case(LASER_ENVELOPE_LINEAR_FALL)
                envelope_name = 'LASER_ENVELOPE_LINEAR_FALL'
                laser_envelope = max(0._8,1.-tlaser/tpulse)

            case default
                envelope_name = 'LASER_ENVELOPE_OFF(undefined)'
                laser_envelope = 0.
        end select

    end function laser_envelope


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Laser propagation according to beam_config
    !>
    !> Calculates time-dependend laser focus position
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine laser_focus(ifocus, tpulse, tlaser, vosc, focus, focus_name)
        implicit none
        integer, intent(in) :: ifocus
        real*8, intent(in) :: tpulse, tlaser, vosc
        real*8, intent(inout) :: focus(3)
        character(30), intent(out) :: focus_name


        select case(ifocus)
            case(LASER_FOCUS_STATIC)
                focus_name = 'LASER_FOCUS_STATIC'
                ! leave focus unchanged

            case(LASER_FOCUS_FOLLOW_XCRIT)
                focus_name = 'LASER_FOCUS_FOLLOW_XCRIT'
                ! TODO: calculate x_crit
                focus(1) = x_crit  ! laser tracks n_c

            case default
                focus_name = 'LASER_FOCUS_STATIC(undefined)'
              ! leave focus unchanged
        end select


    !           case(44,54)  ! Helmholtz solver for vector potential
    !             ! Factor-in pulse shape
    !             amplitude = sqrt(I_laser)  ! amplitude of incoming wave
    !             call density_helmholtz
    !             call em_helmholtz(my_rank,itime,nxh,dxh,theta_inc,amplitude,omega,rho_helm,Az_helm)
    !             Azr(1:nxh) = Real(Az_helm(1:nxh)*exp(yi*pha))
    !             ampl_max = maxval(Azr)
    !             ampl_min = minval(Azr)
    !             if (ampl_max.lt.abs(ampl_min)) ampl_max=-ampl_min
    !
    !           case(64)  ! Helmholtz solver for vector potential circ pol
    !              ! Factor-in pulse shape
    !             amplitude = sqrt(I_laser)  ! amplitude of incoming wave
    !             call density_helmholtz
    !             call em_helmholtz(my_rank,itime,nxh,dxh,theta_inc,amplitude,omega,rho_helm,Az_helm)
    !             Azr(1:nxh) = Real(Az_helm(1:nxh))
    !             ampl_max = maxval(Azr)
    !             ampl_min = minval(Azr)
    !
    !           case(05)  ! propagating fpond
    !             !  Trigger rezoning if laser rezone_frac of the way through plasma
    !             ! - Only works after restart at present
    !             if (restart .and. beam_config ==5 .and. focus(1) >= window_min + x_plasma*rezone_frac) then
    !               if (my_rank==0) then
    !                   write (*,*) 'REZONE'
    !                   write (*,*) 'REZONING CURRENTLY NOT IMPLEMENTED'
    !                  !           read (*,*) go
    !               endif
    !
    !               ! TODO: call rezone
    !               !        window_min = window_min + dt
    !             else
    !               focus(1) = focus(1) + dt  ! propagate forward by c*dt - can include v_g here
    !               propag_laser=propag_laser + dt
    !             endif

    end subroutine



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Laser propagation according to beam_config
    !>
    !> Calculates time-dependend laser intensity and corrects laser focus position
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine laser_update()

        use physvars
        use module_units
        implicit none

        if (beam_config_in == 0) return

        t_laser = t_laser + dt

        !  Laser pulse envelope
        !  =====================
        E_laser = E0 * laser_envelope(beam_envelope, t_pulse, t_laser, vosc, config_names(2))
        I_Laser = unit_epsilon0*unit_c/2.*E_laser**2

        !  Laser focal position and rezoning
        !  ==================================
        call laser_focus(beam_focus, t_pulse, t_laser, vosc, focus, config_names(1))

        !    Deposited laser energy
        !    ==================================
        ! TODO: do something with this
        laser_energy: select case(modulo(beam_config_in, 10))
            case(4)
                elaser = 3./8.*omega**2*sigma**2*vosc**2*t_laser
            case(5)
                elaser = 3./8.*omega**2*sigma**2*vosc**2*t_pulse
            case default
                elaser = 0.
        end select laser_energy

    end subroutine laser_update



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Calculate forces due to external fields
    !> eg: laser, stationary B, E fields depending on beam configuration
    !>
    !> @param[in] p_start  minimum particle number
    !> @param[in] p_finish maximum particle number
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine force_laser(p_start, p_finish)

        use physvars
        implicit none
		  include 'mpif.h'
        real*8 :: E_pon(3), B_em(3), Phipon

        integer, intent(in) :: p_start,p_finish  ! min, max particle nos.
        integer :: p

        ! Include force from laser/external field on electrons
        do p = p_start, p_finish
            call force_laser_at(particles(p)%x, particles(p)%data%v(1), E_pon, B_em, Phipon)

            !if (p==p_start) then
            !  write(17,'(9f8.4)') [Ex(p), Ey(p), Ez(p)], E_pon, [Ex(p), Ey(p), Ez(p)] + E_pon
            !endif

            !  ions assumed not to feel laser, so zero fields
            if ((beam_model == LASER_MODEL_UNIFORM_SINUSOID_ONLYELECTRONS) .and. (particles(p)%data%q > 0.)) then
                E_pon = [ 0., 0., 0. ]
                B_em  = [ 0., 0., 0. ]
            endif

            fpon_max = max(fpon_max, abs(E_pon(1)))
            ! Add external fields to new particle field
            particles(p)%results%pot     = particles(p)%results%pot     + Phipon
            particles(p)%results%e(1:3)  = particles(p)%results%E(1:3)  + E_pon

        end do

    end subroutine force_laser


    subroutine force_laser_at(pos, ux, E_pon, B_em, Phipon)
        use physvars, only : plasma_centre, zl
        implicit none
        real*8, intent(out) :: E_pon(3), B_em(3), Phipon
        real*8, intent(in) :: pos(3), ux

        real*8 :: dc(3)  ! positions relative to centre of plasma
        real*8 :: uxd ! x-momentum
        real*8 :: df(3)  ! position relative to laser focus
        real*8 :: ez_em, az_em, rt

        dxh = (xh_end-xh_start)/nxh  ! HH grid spacing

        df = pos - focus
        dc = pos - plasma_centre
        rt = sqrt(dc(1)**2+dc(2)**2)

        select case(beam_model)

            case(LASER_MODEL_UNIFORM_SINUSOID, LASER_MODEL_UNIFORM_SINUSOID_ONLYELECTRONS)
                select case(beam_polarization)
                    case(LASER_POLARIZATION_X)
                        E_pon  = [ E_laser*sin(omega*t_laser), 0._8, 0._8 ]
                        B_em   = [ 0._8, 0._8, 0._8]
                        Phipon = 0. ! vosc**2/(4.*q) might disturb energy computation
                    case(LASER_POLARIZATION_Y)
                        E_pon = [ 0._8, E_laser*sin(omega*t_laser), 0._8 ]
                        B_em  = [ 0._8, 0._8, 0._8]
                        Phipon = 0.
                    case(LASER_POLARIZATION_Z)
                        E_pon = [ 0._8, 0._8, E_laser*sin(omega*t_laser) ]
                        B_em  = [ 0._8, 0._8, 0._8]
                        Phipon = 0.
                end select


            !		        case(04)  ! standing wave fpond, sin2 pulse
            !		           call fpond( tlaser, tpulse, sigma, vosc, omega, rho_upper, &
            !		                df(1), df(2), df(3), E_pon(1), E_pon(2), E_pon(3), Phipon)
            !
            !                   B_em  = [ 0., 0., 0.]
            !
            !		        case(14)  ! standing wave fpond, linear rise-time, fields reduced
            !		           call fpond_lin( tlaser, tpulse,sigma,vosc,omega,rho_upper, &
            !		                df(1), df(2), df(3), E_pon(1), E_pon(2), E_pon(3), Phipon)
            !
            !                   E_pon(2:3) = E_pon(2:3) / 50.
            !                   B_em  = [ 0., 0., 0.]
            !
            !		        case(24)  ! oblique incidence standing wave, s-pol
            !		           call emobliq( tlaser, tpulse,sigma,vosc,omega,theta_inc,rho_upper, &
            !		                df(1), df(2), df(3), E_pon(1), E_pon(2), E_pon(3), Phipon, ez_em, B_em(1), B_em(2))
            !		           B_em(3) = 0.
            !
            !		        case(34) ! fpond with fully EM components Gaussian spot, sin^2 pulse
            !		           call empond(tlaser,tpulse,sigma,vosc,omega,&
            !		                df(1), df(2), df(3),E_pon(3),B_em(2),B_em(1),az_em,Phipon)
            !
            !		           E_pon(1:2) = [ 0., 0. ]
            !		           B_em(3)    = 0.
            !
            !		        case(44,54)  ! fpond derived from Az_helm; both linear & sin2 pulse forms
            !		           df(1) = x
            !		           uxd  = ux
            !		           call fpond_helm( tlaser, tpulse,sigma,vosc,omega, &
            !		                df(1), df(2), df(3), uxd, Az_helm, nxh, xh_start, xh_end, dxh, focus(1), &
            !		                E_pon(1), E_pon(2), E_pon(3), Phipon)
            !
            !                   B_em  = [ 0., 0., 0.]
            !
            !		        case(64)  ! fpond derived from Az_helm, c-pol light
            !		           df(1) = x
            !		           uxd  = ux
            !		           call fpond_helmc( tlaser, tpulse,sigma,vosc,omega, &
            !		                df(1), df(2), df(3), uxd,Az_helm,nxh,xh_start, xh_end, dxh, focus(1), &
            !		                E_pon(1), E_pon(2), E_pon(3), Phipon)
            !
            !                   B_em  = [ 0., 0., 0.]
            !
            !                case(94)  ! standing wave fpond with transverse fields artificially reduced
            !                   call fpond( tlaser, tpulse,sigma,vosc,omega,rho_upper, &
            !                        df(1), df(2), df(3), E_pon(1), E_pon(2), E_pon(3), Phipon)
            !
            !                   E_pon(2:3) = E_pon(2:3) / 10.
            !                   B_em  = [ 0., 0., 0.]
            !
            !		        case(05)  ! propagating fpond
            !		           call laser_bullet( tlaser, focus(1), tpulse,sigma,vosc,omega, &
            !		                 df(1), df(2), df(3), E_pon(1), E_pon(2), E_pon(3), Phipon)
            !
            !                   B_em  = [ 0., 0., 0.]
            !
            !		        case(06) ! plane wave with Gaussian spot, sin^2 pulse
            !		           call emplane(tlaser,tpulse,sigma,vosc,omega,&
            !		                 df(1), df(2), df(3),E_pon(3),B_em(2),B_em(1),az_em,Phipon)
            !
            !                   E_pon(1:2) = [ 0., 0. ]
            !                   B_em(3)    = 0.
            !
            !		        case(16) ! plane wave with Gaussian spot, linear rise-time
            !		           call emplane_lin(tlaser,tpulse,sigma,vosc,omega,&
            !		                 df(1), df(2), df(3),E_pon(3),B_em(2),B_em(1),az_em,Phipon)
            !
            !                   E_pon(1:2) = [ 0., 0. ]
            !                   B_em(3)    = 0.
            !
            !                case(07)  ! Constant B in z-direction - either charge
            !                   E_pon = [ 0., 0., 0. ]
            !                   B_em  = [ 0., 0., vosc]
            !
            !                case(17)  !  Z-pinch: circular B in x,y
            !                   E_pon = [ 0., 0., 0. ]
            !                   B_em  = [ -vosc*dc(2)/rt, vosc*dc(1)/rt, 0._8]
            !
            !                case(27)  !  Tokamak: circular B in theta + const Bz
            !                   E_pon = [ 0., 0., 0. ]
            !                   B_em  = [ -vosc*dc(2)/rt, vosc*dc(1)/rt, 0._8]
            !
            !                case(37)  !  Mirror in Bz
            !                   E_pon = [ 0., 0., 0. ]
            !                   B_em  = [ 0., 0., 0. ]
            !
            !                   if (dc(3)>-zl/2. .and. dc(3)<zl/2.) then
            !                     B_em(1) = -vosc/5.*dc(1)/rt*(2*dc(3)/zl)
            !                     B_em(2) = -vosc/5.*dc(2)/rt*(2*dc(3)/zl)
            !                     B_em(3) =  vosc   *         (2*dc(3)/zl)**2 + vosc/5.
            !                   endif

            case default  ! no laser
                E_pon = [ 0., 0., 0. ]
                B_em  = [ 0., 0., 0. ]

        end select

    end subroutine force_laser_at



    ! ==================================================================
    !
    !                        PONDEROMOTIVE FORCE
    !
    !  Compute relativistic fpond for standing wave field
    !
    ! ==================================================================

    subroutine fpond_gauss(t,tpulse,sigma0,vosc,omega,xd,yd,zd,epon_x,epon_y,epon_z,phipon)


        real, intent(in) :: t ! time
        real, intent(in) :: tpulse ! pulse duration
        real, intent(in) :: vosc ! quiver strength
        real, intent(in) :: sigma0 ! pulse width (1/e)
        real, intent(in) :: omega ! laser frequency
        real*8, intent(in) :: xd,yd,zd ! position to evaluate force; xd is away from target

        real*8, intent(out) :: phipon, epon_x, epon_y, epon_z ! pond. potential and fields

        real*8 :: xf, yf, zf, Tpon, Rpon, Xpon, gamma, atten

        !  linear rise
        Tpon = 2*vosc**2*min(1.,t/tpulse) * (sin(omega*t))**2

        sigma = sigma0*sqrt(1.+abs(xd)**2/4./sigma0**2) ! take Rayleigh length 2*sigma0

        Rpon = exp((-yd**2-zd**2)/sigma**2)

        !  Standing wave in vacuum; evanescent in overdense plasma
        if (xd.ge.0) then
            xf = sin(2*omega*xd)
            Xpon = cos(omega*xd)**2
        else
            xf = -2*exp(2*xd)
            Xpon = exp(2*xd)
        endif

        yf = -yd/sigma**2     ! suppress radial fpond for scale-model
        zf = -zd/sigma**2
        atten = sigma0**2/sigma**2
        phipon = Tpon*Xpon*Rpon*atten  ! include attenuation factor

        gamma = sqrt(1.+abs(phipon)/2.)

        Epon_x = Tpon*Rpon/4/gamma*xf*atten
        Epon_y = phipon/2/gamma*yf
        Epon_z = phipon/2/gamma*zf
    ! Epon_y = 0.
    ! Epon_z = 0.

    end subroutine fpond_gauss



    ! ==================================================================
    !
    !                        PONDEROMOTIVE FORCE
    !
    !  Compute ponderomotive fields from Helmholtz solution for EM standing wave
    !
    ! ==================================================================

    subroutine fpond_helm(t,tpulse,sigma_in,vosc,omega, &
    x,y,z,ux,Az,nxh,xh_start,xh_end,dxh,x_crit, &
    epon_x,epon_y,epon_z,phipon)

        !   call fpond_helm( tlaser, tpulse,sigma,vosc,omega, &
        !                   xd,yd,zd,uxd,Az_helm,nxh,xh_start, xh_end, dxh, focus(1), &
        !           epon_x,epon_y,epon_z,phipon)

        use module_units
        implicit none
        real, intent(in) :: t     !< time
        real, intent(in) :: tpulse !< pulse duration/rise time
        real, intent(in) :: vosc !< quiver strength
        real*8, intent(in) :: sigma_in !< pulse width (1/e)
        real, intent(in) :: omega !< laser frequency
        real*8, intent(in) :: x,y,z !< position to evaluate force;
                                  !< x absolute; y,z relative to laser axis
        real*8, intent(in) :: ux !< forward momentum for gamma factor
        integer, intent(in) :: nxh !< # 1D Helmholtz grid points
        real, intent(in) :: xh_start !< Start point of HH grid
        real, intent(in) :: xh_end !< End point of HH grid
        real, intent(in) :: dxh !< HH grid spacing
        real, intent(in) :: x_crit  !< target surface
        complex*16, intent(in) :: Az(0:nxh+1) !< Vector pot from Helmholtz solution
        real*8, intent(out) :: phipon, epon_x, epon_y, epon_z ! pond. potential and fields

        real*8 ::  yf, zf, xh, Rpon, gamma, sigma, atten, theta
        real*8 :: Azr_0, Azr_1, Azr_2, Azr_3  ! Vector pot. at control points
        real*8 :: epon1, epon2 ! pond force at control points
        real*8 :: ayi, azi ! vec. pot at particle
        real*8 :: epxi, epyi, epzi ! pond field at particle
        real*8 :: r, sigma0
        real*8 :: xa, b1, b2
        real :: f2d ! 2D switch
        real :: pha  ! Temporal phase
        real*8 :: Z_R  ! Rayleigh length
        complex :: yi = (0.,1.)
        integer ::  i1, i2

        ! fast oscillations
        pha = omega*t

        if (sigma_in <0) then
            sigma0=-sigma_in
            f2d = 0.  ! switch off 2d field components
        else
            sigma0=sigma_in
            f2d=1.
        endif

        ! Normalised Rayleigh length k_p Z_R = omega_0/omega_p * (k_p sigma0)^2
        Z_R = omega*sigma0**2


        if (x.le.x_crit) then
            sigma = sigma0*sqrt(1.+abs(x-x_crit)**2/Z_R**2)  ! Vacuum spot size for Gaussian beam
        else
            sigma = sigma0  ! Don't expand spot inside target
        endif

        r = sqrt(y**2+z**2)
        atten = sigma0**2/sigma**2  ! attenuation factor from focal cone

        ! Radial gradients
        theta = pi*r/4./sigma  ! sigma is HWHM of sin^2 laser spot

        if (x.lt.xh_start .or. x.gt.xh_end .or. r >= 2*sigma) then
            yf = 0.
            zf = 0.
            Rpon = 0.

        else if (r < 2*sigma .and. r/= 0.) then
            yf = -pi/4./sigma*y/r*sin(2*theta)
            zf = -pi/4./sigma*z/r*sin(2*theta)
            Rpon = cos(theta)**2

        else
            ! special case on axis
            Rpon = cos(theta)**2
            yf = 0.
            zf = 0.
        endif



        xh = x-xh_start  ! particle coord on HH grid
        xa = xh/dxh     ! reduced coord
        if (xa.ge.1 .and. xa.le.nxh) then
            ! Only compute force for particles inside HH grid
            ! leave one-point buffer at either end for difference dA/dx
            i1 = max(int(xa),1)  ! lower NGP
            i2 = min(i1+1,nxh)       ! upper NGP
            b2=xa-i1  ! linear weights  W_j = 1-|x_i-x_j|
            b1=1.-b2

            ! Derive fpond from Az - need gradient at both reference points.
            ! include temporal phase factor
            Azr_0 = Real(Az(i1-1)*cexp(yi*pha))
            Azr_1 = Real(Az(i1)*cexp(yi*pha))
            Azr_2 = Real(Az(i2)*cexp(yi*pha))
            Azr_3 = Real(Az(i2+1)*cexp(yi*pha))

            ! pond force at reference points i1, i2
            ! - 2nd order difference without gamma factor, as in emfield
            Epon1 = .5*Azr_1*( Azr_2 - Azr_0 )/dxh
            Epon2 = .5*Azr_2*( Azr_3 - Azr_1 )/dxh

            !        ayi = b1*ay(i1) + b2*ay(i2)
            ayi = 0.
            azi = b1*azr_1 + b2*azr_2
            gamma = sqrt(1. + (ux**2 + ayi**2 + azi**2)/unit_c2)

            ! pond field at particle including total particle gamma
            epxi = (b1*epon1 + b2*epon2)/gamma
            epyi = azi**2/gamma
            epzi = azi**2/gamma

            ! potential and fields, correcting for radial dep. and attenuation factor

            phipon = gamma*Rpon*atten

            Epon_x = epxi*Rpon*atten
            Epon_y = epyi*yf*atten
            Epon_z = epzi*zf*atten
        !     if (i1==1) write(*,'(a)') 'x,abs(Az(i1)),  azi,  Epon, gamma '
        !     write(*,'(i5,5(f12.4))') i1,x,abs(Az(i1)),azi,Epxi,gamma
        else
            phipon=0.
            Epon_x=0.
            Epon_y=0.
            Epon_z=0.
        endif

    end subroutine fpond_helm



    ! ==================================================================
    !
    !                        PONDEROMOTIVE FORCE
    !
    !  Compute ponderomotive fields from Helmholtz solution for EM standing wave
    !  c-pol light
    !
    ! ==================================================================

    subroutine fpond_helmc(t,tpulse,sigma_in,vosc,omega, &
    x,y,z,ux,Az,nxh,xh_start,xh_end,dxh,x_crit, &
    epon_x,epon_y,epon_z,phipon)

        !   call fpond_helmc( tlaser, tpulse,sigma,vosc,omega, &
        !                   xd,yd,zd,uxd,Az_helm,nxh,xh_start, xh_end, dxh, focus(1), &
        !           epon_x,epon_y,epon_z,phipon)

        use module_units
        implicit none
        real, intent(in) :: t     !< time
        real, intent(in) :: tpulse !< pulse duration/rise time
        real, intent(in) :: vosc !< quiver strength
        real*8, intent(in) :: sigma_in !< pulse width (1/e)
        real, intent(in) :: omega !< laser frequency
        real*8, intent(in) :: x,y,z !< position to evaluate force;
                                  !< x absolute; y,z relative to laser axis
        real*8, intent(in) :: ux !< forward momentum for gamma factor
        integer, intent(in) :: nxh !< # 1D Helmholtz grid points
        real, intent(in) :: xh_start !< Start point of HH grid
        real, intent(in) :: xh_end !< End point of HH grid
        real, intent(in) :: dxh !< HH grid spacing
        real, intent(in) :: x_crit  !< target surface
        complex*16, intent(in) :: Az(0:nxh+1) !< Vector pot from Helmholtz solution
        real*8, intent(out) :: phipon, epon_x, epon_y, epon_z ! pond. potential and fields

        real*8 ::  yf, zf, xh, Rpon, gamma, sigma, atten, theta
        real*8 :: Azr_0, Azr_1, Azr_2, Azr_3  ! Vector pot. at control points
        real*8 :: epon1, epon2 ! pond force at control points
        real*8 :: ayi, azi ! vec. pot at particle
        real*8 :: epxi, epyi, epzi ! pond field at particle
        real*8 :: r, sigma0
        real*8 :: xa, b1, b2
        real :: f2d ! 2D switch
        real*8 :: Z_R  ! Rayleigh length
        integer ::  i1, i2


        if (sigma_in <0) then
            sigma0=-sigma_in
            f2d = 0.  ! switch off 2d field components
        else
            sigma0=sigma_in
            f2d=1.
        endif

        ! Normalised Rayleigh length k_p Z_R = omega_0/omega_p * (k_p sigma0)^2
        Z_R = omega*sigma0**2


        if (x.le.x_crit) then
            sigma = sigma0*sqrt(1.+abs(x-x_crit)**2/Z_R**2)  ! Vacuum spot size for Gaussian beam
        else
            sigma = sigma0  ! Don't expand spot inside target
        endif

        r = sqrt(y**2+z**2)
        atten = sigma0**2/sigma**2  ! attenuation factor from focal cone

        ! Radial gradients
        theta = pi*r/4./sigma  ! sigma is HWHM of sin^2 laser spot

        if (x.lt.xh_start .or. x.gt.xh_end .or. r >= 2*sigma) then
            yf = 0.
            zf = 0.
            Rpon = 0.

        else if (r < 2*sigma .and. r/= 0.) then
            yf = -pi/4./sigma*y/r*sin(2*theta)
            zf = -pi/4./sigma*z/r*sin(2*theta)
            Rpon = cos(theta)**2

        else
            ! special case on axis
            Rpon = cos(theta)**2
            yf = 0.
            zf = 0.
        endif



        xh = x-xh_start  ! particle coord on HH grid
        xa = xh/dxh     ! reduced coord
        if (xa.ge.1 .and. xa.le.nxh) then
            ! Only compute force for particles inside HH grid
            ! leave one-point buffer at either end for difference dA/dx
            i1 = max(int(xa),1)  ! lower NGP
            i2 = min(i1+1,nxh)       ! upper NGP
            b2=xa-i1  ! linear weights  W_j = 1-|x_i-x_j|
            b1=1.-b2

            ! Derive fpond from Az - need gradient at both reference points.
            Azr_0 = Real(Az(i1-1))
            Azr_1 = Real(Az(i1))
            Azr_2 = Real(Az(i2))
            Azr_3 = Real(Az(i2+1))

            ! pond force at reference points i1, i2
            ! - 2nd order difference without gamma factor, as in emfield
            Epon1 = .5*Azr_1*( Azr_2 - Azr_0 )/dxh
            Epon2 = .5*Azr_2*( Azr_3 - Azr_1 )/dxh

            !        ayi = b1*ay(i1) + b2*ay(i2)
            ayi = 0.
            azi = b1*azr_1 + b2*azr_2
            gamma = sqrt(1. + (ux**2 + ayi**2 + azi**2)/unit_c2)

            ! pond field at particle including total particle gamma
            epxi = (b1*epon1 + b2*epon2)/gamma
            epyi = azi**2/gamma
            epzi = azi**2/gamma

            ! potential and fields, correcting for radial dep. and attenuation factor

            phipon = gamma*Rpon*atten

            Epon_x = epxi*Rpon*atten
            Epon_y = epyi*yf*atten
            Epon_z = epzi*zf*atten
        !     if (i1==1) write(*,'(a)') 'x,abs(Az(i1)),  azi,  Epon, gamma '
        !     write(*,'(i5,5(f12.4))') i1,x,abs(Az(i1)),azi,Epxi,gamma
        else
            phipon=0.
            Epon_x=0.
            Epon_y=0.
            Epon_z=0.
        endif


    end subroutine fpond_helmc




    ! ==================================================================
    !
    !                        PONDEROMOTIVE FORCE
    !
    !  Compute relativistic fpond for standing wave field, linear rise-time
    !
    ! ==================================================================

    subroutine fpond_lin(t,tpulse,sigma_in,vosc,omega,rho_upper,x,y,z,epon_x,epon_y,epon_z,phipon)

        implicit none
        real, intent(in) :: t ! time
        real, intent(in) :: tpulse ! pulse duration/rise time
        real, intent(in) :: vosc ! quiver strength
        real*8, intent(in) :: sigma_in ! pulse width (1/e)
        real, intent(in) :: omega ! laser frequency
        real*8, intent(in) :: x,y,z ! position to evaluate force; x is distance into target from surface (x_c)
        real, intent(in) :: rho_upper
        real*8, intent(out) :: phipon, epon_x, epon_y, epon_z ! pond. potential and fields

        real*8 :: r, xf, yf, zf, Tpon, Rpon, Xpon, Ypon, Zpon, intensity, gamma, sigma, atten, theta, phi, chi
        real :: pi=3.141592654, a02
        real*8 :: rho0_up, wp_r, ls_r, gamma_s, f2d, sigma0
        real*8 :: Xprop
        real*8 :: Z_R  ! Rayleigh length

         ! fast oscillations
        Tpon = (sin(omega*t))**2

        ! intensity envelope
        a02 = vosc**2
        intensity = a02*min(1.,t/tpulse)

        rho0_up = rho_upper*omega**2   ! Effective density of shelf above xc normalised to rho0

        if (sigma_in <0) then
            sigma0=-sigma_in
            f2d = 0.  ! switch off 2d field components
        else
            sigma0=sigma_in
            f2d=1.
        endif

        !  Standing wave in vacuum; evanescent in overdense plasma
        !  Use standard solution of Helmholtz equation for step profile

        ! Normalised Rayleigh length k_p Z_R = omega_0/omega_p * (k_p sigma0)^2
        Z_R = omega*sigma0**2

        gamma_s = sqrt(1.+4*intensity/rho_upper)  ! gamma factor for EM solution at x=xc
        wp_r = sqrt(rho0_up/gamma_s)  ! effective plasma frequency of upper shelf
        ls_r = 1./wp_r   ! rel. skin depth

        phi = atan(-omega/wp_r) ! Interface phase factor given by tan(phi) = -k * l_s
        chi = omega*x + phi  ! Vacuum phase

        if (x.le.0) then
            xf = omega*sin(2*chi)
            Xpon = sin(chi)    ! laser 'Ez'
            sigma = sigma0*sqrt(1.+abs(x)**2/Z_R**2) !  Vacuum spot size for Gaussian beam
            Xprop = 1.
        else
            xf = -2/ls_r*sin(phi)**2*exp(-2*x/ls_r)
            Xpon = sin(phi)*exp(-x/ls_r)    ! laser Ez inside
            sigma = sigma0  ! Don't expand spot inside target
            Xprop = exp(-5*x/ls_r)
        endif

        r = sqrt(y**2+z**2)

        theta = pi*r/4./sigma  ! sigma is HWHM of sin^2 laser spot
        !  Rpon = exp((-yd**2-zd**2)/sigma**2)  ! Gaussian

        if (r < 2*sigma .and. r/= 0.) then
            yf = -pi/4./sigma*y/r*sin(2*theta)
            zf = -pi/4./sigma*z/r*sin(2*theta)
            Rpon = cos(theta)**2
        else if (r >= 2*sigma) then
            yf = 0.
            zf = 0.
            Rpon = 0.
        else
            Rpon = cos(theta)**2
            yf = 0.
            zf = 0.
        endif



        atten = sigma0**2/sigma**2
        !  phipon = 4*Xpon**2*Rpon*atten*intensity  ! intensity, including attenuation factor
        phipon = Xprop*Rpon*atten*intensity  ! intensity for 'demo' display, including attenuation factor

        gamma = sqrt(1.+abs(phipon*Tpon))  ! relativistic phi_pond

        Epon_x = 2*intensity*Tpon*Rpon/gamma*xf*atten
        Epon_y = f2d*2*intensity*Tpon*Xpon**2/gamma*yf*atten
        Epon_z = f2d*2*intensity*Tpon*Xpon**2/gamma*zf*atten
    !  Epon_z = 0.

    end subroutine fpond_lin



    ! ==================================================================
    !
    !                        PONDEROMOTIVE FORCE
    !
    !  Compute relativistic fpond for standing wave field
    !
    ! ==================================================================

    subroutine fpond_sin2(t,tpulse,sigma0,vosc,omega,xd,yd,zd,epon_x,epon_y,epon_z,phipon)


        real, intent(in) :: t ! time
        real, intent(in) :: tpulse ! pulse duration
        real, intent(in) :: vosc ! quiver strength
        real, intent(in) :: sigma0 ! pulse width (1/e)
        real, intent(in) :: omega ! laser frequency
        real*8, intent(in) :: xd,yd,zd ! position to evaluate force; xd is away from target

        real*8, intent(out) :: phipon, epon_x, epon_y, epon_z ! pond. potential and fields

        real*8 :: xf, yf, zf, Tpon, Xpon, Ypon, Zpon, gamma, sigma, atten

        !  linear rise
        Tpon = 2*vosc**2*min(1.,t/tpulse) * (sin(omega*t))**2

        sigma = sigma0*sqrt(1.+abs(xd)**2/4./sigma0**2) ! take Rayleigh length 2*sigma0

        !  Rpon = exp((-yd**2-zd**2)/sigma**2)  ! Gaussian


        !  Standing wave in vacuum; evanescent in overdense plasma
        if (xd.ge.0) then
            xf = sin(2*omega*xd)
            Xpon = cos(omega*xd)**2
        else
            xf = -2*exp(2*xd)
            Xpon = exp(2*xd)
        endif

        if (abs(yd)<2*sigma) then
            yf = sin(pi*yd/2./sigma)
            Ypon = cos(pi*yd/4./sigma)**2
        else
            yf = 0.
            Ypon = 0.
        endif

        if (abs(zd)<2*sigma) then
            zf = sin(pi*zd/2./sigma)
            Zpon = cos(pi*yd/4./sigma)**2
        else
            zf = 0.
            Zpon = 0.
        endif

        atten = sigma0**2/sigma**2
        phipon = Tpon*Xpon*Ypon*Zpon*atten  ! include attenuation factor

        gamma = sqrt(1.+abs(phipon)/2.)

        Epon_x = Tpon*Ypon*Zpon/4/gamma*xf*atten
        Epon_y = Tpon*Xpon/16./gamma*yf
        Epon_z = Tpon*Xpon/16./gamma*zf


    end subroutine fpond_sin2



    ! ==================================================================
    !
    !                        PONDEROMOTIVE FORCE
    !
    !  Compute relativistic fpond for standing wave field
    !
    ! ==================================================================

    subroutine fpond(t,tpulse,sigma_in,vosc,omega,rho_upper,x,y,z,epon_x,epon_y,epon_z,phipon)

        implicit none
        real, intent(in) :: t ! time
        real, intent(in) :: tpulse ! pulse duration/rise time
        real, intent(in) :: vosc ! quiver strength
        real*8, intent(in) :: sigma_in ! pulse width (1/e)
        real, intent(in) :: omega ! laser frequency
        real*8, intent(in) :: x,y,z ! position to evaluate force; x is distance into target from surface (x_c)
        real, intent(in) :: rho_upper
        real*8, intent(out) :: phipon, epon_x, epon_y, epon_z ! pond. potential and fields

        real*8 :: r, xf, yf, zf, Tpon, Rpon, Xpon, intensity, gamma, atten, theta, gamma_s, wp_r, phi, chi, ls_r, sigma0
        real :: pi=3.141592654, a02
        real :: rho0_up, f2d
        real*8 :: Z_R  ! Rayleigh length
        real*8 :: sigma

         ! fast oscillations
        Tpon = (sin(omega*t))**2

        ! intensity envelope
        a02 = vosc**2
        if (t <= 2*tpulse) then
            intensity = a02*max(0.,sin(pi*t/2./tpulse)**2)
        else
            intensity = 0.
        endif

        rho0_up = rho_upper*omega**2   ! Effective density of shelf above xc normalised to rho0

        if (sigma_in <0) then
            sigma0=-sigma_in
            f2d = 0.  ! switch off 2d field components
        else
            sigma0=sigma_in
            f2d=1.
        endif

        !  Standing wave in vacuum; evanescent in overdense plasma
        !  Use standard solution of Helmholtz equation for step profile

        ! Normalised Rayleigh length k_p Z_R = omega_0/omega_p * (k_p sigma0)^2
        Z_R = omega*sigma0**2

        gamma_s = sqrt(1.+4*intensity/rho_upper)  ! gamma factor for EM solution at x=xc
        wp_r = sqrt(rho0_up/gamma_s)  ! effective plasma frequency of upper shelf
        ls_r = 1./wp_r   ! rel. skin depth

        phi = atan(-omega/wp_r) ! Interface phase factor given by tan(phi) = -k * l_s
        chi = omega*x + phi  ! Vacuum phase

        if (x.le.0) then
            xf = omega*sin(2*chi)
            Xpon = sin(chi)    ! laser 'Ez'
            sigma = sigma0*sqrt(1.+abs(x)**2/Z_R**2)  ! Vacuum spot size for Gaussian beam

        else
            xf = -2/ls_r*sin(phi)**2*exp(-2*x/ls_r)
            Xpon = sin(phi)*exp(-x/ls_r)    ! laser Ez inside
            sigma = sigma0  ! Don't expand spot inside target
        endif

        r = sqrt(y**2+z**2)

        theta = pi*r/4./sigma  ! sigma is HWHM of sin^2 laser spot
        !  Rpon = exp((-yd**2-zd**2)/sigma**2)  ! Gaussian

        if (r < 2*sigma .and. r/= 0.) then
            yf = -pi/4./sigma*y/r*sin(2*theta)
            zf = -pi/4./sigma*z/r*sin(2*theta)
            Rpon = cos(theta)**2
        else if (r >= 2*sigma) then
            yf = 0.
            zf = 0.
            Rpon = 0.
        else
            Rpon = cos(theta)**2
            yf = 0.
            zf = 0.
        endif



        atten = sigma0**2/sigma**2
        phipon = 4*Xpon**2*Rpon*atten*intensity  ! intensity, including attenuation factor

        gamma = sqrt(1.+abs(phipon*Tpon))  ! relativistic phi_pond

        Epon_x = 2*intensity*Tpon*Rpon/gamma*xf*atten
        Epon_y = f2d*2*intensity*Tpon*Xpon**2/gamma*yf*atten
        Epon_z = f2d*2*intensity*Tpon*Xpon**2/gamma*zf*atten
    !  Epon_z = 0.

    end subroutine fpond



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> write time history of laser params
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine laser_hist()

        use physvars

        implicit none

        if (my_rank.eq.0) then
            ! Write out to energy.dat file
            if (itime.eq.1)  write(71,'(a)') '! time  a_L fpond xc'
            write (71,'(f12.5,2(1pe12.3))') t_laser, ampl_max, fpon_max
            write ( *,'(f12.5,2(1pe12.3))') t_laser, ampl_max, fpon_max
        endif
    end subroutine laser_hist



    ! ==================================================================
    !
    !                        TRAVELLING PONDEROMOTIVE FORCE
    !
    !  Compute relativistic fpond for propagating laser field
    !
    ! ==================================================================

    subroutine laser_bullet(t,x0,tpulse,sigma,vosc,omega,x,y,z,epon_x,epon_y,epon_z,phipon)

        implicit none
        real, intent(in) :: t ! time laser on
        real, intent(in) :: x0 ! pulse centre
        real, intent(in) :: tpulse ! pulse duration
        real, intent(in) :: vosc ! quiver strength
        real*8, intent(in) :: sigma ! pulse width (1/e)
        real, intent(in) :: omega ! laser frequency
        real*8, intent(in) :: x,y,z ! position to evaluate force; distance from laser centre (x0,0,0)

        real*8, intent(out) :: phipon, epon_x, epon_y, epon_z ! pond. potential and fields

        real*8 :: xf, yf, zf, Rpon, Xpon, Ypon, Zpon, gamma, phi
        real :: pi=3.141592654, a02

        !  sin^2 time envelope

        phi = pi*x/tpulse/2.
        Rpon = exp((-y**2-z**2)/sigma**2)  ! Gaussian radial profile

        if (x.ge.-tpulse .and. x.lt.tpulse) then
            xf = -pi/4.*sin(2*phi)
            Xpon = cos(phi)**2
        else
            xf = 0.
            Xpon = 0.
        endif

        yf = -2*y/sigma**2
        zf = -2*z/sigma**2

        a02 = vosc**2
        phipon = a02*Xpon*Rpon  ! intensity

        gamma = sqrt(1.+abs(phipon)/2.)  ! relativistic phi_pond

        Epon_x = a02*Rpon/gamma*xf
        Epon_y = a02*Xpon*Rpon/gamma*yf
        Epon_z = a02*Xpon*Rpon/gamma*zf

    end subroutine laser_bullet



    !>     ==================================
    !>
    !>     Helmholtz solver for electromagnetic fields
    !>
    !>     Based on SPLIM project helmholtz.f90
    !>     ==================================


    subroutine em_helmholtz(me,itime,n,dx,theta,a0,w0,rhoe,Az)
        use module_io
        implicit none
        integer, intent(in) :: n !< number of mesh points for 1D Helmholtz fields
        integer, intent(in) :: itime  !< current timestep
        integer, intent(in) :: me  !< current rank
        real, intent(in) :: theta  !< angle of incidence
        real*8, intent(in) ::  a0    !< laser amplitude vosc/c(t)
        real, intent(in) ::  w0    !< laser frequency (w/wp)
        real, intent(in) ::  dx    !< mesh spacing (c/wp)
        real*4, intent(in) :: rhoe(0:n+1)  !< cycle-averaged electron density
        !  real*8, intent(out) :: Ezr(0:n+1), Byr(0:n+1), Azr(0:n+1), epond(0:n+1)
        complex*16, dimension(n) :: alpha,beta,gamma,y !< trisolver work arrays
        complex*16, dimension(0:n+1) :: Az,Ao    !< Vector potential
        complex*16, dimension(n) :: eps  !< permittivity
        real*8 :: rgam(0:n+1)  !<  relativistic gamma
        real*8 :: err(0:n+1)   !<  error check
        complex :: yi, carg  !< complex args
        integer :: i,j,n_iter
        real :: pi, s2th, cth !< constants
        real*8 :: ncrit, errmax  !< critical density

        real*8 :: g0, nu_eff
        integer :: iplas, itav

        itav=1
        n_iter = 3
        ncrit = w0**2  ! critical density in norm. units
        nu_eff=0.2  ! Effective collision rate
        pi = asin(1.0)*2
        yi = (0.,1.)
        err=0.
        s2th=sin(pi/180*theta)**2
        cth = cos(pi/180*theta)


        ! Use gamma from previous step to speed up iteration
        ao=az
        rgam(0)=1.
        rgam(n+1) = 1.
        do i=1,n
            rgam(i) = sqrt(1 + 0.5*abs(ao(i))**2)
            az(i)=(0.,0.)
        end do


        Az(0)=(0.,0.)
        Az(n+1)=(0.,0.)

        if (a0==0) return

        do j=1,n_iter
            do i=1,n
                !  coefficients as for s-pol light
                ! rhoe normalized to nc defined on own 1D grid along laser axis
                eps(i) = 1.-rhoe(i)/ncrit/rgam(i)/(1.+yi*nu_eff)
                y(i)=(0.,0.)
                alpha(i)=1
                beta(i)=-2 + w0**2*dx**2*(eps(i)-s2th)
                gamma(i)=1
            end do

            !  BCs - note additional factor w0=k0 for phase factors
            y(1) = 2*yi*a0*sin(w0*dx*cth)
            carg = yi*w0*dx*cth
            beta(1) = beta(1) + cexp(carg)

            call trisolve(alpha,beta,gamma,y,Az(1:n),n-1,n)

            ! relativistic factor - average old and new potentials
            errmax=0.
            do i=1,n
                rgam(i) = sqrt(1 + 0.5*abs(az(i))**2)
                err(i) = sqrt(abs(az(i)**2-ao(i)**2)/4/a0**2)
            end do

            ! BCs for next iterate (Laplacian of gamma)
            rgam(0) = 2*rgam(1) - rgam(2)
            rgam(n+1) = rgam(n)

            errmax = maxval(err(1:n))
            Ao = Az  ! Store previous iterate

        end do

        iplas = n/2
        if (me==0) write(*,'(i6,2f12.3)') itime,rhoe(iplas),abs(az(iplas))
        if (itime .eq. itav .and. me==0) then
            write (*,'(a20,i2,a10,f12.5,a10,f12.3)') 'Iterate ',j,' error=',errmax,' amplitude',a0
            g0 = sqrt(1+a0**2/2)
            open (file_tempfile_1,file='a_error.dat')
            write(file_tempfile_1,'(a)') '! x, rho, eps, az/a0, gam/g0, err'
            write(file_tempfile_1,'(6(1pe12.3))') (dx*i,rhoe(i),eps(i),abs(az(i))/a0,rgam(i)/g0,err(i),i=1,n)
            close(file_tempfile_1)

        endif


        ! Bcs
        Az(0) = 2*Az(1) - Az(2)  ! gradient continuous
        Az(n+1) = Az(n) ! zero in solid


    end subroutine em_helmholtz



    ! -------------------------------------------------
    !
    !   Tridiagonal matrix solver.
    !
    !   Solves equation system:
    !
    !        alpha_i x_(i-1) + beta_i x_i + gamma_i x_(i+1) = y_i
    !
    ! ------------------------------------------------

    subroutine trisolve(alpha,beta,gamma,y,x,n,nmax)

        implicit none
        integer, intent(in) :: n,nmax
        complex*16, dimension(nmax) :: alpha, beta, gamma,y,x, q
        integer :: i

        q(1) = beta(1)
        x(1) = y(1)/q(1)

        !  forward elimination

        do i = 2,n
            q(i) = beta(i) - alpha(i)*gamma(i-1)/q(i-1)
            x(i) = (y(i) - alpha(i)*x(i-1))/q(i)
        end do

        !  back substitution

        do i=n-1,1,-1
            x(i) = x(i) - gamma(i)/q(i)*x(i+1)
        end do
    end subroutine trisolve



    !     ==================================
    !
    !     Electromagnetic fields
    !
    !     - plane wave with finite rise-time and Gaussian spot
    !  returns laser fields at particle position x,y,z
    !
    !
    !     ==================================
    subroutine emplane(t,tpulse,sigma0,a0,w0,x,y,z,ez,by,bx,az,phipon)

        implicit none
        real*8, intent(in) :: t ! time
        real*8, intent(in) :: tpulse ! pulse duration or rise-time
        real*8, intent(in) :: a0 ! quiver strength
        real*8, intent(in) :: sigma0 ! pulse width (1/e)
        real*8, intent(in) :: w0 ! laser frequency
        real*8, intent(in) :: x,y,z ! position to evaluate force; x is distance into target from surface (x_c)

        real*8, intent(out) :: phipon, ez, by, bx, az ! pond. potential and fields

        real*8 :: tenv, theta_r, gamma0, k0
        real :: pi=3.141592654
        real*8 :: r2, phase, earg

        !     linear rise

        gamma0 = sqrt(1 + a0**2/2)  ! Sarachik & Schappert gamma

        k0=w0
        if (t <= 2*tpulse) then
            tenv = max(0._8,sin(pi*t/2./tpulse)**2)
        else
            tenv= 0.
        endif

        !  tenv = 1.
        phase = w0*t - k0*x
        r2 = y**2+z**2


        if (sigma0.gt.0) then
            !  pulse envelope - Gaussian, 1/e width = sigma
            earg = min(20._8,r2/2/sigma0**2)
            theta_r = exp(-earg)
        else
            !  constant amplitude
            theta_r = 1.
        endif

        ! reconstruct EM fields (s-pol)

        Az = a0*tenv*theta_r*sin(phase)
        Ez = -w0*a0*tenv*theta_r*cos(phase)      ! Ez = -dAz/dt
        By = k0*a0*tenv*theta_r*cos(phase)       ! By = -dAz/dx
        Bx = -y/sigma0**2*a0*tenv*theta_r*sin(phase) ! Bx = dAz/dy


        phipon = Az**2

    end subroutine emplane



    !     ==================================
    !
    !     Electromagnetic fields
    !
    !     - ponderomotive laser model for step-profile
    !  returns laser fields at particle position x,y,z relative to critical surface
    !
    !
    !     ==================================
    subroutine empond(t,tpulse,sigma0,a0,w0,x,y,z,ez,by,bx,az,phipon)

        implicit none
        real*8, intent(in) :: t ! time
        real*8, intent(in) :: tpulse ! pulse duration or rise-time
        real*8, intent(in) :: a0 ! quiver strength
        real*8, intent(in) :: sigma0 ! pulse width (1/e)
        real*8, intent(in) :: w0 ! laser frequency
        real*8, intent(in) :: x,y,z ! position to evaluate force; x is distance into target from surface (x_c)

        real*8, intent(out) :: phipon, ez, by, bx, az ! pond. potential and fields

        real*8 :: tenv, Rpon, dRpon, phase, k0, nonc, gamma_c, wp_r, lskin, phi
        real :: pi=3.141592654
        real*8 :: Z_R  ! Rayleigh length
        real*8 :: r, f_helm, g_helm, chi, atten, theta
        real*8 :: sigma

        !     linear rise

        phase = w0*t
        k0=w0
        nonc = 1./w0**2 ! density normalised to nc
        Tenv = min(1._8,t/tpulse) ! time envelope


        !     Standing wave in vacuum; evanescent in overdense plasma
        !    Use standard solution of Helmholtz equation for step profile

        gamma_c = sqrt(1.+4*a0**2/nonc) ! gamma at surface
        wp_r = 1./sqrt(gamma_c)
        lskin = 1./wp_r   ! Rel. skin depth in EM units

        ! Normalised Rayleigh length k_p Z_R = omega_0/omega_p * (k_p sigma0)^2
        Z_R = w0*sigma0**2

        !   Phase factor given by tan(phi) = -k0 * l_s = k0 * c/wp
        phi = atan(-1./wp_r)
        r = sqrt(y**2+z**2)


        if (x <= 0) then
            !     vacuum solution
            chi = k0*x + phi      ! Vacuum phase
            f_helm = sin(chi)     ! laser 'Ez'
            g_helm = cos(chi)     ! laser 'By'
            sigma = sigma0*sqrt(1.+abs(x)**2/Z_R**2)  ! Vacuum spot size for Gaussian beam

        else
            !   evanescent wave - need Sudan solution here
            f_helm = sin(phi)*exp(-x/lskin)        ! laser Ez inside
            g_helm = cos(phi)*exp(-x/lskin) ! laser By inside
            sigma = sigma0  ! Don't expand spot inside target
        endif


        !  Radial envelope

        theta = pi*r/4./sigma
        atten = sigma0/sigma

        if (r <= 2.*sigma .and. r/=0. ) then
            Rpon = atten*cos(theta)
            dRpon = -pi*y/4./sigma/r*atten*sin(theta)

        else
            Rpon = 0.
            dRpon = 0.
        endif

        ! reconstruct EM fields (s-pol)

        az = 2*a0     *           f_helm * sin(phase) * Rpon * tenv

        ez = -2*w0*a0 *           f_helm * cos(phase) * Rpon * tenv ! Ez = -dAz/dt
        by = -2*k0*a0 *           g_helm * sin(phase) * Rpon * tenv ! By = -dAz/dx
        bx =  2*a0 *              f_helm * sin(phase) * dRpon* tenv ! Bx = dAz/dy

        phipon = (2*a0*f_helm*tenv*Rpon)**2

    end subroutine empond




    ! ==================================================================
    !
    !              EMOBLIQ
    !
    !  Compute relativistic fields for oblique-incidence standing wave field
    !       s-light
    !
    ! ==================================================================

    subroutine emobliq(t,tpulse,sigma_in,vosc,omega,theta,rho_upper,x,y,z,epon_x,epon_y,epon_z,phipon,Ez,Bx,By)

        implicit none
        real, intent(in) :: t ! time
        real, intent(in) :: tpulse ! pulse duration/rise time
        real, intent(in) :: vosc ! quiver strength
        real*8, intent(in) :: sigma_in ! pulse width (1/e)
        real, intent(in) :: omega ! laser frequency
        real, intent(in) :: theta ! angle of incidence

        real*8, intent(in) :: x,y,z ! position to evaluate force; x is distance into target from surface (x_c)
        real, intent(in) :: rho_upper
        real*8, intent(out) :: phipon, epon_x, epon_y, epon_z ! pond. potential and fields
        real*8, intent(out) :: Ez, By, Bx  ! laser fields

        real*8 :: xf, yf, zf, Rpon, Xpon, intensity, gamma, atten, alpha, eps, dXpon, sigma0
        real :: pi=3.141592654, a02
        real :: rho0_up, f2d
        real :: kx, ky, k0, thetar
        real*8 :: r, tphase, chi_s, Tpon, gamma_s, wp_r, ls_r, phi_s
        real*8 :: sigma

        ! phase factors (normalised to nominal wp)
        thetar=pi*theta/180.
        k0 = omega
        kx = k0*cos(thetar)
        ky = k0*sin(thetar)
        tphase = omega*t - ky*y  ! time phase

        Tpon = sin(tphase)**2

        ! intensity envelope
        a02 = vosc**2
        if (t <= 2*tpulse) then
            !    intensity = a02*max(0.,sin(pi*t/2./tpulse)**2)
            intensity = a02
        else
            intensity = 0.
        endif

        rho0_up = rho_upper*omega**2   ! Effective density of shelf above xc normalised to rho0

        if (sigma_in <0) then
            sigma0=-sigma_in
            f2d = 0.  ! switch off 2d field components
        else
            sigma0=sigma_in
            f2d=1.
        endif


        !  Standing wave in vacuum; evanescent in overdense plasma
        !  Use standard solution of Helmholtz equation for step profile

        gamma_s = sqrt(1.+4*intensity/rho_upper)  ! gamma factor for EM solution at x=xc
        wp_r = sqrt(rho0_up/gamma_s)  ! effective plasma frequency of upper shelf
        !  wp_r=1.
        ls_r = 1./wp_r   ! rel. skin depth

        phi_s = atan(-kx*ls_r) ! Interface phase factor given by tan(phi) = -kx * l_r
        chi_s = kx*x + phi_s  ! Vacuum phase (s-pol)

        if (x.le.0) then
            xf = omega*sin(2*chi_s)
            Xpon = sin(chi_s)    ! laser phase
            dXpon = kx*cos(chi_s) ! derivative
            sigma = sigma0*sqrt(1.+abs(x)**2/10./sigma0**2) ! take Rayleigh length 4*sigma0
            eps = 1.  ! refractive index
        else
            xf = -2/ls_r*sin(phi_s)**2*exp(-2*x/ls_r)
            Xpon = sin(phi_s)*exp(-x/ls_r)    ! laser phase inside
            dXpon = -sin(phi_s)/ls_r*exp(-x/ls_r) ! derivative
            sigma = sigma0  ! Don't expand spot inside target
            eps = 1.-(wp_r/omega)**2  ! refractive index inside target
        endif

        r = sqrt(y**2+z**2)

        alpha = pi*r/4./sigma  ! sigma is HWHM of sin^2 laser spot
        !  Rpon = exp((-yd**2-zd**2)/sigma**2)  ! Gaussian

        if (r < 2*sigma .and. r/= 0.) then
            yf = -pi/4./sigma*y/r*sin(2*alpha)
            zf = -pi/4./sigma*z/r*sin(2*alpha)
            Rpon = cos(alpha)**2
        else if (r >= 2*sigma) then
            yf = 0.
            zf = 0.
            Rpon = 0.
        else
            Rpon = cos(alpha)**2
            yf = 0.
            zf = 0.
        endif



        atten = sigma0**2/sigma**2
        phipon = 4*Xpon**2*Rpon*atten*intensity*Tpon  ! intensity, including attenuation factor

        gamma = sqrt(1.+abs(phipon*Tpon))  ! relativistic phi_pond

        Ez = 2*vosc*cos(tphase)*Xpon*sqrt(Rpon*atten)
        Bx = 2*ky*vosc/omega*cos(tphase)*Xpon*sqrt(Rpon*atten)
        By = 2*vosc/omega*sin(tphase)*dXpon*sqrt(Rpon*atten)

        Epon_x = 2*intensity*Tpon*Rpon/gamma*xf*atten
        Epon_y = f2d*2*intensity*Tpon*Xpon**2/gamma*yf*atten
        Epon_z = f2d*2*intensity*Tpon*Xpon**2/gamma*zf*atten
    !  Epon_z = 0.

    end subroutine emobliq



    !     ==================================
    !
    !     Electromagnetic fields
    !
    !     - plane wave with finite rise-time and Gaussian spot
    !  returns laser fields at particle position x,y,z
    !
    !
    !     ==================================


    subroutine emplane_lin(t,tpulse,sigma0,a0,w0,x,y,z,ez,by,bx,az,phipon)

        implicit none
        real, intent(in) :: t ! time
        real, intent(in) :: tpulse ! pulse duration or rise-time
        real, intent(in) :: a0 ! quiver strength
        real*8, intent(in) :: sigma0 ! pulse width (1/e)
        real, intent(in) :: w0 ! laser frequency
        real*8, intent(in) :: x,y,z ! position to evaluate force; x is distance into target from surface (x_c)

        real*8, intent(out) :: phipon, ez, by, bx, az ! pond. potential and fields

        real :: tenv, gamma0
        real*8 :: earg
        real :: k0
        real*8 :: phase, r2, theta_r

        !     linear rise

        gamma0 = sqrt(1 + a0**2/2)  ! Sarachik & Schappert gamma

        k0=w0
        tenv = min(1.,t/tpulse)
        !  tenv = 1.
        phase = w0*t - k0*x
        r2 = y**2+z**2


        if (sigma0.gt.0) then
            !  pulse envelope - Gaussian, 1/e width = sigma
            earg = min(20._8,r2/2/sigma0**2)
            theta_r = exp(-earg)
        else
            !  constant amplitude
            theta_r = 1.
        endif

        ! reconstruct EM fields (s-pol)

        Az = a0*tenv*theta_r*sin(phase)
        Ez = -w0*a0*tenv*theta_r*cos(phase)      ! Ez = -dAz/dt
        By = k0*a0*tenv*theta_r*cos(phase)       ! By = -dAz/dx
        Bx = -y/sigma0**2*a0*tenv*theta_r*sin(phase) ! Bx = dAz/dy


        phipon = Az**2

    end subroutine emplane_lin





















end module module_laser
