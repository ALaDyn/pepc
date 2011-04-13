!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>  Encapsulates anything that is concenred with laser setup and laser-particle interaction
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_laser
      implicit none
      save
      private

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer, public :: beam_config_in = 0 !< Particle or laser beam switch including variations
      integer, public :: beam_config = 0 !< Reduced switch for particle or laser beam

      real, public :: omega     =  0.1    !< frequency
      real, public :: omega_wpl =  0.0    !< frequency omega in wpl_e
      real, public :: omega_hz  =  0.0    !< frequency omega in Hz
      real, public :: lambda    =  1.0    !< laser wavelength
      real, public :: lambda_nm =  1.0    !< laser wavelength in nm
      real*8, public :: rhocrit_nm3 = 0.      !< critical electron density in electrons per nm^3
      real*8, public :: I0_Wpercm2 = 0.   !< initial intensity in W/cm^2
      real*8, public :: E0 = 0.   !< laser field strength amplitude
      real*8, public :: vosc      =  0.1    !< pump strength
      real*8, public :: wpl_e, wpl_i !< electron and ion plasma frequency

      real*8, public :: I_laser           !< Laser intensity (= amplitude**2), updated by each call to laser()
      real, public :: x_offset  =  0.     !< coordinate offset
      real, public :: z_offset  =  0.     !< coordinate offset
      real, public :: focus(3)  = [0., 0., 0.] !< centre of focal spot
      real, public :: tpulse    = 10.     !< pulse duration (in units of 1/omega_p)
      real*8, public :: sigma   =  1.0    !< 1/e pulse width (c/omega_p)
      real, public :: theta_inc =  0.     !< angle of incidence
      real, public :: rho_track =  1.5    !< tracking density for x_crit (/nc)
      real, public :: tlaser    = 0.      !< run time after laser switched on (1/omega_p)
      real*8, public :: elaser        !< deposited laser energy
      real, public :: propag_laser  !< distance travelled by laser after rezoning
      real, public :: intensity     !< normalised intensity = 0.5*vosc^2*omega^2
      real, public :: window_min    !< start of wakefield plasma
      real, public :: rezone_frac=0.75     !< Fraction of box to cross before rezoning switched on
      real, public :: glue_radius=2.e6 !< multiple of box size to catch escaping particles at
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
      real, public :: x_crit         !< critical surface
      real, public :: theta_beam = 0. !< beam elevation angle


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      public setup_laser
      public laser
      public force_laser
      public laser_hist


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !real, parameter :: pi=3.141592654


      contains

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Setup all local module variables and data that depends on them
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine setup_laser()
          use physvars
          implicit none

          beam_config=mod(beam_config_in,10)  ! derived config from s,p variations
          navcycle   = 2.*pi/dt/omega  ! # timesteps in a laser cycle

        end subroutine setup_laser

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!>
		!> Laser propagation according to beam_config
		!>
		!> Calculates time-dependend laser intensity and corrects laser focus position
		!>
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		subroutine laser()

		  use physvars
		  use module_units
		  implicit none

		  integer :: ifile

          if (my_rank==0) then
            write (*,'(/"-- LASER --"/a20,i8," ",a)') ' Beam config: ',beam_config,beam_configs(beam_config)
            write (*,'(a20,i8)') ' type: ',beam_config_in
          end if

          if (beam_config_in == 0) return

          tlaser = tlaser + dt

		  !  Laser pulse envelope
		  !  =====================
		  laser_model: select case(beam_config_in)

		     case(4,34,44,64,6)     ! sin^2 standing wave
		        if (tlaser<2*tpulse)  then
                  I_laser = vosc**2*max(0.,sin(pi*tlaser/2./tpulse)**2)
		        else
		          I_laser=0.
		        endif

		     case(14,94,54,16)  ! standing wave, linear rise
                  I_laser = vosc**2*min(1.,tlaser/tpulse)

		     case(3)        ! constant
		          I_laser = vosc**2

		     case default
		          I_laser = 0.

		  end select laser_model


          ! old fpond model
		  if (itime>0) focus(1) = x_crit  ! laser tracks n_c

          !    Deposited laser energy
          !    ==================================
          laser_energy: select case(beam_config)
          case(4)
             elaser = 3./8.*omega**2*sigma**2*vosc**2*tlaser
          case(5)
             elaser = 3./8.*omega**2*sigma**2*vosc**2*tpulse
          case default
             elaser = 0.
          end select laser_energy


          !    Laser parameter output
          !    ==================================
	      if (my_rank==0) then
	        do ifile = 6,15,9
	           if (beam_config.ne.0)  then
	              write(ifile,'(/"-- LASER --"/(a20,f8.2,a4,f12.6,a4)/6(a20,f9.3/))') &
	                    'tlaser =',          tlaser,'   (',tlaser*unit_t0_in_fs,' fs)' &
	                   ,'amplitude =',       sqrt(I_laser) &
	                   ,'x_crit =',          x_crit &
	                   ,'spot size =',       sigma &
	                   ,'theta =',           theta_beam &
                       ,'elapsed =',         tlaser &
                       ,'steps per cycle =', navcycle
	           else if (beam_config==5 .or. beam_config==6) then
	               write(ifile,'(5(a,f8.2/))') &
	                    'Laser amplitude =',sqrt(I_laser) &
	                  , 'Pulse length',tpulse &
	                  , 'Pulse width', sigma &
	                  , 'Focal position', focus(1)
	           endif
	        end do
	      endif

		end subroutine laser



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

		  integer, intent(in) :: p_start,p_finish  ! min, max particle nos.
		  integer :: p
		  real*8 :: dc(3)  ! positions relative to centre of plasma centre
		  real*8 :: df(3)  ! position relative to laser focus
		  real*8 :: rt
		  real*8 :: E_pon(3)

		  dxh = (xh_end-xh_start)/nxh  ! HH grid spacing

		  if (itime>0 .and. beam_config==4) focus(1) = x_crit  ! laser tracks n_c

		  ! Include force from laser/external field on electrons - ES scheme
	      do p = p_start, p_finish

             df = [ x(p), y(p), z(p) ] - focus
             dc = [ x(p), y(p), z(p) ] - plasma_centre
             rt = sqrt(dc(1)**2+dc(2)**2)

		     laser_model: select case(beam_config_in)

		        case(03)  ! Uniform sinusoid in x (s-pol)
		           E_pon = [ m(p)*vosc*omega*sin(omega*tlaser), 0._8, 0._8 ]

		        case(23)  ! Uniform sinusoid in z (s-pol)
                   E_pon = [ 0._8, 0._8, vosc*omega*sin(omega*tlaser) ]

		        case(04)  ! standing wave fpond, sin2 pulse
                   E_pon = sqrt(I_Laser) * [ 0._8, 0._8, vosc*omega*sin(omega*tlaser) ]

		        case default  ! no laser
                   E_pon = [ 0., 0., 0. ]

		        end select laser_model

                fpon_max = max(fpon_max, abs(E_pon(1)))
  		        ! Add external fields to new particle field
		        Ex(p) = Ex(p) + E_pon(1)
		        Ey(p) = Ey(p) + E_pon(2)
		        Ez(p) = Ez(p) + E_pon(3)
		  end do

		end subroutine force_laser


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
		    write (71,'(f12.5,2(1pe12.3))') tlaser, ampl_max, fpon_max
		  write (*,'(f12.5,2(1pe12.3))') tlaser, ampl_max, fpon_max
		endif
		end subroutine laser_hist




end module module_laser
