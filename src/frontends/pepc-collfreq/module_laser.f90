! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2014 Juelich Supercomputing Centre, 
!                         Forschungszentrum Juelich GmbH,
!                         Germany
! 
! PEPC is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! PEPC is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
! 
! You should have received a copy of the GNU Lesser General Public License
! along with PEPC.  If not, see <http://www.gnu.org/licenses/>.
!

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
    real*8, public :: omega_hz  =  0.0    !< frequency omega in Hz
    real*8, public :: lambda    =  0.0    !< laser wavelength
    real*8, public :: lambda_nm =  0.0    !< laser wavelength in nm
    real*8, public :: rhocrit_nm3 = 0.      !< critical electron density in electrons per nm^3
    real*8, public :: I0_Wpercm2 = 0.   !< initial intensity in W/cm^2
    real*8, public :: E0 = 0.   !< laser field strength amplitude
    real*8, public :: vosc      =  0.1    !< pump strength
    real*8, public :: vosc_vte  =  0.0    !< relative pump strength

    real*8, public :: I_laser           !< Laser intensity (= amplitude**2), updated by each call to laser_update()
    real*8, public :: E_laser           !< Laser amplitude (= intensity**1/2), updated by each call to laser_update()
    real*8, public :: focus(3)  = [0., 0., 0.] !< centre of focal spot
    real*8, public :: t_pulse      = 0.     !< pulse duration (in simulation time units)
    real*8, public :: t_pulse_fs   = 10.     !< pulse duration (in fs)
    real*8, public :: sigma   =  1.0    !< 1/e pulse width (c/omega_p)
    real, public :: theta_inc =  0.     !< angle of incidence
    real, public :: rho_track =  1.5    !< tracking density for x_crit (/nc)
    real*8, public :: t_laser   = 0.      !< run time after laser switched on (in simulation time units)
    real*8, public :: elaser        !< deposited laser energy
    real, public :: intensity     !< normalised intensity = 0.5*vosc^2*omega^2
    real*8, public :: fpon_max !< max amplitudes
    real*8, public :: navcycle     ! # timesteps in a laser cycle

    character*7, public :: beam_configs(0:9)=(/ &
    'off    ','beam   ','i-beam ','laser-u','ES pond','LWFA   ', &
    'EMplane','EM pond',' dust  ','       ' /)

    ! declarations from pepcb physvars module
    real*8, public :: x_crit         !< critical surface


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    public laser_setup
    public laser_update
    public force_laser
    public force_laser_at


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
    subroutine force_laser(particles)
      use module_pepc_types
      use physvars
      implicit none
      include 'mpif.h'
      type(t_particle), intent(inout) :: particles(:)
      real*8 :: E_pon(3), B_em(3), Phipon

      integer(kind_particle) :: p

      ! Include force from laser/external field on electrons
      do p = 1,size(particles,kind=kind(p))
          call force_laser_at(particles(p)%x, E_pon, B_em, Phipon)

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


    subroutine force_laser_at(pos, E_pon, B_em, Phipon)
        implicit none
        real*8, intent(out) :: E_pon(3), B_em(3), Phipon
        real*8, intent(in) :: pos(3)

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

            case default  ! no laser
                E_pon = [ 0., 0., 0. ]
                B_em  = [ 0., 0., 0. ]

        end select

    end subroutine force_laser_at

end module module_laser
