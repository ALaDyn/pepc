! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2016 Juelich Supercomputing Centre, 
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
!>
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_param_dump
      implicit none
      private

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      public PrintPhysicalParameters
      public PrintPeriodicityParameters
      public PrintLaserParameters
      public PrintEnergies
      public WriteParameter
      public WriteTopline
      public WriteHeader

      interface WriteTopline
        module procedure WriteTopline1
        module procedure WriteTopline2
      end interface

      interface WriteHeader
        module procedure WriteHeader1
        module procedure WriteHeader2
      end interface

      interface WriteParameter
        module procedure WriteParameterReal
        module procedure WriteParameterInt
        module procedure WriteParameterLong
        module procedure WriteParameterRealSingle
        module procedure WriteParameterIntSingle
        module procedure WriteParameterLongSingle
        module procedure WriteParameterCharSingle
        module procedure WriteParameterLogicalSingle
        module procedure WriteParameterRealCoord
        module procedure WriteParameterLogicalCoord
      end interface

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  contains

    subroutine  WriteTopline1(ifile,name)
      implicit none
      integer, intent(in) :: ifile
      character(*), intent(in) :: name
      character(62) :: topline
      integer :: lname

      write(topline,'(62("-"))')
      lname  = len(name)
      if (lname > 0) then
        topline(30:30)              = " "
        topline(31:31+lname-1)      = name
        topline(31+lname:31+lname)  = " "
      endif

      write(ifile,'(a)') topline
    end subroutine

    subroutine  WriteTopline2(ifile,name1,name2)
      implicit none
      integer, intent(in) :: ifile
      character(*), intent(in) :: name1, name2
      character(92) :: topline
      integer :: lname

      write(topline,'(92("-"))')

      lname  = len(name1)
      if (lname > 0) then
        topline(30:30)              = " "
        topline(31:31+lname-1)      = name1
        topline(31+lname:31+lname)  = " "
      endif

      lname  = len(name2)
      if (lname > 0) then
        topline(60:60)              = " "
        topline(61:61+lname-1)      = name2
        topline(61+lname:61+lname)  = " "
      endif

      write(ifile,'(a)') topline
    end subroutine



    subroutine  WriteHeader1(ifile,name)
      implicit none
      integer, intent(in) :: ifile
      character(*), intent(in) :: name

      write(ifile,'(a1,a31,2(a20,a10))') "| ", "|", name, "|"

    end subroutine

    subroutine  WriteHeader2(ifile,name1,name2)
      implicit none
      integer, intent(in) :: ifile
      character(*), intent(in) :: name1, name2

      write(ifile,'(a1,a31,2(a20,a10))') "| ", "|", name1, "|", name2,"|"

    end subroutine



    subroutine WriteParameterReal(ifile,name, eval, ival)
      implicit none
      integer, intent(in) :: ifile
      real*8, intent(in) :: eval, ival
      character(*), intent(in) :: name

      write(ifile,'("| ", a28," |", 2(g20.6, "         |"))') name, eval, ival

    end subroutine WriteParameterReal

    subroutine WriteParameterInt(ifile,name, eval, ival)
      implicit none
      integer, intent(in) :: ifile
      integer, intent(in) :: eval, ival
      character(*), intent(in) :: name

      write(ifile,'("| ",a28," |", 2(i20, "         |"))') name, eval, ival

    end subroutine WriteParameterInt

    subroutine WriteParameterLong(ifile,name, eval, ival)
      implicit none
      integer, intent(in) :: ifile
      integer*8, intent(in) :: eval, ival
      character(*), intent(in) :: name

      write(ifile,'("| ",a28," |", 2(i20, "         |"))') name, eval, ival

    end subroutine WriteParameterLong

    subroutine WriteParameterRealSingle(ifile,name, val)
      implicit none
      integer, intent(in) :: ifile
      real*8, intent(in) :: val
      character(*), intent(in) :: name

      write(ifile,'("| ", a28," |", 1(g20.6, "         |"))') name, val

    end subroutine WriteParameterRealSingle

    subroutine WriteParameterIntSingle(ifile,name, val)
      implicit none
      integer, intent(in) :: ifile
      integer, intent(in) :: val
      character(*), intent(in) :: name

      write(ifile,'("| ",a28," |", 1(i20, "         |"))') name, val

    end subroutine WriteParameterIntSingle

    subroutine WriteParameterLongSingle(ifile,name, val)
      implicit none
      integer, intent(in) :: ifile
      integer*8, intent(in) :: val
      character(*), intent(in) :: name

      write(ifile,'("| ",a28," |", 1(i20, "         |"))') name, val

    end subroutine WriteParameterLongSingle

    subroutine WriteParameterCharSingle(ifile,name, val)
      implicit none
      integer, intent(in) :: ifile
      character(*), intent(in) :: val
      character(*), intent(in) :: name

      write(ifile,'("| ",a28," | ", 1(a27, " |"))') name, val

    end subroutine WriteParameterCharSingle

    subroutine WriteParameterRealCoord(ifile,name, val, showlength)
      implicit none
      integer, intent(in) :: ifile
      real*8, intent(in) :: val(3)
      logical, intent(in), optional :: showlength
      character(*), intent(in) :: name

      logical :: showlen

      showlen = .true.
      if (present(showlength)) showlen = showlength

                   write(ifile,'("| ", a28," |       [ ", 3(g13.6,x), "]       |")') " "//name//" ", val(1:3)
      if (showlen) write(ifile,'("| ", a28," |         ", 28x,g13.6,x, "        |")') "|"//name//"|", sqrt(dot_product(val, val))

    end subroutine WriteParameterRealCoord

    subroutine WriteParameterLogicalCoord(ifile,name, val)
      implicit none
      integer, intent(in) :: ifile
      logical, intent(in) :: val(3)
      character(*), intent(in) :: name
      character :: text(3)

      where (val)
        text = 'T'
      elsewhere
        text = 'F'
      end where

      write(ifile,'("| ", a28," |       [ ", 3(12x,a1,x), "]       |")') " "//name//" ", text(1:3)

    end subroutine WriteParameterLogicalCoord


    subroutine WriteParameterLogicalSingle(ifile,name, val)
      implicit none
      integer, intent(in) :: ifile
      logical, intent(in) :: val
      character(*), intent(in) :: name
      character :: text

      if (val) then
        text = 'T'
      else
        text = 'F'
      endif

      write(ifile,'("| ",a28," | ", 1(a27, " |"))') name, text

    end subroutine WriteParameterLogicalSingle



    subroutine PrintPhysicalParameters(ifile)
      use physvars
      use module_laser
      use module_units
      use module_interaction_specific
      implicit none
      integer, intent(in) :: ifile
! xl, yl, zl, plasma_centre
!  call laser_setup()

      call WriteTopline(  ifile, "PARAMETERS", "")
      call WriteHeader(   ifile, "electrons", "ions")
      call WriteTopline(  ifile, "","")
      call WriteParameter(ifile, "number", ne, ni)
      call WriteParameter(ifile, "charge", qe, qi)
      call WriteParameter(ifile, "mass", mass_e, mass_i)
      call WriteParameter(ifile, "number density per nm^3", rhoe_nm3,       rhoi_nm3)
      call WriteParameter(ifile, "number density per cm^3", rhoe_nm3*1.e21, rhoi_nm3*1.e21)
      call WriteParameter(ifile, "av. distance (sim. units)", a_ee, a_ii)
      call WriteParameter(ifile, "Temp (K)",  Te_K,  Ti_K)
      call WriteParameter(ifile, "Temp (eV)", Te_eV, Ti_eV)
      call WriteParameter(ifile, "Temp (Ry)", Te,    Ti)
      call WriteParameter(ifile, "intial Temp (eV)", Te_initial_eV,    Ti_initial_eV)
      call WriteParameter(ifile, "Theta (degen. param.)", &
                                   unit_kB*Te / (unit_hbar**2/(2.*unit_me)                 * (3.*pi**2*rhoe_nm3*unit_abohr_in_nm**3.)**(2./3.)), &
                                   unit_kB*Ti / (unit_hbar**2/(2.*unit_mp_over_me*unit_me) * (3.*pi**2*rhoi_nm3*unit_abohr_in_nm**3.)**(2./3.))  )
      call WriteParameter(ifile, "Gamma (coupl. param.)", &
                                   qe**2 / (4*pi*unit_epsilon0*a_ee) / (unit_kB*Te), &
                                   qi**2 / (4*pi*unit_epsilon0*a_ii) / (unit_kB*Ti)  )
      call WriteParameter(ifile, "therm. velocity", vte, vti)
      call WriteParameter(ifile, "plasma frequency", wpl_e, wpl_i)
      call WriteParameter(ifile, "plasma frequency (fs^-1)", wpl_e/unit_t0_in_fs, wpl_i/unit_t0_in_fs)
      call WriteParameter(ifile, "Debye length", lambdaD_e, lambdaD_i)
      call WriteParameter(ifile, "ion sphere radius", 0._8, a_i)
      call WriteTopline(  ifile, "", "")
      call WriteParameter(ifile, "npart_total", npart_total)
      call WriteParameter(ifile, "special setup (ispecial)", ispecial)
      call WriteParameter(ifile, "workflow setup", workflow_setup)
      call WriteParameter(ifile, "dt", dt)
      call WriteParameter(ifile, "plasma-frequency: max dt", maxdt(MAXDT_OMEGA_PLASMA))
      call WriteParameter(ifile, "position update (vte) : max dt", maxdt(MAXDT_POSUPDATE_VTE))
      call WriteParameter(ifile, "velocity update (vte) : max dt", maxdt(MAXDT_VELUPDATE_VTE))
      call WriteParameter(ifile, "laser-frequency: max dt", maxdt(MAXDT_OMEGA_LASER))
      call WriteParameter(ifile, "position update (vosc): max dt", maxdt(MAXDT_POSUPDATE_VTO))
      call WriteParameter(ifile, "velocity update (vosc): max dt", maxdt(MAXDT_VELUPDATE_VTO))
      call WriteParameter(ifile, "dt (fs)", dt*unit_t0_in_fs)
      call WriteParameter(ifile, "nt", nt)
      call WriteParameter(ifile, "simulation length (fs)", 1._8*nt*dt*unit_t0_in_fs)
      call WriteParameter(ifile, "A_ion", Aion)
      call WriteParameter(ifile, "Z_ion", Zion)
      call WriteParameter(ifile, "Vplas", Vplas)
      call WriteParameter(ifile, "x_plasma", x_plasma)
      call WriteParameter(ifile, "y_plasma", y_plasma)
      call WriteParameter(ifile, "z_plasma", z_plasma)
      call WriteParameter(ifile, "Gamma", physGamma)
      call WriteParameter(ifile, "force constant", 1._8*force_const)
      call WriteParameter(ifile, "eps", 1._8*eps)
      call WriteParameter(ifile, "V(r_ion=0) (eV)", V0_eV)
      call WriteTopline(  ifile, "")
      call WriteParameter(ifile, "theta", 1._8*sqrt(theta2))
      call WriteParameter(ifile, "eps(calc_force)", 1._8*sqrt(eps2))
      call WriteParameter(ifile, "mac_select", mac_select)
      call WriteParameter(ifile, "force_law", force_law)
      call WriteTopline(  ifile, "")
      call WriteParameter(ifile, "beam_config_in", beam_config_in)
      call WriteParameter(ifile, "omega / omega_pl", 1._8*omega_wpl)
      call WriteParameter(ifile, "omega (Hz)", omega_hz)
      call WriteParameter(ifile, "T_laser (fs)", 2._8*pi/omega * unit_t0_in_fs)
      call WriteParameter(ifile, "timesteps per laser cycle (navcycle)", navcycle)
      call WriteParameter(ifile, "lambda", lambda)
      call WriteParameter(ifile, "lambda (nm)", lambda_nm)
      call WriteParameter(ifile, "crit. e-density (1/nm^3)", rhocrit_nm3)
      call WriteParameter(ifile, "laser intensity (W/cm^2)", I0_Wpercm2)
      call WriteParameter(ifile, "laser field strength", E0)
      call WriteParameter(ifile, "laser field strength (V/cm)", E0*unit_E0_in_Vpercm)
      call WriteParameter(ifile, "laser pulse length (fs)", t_pulse_fs)
      call WriteParameter(ifile, "vosc (abohr/t0)", 1._8*vosc)
      call WriteParameter(ifile, "vosc / vtherm_e", vosc/vte)
      call WriteTopline(  ifile, "")
      call WriteParameter(ifile, "rngseed", rngseed)
      !call WriteTopline(  ifile, "")

       call PrintPeriodicityParameters(ifile)

   end subroutine PrintPhysicalParameters


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Outputs time-dependent laser parameters
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine PrintLaserParameters()
      use physvars, only : my_rank
      use module_io
      use module_laser
      use module_units
      use physvars, only : vte
      implicit none
      integer :: ifile

      if (my_rank .ne. 0 .or. beam_config_in.eq.0) return

      do ifile = file_stdout,file_pepc_out,file_pepc_out-file_stdout
        call WriteTopline(  ifile, "LASER")
        call WriteParameter(ifile, "beam_config_in",    beam_config_in)
        call WriteParameter(ifile, "t_laser",           t_laser)
        call WriteParameter(ifile, "t_laser (fs)",      t_laser*unit_t0_in_fs)
        call WriteParameter(ifile, "sin(phase_laser)",  sin(omega*t_laser))
        call WriteParameter(ifile, "Laser amplitude (max)", E_laser)
        call WriteParameter(ifile, "Laser amplitude (cur)", E_laser * sin(omega*t_laser))
        call WriteParameter(ifile, "Laser intensity", I_laser)
        call WriteParameter(ifile, "vosc (abohr/t0)", 1._8*vosc)
        call WriteParameter(ifile, "vosc / vtherm_e", vosc/vte)
        call WriteParameter(ifile, "x_crit",            x_crit)
        call WriteParameter(ifile, "spot size",         sigma)
        call WriteParameter(ifile, "steps per cycle",   navcycle)
        call WriteParameter(ifile, "Pulse envelope",    beam_envelope)
        call WriteParameter(ifile, "Pulse length",      t_pulse)
        call WriteParameter(ifile, "Pulse length (fs)", t_pulse_fs)
        call WriteParameter(ifile, "Focal width",       sigma)
        call WriteParameter(ifile, "Focal position",    focus(1))
        call WriteParameter(ifile, "Beam focus",        config_names(1))
        call WriteParameter(ifile, "Beam envelope",     config_names(2))
        call WriteParameter(ifile, "Beam model",        config_names(3))
        call WriteParameter(ifile, "Beam polarization", config_names(4))
        call WriteTopline(  ifile, "")

      end do

    end subroutine


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Output energies in pretty format
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine PrintEnergies(ifile, epot, ekini, ekine, etot, tempe, tempi, totalmomentum)
      use module_units
      use module_fmm_framework
      implicit none
      integer, intent(in) :: ifile
      real*8, intent(in) :: epot, ekini, ekine, etot, tempe, tempi, totalmomentum(3)

      call WriteTopline(  ifile, "ENERGIES", "")
      call WriteHeader(   ifile, "Ryd", "eV")
      call WriteParameter(ifile, "P.E.   (el+ions)   ", epot, epot*unit_Ryd_in_eV)
      call WriteParameter(ifile, "P.E.   (near field)", potnearfield, potnearfield*unit_Ryd_in_eV)
      call WriteParameter(ifile, "P.E.   (far field) ", potfarfield, potfarfield*unit_Ryd_in_eV)
      call WriteParameter(ifile, "K.E.   (electrons) ", ekine, ekine*unit_Ryd_in_eV)
      call WriteParameter(ifile, "K.E.   (ions)      ", ekini, ekini*unit_Ryd_in_eV)
      call WriteParameter(ifile, "Energy (total)     ", etot, etot*unit_Ryd_in_eV)
      call WriteTopline(  ifile, "", "")
      call WriteHeader(   ifile, "with drift", "without drift")
      call WriteParameter(ifile, "Te (eV)", unit_Ryd_in_eV*ekine/unit_kB*2./3., unit_Ryd_in_eV*tempe)
      call WriteParameter(ifile, "Ti (eV)", unit_Ryd_in_eV*ekini/unit_kB*2./3., unit_Ryd_in_eV*tempi)
      call WriteTopline(  ifile, "", "")
      call WriteParameter(ifile, "Total Momentum", totalmomentum(1:3))
      call WriteTopline(  ifile, "", "")

   end subroutine PrintEnergies



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Output periodicity parameters in pretty format
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine PrintPeriodicityParameters(ifile)
      use physvars
      use module_units
      use module_fmm_framework, only : fmm_extrinsic_correction
      use module_mirror_boxes
      implicit none
      integer, intent(in) :: ifile

      call WriteTopline(  ifile, "", "")
      call WriteParameter(ifile, "periodicity", periodicity)
      call WriteParameter(ifile, "t_lattice_1", t_lattice_1, .false.)
      call WriteParameter(ifile, "t_lattice_2", t_lattice_2, .false.)
      call WriteParameter(ifile, "t_lattice_3", t_lattice_3, .false.)
      call WriteParameter(ifile, "LatticeCenter", LatticeCenter, .false.)
      call WriteParameter(ifile, "LatticeOrigin", LatticeOrigin, .false.)
      call WriteParameter(ifile, "spatial_interaction_cutoff", spatial_interaction_cutoff, .false.)
      call WriteTopline(  ifile, "", "")
      call WriteParameter(ifile, "periodicity_nearest_image", periodicity_nearest_image)
      call WriteParameter(ifile, "mirror_box_layers", mirror_box_layers)
      call WriteParameter(ifile, "num_neighbour_boxes", num_neighbour_boxes)
      call WriteParameter(ifile, "fmm_extrinsic_correction", fmm_extrinsic_correction)
      call WriteTopline(  ifile, "")

   end subroutine PrintPeriodicityParameters

end module module_param_dump
