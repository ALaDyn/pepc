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

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  contains

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




    subroutine PrintPhysicalParameters(ifile)
      use physvars
      use module_laser
      use module_units
      implicit none
      integer, intent(in) :: ifile
! xl, yl, zl, plasma_centre
!  call setup_laser()

      write(ifile,'(92("-"))' )
      write(ifile,'(a1,a31,2(a20,a10))') "| ", "|", "electrons", "|", "ions","|"
      write(ifile,'(92("-"))' )
      call WriteParameterInt(ifile, "number", ne, ni)
      call WriteParameterReal(ifile, "charge", qe, qi)
      call WriteParameterReal(ifile, "mass", mass_e, mass_i)
      call WriteParameterReal(ifile, "number density per nm^3", rhoe_nm3, rhoi_nm3)
      call WriteParameterReal(ifile, "av. distance", a_ee, a_ii)

      call WriteParameterReal(ifile, "Temp (K)",  Te_K,  Ti_K)
      call WriteParameterReal(ifile, "Temp (eV)", Te_eV, Ti_eV)
      call WriteParameterReal(ifile, "Temp (Ry)", Te,    Ti)
      call WriteParameterReal(ifile, "therm. velocity", vte, vti)
      call WriteParameterReal(ifile, "plasma frequency", wpl_e, wpl_i)
      call WriteParameterReal(ifile, "Debye length", lambdaD_e, lambdaD_i)
      call WriteParameterReal(ifile, "ion sphere radius", 0._8, a_i)
      write(ifile,'(92("-"))' )
      call WriteParameterIntSingle(ifile, "npart_total", npart_total)
      call WriteParameterIntSingle(ifile, "special setup (ispecial)", ispecial)
      call WriteParameterRealSingle(ifile, "dt", dt)
      call WriteParameterRealSingle(ifile, "laser-frequency: max dt", maxdt(1))
      call WriteParameterRealSingle(ifile, "plasma-frequency: max dt", maxdt(2))
      call WriteParameterRealSingle(ifile, "position update: max dt", maxdt(3))
      call WriteParameterRealSingle(ifile, "velocity update: max dt", maxdt(4))
      call WriteParameterRealSingle(ifile, "dt (fs)", dt*unit_t0_in_fs)
      call WriteParameterIntSingle(ifile, "nt", nt)
      call WriteParameterRealSingle(ifile, "simulation length (fs)", 1._8*nt*dt*unit_t0_in_fs)
      call WriteParameterIntSingle(ifile, "A_ion", Aion)
      call WriteParameterIntSingle(ifile, "Z_ion", Zion)
      call WriteParameterRealSingle(ifile, "Vplas", Vplas)
      call WriteParameterRealSingle(ifile, "x_plasma", x_plasma)
      call WriteParameterRealSingle(ifile, "y_plasma", y_plasma)
      call WriteParameterRealSingle(ifile, "z_plasma", z_plasma)
      call WriteParameterRealSingle(ifile, "r_sphere", r_sphere)
      call WriteParameterRealSingle(ifile, "Gamma", physGamma)
      call WriteParameterRealSingle(ifile, "force constant", 1._8*force_const)
      call WriteParameterRealSingle(ifile, "eps", 1._8*eps)
      call WriteParameterRealSingle(ifile, "theta", 1._8*theta)
      write(ifile,'(62("-"))' )
      call WriteParameterIntSingle(ifile, "beam_config_in", beam_config_in)
      call WriteParameterIntSingle(ifile, "beam_config", beam_config)
      call WriteParameterRealSingle(ifile, "omega / omega_pl", 1._8*omega_wpl)
      call WriteParameterRealSingle(ifile, "omega (Hz)", 1._8*omega_hz)
      call WriteParameterRealSingle(ifile, "T_laser (fs)", 2._8*pi/omega * unit_t0_in_fs)
      call WriteParameterRealSingle(ifile, "timesteps per laser cycle (navcycle)", navcycle)
      call WriteParameterRealSingle(ifile, "lambda", 1._8*lambda)
      call WriteParameterRealSingle(ifile, "lambda (nm)", 1._8*lambda_nm)
      call WriteParameterRealSingle(ifile, "crit. e-density (1/nm^3)", rhocrit_nm3)
      call WriteParameterRealSingle(ifile, "laser intensity (W/cm^2)", I0_Wpercm2)
      call WriteParameterRealSingle(ifile, "laser field strength (V/cm)", E0*unit_E0_in_Vpercm)

      call WriteParameterRealSingle(ifile, "vosc (abohr/t0)", 1._8*vosc)
      call WriteParameterRealSingle(ifile, "vosc / vtherm_e", vosc/vte)
      write(ifile,'(62("-"))' )

   end subroutine PrintPhysicalParameters


end module module_param_dump
