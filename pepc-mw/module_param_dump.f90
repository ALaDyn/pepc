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

      interface WriteParameter
        module procedure WriteParameterReal
        module procedure WriteParameterInt
        module procedure WriteParameterRealSingle
        module procedure WriteParameterIntSingle
      end interface

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
      call WriteParameter(ifile, "number", ne, ni)
      call WriteParameter(ifile, "charge", qe, qi)
      call WriteParameter(ifile, "mass", mass_e, mass_i)
      call WriteParameter(ifile, "number density per nm^3", rhoe_nm3, rhoi_nm3)
      call WriteParameter(ifile, "av. distance", a_ee, a_ii)

      call WriteParameter(ifile, "Temp (K)",  Te_K,  Ti_K)
      call WriteParameter(ifile, "Temp (eV)", Te_eV, Ti_eV)
      call WriteParameter(ifile, "Temp (Ry)", Te,    Ti)
      call WriteParameter(ifile, "therm. velocity", vte, vti)
      call WriteParameter(ifile, "plasma frequency", wpl_e, wpl_i)
      call WriteParameter(ifile, "Debye length", lambdaD_e, lambdaD_i)
      call WriteParameter(ifile, "ion sphere radius", 0._8, a_i)
      write(ifile,'(92("-"))' )
      call WriteParameter(ifile, "npart_total", npart_total)
      call WriteParameter(ifile, "special setup (ispecial)", ispecial)
      call WriteParameter(ifile, "dt", dt)
      call WriteParameter(ifile, "laser-frequency: max dt", maxdt(1))
      call WriteParameter(ifile, "plasma-frequency: max dt", maxdt(2))
      call WriteParameter(ifile, "position update: max dt", maxdt(3))
      call WriteParameter(ifile, "velocity update: max dt", maxdt(4))
      call WriteParameter(ifile, "dt (fs)", dt*unit_t0_in_fs)
      call WriteParameter(ifile, "nt", nt)
      call WriteParameter(ifile, "simulation length (fs)", 1._8*nt*dt*unit_t0_in_fs)
      call WriteParameter(ifile, "A_ion", Aion)
      call WriteParameter(ifile, "Z_ion", Zion)
      call WriteParameter(ifile, "Vplas", Vplas)
      call WriteParameter(ifile, "x_plasma", x_plasma)
      call WriteParameter(ifile, "y_plasma", y_plasma)
      call WriteParameter(ifile, "z_plasma", z_plasma)
      call WriteParameter(ifile, "r_sphere", r_sphere)
      call WriteParameter(ifile, "Gamma", physGamma)
      call WriteParameter(ifile, "force constant", 1._8*force_const)
      call WriteParameter(ifile, "eps", 1._8*eps)
      call WriteParameter(ifile, "theta", 1._8*theta)
      write(ifile,'(62("-"))' )
      call WriteParameter(ifile, "beam_config_in", beam_config_in)
      call WriteParameter(ifile, "beam_config", beam_config)
      call WriteParameter(ifile, "omega / omega_pl", 1._8*omega_wpl)
      call WriteParameter(ifile, "omega (Hz)", 1._8*omega_hz)
      call WriteParameter(ifile, "T_laser (fs)", 2._8*pi/omega * unit_t0_in_fs)
      call WriteParameter(ifile, "timesteps per laser cycle (navcycle)", navcycle)
      call WriteParameter(ifile, "lambda", 1._8*lambda)
      call WriteParameter(ifile, "lambda (nm)", 1._8*lambda_nm)
      call WriteParameter(ifile, "crit. e-density (1/nm^3)", rhocrit_nm3)
      call WriteParameter(ifile, "laser intensity (W/cm^2)", I0_Wpercm2)
      call WriteParameter(ifile, "laser field strength (V/cm)", E0*unit_E0_in_Vpercm)

      call WriteParameter(ifile, "vosc (abohr/t0)", 1._8*vosc)
      call WriteParameter(ifile, "vosc / vtherm_e", vosc/vte)
      write(ifile,'(62("-"))' )

   end subroutine PrintPhysicalParameters


end module module_param_dump