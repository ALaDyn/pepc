module module_namelist
  use physvars
  use module_pepc
  use module_fmm_framework
  use module_mirror_boxes
  use module_icosahedron
  use module_laser
  use module_pusher
  use module_workflow
  use module_units
  use module_fields
  use module_interaction_specific, only : eps2

  namelist /pepcmw/ &
       ne,  eps, V0_eV, nt, dt, idump, itime_in, idump_vtk, idump_checkpoint, idump_binary, treediags, directforce, & ! fundamental stuff
       ispecial, rhoe_nm3, Zion, Aion, Te_eV, Ti_eV, Te_K, Ti_K, rngseed, &   ! experimental setup
       workflow_setup, &                                             ! workflow
       integrator_scheme, enable_drift_elimination, &                ! pusher configuration
       beam_config_in, vosc,omega, sigma, t_pulse_fs, theta_inc, rho_track, omega_wpl, I0_Wpercm2, lambda_nm, t_laser, & ! laser config
       t_lattice_1, t_lattice_2, t_lattice_3, periodicity, do_extrinsic_correction, &            ! periodicity config
       field_dump_ncells, ngx, ngy, ngz                              ! diagnostics config

  contains

    subroutine write_frontend_parameters_to_file(filename)
      implicit none
      character(*), intent(in) :: filename
      integer, parameter :: filehandle = 91

      if (my_rank ==0) then
        open(filehandle,file=trim(filename),STATUS='UNKNOWN', POSITION = 'APPEND')
        write(filehandle, nml=pepcmw)
        close(filehandle)
      endif

    end subroutine


    subroutine read_frontend_parameters_from_file(filename)
      implicit none
      character(*), intent(in) :: filename
      integer, parameter :: filehandle = 91

      if(my_rank .eq. 0) write(*,*) "reading parameter file, section pepcmw: ", filename
      open(filehandle,file=trim(filename),action='read')
      read(filehandle, nml=pepcmw)
      close(filehandle)

    end subroutine


end module
