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

module module_namelist
  use physvars
  use module_pepc
  use module_laser
  use module_pusher
  use module_units
  use module_interaction_specific, only : eps2

  namelist /pepccollfreq/ &
       ne,  eps, forceconst, V0_eV, nt, dt, idump, itime_in, idump_vtk, idump_checkpoint, idump_binary, treediags, directforce, & ! fundamental stuff
       ispecial, rhoe_nm3, Zion, Aion, Te_eV, Ti_eV, Te_K, Ti_K, Te_initial_eV, Ti_initial_eV, rngseed, periodicity_nearest_image, &   ! experimental setup
       workflow_setup, &                                             ! workflow
       integrator_scheme, enable_drift_elimination, nose_hoover_Q_e, nose_hoover_Q_i, tau_temp_relaxation, &                ! pusher configuration
       beam_config_in, vosc,omega, sigma, t_pulse_fs, theta_inc, rho_track, omega_wpl, I0_Wpercm2, lambda_nm, t_laser, vosc_vte ! laser config

  contains

    subroutine write_frontend_parameters_to_file(filename)
      implicit none
      character(*), intent(in) :: filename
      integer, parameter :: filehandle = 91

      if (my_rank ==0) then
        open(filehandle,file=trim(filename),STATUS='UNKNOWN', POSITION = 'APPEND')
        write(filehandle, nml=pepccollfreq)
        close(filehandle)
      endif

    end subroutine


    subroutine read_frontend_parameters_from_file(filename)
      implicit none
      character(*), intent(in) :: filename
      integer, parameter :: filehandle = 91

      if(my_rank .eq. 0) write(*,*) "reading parameter file, section pepccollfreq: ", filename
      open(filehandle,file=trim(filename),action='read')
      read(filehandle, nml=pepccollfreq)
      close(filehandle)

    end subroutine


end module
