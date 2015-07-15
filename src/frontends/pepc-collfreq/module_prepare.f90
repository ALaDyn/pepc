! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2015 Juelich Supercomputing Centre, 
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

module module_prepare

contains

subroutine frontend_prepare()
  use physvars
  use module_pepc
  use module_fmm_framework
  use module_mirror_boxes
  use module_laser
  use module_pusher
  use module_units
  use module_interaction_specific, only : eps2, include_far_field_if_periodic
  use module_namelist
  use module_io
  use module_interaction_specific, only : force_law
  implicit none

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!  derived parameters (physics)  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  qe     = unit_qe
  mass_e = unit_me
  qi     = unit_qp*Zion
  ni     = ne/Zion ! we assume overall charge neutrality (Zion = ion charge number)
  mass_i = unit_mp*Aion ! Aion = ion mass number
  npart_total = ni+ne

  Vplas    =  ne/rhoe_nm3 / unit_abohr_in_nm**3. ! adjust simulation volume to fit requested electron density while keeping particle number constant
  rhoi_nm3 =  ni/Vplas / unit_abohr_in_nm**3.
  x_plasma = (ne/rhoe_nm3)**(1./3.) / unit_abohr_in_nm ! (assume cubic volume)

  y_plasma = x_plasma
  z_plasma = x_plasma
  xl       = x_plasma
  yl       = y_plasma
  zl       = z_plasma

  a_ee = (Vplas/ne)**(1./3.)
  a_ii = (Vplas/ni)**(1./3.)

  if (Te_K > 0.) Te_eV = unit_kB_in_eVperK * Te_K
  if (Ti_K > 0.) Ti_eV = unit_kB_in_eVperK * Ti_K
  Te    = Te_eV / unit_Ryd_in_eV
  Ti    = Ti_eV / unit_Ryd_in_eV
  Te_K  = Te_eV / unit_kB_in_eVperK
  Ti_K  = Ti_eV / unit_kB_in_eVperK

  tempe = Te
  tempi = Ti

  vte = sqrt(3*unit_kB*Te/mass_e)
  vti = sqrt(3*unit_kB*Ti/mass_i)

  force_const = forceconst/(unit_4piepsilon0)

  wpl_e = sqrt( (ne/Vplas * qe *qe) / (unit_epsilon0 * mass_e) )
  wpl_i = sqrt( (ni/Vplas * qi *qi) / (unit_epsilon0 * mass_i) )

  lambdaD_e = sqrt( (unit_epsilon0*unit_kB*Te)/(qe*qe) * Vplas/ne )
  lambdaD_i = sqrt( (unit_epsilon0*unit_kB*Ti)/(qi*qi) * Vplas/ni )

  t_lattice_1 = [1.0, 0.0, 0.0]*x_plasma
  t_lattice_2 = [0.0, 1.0, 0.0]*y_plasma
  t_lattice_3 = [0.0, 0.0, 1.0]*z_plasma
  fmm_extrinsic_correction = FMM_EXTRINSIC_CORRECTION_REDLACK
  periodicity = [.true., .true., .true.]
  
  LatticeOrigin = LatticeOrigin*[x_plasma, y_plasma, z_plasma]

  if (periodicity_nearest_image) then
    !include_far_field_if_periodic = .false.
    spatial_interaction_cutoff = [x_plasma, y_plasma, z_plasma] * mirror_box_layers
  endif

  a_i       = (4.*pi/3. * ni/Vplas)**(-1./3.)
  physGamma = (qi*qi) / (a_i * unit_kB*Te)

  !if (lambdaD_e > 0.) eps = eps * lambdaD_e
  eps = eps * a_ii

  if (V0_eV .ne. 0.)  eps = force_const * qe*qi / (V0_eV / unit_Ryd_in_eV)

  V0_eV = (force_const * qe*qi / eps) * unit_Ryd_in_eV

  eps2 = eps*eps

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!  parameters (laser)            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (omega_wpl > 0.) omega = omega_wpl * wpl_e
  if (lambda_nm > 0.) omega = unit_c / (lambda_nm / unit_abohr_in_nm)
  omega_wpl = omega / wpl_e
  omega_hz  = omega / unit_t0_in_s
  lambda    = unit_c / omega
  lambda_nm = lambda * unit_abohr_in_nm
  rhocrit_nm3 = omega*omega*unit_epsilon0*mass_e/(qe*qe) / unit_abohr_in_nm**3.

  if (I0_Wpercm2 > 0.) then
    E0   =  sqrt( 2./(unit_epsilon0*unit_c) * (I0_Wpercm2 / unit_P0_in_W * (100*unit_abohr_in_m)**2) )
    vosc = (abs(qe)*E0)/(mass_e*omega)
  else 
    if (vosc_vte > 0.) then
      vosc = vosc_vte * vte
    endif
  endif
  
  E0         = vosc*mass_e*omega/abs(qe)
  I0_Wpercm2 = (unit_epsilon0 * unit_c * E0**2 / 2. ) * unit_P0_in_W / (100*unit_abohr_in_m)**2

  call laser_setup()


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!  parameters (simulation generic)   !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call adjust_maxdt()
  
  trun = itime * dt
  
  tau_temp_relaxation = tau_temp_relaxation*dt
  if (nose_hoover_Q_e < 0) nose_Hoover_Q_e = tau_temp_relaxation*tau_temp_relaxation*3./2.*ne*unit_kB*Te
  if (nose_hoover_Q_i < 0) nose_Hoover_Q_i = tau_temp_relaxation*tau_temp_relaxation*3./2.*ni*unit_kB*Ti


end subroutine frontend_prepare


subroutine adjust_maxdt()
  use physvars
  use module_laser
  use module_units
  use module_interaction_specific, only : eps2, include_far_field_if_periodic, force_law
  use module_io
  implicit none

  real*8 :: maxF
  real*8 :: vte_init
  integer :: ifile

  maxdt(:) = huge(maxdt(1))
  
  maxdt(MAXDT_OMEGA_LASER) = 2./max(omega,1.e-10_8) * 1./19.
  maxdt(MAXDT_OMEGA_PLASMA) = 2./max(wpl_e,1.e-10_8) * 1./19.
  
  vte_init = max(vte, sqrt(3.*unit_kB*Te_initial_eV/unit_Ryd_in_eV/mass_e))
  
  tend = nt*dt
  dt_ori = dt

  if (force_law == 5) then
    ! in r->0 Kelbg is like Plummer with eps = lambda
    eps = unit_hbar / sqrt(mass_e*unit_kB*Te)
  endif
  
  if (eps > 0.) then
    maxF = abs(qe*qi/unit_4piepsilon0 / (eps*eps))
  else
    maxF = epsilon(maxF)
  endif


  if (vte_init>0.) then
    maxdt(MAXDT_POSUPDATE_VTE) = a_ee / vte_init          / 10.
    maxdt(MAXDT_VELUPDATE_VTE) = mass_e * vte_init / maxF !/ 10. ! here we assume that particles never approach each other too close, so that maxF = maxF/10.
  endif
  
  if (vosc>0.) then
    maxdt(MAXDT_POSUPDATE_VTO) = a_ee / vosc          / 10.
    maxdt(MAXDT_VELUPDATE_VTO) = mass_e * vosc / maxF !/ 10. ! here we assume that particles never approach each other too close
  endif

  if (any(maxdt < dt)) then
    if (my_rank == 0) then
      do ifile = file_stdout,file_pepc_out,file_pepc_out-file_stdout
        write(ifile,*) "!!!!!!!! WARNING: timestep dt is too large: dt =", dt, "   maxdt = ", maxdt
        write(ifile,*) "Adjusting parameters appropriately to ensure sufficient resolution while keeping simulation length"
      end do
    end if

    nt = int(tend / minval(maxdt))
    dt = minval(maxdt)
  endif
  
end subroutine

end module


