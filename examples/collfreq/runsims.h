&pepccollfreq
 rngseed=13

! number of electrons
 ne = #NE#

 ispecial = 1
! ispecial = -1: reload particle positions from mpiio-timestamp #itime_in
 itime_in = 5

! number of timesteps
 nt = 300
 dt = 2.0

 beam_config_in = 0121
 vosc_vte        = #VOVTH#
 omega_wpl      = #WWPL#

 Te_eV    = #TE#
 Ti_eV    = #TE#
 rhoe_nm3 = #NENM3# ! 1.0e21cm^-3 = 1.0nm^-3
 Zion     =  1
 Aion     =  1
 eps      =  0.0 ! in units of the Debye length ! is set automatically anyway

workflow_setup = 5 ! [PRE 71, 056408 (2005)] P. Hilse et al: "Collisional absorption of dense plasmas in strong laser fields: quantum statistical results and simulation.", additionally fixing v0/vtherm
tau_temp_relaxation = 5.0
/

&calc_force_coulomb
  ! 3D Kelbg
  force_law  = 5
  ! BH-mac
  mac_select = 0
  ! theta = 0.4
  theta2     = 0.16
/

&libpepc
 debug_level = 0
 num_threads = 16
 interaction_list_length_factor = 1
/

&walk_para_pthreads
 max_particles_per_thread = 2000
/

