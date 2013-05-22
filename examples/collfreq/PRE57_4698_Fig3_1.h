&pepccollfreq

! setup follows one datapoint in Figure 3 in [PRE 57, 4698]

! number of particles, here electrons
 ne = 5000

! initial particle distribution
 ispecial = 1
! ispecial = -1: reload particle positions from mpiio-timestamp #itime_in
 itime_in = 5

! number of timesteps
 nt = 750
 dt = 2.0

! only perform nearest-image-periodicity
!  periodicity_nearest_image = .true.

 beam_config_in = 0121
! I0_Wpercm2     = 4.3E16
vosc_vte        = 0.2
! lambda_nm     = 436.0 ! is automatically set since omega_wpl is given
! t_pulse_fs    = 100.0 ! does not apply here
 omega_wpl      = 3.

 Te_eV    = 1.0
 Ti_K     = 1.000
 rhoe_nm3 = 10.0 ! 1.0e21cm^-3 = 1.0nm^-3
 Zion     =  1
 Aion     =  1
 eps      =  0.0 ! in units of the Debye length ! is set automatically anyway
! V0_eV    = -5.1 ! eps is set explicitly here

!  Available ensemble modes
! pure ES, NVT ensembles
!      1 = NVE - total energy conserved
!      2 = NVT - global Te, Ti conserved
!      3 = global NVT electron Te conserved; ions frozen
!      4 = local NVT: each PE keeps Te clamped; ions frozen
!      5 = local NVT, ions only; electrons added at end of run
! full EM pusher (all E, B components)
!      6 = NVE - total energy conserved
! nonrelativistic push
!      7 = NVE - total energy conserved
! integrator_scheme = 1              ! automatically set by workflow setup
! enable_drift_elimination = .true. ! automatically set by workflow setup

workflow_setup = 5 ! [PRE 71, 056408 (2005)] P. Hilse et al: "Collisional absorption of dense plasmas in strong laser fields: quantum statistical results and simulation.", additionally fixing v0/vtherm

! determies the particle dump interval - 0: never write anything, n: each n-th step, plus first and last step
 idump = 0
! dito for vtk, binary and checkpoint-output
 idump_vtk = 0
 idump_binary = 0
 idump_checkpoint = 1000
! activate output of branch nodes and space-filling curve as vtk-file
!  treediags = .true.

/

&calc_force_coulomb
  ! 3D Kelbg
  force_law  = 5
  ! BH-mac
  mac_select = 0
  ! theta = 0.3
  theta2     = 0.09
/

&libpepc
 num_threads = 60
 debug_level = 0
 interaction_list_length_factor = 1
/

&walk_para_pthreads
 max_particles_per_thread = 250
/

