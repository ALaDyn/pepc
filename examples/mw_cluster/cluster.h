&pepcmw

 treediags = .true.

! number of particles, here electrons
 ne = 310

! initial particle distribution
! 1 homogen, 2: one sphere, 3: two spheres, 4: Plummer (core cut)
 ispecial = 12
! ispecial = -1: reload particle positions from mpiio-timestamp #itime_in
itime_in = 5

! number of timesteps
 nt = 35
 dt = 2.0

! fmm-periodicity framework
! lattice basis vectors
  t_lattice_1 = 1.0   0.0   0.0
  t_lattice_2 = 0.0   1.0   0.0
  t_lattice_3 = 0.0   0.0   1.0
! periodicity in x-, y-, and z-direction
!  periodicity = .true.  .true.  .true.
! extrinsic-to-intrinsic correction
  do_extrinsic_correction = .false.

 beam_config_in = 0221
 I0_Wpercm2     = 1.0E13
 lambda_nm      = 436.0
 t_pulse_fs     = 0.25
 !omega_wpl      = 3.

 Te_eV    =  0.0
 Ti_eV    =  0.0
 rhoe_nm3 = 35.0
 Zion     =  1
 Aion     = 22
 V0_eV    = -5.1

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
integrator_scheme = 1
enable_drift_elimination = .true.

workflow_setup = 3

! determies the particle dump interval
! 0: never write anything
! n: each n-th step, plus first and last step
 idump = 0
! dito for vtk, binary and checkpoint-output
 idump_vtk = 0
 idump_binary = 0
 idump_checkpoint = 0
/

&calc_force_coulomb
  ! 3D coulomb
  force_law  = 3
  ! BH-mac
  mac_select = 0
  ! theta = 0.3
  theta2     = 0.09
/

&libpepc
 debug_level = 0
 interaction_list_length_factor = 2
/

&walk_para_pthreads
 num_walk_threads         = 4
 max_particles_per_thread = 2000
/

