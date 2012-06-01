&pepcmw
rngseed=00001


 treediags = .false.

! number of particles, here electrons
 ne = 4

! initial particle distribution
! 1 homogen, 2: one sphere, 3: two spheres, 4: Plummer (core cut)
 ispecial = 7
! ispecial = -1: reload particle positions from mpiio-timestamp #itime_in
itime_in = 5

! number of timesteps
 nt = 1
 dt = 2.0

! fmm-periodicity framework
! lattice basis vectors
  t_lattice_1 = 1.0   0.0   0.0
  t_lattice_2 = 0.0   1.0   0.0
  t_lattice_3 = 0.0   0.0   1.0
! periodicity in x-, y-, and z-direction
! 3D
  periodicity = .true.  .true.  .true.
! 2D
!   periodicity = .false. .true.  .true.
! 1D
!  periodicity = .false.  .false. .true.
! extrinsic-to-intrinsic correction
  do_extrinsic_correction = .false.

 beam_config_in = 0

 Te_eV    =  0.0
 Ti_eV    =  0.0
 rhoe_nm3 =  1.0
 Zion     =  1
 Aion     =  1
 eps      = 0.0

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

workflow_setup = 0

! determies the particle dump interval
! 0: never write anything
! n: each n-th step, plus first and last step
 idump = 1
! dito for vtk, binary and checkpoint-output
 idump_vtk = 1
 idump_binary = 0
 idump_checkpoint = 50
/

&calc_force_coulomb
  ! 3D coulomb
  force_law  = 3
  ! BH-mac
  mac_select = 0
  ! theta = 0.3
  theta2     = 0.0
/

&libpepc
 debug_level = 2049
 np_mult = -300
/

&walk_para_pthreads
 num_walk_threads         = 4
 max_particles_per_thread = 2000
 !interaction_list_length_factor = 8
/

