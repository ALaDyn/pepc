 &pepcdata
 nmerge = 1
 !perf_anal=.true.
 domain_debug=.false.
 build_debug=.false.
 tree_debug=.true.
 dump_tree=.true.
 load_balance=.false.
 walk_balance=.true.
 !walk_debug=.true.
! particles
  ne = 300
  ni = 300 

 plasma_config = 1
 ! initial_config = 1   ! sphere
 !  initial_config = 2         ! random disc
  !  initial_config=3   ! wire
!   initial_config = 0         ! rectangular slab
  !  initial_config = 10     ! read from parts_all.in
! initial_config = 4 ! ellipsoid
! initial_config = 5 ! wedge
! initial_config = 6 ! hemisphere
! initial_config = 7 ! hollow sphere
! target_geometry = 1 ! sphere
 target_geometry = 0 ! slab
! target_geometry = 8 ! hollow semisphere


! physics stuff

  theta = 0.5
  Te_keV = 0.1 ! Temperatures in keV
  Ti_keV =0.5 
  mass_ratio = 100.
  q_factor = 1.
  coulomb = .true.
  lenjones = .false.
  bond_const = 2.e-5
  r_sphere = 3
  x_plasma = 1    ! plasma disc thickness/ wire length
  y_plasma = 1.     ! plasma width (slab target)
  z_plasma = .01     ! plasma width (slab target)
  xl = 6  ! graphics box size
  yl =6 
  zl = 1


! beam
  !  beam_config = 1  ! fixed beam, initialised at start
 ! beam_config = 2  ! user-controlled, real-time particle source
   beam_config = 0 ! beam off
!  beam_config=4  ! laser fpond
 

  r_beam = 0.05
  u_beam = 0.2
  theta_beam = 0.0
  phi_beam = 0.0
  x_beam = .04
  start_beam = -0.1
  mass_beam = 5.
  rho_beam = -1.

  np_beam = 0 ! initial # beam particles/ dt

  vosc = 6.0
  omega = 0.5
  sigma = 6.
  tpulse = 20.
  lambda = 1.0   ! Wavelength in microns

  ! control
  nt =200
  dt = 1
  eps = 0.5
 restart = .false.
  vis_on = .false.
 ivis =  5
 ivis_fields = 500
  mc_init = .false.
  mc_steps = 1000
  idump = 20000
  iprot=10
  itrack=300
  particle_bcs = 2
  scheme = 1 /

