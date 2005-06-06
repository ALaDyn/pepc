 &pepcdata
 nmerge = 1
 !perf_anal=.true.
 domain_debug=.false.
 load_balance=.true.
 walk_balance=.true.
 walk_summary=.false.
 db_level = 2
! particles
  ne = 500
  ni = 500 

 plasma_config = 1  ! set up plasma target
 ! initial_config = 1   ! sphere
 !   initial_config=7   ! hollow sphere
   target_geometry = 2         ! random disc
 !   initial_config=3   ! wire
 !  initial_config = 0         ! rectangular slab
  !  initial_config = 10     ! read from parts_all.in

! physics stuff

  theta = 0.5
  Te_keV = 0.5 ! Temperatures in keV
  Ti_keV =0. 
  mass_ratio = 2000.
  q_factor = 1.
  coulomb = .true.
  lenjones = .false.
  bond_const = 2.e-3
  r_sphere = 4
  x_plasma = 1.    ! plasma disc thickness/ wire length
  y_plasma = 2.     ! plasma width (slab target)
  z_plasma = 2.     ! plasma width (slab target)
  xl = 3  ! graphics box size
  yl =8 
  zl =8 


! beam
  !  beam_config = 1  ! fixed beam, initialised at start
 ! beam_config = 2  ! user-controlled, real-time particle source
   beam_config_in = 0 ! beam off
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
  nt =400
  dt = 0.5
  eps = 2.5
 restart = .false.
  vis_on = .true.
 ivis = 2
 ivis_domains = 5000 
 ivis_fields = 5000
  mc_init = .false.
  mc_steps = 1000
  idump = 4000
  iprot=20
  itrack=300
  particle_bcs = 1
  scheme = 1 /

