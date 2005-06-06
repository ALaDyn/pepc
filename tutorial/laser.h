 &pepcdata
 nmerge = 1
 !perf_anal=.true.
 domain_debug=.false.
 load_balance=.true.
 walk_balance=.true.
 !walk_debug=.true.
! particles
  ne = 50
  ni = 50 

 plasma_config = 1  ! set up plasma target
 ! initial_config = 1   ! sphere
 !   initial_config=7   ! hollow sphere
   target_geometry = 0         ! random disc
 !   initial_config=3   ! wire
 !  initial_config = 0         ! rectangular slab
  !  initial_config = 10     ! read from parts_all.in

! physics stuff

  theta = 0.5
  Te_keV = 0.5 ! Temperatures in keV
  Ti_keV =0. 
  mass_ratio = 2000.
  q_factor = 1.
  coulomb = .false.
  lenjones = .false.
  bond_const = 2.e-3
  r_sphere = 4
  x_plasma = 1    ! plasma disc thickness/ wire length
  y_plasma = 100.     ! plasma width (slab target)
  z_plasma = 100.     ! plasma width (slab target)
  xl = 37.7 ! graphics box size
  yl =100 
  zl =100 
 ngx=120
 ngy=80
 ngz=80

! beam
  !  beam_config = 1  ! fixed beam, initialised at start
 ! beam_config = 2  ! user-controlled, real-time particle source
!   beam_config = 0 ! beam off
  beam_config_in=14  ! oblique incidence laser 
 

  r_beam = 0.05
  u_beam = 0.2
  phi_beam = 0.0
  x_beam = .04
  start_beam = -0.1
  mass_beam = 5.
  rho_beam = -1.

  np_beam = 0 ! initial # beam particles/ dt

  vosc = 2.0
  omega = .5
  theta_beam = 0.0
  sigma = 20.
  tpulse = 2000.
  lambda = 1.0   ! Wavelength in microns

  ! control
  nt =4000
  dt = 0.5
  eps = 25.
 restart = .false.
  vis_on = .true.
 steering = .true.
 ivis = 2 
 ivis_fields = 2
  mc_init = .false.
  mc_steps = 1000
  idump = 4000
  iprot=20
  itrack=1
  particle_bcs = 2
  scheme = 1 /

