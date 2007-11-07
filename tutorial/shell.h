 &pepcdata
 ncpu_merge = 1
 debug_level = 2
 debug_tree = 0
 balance=0 ! load balance switch
 mac=0
! particles
  ne = 3000
  ni = 3000 

 launch = .true.
 plasma_config = 1  ! set up plasma target
 target_geometry =11   ! hollow sphere
 ramp = .false.  ! add ramp
 lolam = 0.2   ! L/lambda
 rho_min = 0.1  ! min density

! physics stuff

  theta = 0.6
  Te_keV = 0.5 ! Temperatures in keV
  Ti_keV =0. 
  mass_ratio = 100.
  q_factor = 1.
  coulomb = .true.
  lenjones = .false.
  bond_const = 2.e-3
  r_sphere = 2
  x_plasma = 2    ! plasma disc thickness
  y_plasma = 2.     ! plasma width (slab target)
  z_plasma = 2.     ! plasma width (slab target) / wire length
  domain_cut = 6.
  xl = 12  ! graphics box size
  yl =4 
  zl =4
  ngx=90
  ngy=30
  ngz=30
  target_dup = .false.
  displace = 0.,0.,0.

! beam
    beam_config_in = 0  ! fixed beam, initialised at start
  beam_config_in = 6 ! laser fpond with reduced transverse fields 
 

  r_beam = 0.05
  u_beam = 0.2
  theta_beam = 0.0
  phi_beam = 0.0
  x_beam = .04
  start_beam = -0.1
  mass_beam = 5.
  rho_beam = -1.

  np_beam = 0 ! initial # beam particles/ dt

  vosc = 0.3
  omega = 0.5
  sigma = 1.
  tpulse = 30.
  lambda = 1.0   ! Wavelength in microns

  ! control
  nt =1000
  dt = 0.2
  eps = 2.
 restart = .false.
  vis_on = .true.
 steering=.true.
 ivis =  5
 ivis_fields = 5
  mc_init = .false.
  mc_steps = 1000
  idump = 4000
  iprot=10
  itrack=300
  particle_bcs = 1
  scheme = 1 /

