 &pepcdata
 ncpu_merge = 1
 debug_level = 2
 debug_tree = 0
 balance=1 ! load balance switch
 mac=0
! particles
  ne = 2000
  ni = 2000 

 plasma_config = 1  ! set up plasma target
 ! target_geometry =7   ! hollow sphere
   target_geometry = 3         !  wire

! physics stuff

  theta = 0.6
  Te_keV = 0.5 ! Temperatures in keV
  Ti_keV =0. 
  mass_ratio = 2000.
  q_factor = 1.
  coulomb = .true.
  lenjones = .false.
  bond_const = 2.e-3
  r_sphere = 3
  x_plasma = 3.    ! plasma disc thickness
  y_plasma = 4.     ! plasma width (slab target)
  z_plasma = 10.     ! plasma width (slab target) / wire length
  xl = 1  ! graphics box size
  yl =1 
  zl =8 
  target_dup = .false.
  displace = 10.,10.,0.

! beam
  !  beam_config = 1  ! fixed beam, initialised at start
 ! beam_config = 2  ! user-controlled, real-time particle source
   beam_config_in = 4 ! laser fpond 
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

  vosc = 3.0
  omega = 0.5
  sigma = 3.
  tpulse = 100.
  lambda = 1.0   ! Wavelength in microns

  ! control
  nt =50
  dt = 0.2
  eps = 3.
 restart = .false.
  vis_on = .true.
 steering=.false.
 ivis =  5
 ivis_fields = 5000
 ivis_domains = 5000
  mc_init = .false.
  mc_steps = 1000
  idump = 4000
  iprot=1
  itrack=300
  particle_bcs = 1
  scheme = 1 /

