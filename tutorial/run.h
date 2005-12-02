 &pepcdata
! particles
  ne = 5000
  ni = 5000 

 plasma_config = 1  ! set up plasma target
 target_geometry = 0         ! slab 

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
  x_plasma = 2    ! plasma disc thickness/ wire length
  y_plasma = 5.     ! plasma width (slab target)
  z_plasma = 5.     ! plasma width (slab target)
  xl = 8 ! graphics box size
  yl =6 
  zl =6 
 ngx=60
 ngy=40
 ngz=40

! beam
  beam_config_in=14  !  laser 
 

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
  sigma = 1.
  tpulse = 2000.
  lambda = 1.0   ! Wavelength in microns

  ! control
  nt =4000
  dt = 0.5
  eps = 25.
 restart = .false.
  vis_on = .true.
 steering = .false.
 ivis = 10 
 ivis_fields = 10
  mc_init = .false.
  mc_steps = 1000
  idump = 4000
  iprot=20
  itrack=1
  particle_bcs = 2
  scheme = 1 /

