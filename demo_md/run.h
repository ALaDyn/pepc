&pepcdata
 ncpu_merge = 1
 restart = .false.
 debug_level = 4  ! debug switches
 debug_tree = 4
 mac=0            ! asynchron. walk
 walk_scheme=0
 balance=1
 ifreeze=1
 np_mult=-10      ! particle array adjust, neg. -> fix number -15 * 10000
 fetch_mult=3      ! multipole swap adjust
! particles
  ne = 0
  ni = 2000000

 plasma_config = 1  ! set up plasma target
 target_geometry =1   ! sphere

! physics stuff

  theta = 0.5
  Te_keV = 0.5 ! Temperatures in keV
  Ti_keV =0.
  mass_ratio = 2000.
  q_factor = 1.
  coulomb = .true.
  lenjones = .false.
  bond_const = 2.e-3
  force_const = 1.0
  r_sphere = 1.
  x_plasma = .1    ! plasma disc thickness/ wire length
  y_plasma = 2.     ! plasma width (slab target)
  z_plasma = 2.     ! plasma width (slab target)
  xl = 2  ! graphics box size
  yl =2
  zl =2
  ngx=2
  ngy=2
  ngz=2

! beam
   beam_config_in = 0 ! beam off

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
  nt = 2
  dt = 0.2
  eps = .2
  vis_on = .false.
 ivis = 2
 ivis_fields = 5000
 ivis_domains = 5000
  mc_init = .false.
  mc_steps = 1000
  idump =4000
  iprot=1
  itrack=300
  particle_bcs = 1
  scheme = 1 /
