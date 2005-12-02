! Plasma sphere

 &pepcdata
 ncpu_merge = 1
 debug_level = 2
 debug_tree = 0
 mac=0
! particles
  ne = 4000
  ni = 4000 

 plasma_config = 1  ! set up plasma target
 target_geometry = 1   ! sphere

! physics stuff

  theta = 0.5
  Te_keV = 0.5 ! Temperatures in keV
  Ti_keV =0. 
  mass_ratio = 2000.
  q_factor = 1.
  coulomb = .true.
  lenjones = .false.
  bond_const = 2.e-3
  r_sphere = 2. 
  x_plasma = .1    ! plasma disc thickness/ wire length
  y_plasma = 2.     ! plasma width (slab target)
  z_plasma = 2.     ! plasma width (slab target)
  xl = 4  ! graphics box size
  yl =4 
  zl =4 


! beam
   beam_config_in = 0 ! beam off
 
  np_beam = 0 ! initial # beam particles/ dt

  vosc = 6.0
  omega = 0.5
  sigma = 6.
  tpulse = 20.
  lambda = 1.0   ! Wavelength in microns

  ! control
  nt =200
  dt = 0.3
  eps = 1.
 restart = .false.
  vis_on = .true.
 ivis = 5 
 ivis_fields = 5
 ivis_domains = 5000
  mc_init = .false.
  mc_steps = 1000
  idump = 4000
  iprot=1
  itrack=300
  particle_bcs = 1
  scheme = 1 /

