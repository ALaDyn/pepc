 &pepcdata
 debug_level = 2
 debug_tree = 2 
! particles
  ne = 300
  ni = 300 

 plasma_config = 1
 target_geometry = 0 ! slab


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
   beam_config_in = 0 ! beam off
 
  vosc = 6.0
  omega = 0.5
  sigma = 6.
  tpulse = 20.
  lambda = 1.0   ! Wavelength in microns

  ! control
  nt =1
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

