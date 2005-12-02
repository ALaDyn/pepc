 &pepcdata
! particles
  ne = 12000
  ni = 12000 

target_geometry = 1   ! hollow sphere
plasma_config=1

! physics stuff

  theta = 0.5
  Te_keV = 2 ! Temperatures in keV
  Ti_keV =0. 
  mass_ratio = 2000.
  q_factor = 1.
  coulomb = .true.
  lenjones = .false.
  bond_const = 2.e-3
  r_sphere = 1
  x_plasma = 1.    ! plasma disc thickness/ wire length
  y_plasma = 2.     ! plasma width (slab target)
  z_plasma = 2.     ! plasma width (slab target)
  xl = 8  ! graphics box size
  yl =8 
  zl =8 


! beam
   beam_config_in = 0 ! beam off
 
  vosc = 6.0
  omega = 0.5
  sigma = 6.
  tpulse = 20.
  lambda = 1.0   ! Wavelength in microns

  ! control
  nt =400
  dt = 0.5
  eps = 2.
 restart = .false.
  vis_on = .true.
 ivis = 5 
 ivis_fields = 50
  mc_init = .false.
  mc_steps = 1000
  idump = 4000
  iprot=20
  itrack=300
  particle_bcs = 1
  scheme = 1 /

