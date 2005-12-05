 &pepcdata
! nmerge = 1
 mac=0
! particles
  ne = 4200
  ni = 4200 

! set up plasma target
 plasma_config = 1  
! sphere
 target_geometry = 1   

! physics stuff

  theta = 0.7
! Temperatures in keV
  Te_keV = 0.5 
  Ti_keV =0. 
  mass_ratio = 500.
  q_factor = 1.
  coulomb = .true.
  lenjones = .false.
  bond_const = 2.e-3
  r_sphere = 4
! plasma disc thickness/ wire length
  x_plasma = 1.    
! plasma width (slab target)
  y_plasma = 2.     
! plasma width (slab target)
  z_plasma = 2.     
  xl = 5  ! graphics box size
  yl =5 
  zl =5 


! beam
! uniform, sinusoid
 beam_config_in=3  
 
! initial # beam particles/ dt
  np_beam = 0 

! laser pump strength
  vosc = 0.1
  omega = 0.5
  sigma = 6.
  tpulse = 20.
! Wavelength in microns
  lambda = 1.0   

  ! control
  nt =500
  dt = 0.3
  eps = 2.5
 restart = .false.
  vis_on = .false.
 steering = .false.
 ivis = 5
 ivis_fields = 5000
  mc_init = .false.
  mc_steps = 1000
  idump = 4000
  iprot=2
  itrack=300
  particle_bcs = 1
  scheme = 1 /

