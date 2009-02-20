! Plasma sphere

 &pepcdata
 np_mult=10
 fetch_mult=2
 ncpu_merge = 1
 debug_level = 2
 debug_tree = 2
 mac=0
! particles
  ne = 35000
  ni = 35000 

! set up plasma target
 plasma_config = 1  
! sphere
 target_geometry = 1   

! physics stuff

  theta = 0.5
! Temperatures in keV
  Te_keV = 1.0 
  Ti_keV =0.1 
  mass_ratio = 2000.
  q_factor = 1.
  coulomb = .true.
  lenjones = .false.
  bond_const = 2.e-3
  r_sphere = 1. 
! plasma disc thickness/ wire length
  x_plasma = 1    
! plasma width (slab target)
  y_plasma = 1.     
! plasma width (slab target)
  z_plasma = 1.     
! graphics box size
  xl = 4  
  yl =4 
  zl =4 
 ngx=50
 nxh=50
 ngav=50
! beam
   beam_config_in = 0 
 

  vosc = 6.0
  omega = 0.5
  sigma = 6.
  tpulse = 20.
  lambda = 1.0   

  ! control
  nt =10
  dt = 0.5
  eps = 0.5
 restart = .false.
  vis_on = .false.
 ivis = 5 
 ivis_fields = 5
 ivis_domains = 5000
  mc_init = .false.
  mc_steps = 1000
  idump = -4000
  iprot=1
  itrack=300
  particle_bcs = 1
  scheme = 1 /

