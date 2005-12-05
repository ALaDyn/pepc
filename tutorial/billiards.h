!  Input deck for billiards example
!  Forces switched off; reflective container walls

&pepcdata
! particles
  ne =0
  ni = 200 

 plasma_config = 1
! slab
! target_geometry = 0 
! sphere
! target_geometry = 1   
! random disc
! target_geometry = 2  
! wire
! target_geometry = 3  
! ellipsoid
! target_geometry = 4  
! wedge
! target_geometry = 5  
! hemisphere
! target_geometry = 6  
! hollow sphere
! target_geometry = 7  
! hollow hemisphere
target_geometry = 8  


! physics stuff

  theta = 0.5
! Temperatures in keV
  Te_keV = 0.1 
  Ti_keV =0.5 
  mass_ratio = 100.
  q_factor = 1.
  coulomb = .false.
  lenjones = .false.
  bond_const = 2.e-5
  r_sphere = 5
! plasma disc thickness/ wire length
  x_plasma = 1    
! plasma width (slab target)
  y_plasma = 5.     
! plasma width (slab target)
  z_plasma = 3.    
! graphics box size
  xl = 5  
  yl =5 
  zl =3 


! beam (off)
   beam_config_in = 0 

  vosc = 6.0
  omega = 0.5
  sigma = 6.
  tpulse = 20.
! Wavelength in microns
  lambda = 1.0   

  ! control
  nt =2000
  dt = 1
  eps = 0.2
 restart = .false.
  vis_on = .true.
 ivis =  2
 ivis_fields = 500
  mc_init = .false.
  mc_steps = 1000
  idump = 20000
  iprot=10
  itrack=300
  particle_bcs = 2
  ncpu_merge = 1
  scheme = 1 /

