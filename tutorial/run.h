 &pepcdata
 ncpu_merge = 1
 debug_level = 2
 debug_tree = 0
 mac=0
! particles
  ne = 0
  ni = 10000 
 nbuf_max=20000
! set up plasma target
 plasma_config = 1  

! sphere
! initial_config = 1   

! hollow sphere
  target_geometry =7 
 
! disc
!  target_geometry = 2       

! hollow hemisphere
 ! target_geometry =8

! wire
!   initial_config=3   

! rectangular slab
!  initial_config = 0        

! physics stuff

  theta = 0.5
 ! Temperatures in keV
  Te_keV = 0.5
  Ti_keV = 1
  mass_ratio = 2000.
  q_factor = 1.
  coulomb = .false.
  lenjones = .true.
  bond_const = 3.e-3
  r_sphere = 10.
! plasma disc thickness/ wire length
  x_plasma = .05   
! plasma width (slab target)
  y_plasma = 2.    
! plasma width (slab target)
  z_plasma = 2.    
! graphics box size
  xl = 2 
  yl =2 
  zl =2 


! beam
 ! fixed beam, initialised at start
 !  beam_config_in = 1 

 ! user-controlled, real-time particle source
 ! beam_config_in = 2 

! beam off
   beam_config_in = 0 

! laser fpond
!  beam_config_in=4 
 
  vosc = 6.0
  omega = 0.5
  sigma = 6.
  tpulse = 20.
  lambda = 1.0  

  ! control
  nt =500
  dt = 0.1
  eps = .3
 restart = .false.
  vis_on = .true.
 ivis = 20 
 ivis_fields = 5000
 ivis_domains = 5000
  mc_init = .false.
  mc_steps = 1000
  idump = 4000
  iprot=20
  itrack=300
  particle_bcs = 2
  scheme = 5 /

