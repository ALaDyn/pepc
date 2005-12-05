 &pepcdata
! particles
  ne =0
  ni = 2000 

 plasma_config=1

! sphere
! initial_config = 1   

! random disc
!  initial_config = 2    

! wire
!   initial_config=3 

! rectangular slab
!  initial_config = 0 

! hollow sphere
    target_geometry=8  


! physics stuff

  theta = 0.6
! Temperatures in keV
  Te_keV = 0.1 
  Ti_keV =0.0 
  mass_ratio = 2000.
  q_factor = 1.
  coulomb = .false.
  lenjones = .true.
  bond_const = 3.e-3
! sphere/wire/disc radius
  r_sphere = 4
! plasma disc thickness/ wire length
  x_plasma = 1.    
! plasma width (slab target)
  y_plasma = 2.     
! plasma width (slab target)
  z_plasma = 2.     
  xl = 10  ! graphics box size
  yl =10 
  zl =10 
 ngx = 31
 ngy = 31
 ngz = 31

! beam
! fixed beam, initialised at start
  !  beam_config = 1  
! user-controlled, real-time particle source
 ! beam_config = 2  
! beam off
   beam_config_in = 0 
! laser fpond
!  beam_config=4  
 

  r_beam = 0.05
  u_beam = 0.2
  theta_beam = 0.0
  phi_beam = 0.0
  x_beam = .04
  start_beam = -0.1
  mass_beam = 5.
  rho_beam = -1.

! initial # beam particles/ dt
  np_beam = 0 

  vosc = 0.001
  omega = 0.5
  sigma = 6.
  tpulse = 20.
  lambda = 1.0   

  ! control
  nt =4000
  dt = 0.5
  eps = 0.2
 restart = .false.
  vis_on = .true.
 ivis = 5
 ivis_fields = 10
  mc_init = .false.
  mc_steps = 1000
  idump = 4000
  iprot=20
  itrack=300
  particle_bcs = 2
  scheme = 5 /

