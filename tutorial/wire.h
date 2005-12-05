 &pepcdata
 ncpu_merge = 1
 debug_level = 2
 debug_tree = 0
 ! load balance switch
 balance=0
 mac=0
! particles
  ne = 4000
  ni = 4000 

! set up plasma target
 plasma_config = 1 

!  wire
   target_geometry = 3       

! physics stuff

  theta = 0.6
  Te_keV = 0.5
  Ti_keV =0. 
  mass_ratio = 100.
  q_factor = 1.
  coulomb = .true.
  lenjones = .false.
  bond_const = 2.e-3
  r_sphere = 4
! plasma disc thickness
  x_plasma = 3.   
! plasma width (slab target)
  y_plasma = 4.    
! plasma width (slab target) / wire length
  z_plasma = 20.     
  domain_cut = 6.
  xl = 1  ! graphics box size
  yl =1 
  zl =20 
  target_dup = .false.
  displace = 10.,10.,0.

! beam

! laser fpond with reduced transverse fields,  linear rise-time
   beam_config_in = 14 

! laser fpond, sin**2 pulse
!  beam_config_in=4 
 

  r_beam = 0.05
  u_beam = 0.2
  theta_beam = 0.0
  phi_beam = 0.0
  x_beam = .04
  start_beam = -0.1
  mass_beam = 5.
  rho_beam = -1.


  vosc = 1.0
  omega = 0.3
  sigma = 2.
  tpulse = 200.
  lambda = 1.0  

  ! control
  nt =500
  dt = 0.2
  eps = 3.
 restart = .false.
  vis_on = .true.
 steering=.false.
 ivis =  5
 ivis_fields = 5000
  mc_init = .false.
  mc_steps = 1000
  idump = 4000
  iprot=10
  itrack=300
  particle_bcs = 1
  scheme = 1 /

