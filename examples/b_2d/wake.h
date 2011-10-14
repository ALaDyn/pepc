! Plasma disc 

 &pepcdata
 np_mult=-20
 fetch_mult=2
 ncpu_merge = 1
 debug_level = 1
 debug_tree = 1
 mac=0
! Force law 2D
 idim=2
 force_law=2
 force_const=.1592
! np_error = 200 ! uncomment to do error test
 
! Choose sorting routine and load balancing                                                                                                                
! 0: no load balacing, 1: load balancing                                                                                                                   
 weighted = 1                                                                                                                                              
! Choose tree build routine                                                                                                                                
! 0: original, 1: optimized     
! choose_build=0
 curve_type=0
 walk_scheme = 0 
 num_walk_threads =2

! particles
  ne = 1000000
  ni = 1000000 

! set up plasma target
 plasma_config = 1  
! slab in xy plane 
 target_geometry = 0   
 velocity_config = 4 ! 2D Maxwellian
! physics stuff

  theta = 0.5
  eps = 0.2
  
! Temperatures in keV
  Te_keV =0. 
  Ti_keV =0. 
  mass_ratio = 9000.
  rho0=1.0
  q_factor = 1.
  coulomb = .true.
  lenjones = .false.
  bond_const = 2.e-3
  r_sphere = 3. 
! plasma disc thickness/ wire length
  x_plasma = 40.   
! plasma width (slab target)
  y_plasma = 20.     
! plasma width (slab target)
  z_plasma = 0.     
! graphics box size
  xl = 50  
  yl =20 
  zl =4 
 ngx=200
 ngy=100
 nxh=50
 ngav=50
! beam
   beam_config_in = 5 

  vosc = 1.0
  omega = 10.0
  sigma = 5.
  tpulse = 1.
  lambda = 0.8   
  x_offset=-25
  
! fmm-periodicity framework
 ! lattice basis vectors                                                                                   
! t_lattice_1 = 1. 0. 0.
! t_lattice_2 = 0. 1. 0.                                                                                     
! t_lattice_3 = 0. 0. 1.                                                                                    
 ! periodicity in x-, y-, and z-direction                                                                   
 ! periodicity = .true.  .true.  .true.                                                                     
! periodicity = .false.,.false.,.false.                                                                     
  ! extrinsic-to-intrinsic correction                                                                                                                        
!  do_extrinsic_correction = .false.              
!  Available ensemble modes                                                                                                                                
! pure ES, NVT ensembles                                                                                   
!      1 = NVE - total energy conserved                                                                     
!      2 = NVT - global Te, Ti conserved                                                                          
!      3 = global NVT electron Te conserved; ions frozen                                                    
!      4 = local NVT: each PE keeps Te clamped; ions frozen                                                          
!      5 = local NVT, ions only; electrons added at end of run                                                                                             
! full EM pusher (all E, B components)                                                                      
!      6 = NVE - total energy conserved                                                             
! nonrelativistic push                                                                                    
!      7 = NVE - total energy conserved                                                                                                                    
 scheme = 1                          

  ! control
  nt =250
  dt = 0.2
 restart = .false.
  vis_on = .false.
 ivis = 5 
 ivis_fields = 25
 ivis_domains = 5000
  mc_init = .false.
  mc_steps = 1000
  idump = 50
  iprot=10
  itrack=300
  particle_bcs = 1 /



