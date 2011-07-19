! Plasma sphere

 &pepcdata
 np_mult=-20
 fetch_mult=2
 ncpu_merge = 1
 debug_level = 1
 debug_tree = 0
 mac=0
! Choose sorting routine and load balancing                                                                                                                
! 0: no load balacing, 1: load balancing                                                                                                                   
 weighted = 1                                                                                                                                              
! Choose tree build routine                                                                                                                                
! 0: original, 1: optimized     
! choose_build=0
 curve_type=0
 walk_scheme = 0 

! particles
  ne = 3500
  ni = 3500 

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
  r_sphere = 2. 
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
  nt =20
  dt = 0.5
  eps = 0.5
 restart = .false.
  vis_on = .false.
 ivis = 5 
 ivis_fields = 5
 ivis_domains = 5000
  mc_init = .false.
  mc_steps = 1000
  idump = 20
  iprot=5
  itrack=300
  particle_bcs = 1 /



