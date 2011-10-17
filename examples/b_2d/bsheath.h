! Plasma slab with perpendicular magnetic field
! Periodic in y-direction 

 &pepcdata
 np_mult=-20
 fetch_mult=2
 ncpu_merge = 1
 debug_level = 1
 debug_tree = 1
 mac=0

! np_error = 200 ! uncomment to do error test
 
! Choose sorting routine and load balancing
! 0: no load balancing, 1: load balancing
 weighted = 1                                                                     
! Choose tree build routine 
! 0: original, 1: optimized     
 choose_build=1
 curve_type=0  ! Morton curve
 walk_scheme = 0 
 num_walk_threads =1

! particles
  ne = 200
  ni = 200 

 plasma_config = 1  ! set up plasma target
 target_geometry = 0   ! slab in xy plane 
 velocity_config = 4 ! 2D Maxwellian
 idim=2  ! ignore z coord
 force_law=2 ! Force law 2D
 force_const=0.1591549  ! 1/2pi
! force_const=0.

  theta = 0.5
  eps = 1.0  ! smoothing parameter in norm units
  
  vte=1.0  ! Choose Debye norms (vte, wp, lambda_De)
  Te_kev=1.0
  Ti_kev=1.0

  mass_ratio = 100.
  q_factor = 1.
  coulomb = .false.
  lenjones = .false.
  bond_const = 2.e-3
  r_sphere = 3. 
  x_plasma = 15    ! plasma disc thickness/ wire length 
  y_plasma = 15.    ! plasma width (slab target) 
  z_plasma = 0.     ! plasma width (slab target)
  xl = 20  ! graphics box size
  yl =20 
  zl =4 
 ngx=100
 ngy=50
 nxh=50
 ngav=50
 
! external field
  beam_config_in = 0  ! uniform Bz 
  vosc = 0.2

! fmm-periodicity framework
! lattice basis vectors
! t_lattice_1 = 1. 0. 0.
! t_lattice_2 = 0. 1. 0.
! t_lattice_3 = 0. 0. 1.
 ! periodicity in x-, y-, and z-direction
 particle_wrap = .false.,.true.,.false.

  ! extrinsic-to-intrinsic correction
!  do_extrinsic_correction = .false. 

 scheme = 8 ! integration scheme: 2v, TE (Ex, Ey, Bz)                          

  ! control
  nt =100
  dt = 0.1
 restart = .false.
  vis_on = .false.
 ivis = 5 
 ivis_fields = 2
 ivis_domains = 5000
  mc_init = .false.
  mc_steps = 1000
  idump = 1
  iprot=10
  itrack=300
  particle_bcs = 3 /



