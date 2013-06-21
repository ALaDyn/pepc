! Capacitor plates - separated ion/electron slabs
! Select with plasma_config=3, velocity_config=0
! Periodic in y-direction 

 &pepcb
 ncpu_merge = 1
 debug_level = 1
 mac=0

! np_error = 200 ! uncomment to do error test
 
! particles
  ne = 25600
  ni = 25600 

! plasma_config = 3  ! set up ion, electron plates, random particles
! target_geometry = 0   ! slab in xy plane 
 plasma_config=2 ! special start
 ispecial=8  ! slab in xy plane; gridded particles
 velocity_config = 0 ! Cold start
 idim=2  ! ignore z coord
 force_law=2 ! Force law 2D
 force_const=0.1591549  ! 1/2pi
! force_const=0.
    
  theta = 0.5
  eps = 0.1  ! smoothing parameter in norm units
  
  vte=1.0  ! Choose Debye norms (vte, wp, lambda_De)
  Te_kev=1.0
  Ti_kev=1.0

  mass_ratio = 100.
  q_factor = 1.
  coulomb = .false.
  lenjones = .false.
  bond_const = 2.e-3
  r_sphere = 3. 
  x_plasma = 1.    ! plasma disc thickness/ wire length 
  y_plasma = 16.    ! plasma width (slab target) 
  z_plasma = 0.     ! plasma width (slab target)
  displace = 10.,0.,0.
  xl = 20  ! graphics box size
  yl = 16 
  zl =4 
 ngx=50
 ngy=50
 nxh=50
 ngav=50
 
! external field
  beam_config_in = 0  ! uniform Bz 
  vosc = 1.0 

! fmm-periodicity framework
! lattice basis vectors - these get scaled to xl,yl,zl in pepc-b setup
! t_lattice_1 = 1. 0. 0.
 t_lattice_2 = 0. 1. 0.
! t_lattice_3 = 0. 0. 1.
 ! periodicity in x-, y-, and z-direction
 periodicity = .false., .true., false.  ! forces periodic in y (1st neighbour box only)
 particle_wrap = .false.,.true.,.false.

  ! extrinsic-to-intrinsic correction
!  do_extrinsic_correction = .false. 

 scheme = 8 ! integration scheme: 2v, non-rel TE (Ex, Ey, Bz)                          

  ! control
  nt =100
  dt = 0.05
 restart = .false.
  vis_on = .false.
 ivis = 5 
 ivis_fields = 1
 ivis_domains = 5000
  mc_init = .false.
  mc_steps = 1000
  idump = 1
  iprot=10
  itrack=300
  particle_bcs = 3 /

&libpepc
      np_mult=-20
      debug_level = 1
      
      ! Choose sorting routine and 
      ! 0: no load balancing, 1: load balancing
      weighted = 1
      curve_type=0  ! Morton curve
      num_threads = 4
      /

&calc_force_coulomb

/

&walk_para_pthreads
      /

