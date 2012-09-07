! Plasma slab with perpendicular magnetic field
! Periodic in y-direction 

 &pepcb

 ncpu_merge = 1
 debug_level = 1 
 mac=0

! np_error = 200 ! uncomment to do error test
    

! particles
  ne = 160000000
  ni = 160000000

 plasma_config = 2  ! special start
 ispecial=7  ! gridded ions; electrons placed within a_ii/10
 target_geometry = 0   ! slab in xy plane 
 velocity_config = 4 ! 2D Maxwellian
 idim=2  ! ignore z coord
 force_law=2 ! Force law 2D
 force_const=0.1591549  ! 1/2pi
! force_const=0.

  theta = 0.55
  eps = 0.2  ! smoothing parameter in norm units
  
  vte=1.0  ! Choose Debye norms (vte, wp, lambda_De)
  Te_kev=1.0
  Ti_kev=1.0    ! Ti/Te=1/4

  mass_ratio = 100.
  x_plasma = 50.    ! plasma disc thickness/ wire length 
  y_plasma = 125.    ! plasma width (slab target) 
  z_plasma = 0.     ! plasma width (slab target)
  xl = 100  ! graphics box size
  yl = 125 
  zl =4
 displace = 50.,0.,0. 
 ngx=400
 ngy=400
 nxh=50
 ngav=50
 
! external field
  beam_config_in = 7  ! uniform Bz 
  vosc = 1.0 ! Omega_e/omega_p 

! fmm-periodicity framework
! lattice basis vectors
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
  nt = 100
  dt = 0.15
 restart = .false.
  vis_on = .false.
 ivis = 5 
 ivis_fields = 100 
 ivis_domains = 5000
  mc_init = .false.
  mc_steps = 1000
  idump = -500
  iprot=100
  itrack=300
  particle_bcs = 4 /


&libpepc

  np_mult=-100
  debug_level = 1

  ! 0: no load balacing, 1: load 
  weighted = 1                    

  ! space filling curve
  curve_type=0

 /

&calc_force_coulomb

/

&walk_para_pthreads
  num_walk_threads = 32
 /
