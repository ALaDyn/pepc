!  Input deck for constant temperature dynamics (scheme=4) 
!  puts electrons in thermal equilibrium

&pepcdata
 ncpu_merge = 1
 debug_level = 2
 debug_tree = 0
 mac=0
! particles
  ne = 10000
  ni = 10000 
 nbuf_max=30000
 plasma_config = 1  ! set up plasma target
 ! initial_config = 1   ! sphere
  target_geometry =7   ! hollow sphere
 !  target_geometry = 3         ! random disc
 !   initial_config=7   ! hollow sphere
 !   initial_config=3   ! wire
 !  initial_config = 0         ! rectangular slab
  !  initial_config = 10     ! read from parts_all.in

! physics stuff

  theta = 0.3
  Te_keV = 0.5 ! Temperatures in keV
  Ti_keV =0. 
  mass_ratio = 2000.
  q_factor = 1.
  coulomb = .true.
  lenjones = .false.
  bond_const = 2.e-3
  r_sphere = 10.
  x_plasma = .1    ! plasma disc thickness/ wire length
  y_plasma = 2.     ! plasma width (slab target)
  z_plasma = 2.     ! plasma width (slab target)
  xl = 2  ! graphics box size
  yl =2 
  zl =2 


! beam
  !  beam_config = 1  ! fixed beam, initialised at start
 ! beam_config = 2  ! user-controlled, real-time particle source
   beam_config_in = 0 ! beam off
!  beam_config=4  ! laser fpond
 
  vosc = 6.0
  omega = 0.5
  sigma = 6.
  tpulse = 20.
  lambda = 1.0   ! Wavelength in microns

  ! control
  nt =100
  dt = 0.2
  eps = 2.
 restart = .false.
  vis_on = .true.
 ivis = 2 
 ivis_fields = 5000
 ivis_domains = 5000
  mc_init = .false.
  mc_steps = 1000
  idump = 4000
  iprot=1
  itrack=300
  particle_bcs = 1
  scheme = 4 /

