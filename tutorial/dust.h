 &pepcdata
 nmerge = 1
 !perf_anal=.true.
 domain_debug=.false.
 !load_balance=.true.
 !walk_balance=.true.
 !walk_debug=.true.
! particles
  ne = 5000
  ni = 5000 

 ! initial_config = 1   ! sphere
 !  initial_config = 2         ! random disc
 !   initial_config=3   ! wire
   initial_config = 0         ! rectangular slab
  !  initial_config = 10     ! read from parts_all.in


! physics stuff

  theta = 0.5
  Te_keV =.1 ! Temperatures in keV
  Ti_keV =0.0 
  mass_ratio = 2000.
  q_factor = 1.

  r_sphere = 1
  x_plasma = 0.25    ! plasma disc thickness/ wire length
  y_plasma = 2    ! plasma width (slab target)
  z_plasma = 2    ! plasma width (slab target)
  xl = 0.3  ! graphics box size
  yl =2 
  zl =2 


! beam
  !  beam_config = 1  ! fixed beam, initialised at start
 ! beam_config = 2  ! user-controlled, real-time particle source
  beam_config = 3
!   beam_config = 0 ! beam off
!  beam_config=4  ! laser fpond
 

  r_beam = 0.3
  u_beam = 0.0
  theta_beam = 1.57
  phi_beam = 0.0
  x_beam = .04
  start_beam = -0.3
  mass_beam = 5000.
  rho_beam = .2

  np_beam = 300 ! initial # beam particles

  vosc = 6.0
  omega = 0.5
  sigma = 6.
  tpulse = 20.
  lambda = 1.0   ! Wavelength in microns

  ! control
  nt =20000
  dt = 0.2
  eps = 0.3
 restart = .false.
  vis_on = .true.
 steering = .true.
 ivis = 10
 ivis_fields = 5000
  mc_init = .false.
  mc_steps = 1000
  idump = 10000
  iprot=5000
  itrack=3000
  scheme = 4 /

