 &pepcdata
 nmerge =-1 
 !perf_anal=.true.
 domain_debug=.false.
 load_balance=.true.
 walk_balance=.true.
 !walk_debug=.true.
! particles
  ne = 10000
  ni = 10000 

!  initial_config = 1   ! sphere
   initial_config = 2         ! random cylinder 
  !  initial_config=3   ! wire
 !  initial_config = 0         ! rectangular slab
  !  initial_config = 10     ! read from parts_all.in


! physics stuff

  theta = 0.7
  Te_keV = 0.1 ! Temperatures in keV
  Ti_keV =0.2 
  mass_ratio = 200.
  q_factor = 1.

  r_sphere = 2.
  x_plasma = 8.   ! plasma disc thickness/ wire length
  y_plasma = 2.     ! plasma width (slab target)
  z_plasma = 2.
  xl = 8 ! graphics box size
  yl =4.
  zl =4 
  ngx = 80
  ngy = 20
  ngz = 20

! beam
  !  beam_config = 1  ! fixed beam, initialised at start
 ! beam_config = 2  ! user-controlled, real-time particle source
 !  beam_config = 0 ! beam off
  beam_config=5  ! laser bullet (propagating fpond) 
 

  r_beam = 0.05
  u_beam = 0.2
  theta_beam = 0.0
  phi_beam = 0.0
  x_beam = .04
  start_beam = -0.1
  mass_beam = 5.
  rho_beam = -1.

  np_beam = 0 ! initial # beam particles/ dt

  vosc = 0.3
  omega = 5.
  sigma = 8
  tpulse = 1.
  x_offset = -3. ! centre of pulse relative to plasma edge
  lambda = 1.0   ! Wavelength in microns

  ! control
  nt =10
  dt = 0.2
  eps = 3.
 restart = .false.
  vis_on = .true.
 steering = .true.
 ivis = 1
 ivis_fields = 1
  idump = 500000
  iprot=10
  idens=50
  ensemble = 1 /


