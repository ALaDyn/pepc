  &pepcdata
 nmerge = -1
 !perf_anal=.true.
 domain_debug=.false.
 load_balance=.true.
 walk_balance=.true.
 !walk_debug=.true.
! particles
  ne = 320000
  ni = 320000 

!  initial_config = 1   ! sphere
   initial_config = 2         ! random cylinder 
  !  initial_config=3   ! wire
 !  initial_config = 0         ! rectangular slab
  !  initial_config = 10     ! read from parts_all.in


! physics stuff

  theta = 0.75
  Te_keV = 0.3 ! Temperatures in keV
  Ti_keV =0.2 
  mass_ratio = 100.
  q_factor = 1.

  r_sphere = 7
  x_plasma = 40.    ! plasma disc thickness/ wire length
  y_plasma = 20.     ! plasma width (slab target)
  z_plasma = 20.
  xl = 45  ! graphics box size
  yl =20 
  zl =20 
  ngx = 100
  ngy = 41
  ngz = 41

! beam
  !  beam_config = 1  ! fixed beam, initialised at start
 ! beam_config = 2  ! user-controlled, real-time particle source
!   beam_config = 0 ! beam off
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

  vosc = 2.0
  omega = 5.
  sigma = 5.
  tpulse = 3.
  x_offset = -3. ! centre of pulse relative to plasma edge
  lambda = 1.0   ! Wavelength in microns

  ! control
  nt =1000
  dt = 0.5
  eps = 3.
 restart = .true.
  vis_on = .true.
 steering = .true.
 ivis = 5
 ivis_fields = 5
  idump = 10000
  iprot=20
  idens=50
  ivis_skip=4
  ensemble = 1 /

