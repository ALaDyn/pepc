 &pepcdata

 domain_debug=.false.
! particles
  ne = 51200 
  ni = 51200 
 

  initial_config = 1   ! sphere
 !  initial_config = 2         ! random disc
  !  initial_config=3   ! wire
 !  initial_config = 0         ! rectangular slab
  !  initial_config = 10     ! read from parts_all.in


! physics stuff

  theta = 0.5
  Te_keV = 5. ! Temperatures in keV
  Ti_keV =0.02 
  mass_ratio = 200.
  q_factor = 1.

  r_sphere = 6.0
  x_plasma = 0.01    ! plasma disc thickness/ wire length
  y_plasma = 20     ! plasma width (slab target)
  xl = 12  ! graphics box size
  yl =12 
  zl =12 


! beam
  !  beam_config = 1  ! fixed beam, initialised at start
 ! beam_config = 2  ! user-controlled, real-time particle source
   beam_config = 0 ! beam off
!  beam_config=4  ! laser fpond
 

  r_beam = 0.05
  u_beam = 0.2
  theta_beam = 0.0
  phi_beam = 0.0
  x_beam = .04
  start_beam = -0.1
  mass_beam = 5.
  rho_beam = -1.

  np_beam = 0 ! initial # beam particles/ dt

  vosc = 6.0
  omega = 0.5
  sigma = 6.
  tpulse = 20.

  ! control
  nt =2 
  dt = 0.2
  eps = 1.
 restart = .false.
  vis_on = .false.
 ivis = 2
  mc_init = .false.
  mc_steps = 1000
  idump = 1
  ensemble = 1 /

