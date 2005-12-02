 &pepcdata
! nmerge = 1
 mac=0
! particles
  ne = 4200
  ni = 4200 

 plasma_config = 1  ! set up plasma target
 target_geometry = 1   ! sphere

! physics stuff

  theta = 0.7
  Te_keV = 0.5 ! Temperatures in keV
  Ti_keV =0. 
  mass_ratio = 500.
  q_factor = 1.
  coulomb = .true.
  lenjones = .false.
  bond_const = 2.e-3
  r_sphere = 4
  x_plasma = 1.    ! plasma disc thickness/ wire length
  y_plasma = 2.     ! plasma width (slab target)
  z_plasma = 2.     ! plasma width (slab target)
  xl = 5  ! graphics box size
  yl =5 
  zl =5 


! beam
 beam_config_in=3  ! uniform, sinusoid
 
  np_beam = 0 ! initial # beam particles/ dt

  vosc = 0.1
  omega = 0.5
  sigma = 6.
  tpulse = 20.
  lambda = 1.0   ! Wavelength in microns

  ! control
  nt =500
  dt = 0.3
  eps = 2.5
 restart = .false.
  vis_on = .false.
 steering = .false.
 ivis = 5
 ivis_fields = 5000
  mc_init = .false.
  mc_steps = 1000
  idump = 4000
  iprot=2
  itrack=300
  particle_bcs = 1
  scheme = 1 /

