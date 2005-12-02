!  Input deck for billiards example
!  Forces switched off; reflective container walls

&pepcdata
! particles
  ne =0
  ni = 200 

 plasma_config = 1
! target_geometry = 0 ! slab
! target_geometry = 1   ! sphere
! target_geometry = 2  ! random disc
! target_geometry = 3  ! wire
! target_geometry = 4  ! ellipsoid
! target_geometry = 5  ! wedge
! target_geometry = 6  ! hemisphere
! target_geometry = 7  ! hollow sphere
target_geometry = 8  ! hollow hemisphere


! physics stuff

  theta = 0.5
  Te_keV = 0.1 ! Temperatures in keV
  Ti_keV =0.5 
  mass_ratio = 100.
  q_factor = 1.
  coulomb = .false.
  lenjones = .false.
  bond_const = 2.e-5
  r_sphere = 5
  x_plasma = 1    ! plasma disc thickness/ wire length
  y_plasma = 5.     ! plasma width (slab target)
  z_plasma = 3.     ! plasma width (slab target)
  xl = 5  ! graphics box size
  yl =5 
  zl =3 


! beam
   beam_config_in = 0 ! beam off

  vosc = 6.0
  omega = 0.5
  sigma = 6.
  tpulse = 20.
  lambda = 1.0   ! Wavelength in microns

  ! control
  nt =2000
  dt = 1
  eps = 0.2
 restart = .false.
  vis_on = .true.
 ivis =  2
 ivis_fields = 500
  mc_init = .false.
  mc_steps = 1000
  idump = 20000
  iprot=10
  itrack=300
  particle_bcs = 2
  ncpu_merge = 1
  scheme = 1 /

