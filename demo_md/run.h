 &pepcdata
  db_level = 0
! particles
  ne = 100 
  ni = 100

 system_config = 1  ! set up plasma target
 !   initial_config=7   ! hollow sphere
   target_geometry = 0         ! random disc
 !   initial_config=3   ! wire
 !  initial_config = 0         ! rectangular slab
  !  initial_config = 10     ! read from parts_all.in

! physics stuff

  theta = 0.6
  mac=0
  Te_keV = 0.5 ! Temperatures in keV
  Ti_keV =0.1 
  mass_ratio = 2000.
  q_factor = 1.
  coulomb = .true.
  lenjones = .false.
  bond_const = 2.e-3
  r_sphere = 4
  x_plasma = 1.    ! plasma disc thickness/ wire length
  y_plasma = 1.     ! plasma width (slab target)
  z_plasma = 1.     ! plasma width (slab target)
  xl = 3  ! graphics box size
  yl =8 
  zl =8 

  ! control
  nt =10
  dt = 0.5
  eps = 0.
  restart = .false.
  vis_on = .false.
 ivis = 2
 ivis_domains = 5000 
 ivis_fields = 5000
  idump = 1
  iprot=1
  particle_bcs = 1
  scheme = 1 /

