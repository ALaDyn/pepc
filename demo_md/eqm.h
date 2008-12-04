&pepcdata
  db_level = 2
  np_mult = -60
  fetch_mult = 4

! particles

! JUMP upper bound for particles per 4 processors
!   ne = 1000000 ! nearly exact, 346700 does not work...
!   ni = 1000000

! BG/P SMP-Mode upper bound for particles per 4 processors
!   ne = 126000 ! 128000 does not work...
!   ni = 126000

! BG/P VN-Mode upper bound for particles per 4 cores (!)
!   ne = 126000 !  approx., 128000 does not work...
!   ni = 126000

! BG/L CO-Mode upper bound for particles per 4 processors
!   ne = 54480 ! exact, 54482 does not work...
!   ni = 54480

! BG/L VN-Mode upper bound for particles per 4 cores (!)
!   ne = 40758  ! exact, 40760 does not work...
!   ni = 40758

 ne = 12800000
 ni = 12800000

!   ne = 8
!   ni = 8
!   ne = 16
!   ni = 16
!   ne = 128
!   ni = 128
!   ne = 256
!   ni = 256
!   ne = 512
!   ni = 512
!   ne = 1024
!   ni = 1024
!   ne = 2048
!   ni = 2048
!   ne = 4096
!   ni = 4096
!   ne = 8192
!   ni = 8192
!   ne = 16384
!   ni = 16384
!   ne = 32768
!   ni = 32768
!  ne = 65536
!  ni = 65536
!  ne = 131072
!  ni = 131072
!  ne = 262144
!  ni = 262144
!  ne = 524288
!  ni = 524288
!  ne = 1048576
!  ni = 1048576
!  ne = 2097152
!  ni = 2097152
!  ne = 4194304
!  ni = 4194304
!  ne = 8388608
!  ni = 8388608
!  ne = 16777216
!  ni = 16777216
!  ne = 33554432
!  ni = 33554432
!  ne = 67108864
!  ni = 67108864

 system_config = 2  ! set up plasma target
 ispecial = 3
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
  xl =1  ! graphics box size
  yl =1
  zl =1

  ! control
  nt = 2
  dt = 0.01
  eps = 1.
  restart = .false.
  vis_on = .false.
 ivis = 2
 ivis_domains = 5000
 ivis_fields = 5000
  idump = 1
  iprot=1
  particle_bcs = 1
  scheme = 1 /
