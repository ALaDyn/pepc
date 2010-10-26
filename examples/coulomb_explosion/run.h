&pepcdata

 db_level = 0
 np_mult = -45
theta=0.5
eps=0.01
ngx=200

! number of particles, here ions
! ne = 200000
! ni = 200000
! ne = 12800000
! ni = 12800000
! ni = 12800000
ne=0
ni=50000

! sphere radius
  r_sphere=1.0

! initial particle distribution
  ! 1 homogen, 2: one sphere, 3: two spheres, 4: Plummer (core cut)
 ispecial = 11

! number of timesteps
 nt = 50
dt=0.03
xl=8
! communication scheme
! 0: point-to-point, 1: collectives
 scheme = 0

! Choose sorting routine and load balancing
! 0: no load balacing, 1: load balacing
 weighted = 1
! 0: Simple sort, 1: pbalsort, 2: sl_sort_part, 3: sl_sort_keys
 choose_sort = 3

! Choose tree build routine
! 0: original, 1: optimized
 choose_build = 0

! determines the particle dump interval
! 0: never write anything
! n: each n-th step, plus first and last step
 idump = 1  /
