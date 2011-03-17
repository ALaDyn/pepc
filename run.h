&pepcdata

 db_level = 0
 np_mult = -45
 num_walk_threads = 3

! number of particles, here ions
 ne = 100000
 ni = 100000

! initial particle distribution
  ! 1 homogen, 2: one sphere, 3: two spheres, 4: Plummer (core cut)
 ispecial = 1

! number of timesteps
 nt = 100
 
! fmm-periodicity framework
! lattice basis vectors
  t_lattice_1 = 1.0   0.0   0.0
  t_lattice_2 = 0.0   1.0   0.0
  t_lattice_3 = 0.0   0.0   1.0
! periodicity in x-, y-, and z-direction
  periodicity = .true.  .true.  .true.
! extrinsic-to-intrinsic correction
  do_extrinsic_correction = .false.

! Choose sorting routine and load balancing
! 0: no load balancing, 1: load balancing
 weighted = 1
! 1: pbalsort, 2: sl_sort_part, 3: sl_sort_keys
 choose_sort = 3

! determinates the particle dump interval
! 0: never write anything
! n: each n-th step, plus first and last step
 idump = 0  /
