&pepcdata

 db_level = 2
 np_mult = -300
 num_walk_threads = 4
 max_particles_per_thread = 2000

! number of particles, here ions
 ne = 1000
 ni = 1000

! initial particle distribution
! 1 homogen, 2: one sphere, 3: two spheres, 4: Plummer (core cut)
 ispecial = 1
! ispecial = -1: reload particle positions from mpiio-timestamp #itime_in
itme_in = 5

! number of timesteps
 nt = 5
 
! fmm-periodicity framework
! lattice basis vectors
  t_lattice_1 = 1.0   0.0   0.0
  t_lattice_2 = 0.0   1.0   0.0
  t_lattice_3 = 0.0   0.0   1.0
! periodicity in x-, y-, and z-direction
!  periodicity = .true.  .true.  .true.
! extrinsic-to-intrinsic correction
  do_extrinsic_correction = .false.

! Choose sorting routine and load balancing
! 0: no load balancing, 1: load balancing
 weighted = 1

! determies the particle dump interval
! 0: never write anything
! n: each n-th step, plus first and last step
 idump = 0
! dito for vtk, binary and checkpoint-output
 idump_vtk = 1
 idump_binary = 0
 idump_checkpoint = 5
  /
