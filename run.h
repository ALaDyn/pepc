&pepcdata

 db_level = 2
 np_mult = -90
! number of particles, here ions
 ni = 400000

! initial particle distribution
  ! 1 homogen, 2: one sphere, 3: two spheres, 4: Plummer (core cut)
 ispecial = 4 

! number of timesteps
 nt = 3
 
! communication scheme
! 0: point-to-point, 1: collectives
 scheme = 0

! determinates the particle dump interval
! 0: never write anything
! n: each n-th step, plus first and last step
 idump = 0  /
