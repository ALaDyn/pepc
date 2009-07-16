&pepcdata

! number of particles, here ions
 ni = 1000

! initial particle distribution
! 1 homogen, 2: one sphere, 3: two spheres
 ispecial = 1 

! number of timesteps
 nt = 10
 
! communication scheme
! 0: point-to-point, 1: collectives
 scheme = 0

! determinates the particle dump interval
! 0: never write anything
! n: each n-th step, plus first and last step
 idump = 2  /
