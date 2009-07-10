&pepcdata

! debug level, in general it should be set to 0
 db_level = 0

! memory allocation parameter, _might/will_ 
! depend on target architecture, in general they 
! should be set to np_mult = -30 and fetch_mult = 3

 np_mult = -30 
 fetch_mult = 3

! number of electrons and ions
 ne = 0
 ni = 1000

! initial particle distribution
! 1 homogen, 2: one sphere, 3: two spheres
 ispecial = 1 

! number of timesteps
 nt = 10
 
! communication scheme
! 0: point-to-point, 1: collectives
 scheme = 1 

! determinates the particle dump interval
! 0: never write anything
! n: each n-th step, plus first and last step
 idump = 2  /
