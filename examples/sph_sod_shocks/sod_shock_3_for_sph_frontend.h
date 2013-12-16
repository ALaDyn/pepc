&pepcsph

! number of particles, here ions
  ni = 1000

! initial particle distribution
! ispecial = 15 for 1D shock problem 2
ispecial = 15


! determies the particle dump interval
! 0: never write anything
! n: each n-th step, plus first and last step
! idump = 0
! dito for vtk, binary and checkpoint-output
! idump_vtk = 1
! idump_binary = 0

! vtk dump every dump_time steps
 dump_time = 100

! checkpoint files every cp_time steps
! cp_time = 1


/



&libpepc

 num_threads = 8

/


&walk_para_pthreads

! max_particles_per_thread = 2000


/


&calc_force_nearestneighbour

/
