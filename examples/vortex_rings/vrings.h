&pepcv

 ispecial = 5

 g = 2
 rmax = 0.4
 r_torus = 1.
 nc = 20
 torus_offset = 1. 0. 0.

 nu = 0.00555

 rem_freq = 0
 eps = 0.2
 thresh = 0.1D-08
 m_h = 0.02581

 n_in = 1000000

 dt = 0.5
 ts = 0.
 te = 500

 dump_time = 1
 cp_time = 10
/

&calc_force_vortex
  force_law  = 22
  mac_select = 0
  theta2     = 0.16
/

&libpepc
 debug_level = 0
 weighted = 0
/

&walk_para_pthreads
 num_walk_threads         = 4
 max_particles_per_thread = 2000
 defer_list_length_factor = 4
/

&calc_force_vortex

  
/

