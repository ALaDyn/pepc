&pepcv

 ispecial = 3

 g = 4.
 nu = 0.

 rem_freq = 1
 thresh = 0.1D-04
 m_h = 0.1

 n_in = 1000

 dt = 0.5
 ts = 0.
 te = 10.

 dump_time = 1
 cp_time = 1
/

&calc_force_vortex
  force_law  = 2
  mac_select = 0
  theta2     = 0.36
/

&libpepc
 debug_level = 99
/

&walk_para_pthreads
 num_walk_threads         = 2
 max_particles_per_thread = 2000
 defer_list_length_factor = 8
/

&calc_force_vortex

  
/

