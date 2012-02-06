&pepcv

 ispecial = 3

 g = 4.
 nu = 0.

 rem_freq = 0
 !eps = 0.1
 eps = 0.4641588834
 thresh = 0.0D-04
 m_h = 0.1

 n_in = 20000

 dt = 0.5
 ts = 0.
 te = 25.

 dump_time = 0
 cp_time = 0
/

&calc_force_vortex
  force_law  = 62
  mac_select = 0
  theta2     = 0.0
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

