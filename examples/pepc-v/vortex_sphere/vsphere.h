&pepcv

 ispecial = 3

 g = 4.
 nu = 0.

 rem_freq = 0
 !eps = 0.1
 eps = 0.4641588834
 thresh = 0.0D-04
 m_h = 0.1

 n_in = 2000

 dt = 0.5
 ts = 0.
 te = 50.

 dump_time = 0
 cp_time = 0
/

&calc_force_vortex
  force_law  = 62
  mac_select = 1
  theta2     = 0.09
/

&libpepc
 debug_level = 0
 weighted = 0
 interaction_list_length_factor = 4
 num_threads = 4
/

&walk_para_pthreads
 max_particles_per_thread = 2000
/

&calc_force_vortex

  
/

