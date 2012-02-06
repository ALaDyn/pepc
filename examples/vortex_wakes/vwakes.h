&pepcv

 ispecial = 5

 g = 2
 rmax = 0.4
 r_torus = 1.
 nc = 120
 torus_offset = 0.8 0. 0.

 nu = 0.000914

 rem_freq = 2
 eps = 0.2
 thresh = 0.1D-05
 m_h = 0.05236

 n_in = 1000000

 dt = 0.5
 ts = 0.
 te = 200

 dump_time = 1
 cp_time = 0
/

&calc_force_vortex
  force_law  = 22
  mac_select = 0
  theta2     = 0.36
/

&libpepc
 debug_level = 0
 weighted = 0
/

&walk_para_pthreads
 num_walk_threads         = 4
 max_particles_per_thread = 2000
 defer_list_length_factor = 8
/

&calc_force_vortex

  
/

