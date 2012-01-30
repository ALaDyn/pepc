&pepcv

 ispecial = 2

 g = 23.8
 rmax = 0.4
 r_torus = 1.
 nc = 15
 nphi = 244
 torus_offset = 2.8 0. 1.

 nu = 0.00555

 rem_freq = 2
 eps = 0.2
 thresh = 0.5D-06
 m_h = 0.02581

 n_in = 1000000

 dt = 0.02
 ts = 0.
 te = 30.

 dump_time = 5
 cp_time = 10
/

&calc_force_vortex
  force_law  = 22
  mac_select = 0
  theta2     = 0.16
/

&libpepc
 debug_level = 0
/

&walk_para_pthreads
 num_walk_threads         = 4
 max_particles_per_thread = 2000
 defer_list_length_factor = 1
/

&calc_force_vortex

  
/

