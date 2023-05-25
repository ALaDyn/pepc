&pepcv

 ispecial = 6

 g = 1.0
 Lref = 1.0
 r_torus = 1.0
 torus_offset = 0. 0. 0.
 nv_on_Lref = 16
 input_itime = 0
 thresh = 2.332E-6

 nu = 0.002 

 n_in = 1000000

 ts = 0.
 te = 100. 

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
 interaction_list_length_factor = 4
 num_threads = 10
/

&walk_para_pthreads
 max_particles_per_thread = 15000
/

&calc_force_vortex

  
/

