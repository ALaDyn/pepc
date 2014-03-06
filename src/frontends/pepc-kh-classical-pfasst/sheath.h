&pepc_nml

! total number of particles
np = 0

! output particles every ... steps (no output if 0)
pdump = 0

! output fields every ... steps (no output if 0)
fdump = 10

! output checkpoints every ... steps (no output if 0)
cdump = 0

! periodicity
mirror_box_layers = 10

do_extrinsic_correction = .false.

! WM_BORIS_SDC   = 1
! WM_BORIS_MLSDC = 2
! WM_BORIS       = 3
! WM_BENEDIKT    = 4
workingmode = 1

/

&time_nml

! time at end of simulation
te = 500

! number of steps from t = 0 to t = te
nsteps = 150

! resume from this step
nresume = 0
/

&physics_nml

! ratio of mi / me
m_ratio = 16.0

! ratio of Ti / Te
T_ratio = 1.0

! static magnetic field
B0 = 2.0

! extent of the plasma slab
l_plasma = 51.2 204.8 0.0

! halfwidth of the shear layer 
shear_halfwidth = 2.0

! strength of the shear (in units of wci x shear_halfwidth)
shear_strength = 1.0

! number of ions in initial configuration
ni = 20000
/

&field_grid_nml

! nx x ny grid points
n = 64 128

! shifted right by 50
offset = 0.0 0.0

! physical space spanned by the grid
extent = 51.2 204.8
/

&calc_force_log2d

mac_select = 0

theta2 = 0.36

eps2 = 0.05

include_far_field_if_periodic = .false.
/

&libpepc
 debug_level = 0
 np_mult = -40

! Choose sorting routine and load balancing
! 0: no load balancing, 1: load balancing
 weighted = 1

! type of space-filling curve, 0=Z-curve, 1=Hilbert-curve
 curve_type = 1
 
 interaction_list_length_factor = 8
 num_threads = 8
/

&walk_para_pthreads
 max_particles_per_thread = 100
/


&pfasst
  niter   = 1
  nlevels = 1
  
  ! level parameters (one entry per level and line, finer levels are more left, coarser levels more right)
  nsweeps = 1 1
  nnodes  = 3 3
  ! true  --> evaluate forces directly: O(N**2) but exact
  ! false --> evaluate forces with PEPC: O(N log N) but multipole approximation
  directforce = .false. .false.
  ! 1st digit: 0 - no external field, 1 - full external field
  ! 2nd digit: 0 - no internal field, 1 - full internal field
  feval_mode = 11 10
  ! multipole acceptance parameter
  theta   = 0.3 0.3
  ! end of level parameters
  
  res_tol = 0.
  
  num_space_instances = 1
  color_space_div     = .true.
  color_time_div      = .true.
  
  ! create fort.601 files with timings
  echo_timings        = .false.
  ! unused
  echo_errors         = .true.
/

