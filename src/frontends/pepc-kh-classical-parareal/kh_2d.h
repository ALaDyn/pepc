&pepc_nml

! total number of particles
np = 0

! output particles every ... steps (no output if 0)
pdump = 0

! output fields every ... steps (no output if 0)
fdump = 0

! output checkpoints every ... steps (no output if 0)
cdump = 1

! periodicity
mirror_box_layers = 7

do_extrinsic_correction = .false.
/

&time_nml

! size of time step
dt = 0.01

! time at beginning
tresume = 0.0

! number of steps from t = 0 to t = te
nsteps = 5

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
l_plasma = 73.0, 292.0, 0.0

! halfwidth of the shear layer
shear_halfwidth = 2.0

! strength of the shear (in units of wci x shear_halfwidth)
shear_strength = 1.0

! number of ions in initial configuration
ni = 250000
/

&field_grid_nml

! nx x ny grid points
n = 128 512

! shifted right by 50
offset = 0.0 0.0

! physical space spanned by the grid
extent = 73.0 292.0
/

&calc_force_log2d

mac_select = 0

theta2 = 0.09

eps2 = 0.04

include_far_field_if_periodic = .false.
/

&libpepc
 debug_level = 0
 np_mult = 8.0

! Choose sorting routine and load balancing
! 0: no load balancing, 1: load balancing
 weighted = 1

! type of space-filling curve, 0=Z-curve, 1=Hilbert-curve
 curve_type = 1

 interaction_list_length_factor = 8
 num_threads = 6
/

&walk_para_pthreads
 max_particles_per_thread = 100
/
