&pepc_nml

! total number of particles, this must conform to:
! np = 2 x n**2 x Lx x Ly / (GCD(Lx, Ly))**2, n = 1, 2, 3, ...
np = 5000000

! output particles every ... steps (no output if 0)
pdump = 0

! output fields every ... steps (no output if 0)
fdump = 1

! output checkpoints every ... steps (no output if 0)
cdump = 1

! BH opening angle
! theta = 0.3

! periodicity
t_lattice_1 = 1.0  0.0  0.0
t_lattice_2 = 0.0  125.0  0.0
t_lattice_3 = 0.0  0.0  1.0

periodicity = .false.  .true.  .false.

mirror_box_layers = 1
/

&time_nml

! time at end of simulation
te = 450.0

! number of steps from t = 0 to t = te
nsteps = 4500

! resume from this step
nresume = 0

/

&physics_nml

! ratio of mi / me
m_ratio = 100.0

! ratio of Ti / Te
T_ratio = 1.0

! static magnetic field
B0 = 1.0

! extent of the plasma slab
l_plasma = 50.0, 125.0, 0.0

/

&field_grid_nml

! nx x ny grid points
n = 512 640

! shifted right by 50
offset = -50.0 0.0

! physical space spanned by the grid
extent = 100.0 125.0

/

&calc_force_coulomb

! 3D coulomb
force_law  = 2

! BH-mac
mac_select = 0

! theta = 0.0
theta2     = 0.09

! eps squared
eps2 = 0.04

include_far_field_if_periodic = .false.
/

&libpepc
 debug_level = 0
 np_mult = -750

! Choose sorting routine and load balancing
! 0: no load balancing, 1: load balancing
 weighted = 1

! type of space-filling curve, 0=Z-curve, 1=Hilbert-curve
 curve_type = 1
 
 interaction_list_length_factor = 8
/

&walk_para_pthreads
 num_walk_threads = 48
 max_particles_per_thread = 100
/

