
module module_physvars
  use treetypes
  use tree_walk_pthreads
  use module_fmm_framework
  implicit none

  real, parameter :: pi=3.141592654


  real, allocatable ::  rhoe_loc(:,:,:), rhoi_loc(:,:,:), & !< field arrays for time-averages
                          ex_ave(:), ey_ave(:,:), ez_ave(:)      !< time-ave field 
  real, allocatable ::  rhoi(:,:,:), rhoe(:,:,:) !< densities
  real, allocatable ::  rhoi2d_loc(:,:), rhoe2d_loc(:,:) !< 2D densities
  real*8, allocatable ::  ex2d(:,:), ey2d(:,:), pot2d(:,:) !< 2D fields
  real, allocatable ::  g_ele(:,:,:), g_ion(:,:,:)  !< particle counts (weights for T)
  real, allocatable ::  ex_loc(:,:,:), ey_loc(:,:,:), ez_loc(:,:,:)  !< E-field 
  real, allocatable ::  bx_loc(:,:,:), by_loc(:,:,:), bz_loc(:,:,:)  ! B-field 
  real, allocatable ::  jxe_loc(:,:,:), jye_loc(:,:,:), jze_loc(:,:,:)  !< elec current
  real, allocatable ::  Te_loc(:,:,:), Ti_loc(:,:,:) !< local temps

  real, allocatable :: xslice(:),  yslice(:),  zslice(:), &     !< Rezoning slice
                      uxslice(:), uyslice(:), uzslice(:), &     !< velocity
                              qslice(:),  mslice(:)    !< charge and mass
  real, allocatable :: vbuffer(:,:), vbuf_local(:,:) !< vis buffers

  real, allocatable :: rho_helm(:)  !< Helmholtz density
  complex, allocatable :: Az_helm(:)   !< Helmholtz vector potential

  integer :: choose_sort, weighted, choose_build
  type(t_calc_force_params) ::cf_par

  !  physics, target data

  integer, parameter :: maxlayers=10
  integer :: ni, ne       !<  # ions, electrons
  integer :: nep, nip     !< # particles/electrons/ions per process

  real :: xl, yl, zl     !< simulation box dimensions
  real :: vte=0., vti=0.       !< electron, ion thermal velocities; set vte=1 to choose Debye unit system
  real :: Te_keV, Ti_keV !< electron, ion emperatures in keV
  real :: T_scale = 1       !< factor for rescaling Te after restart 
  real :: force_const    !< force constant depending on unit system
  real :: eps = 1.       !< potential/force law cutoff
  integer :: force_law = 3 !< force law  0 = no interactions; 2 = 2D Coulomb; 3 = 3D Coulomb
  integer :: idim=3      !< # dimensions (velocity and position updates - indep of force_law)
  real :: bond_const=0.  !< bonding force constant for ion crystal
  real :: mass_ratio=1836   !< ion:electron mass ratio
  real :: mass_proton=1836  !< proton:electron mass ration
  real :: qe, qi         !< electron, ion charge
  real :: mass_e, mass_i   !< electron, ion mass
  real :: r_sphere       !< initial radius of plasma sphere
  real :: x_plasma       !< initial plasma length (slab or disc targets)
  real :: y_plasma       !< initial plasma y-width (slab)
  real :: z_plasma       !< initial plasma z-width (slab)
  real :: plasma_centre(3) !< vector defining centre of plasma target
  real :: displace(3)=0.    !< particle displacement vector for restart or 2nd target
  integer :: n_layer(maxlayers)=0   !< Additional target parameters
  real, dimension(maxlayers) :: rho_layer, Zion_layer, mratio_layer, x_layer, y_layer, z_layer, r_layer

  real :: x_crit         !< critical surface
  real :: x_offset       !< coordinate offset
  real :: z_offset       !< coordinate offset
  real :: rho0           !< electron density (1)
  real :: rho_track      !< tracking density for x_crit (/nc)
  real :: rho_upper=1.0  !< shelf/profile density above x_crit (/nc)
  real :: Vplas          !< plasma volume
  real :: Aplas          !< plasma area
  real :: Qplas          !< plasma charge
  real :: a_ii           !< mean ion spacing
  real :: a_ee           !< mean electron spacing
  real :: r_neighbour    !< nearest-neighbour search radius
  real :: delta_mc       !< step size for MC config
  real :: uthresh        !< velocity (u^2) threshold for vis_parts
  real :: rho_min        !< min density for exponential ramp
  real :: lolam          !< L/lambda density scale-length
  real :: Zion           !< Ionic charge
  real :: fnn            !< Near-neighbour factor
  real*8 :: Ukine          !< Electron kinetic energy
  real*8 :: Ukini          !< Ion kinetic energy
  real*8 :: Umagnetic           !< Magnetic energy
  real*8 :: Ubeam          !< Beam energy
  real :: tpert=0.1	         !< Temperature perturbation in transport test
  real :: kpert=2 	 !< normalised wave number Lx/lambda

  real :: nu_ei !< norm. collision frequency
  real :: sigma_e !< electrical conductivity
  real :: range_hot !< hot electron range
  real :: t_sat !< saturation time for n_hot
  real :: t_foil !< traversal time of hot e front

  ! particle beam stuff
  integer :: np_beam    !< # beam particles
  integer :: nb_pe    !< # beam particles per PE
  real :: x_beam        !< beam length
  real :: r_beam        !< beam radius
  real :: start_beam    !< starting position (x-axis)
  real :: u_beam, theta_beam, phi_beam        !< beam velocity and angles
  real :: rho_beam      !< beam density as fraction of plasma density (=1)
  real :: mass_beam     !< mass of beam particles
  real :: qeb           !< beam particle charge

 ! laser parameters
  real :: tpulse        !< pulse duration (in units of 1/omega_p)
  real :: sigma         !< 1/e pulse width (c/omega_p)
  real :: vosc          !< pump strength    (c)
  real :: omega         !< frequency  (omega_p)
  real :: lambda        !< laser wavelength
  real :: theta_inc     !< angle of incidence
  real :: focus(3)      !< centre of focal spot
  real :: tlaser        !< run time after laser switched on (1/omega_p)
  real :: elaser        !< deposited laser energy
  real :: propag_laser  !< distance travelled by laser after rezoning
  real :: intensity     !< normalised intensity = 0.5*vosc^2*omega^2
  real :: window_min    !< start of wakefield plasma
  real :: rezone_frac=0.75     !< Fraction of box to cross before rezoning switched on
  real :: glue_radius=2.e6 !< multiple of box size to catch escaping particles at
  real :: fpon_max, ampl_max !< max amplitudes

  integer :: nxh=100 !< 1D Helmholtz grid dimension
  real :: dxh !< HH grid spacing
  real :: xh_start=0.
  real :: xh_end=10.  !< Start and end points of Helmholtz grid
  integer :: ngav=100 !< Time-ave grid dimension
  real :: xgav_start=0.
  real :: xgav_end=10.  !< Limits for time-ave grid
  real :: xgav_pos(1:3)=(/0.,1.,2./)  !< Limits for time-ave grid - radial field positions

!  Variables needing 'copy' for tree routines

  integer :: npart_total  !< Total # particles (npart)
  integer :: np_local  !< Local # particles (npp)
  integer :: np_alloc  !< Memory allocation for particle arrays
  real :: np_mult=1.5   !< particle array safety margin
  integer :: fetch_mult=3 !< fetch array factor
  integer :: ipefile=20  !< Output file number for process-dep diagnostics
  logical :: dynamic_memalloc=.false. !< Switch for memory manager

!  Associated MPI stuff

  integer :: my_rank       !< Rank of current task
  integer :: n_cpu   !< # cpus used by program
  integer :: ifile_cpu    !< O/P stream
  integer :: debug_rank=0 !< CPU # for printed debug IO

! Control stuff
  logical :: launch = .true. !< Start simulation immediately after setup
  logical :: vis_on=.true.  !< online visualisation on/off
  logical :: netcdf=.false.  !< netcdf write off
  logical :: mc_init = .false. !< MC initialisation switch
  logical :: restart = .false.  !< Restart switch: config read from parts_all.in
  logical :: coulomb = .true.  !< Compute Coulomb forces
  logical :: bfields = .false.  !< Include magnetic fields
  logical :: bonds = .false. !< Include SHO bonding potential for ions
  logical :: lenjones = .false. !< Include short-range Lennard-Jones potential
  logical :: steering = .false.  !< VISIT steering switch
  logical :: ramp = .false.  !< profile-ramp switch
  logical :: te_perturb = .false.  !< Temperature perturbation switch

  integer :: mc_steps !< # MC steps
  integer :: plasma_config = 1  !< Switch for initial configuration (positions, velocities)
  integer :: target_geometry = 0  !< Geometry for plasma target
  integer :: layer_geometry = 0  !< Geometry for 2nd layer
  integer :: velocity_config = 1  !< Velocity distrib. (Maxw) 
  integer :: foam_geom(1:3) = (/1,1,1/)  !< Foam array dimensions
  integer :: beam_config_in = 0 !< Particle or laser beam switch including variations 
  integer :: beam_config = 0 !< Reduced switch for particle or laser beam 
  integer :: ispecial       !< Switch to select special electron configs 
  integer :: scheme = 1 !< Integrator scheme switch: 2-4= const. Te dynamics, 6=EM
  integer :: particle_bcs = 1 !< Particle BC switch: 1=open, 2=reflective, 3=periodic
  logical :: particle_wrap(3) = (/.false.,.false.,.false./)
  integer :: debug_level =1 !< Debug level for printed O/P
  integer :: debug_tree =0 !< Debug level for tree diagnostics O/P
  integer :: ncpu_merge=1  !< Restart control: -1= split data amoung all CPUs 
  integer :: new_label  !< Rezone parameter
  integer :: proton_label=1  ! Start label of protons (number in nlayer(x))
  integer :: nproton=0    ! # protons
  integer :: nslice            !< # particles in rezoning slice (determined in predef)
  integer :: np_error=0            !< # particles in error test sample
  
  character*4 csubme   !< Character string of data subdirectory 'data/peXXXX'

   real :: dt             !< timestep
   real :: trun           !< total run time including restarts
   real :: convert_fs     !< conversion factor from wp^-1 to fs
   real :: convert_mu     !< conversion factor from c/wp to microns
   real :: convert_keV     !< conversion factor from norm energy to keV/particle
   real :: convert_erg     !< conversion factor from norm energy to ergs
   real :: domain_cut      !< cutoff height for domain boxes
   integer :: nt, itime   !< # timesteps and current timestep
   integer :: itime_start !< restart time stamp
   integer :: idump=100      !< output frequency (timesteps)
   integer :: iprot=1       !< protocoll frequency
   integer :: ivis=5        !< frequency for particle shipping to VISIT
   integer :: ivis_fields=10    !<  frequency for field shipping to VISIT
   integer :: ivis_domains=10    !<  frequency for domain shipping to VISIT
   integer :: vis_select = 1  !< select switch for particles
   integer :: field_select(1:4) = 0. !< field selection switches for vis.
   integer :: itrack=1       !< frequency for computing ion density (tracking)
   integer :: navcycle     !< # timesteps in a laser cycle 
   integer :: ngx=50, ngy=10, ngz=10  !< Plot grid dimensions
   integer :: ncid         !< NetCDF id
   integer :: nbuf_max=10000     !< Max vis buffer size
   integer :: ndom_max=1000     !< Max # domains
   integer :: attrib_max = 22 !< max # attributes for vis buffer
   real :: ops_per_sec = 0.  !< Work load (interactions per sec)
   real :: work_tot = 0.  !< Integrated work load (npart*nlist)

   ! constrain
   real :: constrain_proof !< quality of getting crossing points
   real :: len_tripod      !< length of tripod
   integer :: struct_step  !< skip
   integer :: number_faces !< # faces in target container

 ! tree stuff
  real :: theta=0.5       !< Clumping parameter
  real :: force_tolerance=1.      !< Permitted error in force calculation
  integer :: nint_max = 2500
  integer :: mac = 0  !< MAC (default=BH)
  integer :: walk_scheme = 0  !< Asynch/Collective walk algorithm 
  integer :: balance = 1  !< Load balancing switch
  integer :: ifreeze = 1  !< Tree-freezing control
  integer :: curve_type = 0  !< Morton or Hilbert curve

  character*10 :: plasma_configs(0:50)= (/ &
       'no plasma ','rand-geom ','special   ','i+e slab  ','i+e sphere', &
       '          ','          ','imp sphere','          ','          ', &
       'ion overl ','slab+prot ','microdot  ','doped     ','micro-rear', &
       'slab+foam ','dish+prot ','          ','          ','          ', &
       'tincan    ','can+prot  ','tube+prot ','rev-can   ','          ', &
       '          ','          ','          ','          ','          ', & 
       '          ','          ','microd-2  ','          ','          ', &
       'microd-3  ','          ','          ','          ','          ', & 
       'foam      ','pipes     ','          ','          ','          ', &
       '          ','          ','          ','          ','          ', &
       '          '  /)

  character*10 :: geometries(0:10)= (/ &
       'slab      ', & 
       'sphere    ','disc      ','wire      ','ellipsoid ','wedge     ', &
       'hemisphere','hollow sph','hollow hsp','          ','special   ' /)
  character*7 :: beam_configs(0:9)=(/ &
       'eqm    ','beam   ','i-beam ','laser-u','ES pond','LWFA   ', &
       'EMplane','EM pond',' dust  ','       ' /)
  character*7 :: schemes(1:8)=(/ &
       'const U','Te, Ti ','glob Te','loc  Te','Ti only','Full 3V','non-rel','2V  TE ' /)

!  Redundant variables retained for namelist compatibility
  real :: q_factor=1.0
  logical :: target_dup = .false.  

 namelist /pepcdata/ nep, nip, ne, ni, &
       theta, mac, nint_max, mass_ratio, mass_proton, rho0, &
       plasma_config, target_geometry, velocity_config, ispecial, &
       vte, Te_keV, Ti_keV, T_scale, Zion, &
       r_sphere, x_plasma, y_plasma, z_plasma, delta_mc, &
       n_layer, x_layer, y_layer, z_layer, r_layer, Zion_layer, rho_layer, mratio_layer, layer_geometry, &
       xl, yl, zl, displace, bond_const, fnn, rho_min, lolam, &
       beam_config_in, np_beam, np_error, &
       r_beam, u_beam, theta_beam, phi_beam, x_beam, start_beam, rho_beam, mass_beam, & 
       lambda, sigma, tpulse, vosc, omega, focus, x_offset,  z_offset, &
       nt, dt, mc_steps, idump, ivis, ivis_fields, ivis_domains, iprot, itrack, &
       ncpu_merge, ngx, ngy, ngz, choose_sort, weighted, choose_build, &
       vis_on, steering,  vis_select, mc_init, restart, scheme, particle_bcs, &
       coulomb,  bfields,  bonds, lenjones, target_dup, ramp, &
       debug_level, debug_tree, ncpu_merge, balance, curve_type, ifreeze, walk_scheme, &
       constrain_proof, len_tripod, struct_step, uthresh, domain_cut, glue_radius, &
       debug_rank, np_mult, fetch_mult,nbuf_max, te_perturb, tpert, kpert, &
       q_factor, netcdf, launch, foam_geom, force_tolerance, dynamic_memalloc, &
       eps, force_const, idim, force_law, &
       xh_start, xh_end, nxh, xgav_start, xgav_end, ngav, xgav_pos, &
       num_walk_threads, particle_wrap, periodicity,t_lattice_1,t_lattice_2,t_lattice_3
!       t1, t2, t3, periodicity, do_extrinsic_correction              ! periodicity config
 
end module module_physvars

