!  ================================
!
!         SETUP
!
!   $Revision$
!
!     Initialise physics constants and 
!      simulation variables
!
!  ================================


subroutine setup
  use physvars
  implicit none

  integer :: k, npb_pe, ne_rest, ni_rest

  include 'namelist.h' 

  !  Default input set

  ! switches

  plasma_config = 1  ! plasma target
  target_geometry = 1         ! random sphere


  ! particles
  nep = 0 ! # plasma electrons per PE
  nip = 0
  ne = 100  ! Total # plasma electrons
  ni = 100  ! total # plasma ions
  mc_steps = 10

  xl = 2
  yl = 2
  zl = 2
  
  ! physics stuff
  force_const = 1./3.
  bond_const = 0.1
  rho0 = 1.0
  mac = 0        ! Multipole acceptance criterion (BH by default)
  theta = 0.5
  Te_keV = 1.
  Ti_keV = Te_keV/10.
  mass_ratio = 10.
  Zion = 1.
  uthresh = -1.

  r_sphere = 0.5
  x_plasma = 0.1    ! plasma disc thickness (2) or wire length (3)
  y_plasma = 1.     ! plasma width (slab target)
  z_plasma = 1.     ! plasma width (slab target)
  eps = 0.1
  fnn = 5  ! Neighbour search radius multiplier (x a_ii)
  delta_mc = r_sphere/5.
  displace(1:3) = (/0.,0.,0./)
  ! beam

  !  beam_config = 1  ! fixed beam, initialised at start
  ! beam_config = 2  ! user-controlled, real-time particle source
  beam_config_in = 0 ! beam off
  ! beam_config = 4  ! laser ponderomotive force

  r_beam = 0.15
  u_beam = 0.2
  theta_beam = 0.
  phi_beam = 0.3
  x_beam = .04
  start_beam = 0.4
  rho_beam = -5.
  mass_beam = 3.
  np_beam = 0 ! initial # beam particles/ dt

  ! laser stuff
  sigma = 1.
  tpulse = 10.
  vosc = 0.1
  omega = 0.5
  theta_inc = 0.
  x_offset = 0.
  z_offset = 0.
  tlaser = 0.    ! time since laser switched on
  lambda = 1.0
  rho_track = 1.5
  lolam = 0.1
  rho_min = 0.1

  ! control
  nt = 600
  dt = 0.2
  trun = 0.
  ivis = 1
  ivis_fields = 1
  ivis_domains = 1
  start_step = 0
  itrack = 10
  domain_cut = zl

  ngx = 25   ! Grid size for plots
  ngy = 25
  ngz = 25
  ! constrain
!!$    len_tripod = .001
  constrain_proof = .001
  struct_step = 0

  ! Read actual inputs from namelist file
  open(10,file='run.h')
  read (10,NML=pepcdata)

 ! # particles in primary target component
  npart_total = ni+ne

! Adjust local numbers if total non-multiple of # PEs
  if (my_rank==0) then 
     ne_rest = mod(ne,n_cpu)
     ni_rest = mod(ni,n_cpu)
  else
     ne_rest = 0
     ni_rest = 0
  endif

  np_local = npart_total/n_cpu + ne_rest + ni_rest  ! total # particles on this CPU
  nep = ne/n_cpu + ne_rest ! local # electrons and ions - may be adjusted later
  nip = ni/n_cpu + ni_rest

  new_label = npart_total  ! Rezone label

  if (target_dup .or. scheme==5) np_mult = np_mult*2  ! double up particle array size if multi-target or ion config mode 

  if (.not. restart) npart_total = npart_total + 2*SUM(n_layer) ! Include particles from remaining layers for array setup.

  beam_config=mod(beam_config_in,10)  ! derived config from s,p variations


  lolam = lolam*2.*pi/omega  ! normalise scale-length

end subroutine setup




