!  ================================
!
!         SETUP
!
!   $Revision$
!
!     Initialise constants and 
!      simulation variables
!
!  ================================


subroutine setup
  use physvars
  use treevars
  implicit none

  integer :: k, npb_pe
  real :: Qplas, Aplas



  namelist /pepcdata/ nep, nip, ne, ni, &
       theta, mass_ratio, q_factor, eps, &
       plasma_config, target_geometry, ispecial, &
       Te_keV, Ti_keV, T_scale, &
       r_sphere, x_plasma, y_plasma, z_plasma, delta_mc, &
       xl, yl, zl, displace, bond_const, fnn, rho_min, lolam, &
       beam_config_in, np_beam, idim, &
       r_beam, u_beam, theta_beam, phi_beam, x_beam, start_beam, rho_beam, mass_beam, & 
       lambda, sigma, tpulse, vosc, omega, focus, x_offset,  z_offset, &
       nt, dt, mc_steps, idump, ivis, ivis_fields, ivis_domains, iprot, itrack, nmerge, ngx, ngy, ngz, &
       vis_on, steering, domain_debug,  mc_init, restart, scheme, particle_bcs, &
       load_balance, walk_balance, walk_debug, force_debug, prefetch_debug, &
       dump_tree, perf_anal, coulomb, bonds, lenjones, target_dup, ramp, &
       prefetch, walk_summary, branch_debug, tree_debug, build_debug, &
       db_level, &
       constrain_proof, len_tripod, use_multipoles, struct_step, uthresh, bfield_on

  !  Default input set

  ! switches
  tree_debug = .false.
  domain_debug = .false.
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
  theta = 0.8
  Te_keV = 1.
  Ti_keV = Te_keV/10.
  mass_ratio = 10.
  q_factor = 1.
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
  np_beam = 64 ! initial # beam particles/ dt

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
  itime_start = 0
  itrack = 10
  sumprefetches = 0

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


  ! Derived parameters

  if (nep > 0) then
     ! particles specified per processor in input file
     ne = nep*n_cpu  ! total # electrons
     ni = nip*n_cpu  ! total # ions
  else
     ! total # particles specified in input file 
     nep = ne/n_cpu
     nip = ni/n_cpu
     npp = nep+nip
  endif

  !  npb_pe = np_beam/num_pe
  if (.not. restart) then
     if (nip*n_cpu /= ni .or. nep*n_cpu /= ne ) then
        ne = (nep+1)*n_cpu
        ni = (nip+1)*n_cpu

        if (my_rank==0) then
           write(*,'(//a//)') '*** Warning: number each particle species (ne, ni) must be divisible by # processors ***'
           write(*,'(a,i6)') '*** Resetting to ',ne,ni
        endif
        !     else if (npb_pe*n_cpu /= np_beam) then 
        !        np_beam = (npb_pe+1)*n_cpu
        !        if (my_rank==0) then
        !           write(*,'(//a)') '*** Warning: number of beam particles must be divisible by # processors ***'
        !           write(*,'(a,i6)') '*** Resetting to ',np_beam
        !        endif

     endif
  endif

  npart_total = ni+ne
  npp = npart_total/n_cpu  ! total # particles per processor
  new_label = npart_total  ! Rezone label

  geometry: select case(target_geometry)

  case(0) ! slab
     Vplas = x_plasma * y_plasma * z_plasma  ! plasma volume
     Aplas = x_plasma * y_plasma ! plasma area
     focus = (/xl / 2 + x_offset, yl / 2., zl / 2./) ! Centre of laser focal spot
     plasma_centre =  (/xl / 2., yl / 2., zl / 2./) ! Centre of plasma
     number_faces = 6

  case(1) ! sphere
     Vplas = 4 * pi * r_sphere**3 / 3.
     Aplas = pi*r_sphere**2
     focus = (/xl / 2. - r_sphere, yl / 2., zl / 2./) ! Centre of laser focal spot
     plasma_centre = (/xl / 2., yl / 2., zl / 2./) ! Centre of plasma
     number_faces = 1

  case(2) ! disc
     Vplas = pi * r_sphere**2 * x_plasma
     Aplas = x_plasma*y_plasma
     focus = (/xl / 2. - x_plasma / 2., yl / 2., zl / 2./) ! Centre of laser focal spot
     plasma_centre = (/xl / 2., yl / 2., zl / 2./) ! Centre of plasma        
     number_faces = 3

  case(3) ! wire
     Vplas = pi * r_sphere**2 * z_plasma
     Aplas = pi*r_sphere**2
     focus = (/xl / 2. - r_sphere + x_offset, yl / 2., zl / 2. + z_offset/) ! Centre of laser focal spot
     plasma_centre = (/xl / 2., yl / 2., zl / 2./) ! Centre of plasma
     number_faces = 3

  case(4) ! ellipsoid
     Vplas = 4 * pi * x_plasma * y_plasma * z_plasma / 3.
     Aplas = pi*x_plasma*y_plasma*2
     focus = (/xl / 2. - x_plasma * r_sphere, yl / 2., zl / 2./) ! Centre of laser focal spot
     plasma_centre = (/xl / 2., yl / 2., zl / 2./) 
     number_faces = 1

  case(5) ! wedge
     Vplas = .5 * x_plasma * y_plasma * z_plasma
     Aplas = .5*x_plasma*y_plasma
     focus = (/xl / 2. - x_plasma / 2., yl / 2., zl / 2./)
     plasma_centre = (/xl / 2., yl / 2., zl / 2./)
     number_faces = 5

  case(6) ! hemisphere
     Vplas = 4 * pi * r_sphere**3 / 6.
     Aplas = pi*r_sphere**2/2.
     focus = (/xl / 2. - r_sphere / 2., yl / 2., zl / 2./)
     plasma_centre = (/xl / 2., yl / 2., zl / 2./)
     number_faces = 2

  case(7) ! hollow sphere
     Vplas = (4 * pi / 3.) * (r_sphere**3 - (r_sphere - x_plasma)**3)
     Aplas = pi*(r_sphere**2-(r_sphere-x_plasma)**2)
     focus = (/xl / 2. - r_sphere / 2., yl / 2., zl / 2./)
     plasma_centre = (/xl / 2., yl / 2., zl / 2./)
     number_faces = 2

  case(8) ! hollow hemisphere
     Vplas = (4 * pi / 6.) * (r_sphere**3 - (r_sphere - x_plasma)**3)
     Aplas = pi/2.*(r_sphere**2-(r_sphere-x_plasma)**2)
     focus = (/xl / 2. - r_sphere / 2., yl / 2., zl / 2./)
     plasma_centre = (/xl / 2., yl / 2., zl / 2./)
     number_faces = 3

  end select geometry

  window_min = plasma_centre(1) - x_plasma/2.
  propag_laser=focus(1)

  if (plasma_config==2) then ! Electrons only user-defined config (special_start)
     Vplas = x_plasma * y_plasma * z_plasma
     Aplas = x_plasma * y_plasma
     focus = (/xl /4., yl / 2., zl / 2./) ! Centre of laser focal spot
     plasma_centre =  (/xl / 2., yl / 2., zl / 2./) ! Centre of plasma
     number_faces = 6  
  endif

  vte = sqrt(Te_keV/511.)  ! convert from keV to /c


  if (idim==3) then    ! 3D
     if (ne > 0) then
        qe = -Vplas*rho0/ne
        qi = -qe*q_factor
        mass_e = -qe
        Qplas = abs(qe)*ne
     else
        ! ions only
        qi = Vplas*rho0/ni
        qe = -qi
        mass_e = qi
        Qplas = abs(qi)*ni
     endif
     a_ii = (Vplas/ni)**(1./3.)

  else                ! 2D - use areal density instead of volume
     if (ne > 0) then
        qe = -Aplas*rho0/ne
        qi = -qe*q_factor
        mass_e = -qe
        Qplas = abs(qe)*ne
     else
        ! ions only
        qi = Aplas*rho0/ni
        qe = -qi
        mass_e = qi
        Qplas = abs(qi)*ni
     endif
     a_ii = (Aplas/ni)**(1./2.) ! assume 2D slab mode
  endif

  vti = sqrt(Ti_keV/511./mass_ratio)
  mass_i = mass_e*mass_ratio
  convert_fs = 10.*omega*lambda/(6*pi)     ! convert from wp^-1 to fs
  convert_mu = omega/2./pi*lambda          ! convert from c/wp to microns
  lolam = lolam*2.*pi/omega  ! normalise scale-length
  convert_keV = 2./3./Qplas*511     ! convert from code energy units to keV/particle

  r_neighbour = fnn*a_ii  ! Nearest neighbour search radius

  navcycle = 2*pi/dt/omega  ! # timesteps in a laser cycle

  nu_ei = 1./40./pi*a_ii**3/vte/eps**2  ! collision frequency (fit to Okuda & Birdsall)

  sigma_e = 1./nu_ei   ! Spitzer conductivity

  intensity = 0.2*vosc**2*omega**2  ! normalised laser intensity

  beam_config=mod(beam_config_in,10)  ! derived config from s,p variations

!  if (scheme > 1 .and. beam_config >= 3) then
!     if (my_rank==0) write(*,*) 'Constant-Te mode: turning beam off'
!     beam_config=0
!  endif

  if ( beam_config >=3 .and. beam_config<=6) then
     rho_beam= vosc
     r_beam=sigma
   
  else if (scheme == 5) then
  ! ion crystal eqm mode
     r_beam = a_ii
     u_beam = Ti_keV
     rho_beam = log(bond_const)
  endif


! phys -> tree variables (eventually goes into interface)
  me = my_rank
  num_pe = n_cpu
  ipefile = ifile_cpu
  npart = npart_total

end subroutine setup




