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
  use utils

  implicit none
  integer :: k, npb_pe

 

  namelist /pepcdata/ nep, nip, ne, ni, &
       theta, mass_ratio, q_factor, eps, &
       initial_config, ispecial, &
       Te_keV, Ti_keV, T_scale, &
       r_sphere, x_plasma, y_plasma, z_plasma, delta_mc, &
       xl, yl, zl, displace, bond_const, fnn, rho_min, lolam, &
       beam_config, np_beam, &
       r_beam, u_beam, theta_beam, phi_beam, x_beam, start_beam, rho_beam, mass_beam, & 
       lambda, sigma, tpulse, vosc, omega, focus, x_offset,  z_offset, &
       nt, dt, mc_steps, idump, ivis, ivis_fields, iprot, itrack, nmerge, ngx, ngy, ngz, &
       vis_on, steering, domain_debug,  mc_init, restart, scheme, particle_bcs, &
       load_balance, walk_balance, walk_debug, force_debug, prefetch_debug, &
       dump_tree, perf_anal, coulomb, bonds, lenjones, target_dup, ramp, &
       prefetch, walk_summary, branch_debug, tree_debug, &
       constrain_proof, len_tripod, use_multipoles, struct_step, uthresh

  !  Default input set

  ! switches
  tree_debug = .false.
  domain_debug = .false.
  initial_config = 1         ! random sphere
  ! initial_config = 2         ! random disc
  !    initial_config = 3         ! rectangular slab
  !  initial_config = 10     ! read from parts_all.in

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
  beam_config = 0 ! beam off
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
  itime_start = 0
  itrack = 10

  ngx = 100   ! Grid size for plots
  ngy = 50
  ngz = 50
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
     ne = nep*num_pe  ! total # electrons
     ni = nip*num_pe  ! total # ions
  else
     ! total # particles specified in input file 
     nep = ne/num_pe
     nip = ni/num_pe
     npp = nep+nip
  endif

!  npb_pe = np_beam/num_pe
  if (.not. restart) then
     if (nip*num_pe /= ni .or. nep*num_pe /= ne ) then
        ne = (nep+1)*num_pe
        ni = (nip+1)*num_pe

        if (me==0) then
           write(*,'(//a//)') '*** Warning: number each particle species (ne, ni) must be divisible by # processors ***'
           write(*,'(a,i6)') '*** Resetting to ',ne,ni
        endif
!     else if (npb_pe*num_pe /= np_beam) then 
!        np_beam = (npb_pe+1)*num_pe
!        if (me==0) then
!           write(*,'(//a)') '*** Warning: number of beam particles must be divisible by # processors ***'
!           write(*,'(a,i6)') '*** Resetting to ',np_beam
!        endif

     endif
  endif

  npart = ni+ne
  npp = nep + nip  ! total # particles per processor
  new_label = npart  ! Rezone label

  geometry: select case(initial_config)

    case(0) ! slab
        Vplas = x_plasma * y_plasma * z_plasma
        focus = (/xl / 2., yl / 2., zl / 2./) ! Centre of laser focal spot
        plasma_centre =  (/xl / 2., yl / 2., zl / 2./) ! Centre of plasma
        number_faces = 6
     
    case(1) ! sphere
        Vplas = 4 * pi * r_sphere**3 / 3.
        
        focus = (/xl / 2. - r_sphere, yl / 2., zl / 2./) ! Centre of laser focal spot
        plasma_centre = (/xl / 2., yl / 2., zl / 2./) ! Centre of plasma
        number_faces = 1
        
    case(2) ! disc
        Vplas = pi * r_sphere**2 * x_plasma
        focus = (/xl / 2. - x_plasma / 2., yl / 2., zl / 2./) ! Centre of laser focal spot
        plasma_centre = (/xl / 2., yl / 2., zl / 2./) ! Centre of plasma        
        number_faces = 3

    case(3) ! wire
        Vplas = pi * r_sphere**2 * z_plasma
        focus = (/xl / 2. - r_sphere + x_offset, yl / 2., zl / 2. + z_offset/) ! Centre of laser focal spot
        plasma_centre = (/xl / 2., yl / 2., zl / 2./) ! Centre of plasma
        number_faces = 3
        
    case(4) ! ellipsoid
        Vplas = 4 * pi * x_plasma * y_plasma * z_plasma * r_sphere**3 / 3.
        focus = (/xl / 2. - x_plasma * r_sphere, yl / 2., zl / 2./) ! Centre of laser focal spot
        plasma_centre = (/xl / 2., yl / 2., zl / 2./) 
        number_faces = 1

    case(5) ! wedge
        Vplas = .5 * x_plasma * y_plasma * z_plasma
        focus = (/xl / 2. - x_plasma / 2., yl / 2., zl / 2./)
        plasma_centre = (/xl / 2., yl / 2., zl / 2./)
        number_faces = 5

    case(6) ! hemisphere
        Vplas = 4 * pi * r_sphere**3 / 6.
        focus = (/xl / 2. - r_sphere / 2., yl / 2., zl / 2./)
        plasma_centre = (/xl / 2., yl / 2., zl / 2./)
        number_faces = 2

    case(7) ! hollow sphere
        Vplas = (4 * pi / 3.) * (r_sphere**3 - (r_sphere - x_plasma)**3)
        focus = (/xl / 2. - r_sphere / 2., yl / 2., zl / 2./)
        plasma_centre = (/xl / 2., yl / 2., zl / 2./)
        number_faces = 2

    case(8) ! hollow hemishpere
        Vplas = (4 * pi / 6.) * (r_sphere**3 - (r_sphere - x_plasma)**3)
        focus = (/xl / 2. - r_sphere / 2., yl / 2., zl / 2./)
        plasma_centre = (/xl / 2., yl / 2., zl / 2./)
        number_faces = 3

     case(10) ! Electrons only user-defined config (special_start)
        Vplas = x_plasma * y_plasma * z_plasma
        focus = (/xl /4., yl / 2., zl / 2./) ! Centre of laser focal spot
        plasma_centre =  (/xl / 2., yl / 2., zl / 2./) ! Centre of plasma
        number_faces = 6  
    end select geometry

  vte = sqrt(Te_keV/511.)  ! convert from keV to /c
  if (ne > 0) then
     qe = -Vplas*rho0/ne
     qi = -qe*q_factor
     mass_e = -qe
  else
     ! ions only
     qi = Vplas*rho0/ni
     qe = -qi
     mass_e = qi
  endif

  a_ii = (Vplas/ni)**(1./3.)
  vti = sqrt(Ti_keV/511./mass_ratio)
  mass_i = mass_e*mass_ratio
  convert_fs = 10.*omega*lambda/(6*pi)     ! convert from wp^-1 to fs
  convert_mu = omega/2./pi*lambda          ! convert from c/wp to microns
  lolam = lolam*2.*pi/omega  ! normalise scale-length

  r_neighbour = fnn*a_ii  ! Nearest neighbour search radius

  navcycle = 2*pi/dt/omega  ! # timesteps in a laser cycle

  nu_ei = 1./40./pi*a_ii**3/vte/eps**2  ! collision frequency (fit to Okuda & Birdsall)

  sigma_e = 1./nu_ei   ! Spitzer conductivity

  intensity = 0.2*vosc**2*omega**2  ! normalised laser intensity

!  if (scheme > 1 .and. beam_config >= 4) then
!     if (me==0) write(*,*) 'Constant-Te mode: turning beam off'
!     beam_config=0
!  endif

  if ( beam_config ==4 .or. beam_config==6) then
     rho_beam= vosc
     r_beam=sigma

  else if (scheme == 5) then
  ! ion crystal eqm mode
     r_beam = a_ii
     u_beam = Ti_keV
     rho_beam = log(bond_const)
  endif



end subroutine setup
