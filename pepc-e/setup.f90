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
  use tree_utils
  implicit none
  include 'mpif.h'

  integer :: k, npb_pe
  real :: Qplas, Aplas

  type (particle_p1) :: ship_props_a, get_props_a
  integer, parameter :: nprops_particle=10   ! # particle properties to ship
  integer, dimension(nprops_particle) :: blocklengths, displacements, types
  integer*8 :: send_base, receive_base
  integer :: ierr
  integer*8, dimension(nprops_particle) :: address

 
  namelist /pepcdata/ nep, nip, np_mult, fetch_mult, ne, ni, &
       mac, theta, mass_ratio, q_factor, eps, &
       system_config, target_geometry, ispecial, &
       Te_keV, Ti_keV, T_scale, &
       r_sphere, x_plasma, y_plasma, z_plasma, delta_mc, &
       xl, yl, zl, displace, bond_const, fnn, rho_min, lolam, &
       beam_config_in, np_beam, idim, &
       r_beam, u_beam, theta_beam, phi_beam, x_beam, start_beam, rho_beam, mass_beam, & 
       lambda, sigma, tpulse, vosc, omega, focus, x_offset,  z_offset, &
       nt, dt, mc_steps, idump, ivis, ivis_fields, ivis_domains, iprot, itrack, ngx, ngy, ngz, &
       vis_on, steering,  mc_init, restart, scheme, particle_bcs, &
       coulomb, bonds, lenjones, target_dup, ramp, &
       db_level, &
       constrain_proof, len_tripod, struct_step, uthresh, bfield_on

  !  Default input set

 
  system_config = 1  ! plasma target
  target_geometry = 1         ! random sphere

  db_level=1

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
  force_const = 1.
  bond_const = 0.1
  rho0 = 1.0
  mac = 0        ! Multipole acceptance criterion (BH by default)
  theta = 0.5
  err_f = 0.01   ! force error tolerance
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


  ! control
  nt = 600
  dt = 0.2
  trun = 0.
  ivis = 1
  ivis_fields = 1
  ivis_domains = 1
  itime_start = 0
  itrack = 10

  ngx = 25   ! Grid size for plots
  ngy = 25
  ngz = 25
  ! constrain
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
     np_local = nep+nip
!     npp = nep+nip
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
 

     endif
  endif

  npart_total = ni+ne
  np_local = npart_total/n_cpu  ! initial total # particles per processor
!  npp = npart_total/n_cpu  ! initial total # particles per processor
  nppm = np_local*1.5

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

  if (system_config==2) then ! Electrons only user-defined config (special_start)
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

  ! array allocation

  allocate ( x(nppm), y(nppm), z(nppm), ux(nppm), uy(nppm), uz(nppm), & 
       q(nppm), m(nppm), Ex(nppm), Ey(nppm), Ez(nppm), pot(nppm), pelabel(nppm), number(nppm), work(nppm) )

  allocate (vbuffer(0:attrib_max-1,nbuf_max), vbuf_local(0:attrib_max-1,nbuf_max))


  blocklengths(1:nprops_particle) = 1   

  types(1:8) = MPI_REAL8
  types(9) = MPI_INTEGER8
  types(10) = MPI_INTEGER

!  receive_base=LOC(get_props_a%x)
  call LOCADDRESS( get_props_a%x, receive_base, ierr )  ! Base address for receive buffer
  call LOCADDRESS( ship_props_a%x, send_base, ierr )  ! Base address for send buffer

!  if (me==0) write(*,'(a30,o21)') 'Particle address base:',receive_base
!  call MPI_GET_ADDRESS( get_props_a%x, receive_base, ierr )  ! Base address for receive buffer
!  call MPI_GET_ADDRESS( ship_props_a%x, send_base, ierr )  ! Base address for send buffer

  call LOCADDRESS( ship_props_a%x, address(1), ierr )
  call LOCADDRESS( ship_props_a%y, address(2), ierr )
  call LOCADDRESS( ship_props_a%z, address(3), ierr )
  call LOCADDRESS( ship_props_a%ux, address(4), ierr )
  call LOCADDRESS( ship_props_a%uy, address(5), ierr )
  call LOCADDRESS( ship_props_a%uz, address(6), ierr )
  call LOCADDRESS( ship_props_a%q, address(7), ierr )
  call LOCADDRESS( ship_props_a%m, address(8), ierr )
  call LOCADDRESS( ship_props_a%work, address(9), ierr )
  call LOCADDRESS( ship_props_a%label, address(10), ierr )

  displacements(1:nprops_particle) = address(1:nprops_particle) - send_base  !  Addresses relative to start of particle (receive) data

  call MPI_TYPE_STRUCT( nprops_particle, blocklengths, displacements, types, mpi_type_particle_p1, ierr )   ! Create and commit
  call MPI_TYPE_COMMIT( mpi_type_particle_p1, ierr)

  if (my_rank == 0) write(*,*) "Starting PEPC-E with",n_cpu," Processors, simulating",np_local," Particles on each Processor in",nt,"timesteps..."

end subroutine setup




