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


subroutine setup()
  use physvars
  use tree_utils
  use module_fmm_framework
  implicit none
  include 'mpif.h'

  real :: Qplas, Aplas

  type (particle_p1) :: ship_props_a, get_props_a
  integer, parameter :: nprops_particle=10   ! # particle properties to ship
  integer, dimension(nprops_particle) :: blocklengths, displacements, types
  integer*8 :: send_base, receive_base
  integer :: ierr, npart_tmp
  integer*8, dimension(nprops_particle) :: address

  character(50) :: parameterfile
  integer :: read_param_file

  integer*4 IARGC



  namelist /pepcdata/ nep, nip, np_mult, fetch_mult, ne, ni, &
       mac, theta, mass_ratio, q_factor, eps, &
       system_config, target_geometry, ispecial, choose_sort, weighted, choose_build, &
       Te_keV, Ti_keV, T_scale, &
       r_sphere, x_plasma, y_plasma, z_plasma, delta_mc, &
       xl, yl, zl, displace, bond_const, rho_min, lolam, &
       beam_config_in, np_beam, idim, &
       r_beam, u_beam, theta_beam, phi_beam, x_beam, start_beam, rho_beam, mass_beam, & 
       lambda, sigma, tpulse, vosc, omega, focus, x_offset,  z_offset, &
       nt, dt, mc_steps, idump, ivis, ivis_fields, ivis_domains, iprot, itrack, ngx, ngy, ngz, &
       vis_on, steering,  mc_init, restart, scheme, &
       coulomb, bonds, lenjones, target_dup, ramp, &
       db_level, &
       constrain_proof, len_tripod, struct_step, uthresh, bfield_on, &
       t_lattice_1, t_lattice_2, t_lattice_3, periodicity, do_extrinsic_correction


  !  Default input set
 
  system_config   =   2
  target_geometry =   0

  db_level        =   0

  np_mult         = -45
  fetch_mult      =   3

  ispecial        =   1

  choose_sort     =   3
  weighted        =   1
  choose_build    =   0

  scheme          =   0

  ! particles
  nep = 0    ! # plasma electrons per PE
  nip = 0
  ne  = 0    ! Total # plasma electrons
  ni  = 1000 ! total # plasma ions
  mc_steps = 10

  xl = 1
  yl = 1
  zl = 1

  ! physics stuff
  force_const = 1.
  bond_const  = 2.e-3
  rho0        = 1.0
  mac         = 0
  theta       = 0.6
  err_f       = 0.01
  Te_keV      = 0.5
  Ti_keV      = 0.1
  mass_ratio  = 2000.
  q_factor    = 1.
  uthresh     = -1.

  r_sphere      = 4
  x_plasma      = 1.
  y_plasma      = 1.
  z_plasma      = 1.
  eps           = 0.01
  delta_mc      = r_sphere/5.
  displace(1:3) = (/0.,0.,0./)


  ! control
  nt           = 10
  dt           = 0.01
  trun         = 0.
  ivis         = 2
  ivis_fields  = 5000
  ivis_domains = 5000
  itime_start  = 0
  itrack       = 10

  restart      = .false.
  vis_on       = .false.
  idump        = 0
  iprot        = 1

  ngx = 25   
  ngy = 25
  ngz = 25

  ! constrain
  constrain_proof = .001
  struct_step = 0

  ! rank 0 reads in first command line argument
  read_param_file = 0
  if (my_rank .eq. 0) then
     if( IARGC() .ne. 0 ) then
        call GETARG(1, parameterfile)
        read_param_file = 1
     end if
  end if

  call MPI_BCAST( read_param_file, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )

  ! broadcast file name, read actual inputs from namelist file

  if (read_param_file .eq. 1) then

     call MPI_BCAST( parameterfile, 50, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr )

     if(my_rank .eq. 0) write(*,*) "reading parameter file: ", parameterfile
     open(10,file=parameterfile)
     read(10,NML=pepcdata)
     close(10)
  else
     if(my_rank .eq. 0) write(*,*) "##### using default parameter #####"
  end if

  ! Derived parameters

  if (nep > 0) then
     ! particles specified per processor in input file
     ne = nep*n_cpu  ! total # electrons
     ni = nip*n_cpu  ! total # ions     
  endif

  npart_total = ni+ne

  if ((system_config == 2) .and. (ispecial == 7)) then ! Madelung setup needs ne=ni and (ne+ni)mod 8==0 and (ne+ni)/8==2**sth

    npart_tmp   = nint((npart_total/8)**(1./3.))
    npart_total = (npart_tmp**3)*8

    if (my_rank == 0) write(*,*) "Using 3D-Madelung Setup: Total particle number must be representable by 8*k^3. Setting npart_total =", npart_total

    ne = npart_total/2
    ni = ne
  end if

  ! total # particles specified in input file
  nep = ne/n_cpu
  nip = ni/n_cpu
  if (nep*n_cpu /= ne .and. mod(ne,n_cpu) > my_rank)  nep = ne/n_cpu+1
  if (nip*n_cpu /= ni .and. mod(ni,n_cpu) > my_rank)   nip = ni/n_cpu+1

  np_local = nep+nip

  if (n_cpu.eq.1) then
     nppm=int(1.5*npart_total + 1000)  ! allow for additional ghost particles for field plots
!  else if (np_mult<0) then 
!     nppm = abs(np_mult)*max(npart_total/n_cpu,1000) ! allow 50% fluctuation
  else
     nppm = int(1.5*max(npart_total/n_cpu,1000)) ! allow 50% fluctuation
  end if

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

  displacements(1:nprops_particle) = int(address(1:nprops_particle) - send_base)  !  Addresses relative to start of particle (receive) data

  call MPI_TYPE_STRUCT( nprops_particle, blocklengths, displacements, types, mpi_type_particle_p1, ierr )   ! Create and commit
  call MPI_TYPE_COMMIT( mpi_type_particle_p1, ierr)

  if (my_rank == 0) write(*,*) "Starting PEPC-E with",n_cpu," Processors, simulating",np_local, &
			" Particles on each Processor in",nt,"timesteps..."

end subroutine setup




