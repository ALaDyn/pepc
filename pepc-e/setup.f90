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


subroutine pepc_setup()
  use physvars
  use tree_utils
  use module_fmm_framework
  use tree_walk_utils
  use tree_walk_communicator
  implicit none
  include 'mpif.h'

  real :: Qplas, Aplas

  integer :: ierr, npart_tmp

  character(50) :: parameterfile
  integer :: read_param_file

  integer*4 IARGC



  namelist /pepcdata/ nep, nip, np_mult, ne, ni, num_walk_threads, max_particles_per_thread, &
       mac, theta, mass_ratio, q_factor, eps, curve_type, &
       system_config, target_geometry, ispecial, weighted, &
       Te_keV, Ti_keV, T_scale, &
       r_sphere, x_plasma, y_plasma, z_plasma, &
       xl, yl, zl, &
       idim, &
       nt, dt, idump, &
       db_level, itime_in, idump_vtk, idump_checkpoint, idump_binary, &
       t_lattice_1, t_lattice_2, t_lattice_3, periodicity, do_extrinsic_correction


  !  Default input set
  system_config   =   2
  target_geometry =   0

  db_level        =   0

  np_mult         = -45

  ispecial        =   1

  weighted        =   1

  ! particles
  nep = 0    ! # plasma electrons per PE
  nip = 0
  ne  = 0    ! Total # plasma electrons
  ni  = 1000 ! total # plasma ions

  xl = 1
  yl = 1
  zl = 1

  ! physics stuff
  force_const = 1.
  mac         = 0
  theta       = 0.6
  Te_keV      = 0.5
  Ti_keV      = 0.1
  mass_ratio  = 2000.
  q_factor    = 1.
  rho0        = 1.

  r_sphere      = 4
  x_plasma      = 1.
  y_plasma      = 1.
  z_plasma      = 1.
  eps           = 0.01

  ! control
  nt           = 10
  dt           = 0.01
  trun         = 0.

  idump        = 0
  idump_vtk    = 0
  idump_checkpoint  = 0
  idump_binary = 0

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

  if (system_config == 2) then
    if (ispecial == 7) then ! Madelung setup needs ne=ni and (ne+ni)mod 8==0 and (ne+ni)/8==2**sth
      npart_tmp   = nint((npart_total/8)**(1./3.))
      npart_total = (npart_tmp**3)*8
      if (my_rank == 0) write(*,*) "Using 3D-Madelung Setup: Total particle number must be representable by 8*k^3. Setting npart_total =", npart_total
      ne = npart_total/2
      ni = ne
    end if
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
     plasma_centre =  (/xl / 2., yl / 2., zl / 2./) ! Centre of plasma
     number_faces = 6

  case(1) ! sphere
     Vplas = 4 * pi * r_sphere**3 / 3.
     Aplas = pi*r_sphere**2
     plasma_centre = (/xl / 2., yl / 2., zl / 2./) ! Centre of plasma
     number_faces = 1

  case(2) ! disc
     Vplas = pi * r_sphere**2 * x_plasma
     Aplas = x_plasma*y_plasma
     plasma_centre = (/xl / 2., yl / 2., zl / 2./) ! Centre of plasma        
     number_faces = 3

  case(3) ! wire
     Vplas = pi * r_sphere**2 * z_plasma
     Aplas = pi*r_sphere**2
     plasma_centre = (/xl / 2., yl / 2., zl / 2./) ! Centre of plasma
     number_faces = 3

  case(4) ! ellipsoid
     Vplas = 4 * pi * x_plasma * y_plasma * z_plasma / 3.
     Aplas = pi*x_plasma*y_plasma*2
     plasma_centre = (/xl / 2., yl / 2., zl / 2./) 
     number_faces = 1

  case(5) ! wedge
     Vplas = .5 * x_plasma * y_plasma * z_plasma
     Aplas = .5*x_plasma*y_plasma
     plasma_centre = (/xl / 2., yl / 2., zl / 2./)
     number_faces = 5

  case(6) ! hemisphere
     Vplas = 4 * pi * r_sphere**3 / 6.
     Aplas = pi*r_sphere**2/2.
     plasma_centre = (/xl / 2., yl / 2., zl / 2./)
     number_faces = 2

  case(7) ! hollow sphere
     Vplas = (4 * pi / 3.) * (r_sphere**3 - (r_sphere - x_plasma)**3)
     Aplas = pi*(r_sphere**2-(r_sphere-x_plasma)**2)
     plasma_centre = (/xl / 2., yl / 2., zl / 2./)
     number_faces = 2

  case(8) ! hollow hemisphere
     Vplas = (4 * pi / 6.) * (r_sphere**3 - (r_sphere - x_plasma)**3)
     Aplas = pi/2.*(r_sphere**2-(r_sphere-x_plasma)**2)
     plasma_centre = (/xl / 2., yl / 2., zl / 2./)
     number_faces = 3

  end select geometry

  if (system_config==2) then ! Electrons only user-defined config (special_start)
     Vplas = x_plasma * y_plasma * z_plasma
     Aplas = x_plasma * y_plasma
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
  convert_keV = 2./3./Qplas*511     ! convert from code energy units to keV/particle

  allocate ( x(nppm), y(nppm), z(nppm), ux(nppm), uy(nppm), uz(nppm), & 
       q(nppm), m(nppm), Ex(nppm), Ey(nppm), Ez(nppm), pot(nppm), pelabel(nppm), work(nppm) )


  if (my_rank == 0) then
     write(*,*) "Starting PEPC-E with",n_cpu," Processors, simulating",np_local, &
                         " Particles on each Processor in",nt,"timesteps..."
     write(*,*) "Using",num_walk_threads,"worker-threads and 1 communication thread in treewalk on each processor (i.e. per MPI rank)"
     write(*,*) "Maximum number of particles per work_thread = ", max_particles_per_thread
  end if

end subroutine pepc_setup




