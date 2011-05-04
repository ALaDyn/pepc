!  ================================
!
!         SETUP
!
!   $Revision: 1841 $
!
!     Initialise constants and 
!      simulation variables
!
!  ================================


subroutine setup()
  use physvars
  use tree_utils
  use module_fmm_framework
  use tree_walk_utils
  use tree_walk_communicator
  use module_icosahedron
  use module_laser
  use module_pusher
  use module_workflow
  use module_units
  use module_param_dump
  use module_fields
  implicit none
  include 'mpif.h'

  integer :: ierr, ifile

  character(50) :: parameterfile
  integer :: read_param_file

  integer*4 IARGC


  namelist /pepcdata/ &
       np_mult, num_walk_threads, mac, theta, max_particles_per_thread, &
       weighted, &                                      ! algorithm parameters
       ne,  eps, nt, dt, idump, db_level, itime_in, idump_vtk, idump_checkpoint, idump_binary, & ! fundamental stuff
       ispecial, rhoe_nm3, Zion, Aion, Te_eV, Ti_eV, Te_K, Ti_K, &   ! experimental setup
       workflow_setup, &                                             ! workflow
       integrator_scheme, enable_drift_elimination, &                ! pusher configuration
       beam_config_in, vosc,omega, sigma, tpulse, theta_inc, rho_track, omega_wpl, I0_Wpercm2, & ! laser config
       t_lattice_1, t_lattice_2, t_lattice_3, periodicity, do_extrinsic_correction, &            ! periodicity config
       field_dump_ncells, ngx, ngy, ngz                              ! diagnostics config


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!  default parameters            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  db_level        = 0

  np_mult         = -45

  ispecial        = 1

  weighted        = 1

  ! particles
  nep = 0    ! # plasma electrons per PE
  nip = 0
  ne  = 1000    ! Total # plasma electrons
  ni  = 0 ! total # plasma ions

  xl = 1
  yl = 1
  zl = 1

  ! physics stuff
  force_const = 1.
  mac         = 0
  theta       = 0.6
  Te_eV       = 50
  Ti_eV       = 10

  r_sphere      = 4
  x_plasma      = 1.
  y_plasma      = 1.
  z_plasma      = 1.
  eps           = 1.

  ! control
  nt           = 10
  dt           = 0.01
  trun         = 0.

  idump        = 0
  idump_vtk    = 0
  idump_checkpoint  = 0
  idump_binary = 0

  ngx = 25
  ngy = 25
  ngz = 25


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!  read parameter file           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! rank 0 reads in first command line argument
  read_param_file = 0
  if (my_rank .eq. 0) then
     if( IARGC() .ne. 0 ) then
        call GETARG(1, parameterfile)
        read_param_file = 1
     end if
  end if

  call MPI_BCAST( read_param_file, 1, MPI_INTEGER, 0, MPI_COMM_PEPC, ierr )

  ! broadcast file name, read actual inputs from namelist file

  if (read_param_file .eq. 1) then

     call MPI_BCAST( parameterfile, 50, MPI_CHARACTER, 0, MPI_COMM_PEPC, ierr )

     if(my_rank .eq. 0) write(*,*) "reading parameter file: ", parameterfile
     open(10,file=parameterfile)
     read(10,NML=pepcdata)
     close(10)
  else
     if(my_rank .eq. 0) write(*,*) "##### using default parameter #####"
  end if


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!  derived parameters (physics)  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  qe     = unit_qe
  mass_e = unit_me
  qi     = unit_qp*Zion
  ni     = ne/Zion ! we assume overall charge neutrality (Zion = ion charge number)
  mass_i = unit_mp*Aion ! Aion = ion mass number
  npart_total = ni+ne

  if (ispecial == 7) then ! Madelung setup needs ne=ni and (ne+ni)mod 8==0 and (ne+ni)/8==2**sth
     npart_total = (nint((npart_total/8)**(1./3.))**3)*8
     if (my_rank == 0) write(*,*) "Using 3D-Madelung Setup: Total particle number must be representable by 8*k^3. Setting npart_total =", npart_total
     ne = npart_total/2
     ni = ne
  end if

  if (ispecial == 13) then ! ion lattice setup needs integer number of ions per edge
     ni = nint((ne/Zion)**(1./3.))**3
     ne = ni * Zion
     npart_total = ne + ni
     if (my_rank == 0) write(*,*) "Using Ion Lattice Setup: Number of Ions per edge must be integer, setting ne =", ne, ", ni=", ni
  end if
  
  if (ispecial == 12) then
    npart_total   = get_nextlower_particles(npart_total/2)*2
    if (my_rank == 0) write(*,*) "Using Mackay Icosahedron: Total particle number must be two times a magic cluster number. Setting npart_total =", npart_total
    ne = npart_total/2
    ni = ne
  end if


  Vplas    =  ne/rhoe_nm3 / unit_abohr_in_nm**3. ! adjust simulation volume to fit requested electron density while keeping particle number constant
  rhoi_nm3 =  ni/Vplas / unit_abohr_in_nm**3.
  x_plasma = (ne/rhoe_nm3)**(1./3.) / unit_abohr_in_nm ! (assume cubic volume)
  y_plasma = x_plasma
  z_plasma = x_plasma
  xl       = x_plasma
  yl       = y_plasma
  zl       = z_plasma

  plasma_centre =  (/xl / 2., yl / 2., zl / 2./) ! Centre of plasma
  a_ee = (Vplas/ne)**(1./3.)
  a_ii = (Vplas/ni)**(1./3.)
  r_sphere = ((3.*Vplas)/(4.*pi))**(1./3.)

  if (Te_K > 0.) Te_eV = unit_kB_in_eVperK * Te_K
  if (Ti_K > 0.) Ti_eV = unit_kB_in_eVperK * Ti_K
  Te    = Te_eV / unit_Ryd_in_eV
  Ti    = Ti_eV / unit_Ryd_in_eV
  Te_K  = Te_eV / unit_kB_in_eVperK
  Ti_K  = Ti_eV / unit_kB_in_eVperK

  vte = sqrt(3*unit_kB*Te/mass_e)
  vti = sqrt(3*unit_kB*Ti/mass_i)

  force_const = 1./(unit_4piepsilon0)

  wpl_e = sqrt( (ne/Vplas * qe *qe) / (unit_epsilon0 * mass_e) )
  wpl_i = sqrt( (ni/Vplas * qi *qi) / (unit_epsilon0 * mass_i) )

  lambdaD_e = sqrt( (unit_epsilon0*unit_kB*Te)/(qe*qe) * Vplas/ne )
  lambdaD_i = sqrt( (unit_epsilon0*unit_kB*Ti)/(qi*qi) * Vplas/ni )

  t_lattice_1 = t_lattice_1*x_plasma
  t_lattice_2 = t_lattice_2*y_plasma
  t_lattice_3 = t_lattice_3*z_plasma

  a_i       = (4.*pi/3. * ni/Vplas)**(-1./3.)
  physGamma = (qi*qi) / (a_i * unit_kB*Te)

  eps = eps * lambdaD_e

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!  parameters (laser)            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (omega_wpl > 0.) omega = omega_wpl * wpl_e
  omega_wpl = omega / wpl_e
  omega_hz  = omega / unit_t0_in_s
  lambda    = unit_c / omega
  lambda_nm = lambda * unit_abohr_in_nm
  rhocrit_nm3 = omega*omega*unit_epsilon0*mass_e/(qe*qe) / unit_abohr_in_nm**3.

  if (I0_Wpercm2 > 0.) then
    E0   = unit_Z0 *sqrt(I0_Wpercm2*1.E4 / unit_P0_in_W ) * unit_abohr_in_m
    vosc = (abs(qe)*E0)/(mass_e*omega)
  endif
  E0         = vosc*mass_e*omega/abs(qe)
  I0_Wpercm2 = (E0 / unit_abohr_in_m / unit_Z0)**2. * unit_P0_in_W * 1.E-4

  call setup_laser()


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!  parameters (simulation generic)   !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  maxdt(1) = 2*pi/max(omega,1.e-10) * 1./100.
  maxdt(2) = 2*pi/max(wpl_e,1.e-10_8) * 1./100.
  maxdt(3) = lambdaD_e/10./vte
  maxdt(4) = abs(mass_e/qe * eps*eps / (10.*qe) * vte/10.)

  if (any(maxdt < dt)) then
    if (my_rank == 0) then
      do ifile = 6,24,18
        write(ifile,*) "!!!!!!!! WARNING: timestep dt is too large: dt =", dt, "   maxdt = ", maxdt
        write(ifile,*) "Adjusting parameters appropriately to ensure sufficient resolution while keeping simulation length"
      end do
    end if

    nt = nt * dt / minval(maxdt)
    dt = minval(maxdt)
  endif


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!  parameter output              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (my_rank == 0) then
    call PrintPhysicalParameters(6)
    call PrintPhysicalParameters(24)
  endif


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!  derived parameters (tree code)  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  nep = ne/n_cpu
  nip = ni/n_cpu
  if (nep*n_cpu /= ne .and. mod(ne,n_cpu) > my_rank) then
     nep = ne/n_cpu+1
  end if
  if (nip*n_cpu /= ni .and. mod(ni,n_cpu) > my_rank) then
     nip = ni/n_cpu+1
  end if

  np_local = nep+nip

  if (n_cpu.eq.1) then
     nppm = nint(1.5*npart_total + 1000)  ! allow for additional ghost particles for field plots
  else
     nppm = nint(1.5*max(npart_total/n_cpu,1000)) ! allow 50% fluctuation
  end if


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!  array allocation                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  allocate ( &
       x(nppm), y(nppm), z(nppm), &
       ux(nppm), uy(nppm), uz(nppm), &
       q(nppm), m(nppm), &
       Ex(nppm), Ey(nppm), Ez(nppm), &
       Ax(nppm), Ay(nppm), Az(nppm), &
       Bx(nppm), By(nppm), Bz(nppm), pot(nppm), pelabel(nppm), work(nppm) )

  if (my_rank == 0) then
    write(*,*) "Starting PEPC-MW with",n_cpu," Processors, simulating",np_local, &
			" Particles on each Processor in",nt,"timesteps..."
	write(*,*) "Using",num_walk_threads,"worker-threads in treewalk on each processor (i.e. per MPI rank)"
    write(*,*) "Maximum number of particles per work_thread = ", max_particles_per_thread
  end if


end subroutine setup




