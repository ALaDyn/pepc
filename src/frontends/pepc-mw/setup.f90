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


subroutine pepc_setup()
  use physvars
  use tree_utils
  use module_setup
  use module_fmm_framework
  use module_mirror_boxes
  use module_particle_setup
  use module_icosahedron
  use module_laser
  use module_pusher
  use module_workflow
  use module_units
  use module_fields
  implicit none
  include 'mpif.h'

  integer :: ierr, ifile

  character(255) :: parameterfile
  integer :: read_param_file

  namelist /pepcmw/ &
       mac, theta, &
       ne,  eps, V0_eV, nt, dt, idump, itime_in, idump_vtk, idump_checkpoint, idump_binary, treediags, & ! fundamental stuff
       ispecial, rhoe_nm3, Zion, Aion, Te_eV, Ti_eV, Te_K, Ti_K, rngseed, &   ! experimental setup
       workflow_setup, &                                             ! workflow
       integrator_scheme, enable_drift_elimination, &                ! pusher configuration
       beam_config_in, vosc,omega, sigma, t_pulse_fs, theta_inc, rho_track, omega_wpl, I0_Wpercm2, lambda_nm, & ! laser config
       t_lattice_1, t_lattice_2, t_lattice_3, periodicity, do_extrinsic_correction, &            ! periodicity config
       field_dump_ncells, ngx, ngy, ngz                              ! diagnostics config


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!  default parameters            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ispecial        = 1

  ! particles
  ne  = 200! Total # plasma electrons
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
  nt           = 100
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
  ! read in first command line argument
  call libpepc_get_para_file(read_param_file, parameterfile, my_rank)
  if (read_param_file .eq. 1) then

     if(my_rank .eq. 0) write(*,*) "reading parameter file: ", parameterfile
     open(10,file=parameterfile)
     read(10,NML=pepcmw)
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

  Vplas    =  ne/rhoe_nm3 / unit_abohr_in_nm**3. ! adjust simulation volume to fit requested electron density while keeping particle number constant
  rhoi_nm3 =  ni/Vplas / unit_abohr_in_nm**3.
  x_plasma = (ne/rhoe_nm3)**(1./3.) / unit_abohr_in_nm ! (assume cubic volume)
  y_plasma = x_plasma
  z_plasma = x_plasma
  xl       = x_plasma
  yl       = y_plasma
  zl       = z_plasma

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

  if (V0_eV > 0.) eps = force_const * qe*qi / (V0_eV / unit_Ryd_in_eV)

  V0_eV = (force_const * qi / eps) * unit_Ryd_in_eV

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!  parameters (laser)            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (omega_wpl > 0.) omega = omega_wpl * wpl_e
  if (lambda_nm > 0.) omega = unit_c / (lambda_nm / unit_abohr_in_nm)
  omega_wpl = omega / wpl_e
  omega_hz  = omega / unit_t0_in_s
  lambda    = unit_c / omega
  lambda_nm = lambda * unit_abohr_in_nm
  rhocrit_nm3 = omega*omega*unit_epsilon0*mass_e/(qe*qe) / unit_abohr_in_nm**3.

  if (I0_Wpercm2 > 0.) then
    E0   =  sqrt( 2./(unit_epsilon0*unit_c) * (I0_Wpercm2 / unit_P0_in_W * (100*unit_abohr_in_m)**2) )
    vosc = (abs(qe)*E0)/(mass_e*omega)
  endif
  E0         = vosc*mass_e*omega/abs(qe)
  I0_Wpercm2 = (unit_epsilon0 * unit_c * E0**2 / 2. ) * unit_P0_in_W / (100*unit_abohr_in_m)**2

  call laser_setup()


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!  parameters (simulation generic)   !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  maxdt(1) = 2*pi/max(omega,1.e-10) * 1./100.
  maxdt(2) = 2*pi/max(wpl_e,1.e-10_8) * 1./100.

  if (vte>0.) then
    maxdt(3) = lambdaD_e/10./vte
    maxdt(4) = abs(mass_e/qe * eps*eps / (10.*qe) * vte/10.)
  else
    maxdt(3:4) = 1./epsilon(maxdt(3))
  endif

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
  !!!!!!!!!!!!!!!  derived parameters (tree code)  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call update_particle_numbers(ne + ni)

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
       Bx(nppm), By(nppm), Bz(nppm), pot(nppm), pelabel(nppm), work(nppm), energy(1:3,nppm) )

  if (my_rank == 0) then
    write(*,*) "Starting PEPC-MW with",n_cpu," Processors, simulating",np_local, &
			" Particles on each Processor in",nt,"timesteps..."
    !write(*,*) "Using",num_walk_threads,"worker-threads in treewalk on each processor (i.e. per MPI rank)"
    !write(*,*) "Maximum number of particles per work_thread = ", max_particles_per_thread
  end if


end subroutine pepc_setup




