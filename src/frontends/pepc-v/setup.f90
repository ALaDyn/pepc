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
  use tree_walk_pthreads
  implicit none
  include 'mpif.h'

  integer :: ierr

  character(50) :: parameterfile
  integer :: read_param_file

  integer*4 IARGC



  namelist /pepcdata/ np_mult, n, num_walk_threads, max_particles_per_thread, &
       mac, theta, eps, ispecial, weighted, curve_type, dt, ts, te, db_level, &
       h, m_h, nu, rem_freq, &
       rmax, r_torus, nc, nphi, g, torus_offset


  !  Default input set
  db_level        =   0
  np_mult         = -45
  ispecial        =   1
  weighted        =   1
  curve_type      =   0


  ! particles
  np = 0
  n  = 1000 ! total # vortex particles

  ! physics stuff
  force_const  = 0.25D00/pi  ! 3D prefactor for u and af
  mac          = 0
  theta        = 0.6
  eps          = 0.01      ! this is my smoothing radius, will adapt it in special_start
  h            = 0.
  m_h          = 0.
  rem_freq     = 0
  rmax         = 0.
  r_torus      = 0.
  nc           = 0
  nphi         = 0
  g            = 0
  torus_offset = [0., 0., 0.]

  ! control
  ts           = 0.
  te           = 0.
  dt           = 0.01

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

  init: select case(ispecial)
  case(1,2)                         ! Vortex rings
     rl = rmax/(2*nc+1)
     ns = 1+4*nc*(nc+1)
     np = int(2.0*Ns*Nphi/n_cpu)
     n  = 2*Ns*Nphi
     kernel_c = sqrt(nu*rem_freq*dt)/m_h
  case default
     write(*,*) 'ERROR: need to specify setup via ispecial in your .h file'
     call MPI_ABORT(MPI_COMM_WORLD,ierr)
     stop
  end select init

  trun = ts
  nt = int((te-ts)/dt) ! Number of timesteps
  rk_stages = 2   ! TODO: inflexible RK time integration scheme, hard-wired so far

  if (n_cpu.eq.1) then
     nppm=int(1.5*n + 1000)  ! allow for additional ghost particles for field plots
!  else if (np_mult<0) then 
!     nppm = abs(np_mult)*max(npart_total/n_cpu,1000) ! allow 50% fluctuation
  else
     nppm = int(1.5*max(n/n_cpu,1000)) ! allow 50% fluctuation
  end if

  allocate ( vortex_particles(nppm) )

  if (my_rank == 0) then
     write(*,*) "Starting PEPC-MINI with",n_cpu," Processors, simulating",np, &
                         " Particles on each Processor in",nt,"timesteps..."
     write(*,*) "Using",num_walk_threads,"worker-threads and 1 communication thread in treewalk on each processor (i.e. per MPI rank)"
     write(*,*) "Maximum number of particles per work_thread = ", max_particles_per_thread
  end if

end subroutine pepc_setup




