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

  use physvars, only: &
       particles, &
       theta, &
       q_factor, &
       r_sphere, &
       nt, &
       np_local, &
       npart_total, &
       nppm, &
       ni, &
       ne, &
       nip, &
       nep, &
       trun, &
       n_cpu, &
       my_rank, &
       mac, &
       ispecial, &
       idim, &
       eps, &
       dt, &
       dump_time, &
       cp_time

  use module_mirror_boxes, only: &
       t_lattice_1, &
       t_lattice_2, &
       t_lattice_3, &
       periodicity

  use module_pepc, only: &
       pepc_get_para_file

  use module_interaction_specific_types, only: &
       num_neighbour_particles


  implicit none
  include 'mpif.h'

  integer :: ierr, npart_tmp

  character(255) :: parameterfile
  logical :: read_param_file


  namelist /pepcsph/ ne, ni, &
       mac, theta, q_factor, eps, ispecial, &
       r_sphere, idim, nt, dt, &
       t_lattice_1, t_lattice_2, t_lattice_3, periodicity, &
       dump_time, cp_time, num_neighbour_particles


  !  Default input set
  ispecial        =   6

  ! particles
  nep = 0    ! # plasma electrons per PE
  nip = 0
  ne  = 0    ! Total # plasma electrons
  ni  = 10000 ! total # plasma ions

  ! physics stuff
  mac         = 0
  theta       = 0.6
  q_factor    = 1.

  r_sphere      = 4.
  eps           = 0.01

  ! control
  nt           = 1
  dt           = 0.01
  trun         = 0.

  ! read in first command line argument
  call pepc_get_para_file(read_param_file, parameterfile, my_rank)

  if (read_param_file) then
     if(my_rank .eq. 0) write(*,*) "reading parameter file, section pepcsph: ", parameterfile
     open(10,file=parameterfile)
     read(10,NML=pepcsph)
     close(10)
  else
     if(my_rank .eq. 0) write(*,*) "##### using default parameter #####"
  end if


  write(*,*) 'periodicity:', periodicity
  write(*,*) 't_lattice_1:', t_lattice_1




  ! Derived parameters

  if (nep > 0) then
     ! particles specified per processor in input file
     ne = nep*n_cpu  ! total # electrons
     ni = nip*n_cpu  ! total # ions     
  endif

  npart_total = ni+ne

  if (ispecial == 7) then ! Madelung setup needs ne=ni and (ne+ni)mod 8==0 and (ne+ni)/8==2**sth

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


  allocate ( particles(nppm) )


  if (my_rank == 0) then
     write(*,*) "Starting PEPC-SPH with",n_cpu," Processors, simulating",np_local, &
                         " Particles on each Processor in",nt,"timesteps..."
  end if

end subroutine pepc_setup



