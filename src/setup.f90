!  ================================
!
!         SETUP
!
!     Initialise constants and 
!      simulation variables
!
!  ================================


subroutine setup

  use treevars
  use utils

  implicit none
  integer :: ibig, machinebits, maxleaf, maxtwig,k, npb_pe
  real :: q_factor, x_offset, z_offset
  character*7 :: configs(0:10)= (/   'default', &
       'sphere ','disc   ','wire   ','slab   ','       ', &
       '       ','       ','       ','       ','restart' /)
  character*7 :: beam_configs(0:9)=(/ &
       'eqm    ','beam   ','i-beam ','laser-u','fpond  ', &
       '       ','       ','       ','       ','       ' /)

  namelist /pepcdata/ nep, nip, ne, ni, &
       theta, mass_ratio, q_factor, eps, &
       initial_config, &
       Te_keV, Ti_keV, &
       r_sphere, x_plasma, y_plasma, delta_mc, &
       xl, yl, zl, displace, &
       beam_config, np_beam, &
       r_beam, u_beam, theta_beam, phi_beam, x_beam, start_beam, rho_beam, mass_beam, & 
       lambda, sigma, tpulse, vosc, omega, focus, x_offset,  z_offset, &
       nt, dt, mc_steps, idump, ivis, ivis_fields, nmerge, ngx, ngy, ngz, &
       vis_on, domain_debug,  mc_init, restart, ensemble, &
       load_balance, walk_balance, walk_debug, dump_tree, perf_anal




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
  rho0 = 1.0
  theta = 0.8
  Te_keV = 1.
  Ti_keV = Te_keV/10.
  mass_ratio = 10.
  q_factor = 1.

  r_sphere = 0.5
  x_plasma = 0.1    ! plasma disc thickness (2) or wire length (3)
  y_plasma = 1.     ! plasma width (slab target)
  eps = 0.1
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

  ! control
  nt = 600
  dt = 0.2
  trun = 0.
  ivis = 1
  ivis_fields = 1
  itime_start = 0
  ngx = 100   ! Grid size for plots
  ngy = 50
  ngz = 50
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

  npb_pe = np_beam/num_pe
  if (.not. restart) then
     if (nip*num_pe /= ni .or. nep*num_pe /= ne ) then
        ne = (nep+1)*num_pe
        ni = (nip+1)*num_pe

        if (me==0) then
           write(*,'(//a//)') '*** Warning: number each particle species (ne, ni) must be divisible by # processors ***'
           write(*,'(a,i6)') '*** Resetting to ',ne,ni
        endif
     else if (npb_pe*num_pe /= np_beam) then 
        np_beam = (npb_pe+1)*num_pe
        if (me==0) then
           write(*,'(//a)') '*** Warning: number of beam particles must be divisible by # processors ***'
           write(*,'(a,i6)') '*** Resetting to ',np_beam
        endif

     endif
  endif

  npart = ni+ne
  npp = nep + nip  ! total # particles per processor


  if (initial_config ==1) then
     ! sphere
     Vplas = 4*pi/3.*r_sphere**3

     focus = (/ xl/2.-r_sphere,yl/2.,zl/2. /) ! Centre of laser focal spot
     plasma_centre =  (/ xl/2., yl/2., zl/2. /) ! Centre of plasma

  else if (initial_config ==2) then
     ! disc
     Vplas = pi*r_sphere**2*x_plasma
     focus = (/ xl/2.-x_plasma/2., yl/2., zl/2. /) ! Centre of laser focal spot
     plasma_centre =  (/ xl/2., yl/2., zl/2. /) ! Centre of plasma

  else if (initial_config ==3) then
     ! wire
     Vplas = pi*r_sphere**2*x_plasma
     focus = (/ xl/2.-r_sphere+x_offset, yl/2., zl/2.+z_offset  /) ! Centre of laser focal spot
     plasma_centre =  (/ xl/2., yl/2., zl/2. /) ! Centre of plasma

  else
     ! slab
     Vplas = x_plasma*y_plasma**2
     focus = (/ xl/2.-x_plasma/2., yl/2., zl/2. /) ! Centre of laser focal spot
     plasma_centre =  (/ xl/2., yl/2., zl/2. /) ! Centre of plasma

  endif


  vte = sqrt(Te_keV/511.)  ! convert from keV to /c
  qe = -Vplas*rho0/ne
  qi = -qe*q_factor
  mass_e = -qe
  convert_fs = 10.*omega*lambda/(6*pi)     ! convert from wp^-1 to fs

  if ( ensemble == 2 ) then
     ! special case for velocity-clamp mode: set ion mass = 100*electron mass
     vti = sqrt(Ti_keV/511./100.)
     mass_i = 100*mass_e
     if (beam_config==4) then 
        if (me==0) write(*,*) 'Constant-Te mode: turning laser off'
        beam_config=0
     endif

  else if (ensemble == 3) then
     vti = sqrt(Ti_keV/511./mass_ratio)
     mass_i = mass_e*mass_ratio
     if (beam_config==4) then 
        if (me==0) write(*,*) 'Constant-Te mode: turning laser off'
        beam_config=0
     endif

  else
     ! default is to use actual, user-defined mass ratio
     vti = sqrt(Ti_keV/511./mass_ratio)
     mass_i = mass_e*mass_ratio
  endif

  if ( beam_config ==4 ) then
     rho_beam= vosc
     r_beam=sigma
  endif


  !  npartm = npart + nt*np_beam  ! Max # particles permitted
  npartm = 2*npart  ! allow 50% fluctuation
  nppm = max(npartm/num_pe,1000)
  nshortm = 800    ! Max shortlist length: leave safety factor for nshort_list in FORCES
  
  ! Estimate of interaction list length - Hernquist expression
  if (theta >0 ) then
     nintmax = 2.5*24*log(1.*npartm)/theta**2
  else
     nintmax = npartm
  endif
  max_list_length = 0 ! current max length of all interaction lists

  ! tree stuff

  idim = 3               ! # dimensions (2 or 3)
  nlev = 20                     ! max refinement level
  iplace = 2_8**(idim*nlev)           ! place holder bit
!  nbaddr = 15                  ! # bits for cell address in hash-table
  !  Space for # table and tree arrays
  !  TODO: need good estimate for max # branches

  size_tree = max(2*nintmax+5*nppm,1000)+1
  maxaddress = size_tree
   nbaddr = log(1.*maxaddress)/log(2.) + 1
  maxaddress = 2**nbaddr
  hashconst = maxaddress-1

  free_lo = 1024      ! lowest free address for collision resolution (from 4th level up)



! Some MPI constants
  lastpe = num_pe - 1          ! # of last PE
  me_minus_one = me - 1
  me_plus_one = me + 1

  if (me==0) then
     do ifile = 6,15,9
        write (ifile,*) 'Max address in #-table: ',2**nbaddr-1
        machinebits = bit_size(1_8)    ! # bits in integer variable (hardware) 
        write (ifile,*) 'Machine bit-size = ',machinebits
        write (ifile,*) 'Max permitted particles / PE: ', nppm
        write (ifile,*) 'Max size of interaction lists: ', nintmax
	write (ifile,*) 'Shortlist length: ',nshortm
        write (ifile,*) 'Electrons: ', ne
        write (ifile,*) 'Ions: ', ni
        write (ifile,*) 'Beam particles: ', np_beam
 	write (ifile,*) 'Beam angles ',theta_beam, phi_beam
        write (ifile,*) 'Particles per PE: ', npp
        write (ifile,*) 'Memory needed for lists = ',4*nintmax*nshortm*8/2**20,' Mbytes'
        write (ifile,*) 'Plasma config: ',configs(initial_config)
        write (ifile,*) 'Laser config: ',beam_configs(beam_config)
        write (ifile,'(a,1pe12.3)') 'Plasma volume: ',Vplas
        write (ifile,'(a,1pe12.3)') 'Sphere radius: ',r_sphere
        write (ifile,'(a,1pe12.3)') 'Electron charge: ',qe
        write (ifile,'(a,1pe12.3)') 'Electron mass: ',mass_e
        write (ifile,'(a,1pe12.3)') 'Ion mass: ',mass_i
        write (ifile,'(a,1pe12.3)') 'Te: ',Te_keV
        write (ifile,'(a,1pe12.3)') 'Ti: ',Ti_keV

        write (ifile,'(a,1pe12.3)') 'Cutoff radius: ',eps
        write (ifile,'(/a/a)') 'Switches:','--------'
        write (ifile,'(a20,l3)') 'load balance: ',load_balance
        write (ifile,'(a20,l3)') 'walk balance: ',walk_balance
        write (ifile,'(a20,l3)') 'restart: ',restart

        write (ifile,'(a20,l3)') 'domain debug: ',domain_debug
        write (ifile,'(a20,l3)') 'walk debug: ',walk_debug
        write (ifile,'(a20,l3)') 'dump tree: ',dump_tree
        write (ifile,'(a20,l3)') 'performance anal.: ',perf_anal
        write (ifile,'(a20,l3)') 'visit: ',vis_on

        write (ifile,'(/a)') 'Other inputs:'
	write(ifile,NML=pepcdata)

        !  ibig = 2**63 - 1
        !  write (*,'(i20,b64/o24,z20)') ibig,ibig,ibig,ibig
     end do
  endif


  ! array allocation

  allocate ( x(nppm), y(nppm), z(nppm), ux(nppm), uy(nppm), uz(nppm), & 
       q(nppm), m(nppm), ax(nppm), ay(nppm), az(nppm), pot(nppm), work(nppm), &
       pepid(nppm), pelabel(nppm), pekey(nppm) )    ! Reserve particle array space N/NPE

  allocate ( nterm(nshortm), intlist(nintmax,nshortm), nodelist(nintmax,nshortm) )      ! Space for interaction lists


  allocate ( htable(0:maxaddress), all_addr(0:maxaddress), free_addr(maxaddress), point_free(0:maxaddress), &
       nbranches(num_pe+2), igap(num_pe+3), &
       treekey(size_tree), branch_key(size_tree), branch_owner(size_tree), &
       pebranch(size_tree), &
       requested_keys(size_tree, num_pe), fetched_keys(size_tree, num_pe) )

  all_addr = (/ (k,k=0,maxaddress) /)      ! List of all possible # table addresses
  htable%node = 0
  htable%key = 0
  htable%link = -1
  htable%leaves = 0
  htable%childcode = 0

  ! Allocate memory for tree node properties

  maxtwig = size_tree
  maxleaf = size_tree

  allocate ( first_child(-maxtwig:maxleaf), n_children(-maxtwig:maxleaf), node_level(-maxtwig:maxleaf) )

  allocate ( charge(-maxtwig:maxleaf), &                    ! charge
       abs_charge(-maxtwig:maxleaf), &                ! absolute charge
       xcoc(-maxtwig:maxleaf), ycoc(-maxtwig:maxleaf), zcoc(-maxtwig:maxleaf), &    ! centre of charge 
       xshift(-maxtwig:maxleaf), yshift(-maxtwig:maxleaf), zshift(-maxtwig:maxleaf), &    ! shift vector
       xdip(-maxtwig:maxleaf), ydip(-maxtwig:maxleaf), zdip(-maxtwig:maxleaf), &          ! dipole moment
       xxquad(-maxtwig:maxleaf), yyquad(-maxtwig:maxleaf), zzquad(-maxtwig:maxleaf), &       ! quadrupole moment
       xyquad(-maxtwig:maxleaf), yzquad(-maxtwig:maxleaf), zxquad(-maxtwig:maxleaf) )      !

  allocate ( pack_child(size_tree), get_child(size_tree) )    ! Multipole shipping buffers

  allocate (rhoe(0:ngx+1,0:ngy+1,0:ngz+1), rhoi(0:ngx+1,0:ngy+1,0:ngz+1) )   ! Field arrays


!  MPI stuff

  allocate (send_counts(num_pe+2), send_strides(num_pe+3), recv_counts(num_pe+2), recv_strides(num_pe+3) )  ! buf lengths and strides
  allocate ( stat_pe(MPI_STATUS_SIZE, num_pe), & ! status
                         pe_handle(2*num_pe), &  ! Handles for non-blocking comm  2*num_pe
                         send_key_handle(num_pe), &  ! (num_pe)
                         recv_key_handle(num_pe), &
                         send_child_handle(num_pe), & ! (num_pe)
                         recv_child_handle(num_pe) )



 ! Create new contiguous datatype for shipping particle properties (12 arrays)
 
  blocklengths(1:nprops_particle) = 1   


  types(1:9) = MPI_REAL8
  types(10) = MPI_INTEGER8
  types(11:12) = MPI_INTEGER

  call MPI_ADDRESS( get_props%x, receive_base, ierr )  ! Base address for receive buffer
  call MPI_ADDRESS( ship_props%x, send_base, ierr )  ! Base address for send buffer



  call MPI_ADDRESS( ship_props%x, address(1), ierr )
  call MPI_ADDRESS( ship_props%y, address(2), ierr )
  call MPI_ADDRESS( ship_props%z, address(3), ierr )
  call MPI_ADDRESS( ship_props%ux, address(4), ierr )
  call MPI_ADDRESS( ship_props%uy, address(5), ierr )
  call MPI_ADDRESS( ship_props%uz, address(6), ierr )
  call MPI_ADDRESS( ship_props%q, address(7), ierr )
  call MPI_ADDRESS( ship_props%m, address(8), ierr )
  call MPI_ADDRESS( ship_props%work, address(9), ierr )
  call MPI_ADDRESS( ship_props%key, address(10), ierr )
  call MPI_ADDRESS( ship_props%label, address(11), ierr )
  call MPI_ADDRESS( ship_props%pid, address(12), ierr )

  displacements(1:nprops_particle) = address(1:nprops_particle) - send_base  !  Addresses relative to start of particle (receive) data

  call MPI_TYPE_STRUCT( nprops_particle, blocklengths, displacements, types, mpi_type_particle, ierr )   ! Create and commit
  call MPI_TYPE_COMMIT( mpi_type_particle, ierr)

 ! Create new contiguous datatype for shipping multipole properties (18 arrays)
 
  blocklengths(1:nprops_multipole) = 1   


  types(1) = MPI_INTEGER8
  types(2:3) = MPI_INTEGER
  types(4) = MPI_INTEGER8
  types(5:18) = MPI_REAL8

  call MPI_ADDRESS( node_dummy%key, send_base, ierr )  ! Base address for send buffer



  call MPI_ADDRESS( node_dummy%key, address(1), ierr )
  call MPI_ADDRESS( node_dummy%byte, address(2), ierr )
  call MPI_ADDRESS( node_dummy%leaves, address(3), ierr )
  call MPI_ADDRESS( node_dummy%next, address(4), ierr )
  call MPI_ADDRESS( node_dummy%q, address(5), ierr )
  call MPI_ADDRESS( node_dummy%absq, address(6), ierr )
  call MPI_ADDRESS( node_dummy%xcoc, address(7), ierr )
  call MPI_ADDRESS( node_dummy%ycoc, address(8), ierr )
  call MPI_ADDRESS( node_dummy%zcoc, address(9), ierr )
  call MPI_ADDRESS( node_dummy%xdip, address(10), ierr )
  call MPI_ADDRESS( node_dummy%ydip, address(11), ierr )
  call MPI_ADDRESS( node_dummy%zdip, address(12), ierr )
  call MPI_ADDRESS( node_dummy%xxquad, address(13), ierr )
  call MPI_ADDRESS( node_dummy%yyquad, address(14), ierr )
  call MPI_ADDRESS( node_dummy%zzquad, address(15), ierr )
  call MPI_ADDRESS( node_dummy%xyquad, address(16), ierr )
  call MPI_ADDRESS( node_dummy%yzquad, address(17), ierr )
  call MPI_ADDRESS( node_dummy%zxquad, address(18), ierr )
  
  displacements(1:nprops_multipole) = address(1:nprops_multipole) - send_base   !  Addresses relative to start of particle (receive) data

  call MPI_TYPE_STRUCT( nprops_multipole, blocklengths, displacements, types, mpi_type_multipole, ierr )   ! Create and commit
  call MPI_TYPE_COMMIT( mpi_type_multipole, ierr)
  

end subroutine setup
