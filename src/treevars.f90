!
!   Header module for treemp: global arrays and variables
!
!   $ Revision:   $ 
!   All key variables defined 64 bit (8 byte)
!
 
module treevars

  use my_mpidefs          ! Machine-specific MPI defs (4-byte IBM, 8-byte Cray)

 ! Constants

  real, parameter :: pi=3.141592654
  integer, dimension(0:7) :: bitarr = (/ 0,1,2,3,4,5,6,7 /)    ! Array of bit positions

 ! Hash table datatype
 
  type hash
     integer   :: node          ! Address of particle/pseudoparticle data
     integer*8 :: key           ! Key
     integer   :: link          ! Pointer to next empty address in table in case of collision
     integer   :: leaves        ! # leaves contained within twig (=1 for leaf, npart for root)
     integer   :: childcode     ! Byte code indicating position of children (twig node); particle label (leaf node)
     integer*8 :: next          ! Pointer to next key to examine in tree-walk
     integer   :: owner         ! Node owner (for branches)
  end type hash

  !  hash table

  type (hash), allocatable :: htable(:)
  
  
 ! Data structure for shipping single particles
  
  type particle
     real :: x    ! coords
     real :: y
     real :: z
     real :: ux    ! momenta
     real :: uy
     real :: uz 
     real :: q     ! charge
     real :: m     ! mass
     real :: work  ! work load from force sum
     integer*8 :: key           ! Key
     integer :: label    ! label
     integer :: pid      ! owner
  end type particle

  type (particle) :: ship_props, get_props



 ! Data structure for shipping multiple moments of child nodes
  
  type multipole
     integer*8 :: key     ! key
     integer   :: byte    ! byte code
     integer   :: leaves  ! # leaves contained
     integer*8 :: next    ! next key on walk
     real :: q    	! net charge sum
     real :: absq  	!  absolute charge sum
     real :: xcoc  	! centre of charge
     real :: ycoc    
     real :: zcoc
     real :: xdip  	! dipole moment
     real :: ydip
     real :: zdip
     real :: xxquad  	! quadrupole moment
     real :: yyquad  
     real :: zzquad
     real :: xyquad
     real :: yzquad
     real :: zxquad
  end type multipole

  type (multipole), allocatable :: pack_child(:), get_child(:)
  type (multipole) :: node_dummy

!  Associated MPI stuff

  integer, parameter :: nprops_particle=12, &    ! # particle properties to ship
  			nprops_multipole=18      ! Number of multipole properties to ship
  integer, dimension(nprops_multipole) :: blocklengths, displacements, types, address
  integer :: send_base, receive_base, mpi_type_particle, mpi_type_multipole
  

  !  tree variables

  integer*8, allocatable :: &
                                requested_keys(:,:), &  ! Local multipole nodes required elsewhere
                                fetched_keys(:,:), &  ! Remote nodes fetched during tree walk
                                treekey(:), &       ! keys of all twig and leaf nodes
                                branch_key(:), &    ! keys of branch nodes covering all domains
                                pebranch(:)         ! keys of branch nodes covering local domain


  integer, allocatable :: &
                                nbranches(:), &       ! # branches in local domain
                                igap(:), &    !  stride lengths of local branch arrays
                                branch_owner(:), &    ! owners of branch nodes covering all domains
                                all_addr(:), &  ! List of all possible #table addresses
                                free_addr(:), &    ! List of free #table addresses (for HASHENTRY routine)
                                point_free(:), &   ! Pointer to free address index

                                nreqs_total(:), &    ! total # nodes requested from local PE during tree walk
                                nfetch_total(:)   ! total # non-local nodes fetched during tree walk 

  real, allocatable  :: &                ! Tree node properties:
                               charge(:), &                          ! charge
                               abs_charge(:), &                      ! absolute charge
                               xcoc(:), ycoc(:), zcoc(:), &          ! centre of charge 
                               xshift(:), yshift(:), zshift(:), &    ! shift vector
                               xdip(:), ydip(:), zdip(:), &          ! dipole moment
                               xxquad(:), yyquad(:), zzquad(:), &    ! quadrupole moment
                               xyquad(:), yzquad(:), zxquad(:)       !

  integer*8, allocatable ::    first_child(:)   ! key of first child

  integer, allocatable ::       n_children(:), &     ! # children
                                node_level(:)       ! refinement level
 

  integer*8 ::  hashconst, &   ! hashing constants
	        hashchild=7, &
                iplace         ! value of place holder bit = 2^(2*nlev)

  integer :: idim, &           ! # dimensions (2, 3)
             nlev, &           ! max refinement level
             nbaddr, &         ! # bits in hashing function
             nleaf, &          ! total # leaf nodes in local #table 
             ntwig, &          ! total # twig nodes in local #table
             nleaf_me, &       ! total # leaves in local domain
             ntwig_me, &       ! total # twigs in local domain
             nlist, &          ! # particles/PE + boundary bodies (1 or 2)
             nnodes, &         ! nleaf + ntwig
             nnodes_me, &      ! nleaf_me + ntwig_me
             nnodes_pw, &      ! # nnodes prior to tree-walk
             nleaf_pw, &       ! # leaves prior to walk
             ntwig_pw, &       ! # twigs prior to walk
             nbranch, &        ! min # branch nodes covering local domain
             nbranch_sum, &    ! total # branch nodes covering all domains
             nintmax, &        ! max # terms allowed in interaction list
             max_list_length, & ! current max list length
             maxaddress, &     ! max address allowed in #table
             size_tree, &      ! array space needed for local tree
             free_lo, &        ! min address allowed for resolving collisions
	     tablehigh, &      ! highest current address in #table 
             sum_unused, &     ! # free addresses
             ifile, ipefile, &             ! local O/P stream
             npartm, &         ! absolute max # particles
             nppm, &           ! max # particles/PE
             nshortm, &        ! shortlist length
             iused, &          ! counter for collision resolution array free_addr()
  	     nmerge = 1        ! merge factor for data sets


  real :: xmin, xmax    ! box limits
  real :: ymin, ymax  
  real :: zmin, zmax
  real :: boxsize       ! length of box

 ! Force control
  logical :: load_balance=.true.   ! Balances particles in || sort according to work load
  logical :: walk_balance=.false.  ! Controls balancing of shortlists in walk/force sum

  ! Debugging switches (all off by default)
  logical :: tree_debug=.false.
  logical :: domain_debug = .false.
  logical :: branch_debug=.false.
  logical :: prefetch_debug=.false.
  logical :: walk_debug=.false.
  logical :: force_debug=.false.
  logical :: dump_tree=.false.
  logical :: perf_anal=.false.  ! Performance analysis mode: turns off all diagnostic routines

  !  particle data - dimensions filled in when #PEs known

  integer, allocatable :: pelabel(:)             ! particle label
  real, allocatable :: x(:),  y(:),  z(:), &     ! position
                      ux(:), uy(:), uz(:), &     ! velocity
                              q(:),  m(:), &     ! charge and mass
			ax(:), ay(:), az(:), &      ! accelerations
			pot(:), &	         ! potential
			work(:)                  ! interaction work load 

  integer*8, allocatable ::   pekey(:), &  ! local particle keys                             
                              intlist(:,:) ! interaction key-list

  integer, allocatable ::     pepid(:), & ! owner
                              nterm(:),  nodelist(:,:) ! interaction node-list

  real, allocatable ::  rhoe(:,:,:), rhoi(:,:,:)  ! field arrays

  !  physics data

  integer :: npart,ni,ne
  integer :: npp, nep, nip     ! # particles/electrons/ions per PE
  real :: xl, yl, zl, theta
  real :: vte, vti       ! electron, ion thermal velocities
  real :: Te_keV, Ti_keV ! electron, ion emperatures in keV
  real :: T_scale = 1       ! factor for rescaling Te after restart 
  real :: force_const    ! force constant depending on unit system
  real :: bond_const     ! bonding force constant for ion crystal
  real :: mass_ratio     ! ion:electron mass ratio
  real :: qe, qi         ! electron, ion charge
  real :: mass_e, mass_i   ! electron, ion mass
  real :: r_sphere       ! initial radius of plasma sphere
  real :: x_plasma       ! initial plasma length (slab or disc targets)
  real :: y_plasma       ! initial plasma width (slab)
  real :: plasma_centre(3) ! vector defining centre of plasma target
  real :: x_crit         ! critical surface
  real :: rho0           ! electron density (1)
  real :: Vplas          ! plasma volume
  real :: a_ii           ! mean ion spacing
  real :: r_neighbour    ! nearest-neighbour search radius
  real :: eps            ! potential/force law cutoff
  real :: delta_mc       ! step size for MC config
  real :: displace(3)    ! particle displacement vector for restart (change of view box)

  ! particle beam stuff
  integer :: np_beam    ! # beam particles
  integer :: nb_pe    ! # beam particles per PE
  real :: x_beam        ! beam length
  real :: r_beam        ! beam radius
  real :: start_beam    ! starting position (x-axis)
  real :: u_beam, theta_beam, phi_beam        ! beam velocity and angles
  real :: rho_beam      ! beam density as fraction of plasma density (=1)
  real :: mass_beam     ! mass of beam particles
  real :: qeb           ! beam particle charge

 ! laser parameters
  real :: tpulse        ! pulse duration (in units of 1/omega_p)
  real :: sigma         ! 1/e pulse width (c/omega_p)
  real :: vosc          ! pump strength    (c)
  real :: omega         ! frequency  (omega_p)
  real :: lambda        ! laser wavelength
  real :: focus(3)      ! centre of focal spot
  real :: tlaser        ! run time after laser switched on (1/omega_p)
  
  
! Control stuff
  logical :: vis_on=.true.  ! online visualisation on/off
  logical :: mc_init = .false. ! MC initialisation switch
  logical :: restart = .false.  ! Restart switch: config read from parts_all.in
  logical :: coulomb = .true.  ! Compute Coulomb forces
  logical :: bonds = .false. ! Include SHO bonding potential for ions
  logical :: lenjones = .false. ! Include short-range Lennard-Jones potential
  logical :: steering = .false.  ! VISIT steering switch
  integer :: mc_steps
  integer :: initial_config = 4  ! Switch for initial configuration (positions, velocities)
  integer :: beam_config = 0 ! Switch for particle beam mode
  integer :: ensemble = 1 ! Canonical ensemble switch: 2-4= const. Te dynamics
  integer :: particle_bcs = 1 ! Particle BC switch: 1=open, 2=reflective
  
   real :: dt             ! timestep
   real :: trun           ! total run time including restarts
   real :: convert_fs     ! conversion factor from wp^-1 to fs
   integer :: nt, itime   ! # timesteps and current timestep
   integer :: itime_start ! restart time stamp
   integer :: idump       ! output frequency (timesteps)
   integer :: iprot=1       ! protocoll frequency
   integer :: ivis        ! frequency for particle shipping to VISIT
   integer :: ivis_fields    !  frequency for field shipping to VISIT
   integer :: ngx, ngy, ngz  ! Plot grid dimensions

end module treevars










