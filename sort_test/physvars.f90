module physvars

  real, parameter :: pi=3.141592654

  type particle_p1
     real*8 :: x    ! coords
     real*8 :: y
     real*8 :: z
     real*8 :: ux    ! momenta
     real*8 :: uy
     real*8 :: uz 
     real*8 :: q     ! charge
     real*8 :: m     ! mass
     real*8 :: work  ! work load from force sum
     integer*8 :: label    ! label
  end type particle_p1
  
  integer :: mpi_type_particle_p1


! particle arrays
  real*8, allocatable :: x(:),  y(:),  z(:), &     ! position
                         ux(:), uy(:), uz(:), &     ! velocity
                         q(:),  m(:), &     ! charge and mass
                         Ex(:), Ey(:), Ez(:), &   ! E-field
                         pot(:), &	         ! scalar potential
                         work(:)	         ! work load (interaction list length)

  integer*8, allocatable ::  pelabel(:) ! particle label
    
  real, allocatable ::  rhoe_loc(:,:,:), rhoi_loc(:,:,:)  ! field arrays for time-averages
  real, allocatable ::  rhoi(:,:,:), rhoe(:,:,:)
  real, allocatable ::  ex_loc(:,:,:), ey_loc(:,:,:), ez_loc(:,:,:)  ! E-field 
  real, allocatable ::  bx_loc(:,:,:), by_loc(:,:,:), bz_loc(:,:,:)  ! B-field 
  real, allocatable ::  jxe_loc(:,:,:), jye_loc(:,:,:), jze_loc(:,:,:)  ! elec current

  real*4, allocatable :: vbuffer(:,:), vbuf_local(:,:) !< vis buffers

  real, allocatable :: xslice(:),  yslice(:),  zslice(:), &     ! Rezoning slice
                      uxslice(:), uyslice(:), uzslice(:), &     ! velocity
                              qslice(:),  mslice(:)    ! charge and mass

  !  physics data

  integer*8 :: ni, ne       !  # ions, electrons
  integer :: nep, nip     ! # particles/electrons/ions per PE
  real :: xl, yl, zl      ! box size
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
  real :: y_plasma       ! initial plasma y-width (slab)
  real :: z_plasma       ! initial plasma z-width (slab)
  real :: plasma_centre(3) ! vector defining centre of plasma target
  real :: x_crit         ! critical surface
  real :: x_offset       ! coordinate offset
  real :: z_offset       ! coordinate offset
  real :: rho0           ! electron density (1)
  real :: rho_track      ! tracking density for x_crit (/nc)
  real :: rho_upper      ! shelf/profile density above x_crit (/nc)
  real :: Vplas          ! plasma volume
  real :: a_ii           ! mean ion spacing
  real :: r_neighbour    ! nearest-neighbour search radius
  real :: eps            ! potential/force law cutoff
  real :: err_f          ! user-specified force error tolerance
  real :: delta_mc       ! step size for MC config
  real :: displace(3)    ! particle displacement vector for restart (change of view box)
  real :: uthresh        ! velocity (u^2) threshold for vis_parts
  real :: rho_min        ! min density for exponential ramp
  real :: lolam          ! L/lambda density scale-length
  real :: q_factor       ! Charge factor
  real :: fnn            ! Near-neighbour factor
  real*8 :: Ukine          ! Electron kinetic energy
  real*8 :: Ukini          ! Ion kinetic energy
  real*8 :: Umagnetic           ! Magnetic energy
  real*8 :: Ubeam          ! Beam energy

  real :: nu_ei ! norm. collision frequency
  real :: sigma_e ! electrical conductivity
  real :: range_hot ! hot electron range
  real :: t_sat ! saturation time for n_hot
  real :: t_foil ! traversal time of hot e front

 ! tree stuff
  real :: theta       ! Clumping parameter
  integer :: mac = 0  ! MAC (default=BH)

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
  real :: theta_inc     ! angle of incidence
  real :: focus(3)      ! centre of focal spot
  real :: tlaser        ! run time after laser switched on (1/omega_p)
  real :: elaser        ! deposited laser energy
  real :: propag_laser  ! distance travelled by laser after rezoning
  real :: intensity     ! normalised intensity = 0.5*vosc^2*omega^2
  real :: window_min    ! start of wakefield plasma
  real :: rezone_frac=0.75     ! Fraction of box to cross before rezoning switched on



!  Variables needing 'copy' for tree routines
  integer*8 :: npart_total  ! Total # particles (npart)
!  integer :: npp  ! Total # particles (npart)
  integer :: np_local 
  integer :: nppm  ! Total # particles (npart)
  real :: np_mult=1.5
  integer :: fetch_mult=2

!  Associated MPI stuff

  integer :: my_rank       ! Rank of current task
  integer :: n_cpu   ! # cpus used by program

! Control stuff
  logical :: vis_on=.true.  ! online visualisation on/off
  logical :: mc_init = .false. ! MC initialisation switch
  logical :: restart = .false.  ! Restart switch: config read from parts_all.in
  logical :: coulomb = .true.  ! Compute Coulomb forces
  logical :: bfield_on = .false. ! Include self-induced B-field in force calculation
  logical :: bonds = .false. ! Include SHO bonding potential for ions
  logical :: lenjones = .false. ! Include short-range Lennard-Jones potential
  logical :: steering = .false.  ! VISIT steering switch
  logical :: target_dup = .false. ! Target duplication switch
  logical :: ramp = .false.  ! profile-ramp switch
  integer :: mc_steps
  integer :: system_config = 1  ! Switch for initial configuration (positions, velocities)
  integer :: target_geometry = 0  ! Geometry for plasma target
  integer :: idim=3  ! # dimensions (velocity and position updates)
  integer :: beam_config_in = 0 ! Particle or laser beam switch including variations 
  integer :: beam_config = 0 ! Reduced switch for particle or laser beam 
  integer :: ispecial       ! Switch to select special electron configs 
  integer :: scheme = 1 ! Integrator scheme switch: 2-4= const. Te dynamics, 6=EM
  integer :: particle_bcs = 1 ! Particle BC switch: 1=open, 2=reflective
  integer :: debug_level =0 ! Debug level for printed O/P

   real :: dt             ! timestep
   real :: trun           ! total run time including restarts
   real :: convert_fs     ! conversion factor from wp^-1 to fs
   real :: convert_mu     ! conversion factor from c/wp to microns
   real :: convert_keV     ! conversion factor from norm energy to keV/particle
   integer :: nt, itime   ! # timesteps and current timestep
   integer :: itime_start ! restart time stamp
   integer :: idump       ! output frequency (timesteps)
   integer :: db_level = 1  ! printed o/p debug level
   integer :: iprot=1       ! protocoll frequency
   integer :: ivis=5        ! frequency for particle shipping to VISIT
   integer :: ivis_fields=10    !  frequency for field shipping to VISIT
   integer :: ivis_domains=10    !  frequency for domain shipping to VISIT
   integer :: vis_select = 1  !< select switch for particles
   integer :: nbuf_max=10000     !< Max vis buffer size
   integer :: ndom_max=1000     !< Max # domains
   integer :: attrib_max = 22 !< max # attributes for vis buffer
   real :: ops_per_sec = 0.  !< Work load (interactions per sec)


   integer :: itrack       ! frequency for computing ion density (tracking)
   integer :: navcycle     ! # timesteps in a laser cycle 
   integer :: ngx, ngy, ngz  ! Plot grid dimensions
   integer :: ifile_cpu    ! O/P stream

   ! constrain
   real :: constrain_proof ! quality of getting crossing points
   real :: len_tripod      ! length of tripod
   real :: number_faces    ! # faces to constrain the particles in
   integer :: struct_step

  character*10 :: plasma_configs(0:2)= (/ &
       'no plasma ','plasma    ','special   ' /)

  character*10 :: geometries(0:10)= (/ &
       'slab      ', & 
       'sphere    ','disc      ','wire      ','ellipsoid ','wedge     ', &
       'hemisphere','hollow sph','hollow hsp','          ','special   ' /)
  character*7 :: beam_configs(0:9)=(/ &
       'eqm    ','beam   ','i-beam ','laser-u','ES pond','LWFA   ', &
       'EMplane','EM pond',' dust  ','       ' /)
  character*7 :: schemes(1:6)=(/ &
       'const U','Te, Ti ','glob Te','loc  Te','Ti only','Full 3V' /)

end module physvars




