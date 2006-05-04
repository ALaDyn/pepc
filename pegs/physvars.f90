!
!  Module defining physics constants and arrays
!
module physvars

  real, allocatable ::  rhoe(:,:,:), rhoi(:,:,:), phi_g(:,:,:)    ! field arrays
  real, allocatable ::  Ex_g(:,:,:),   Ey_g(:,:,:),  Ez_g(:,:,:)

  real, allocatable :: xslice(:),  yslice(:),  zslice(:)    ! Rezoning slice
  real, allocatable :: uxslice(:), uyslice(:), uzslice(:)   ! velocity
  real, allocatable :: qslice(:),  mslice(:)                ! charge and mass

  real :: x_star(2), y_star(2), z_star(2)  ! Position of stars
  real :: ux_star(2), uy_star(2), uz_star(2) ! Velocity of stars
  real :: m_star(2) ! Star masses
  real :: pot_star(2) ! potentials

  !  physics data

  real, parameter :: pi = 3.141592654

  real, parameter :: rsun =   6.960e8                 ! solar radius
  real, parameter :: sunmass =1.989e30                ! solar mass
  real, parameter :: au =     1.496e11                ! astronomical unit
  real, parameter :: year =   31557600e0              ! year in seconds
  real, parameter :: gammasi =6.673e-11*year**2       ! gravitational constant in SI units
  real, parameter :: gamma =  gammasi*sunmass/(au**3) ! gravitational constant in AU/YEAR/solarmass
  real, parameter :: sigmasi =5.671e-8*year**3        ! viscosity constant in SI units
  real, parameter :: sigma =  sigmasi/sunmass         ! viscosity constant
  real, parameter :: kbmusi = 1.6629e4*year**2        ! Boltzmann constant in SI units
  real, parameter :: kbmu =   kbmusi*0.5/(2.35*au**2)  ! Boltzmann constant
  real, parameter :: gam  =   5.e0/3.e0               ! eqs

  real :: mdisc          ! disc mass

  integer :: ni,ne, np_beam
  integer :: nep, nip     ! # particles/electrons/ions per PE
  real :: vte, vti       ! electron, ion thermal velocities
  real :: Te_keV, Ti_keV ! electron, ion emperatures in keV
  real :: T_scale = 1       ! factor for rescaling Te after restart 
  real :: bond_const     ! bonding force constant for ion crystal
  real :: mass_ratio     ! ion:electron mass ratio
  real :: qe, qi         ! electron, ion charge
  real :: mass_e, mass_i   ! electron, ion mass
  real :: r_sphere       ! initial radius of plasma sphere
  real :: x_plasma       ! initial plasma length (slab or disc targets)
  real :: y_plasma       ! initial plasma y-width (slab)
  real :: z_plasma       ! initial plasma z-width (slab)
  real :: box_centre(3) ! vector defining centre of plasma target
  real :: x_crit         ! critical surface
  real :: rho0           ! electron density (1)
  real :: Vplas          ! plasma volume
  real :: a_ii           ! mean ion spacing
  real :: r_neighbour    ! nearest-neighbour search radius


 
  real :: xl, yl, zl    ! initial simulation region side lengths
  real ::  theta  ! clumping parameter
  real :: force_const    ! force constant depending on unit system
  real :: eps            ! potential/force law cutoff
  real :: displace(3)    ! particle displacement vector for restart (change of view box)

  ! Control stuff

  integer :: mc_steps
  integer :: initial_config = 4  ! Switch for initial configuration (positions, velocities)
  integer :: beam_config = 0 ! Switch for particle beam mode
  integer :: scheme = 1 ! Canonical ensemble switch: 2-4= const. Te dynamics
  integer :: particle_bcs = 1 ! Particle BC switch: 1=open, 2=reflective

  ! Control

  logical :: balance=.true.   ! Balances particles in || sort according to work load
  logical :: walk_balance=.true.   ! Does detailed balancing of particle groups
  logical :: restart = .false.  ! Restart switch: config read from parts_all.in
  logical :: coulomb = .true.  ! Compute Coulomb forces
  logical :: vis_on=.false.  ! online visualisation on/off
  logical :: steering = .false.  ! VISIT steering switch

  integer :: nt, itime   ! # timesteps and current timestep
  integer :: itime_start ! restart time stamp
  integer :: idump=100   ! output frequency (timesteps)
  integer :: iprot=1       ! protocoll frequency
  integer :: ivis=100        ! frequency for particle shipping to VISIT
  integer :: ivis_fields=100    !  frequency for field shipping to VISIT
  integer :: idens=10       ! frequency for computing densities (tracking) 
  integer :: ngx, ngy, ngz  ! Plot grid dimensions
  integer :: nmerge  ! dataset-merging factor
  integer :: debug_tree = 1
  real :: dt             ! timestep
  real :: trun           ! total run time including restarts

  integer :: my_rank
  integer :: n_cpu
  integer :: np_mult = 1
  integer :: fetch_mult = 1

end module physvars
