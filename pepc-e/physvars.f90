module physvars

  real, parameter :: pi=3.141592654

! particle arrays
  real*8, allocatable :: x(:),  y(:),  z(:), &     ! position
                      ux(:), uy(:), uz(:), &     ! velocity
                              q(:),  m(:), &     ! charge and mass
			Ex(:), Ey(:), Ez(:), &   ! E-field
			pot(:), &	         ! scalar potential
			work(:)	         ! work load (interaction list length)

  integer, allocatable ::  pelabel(:)  ! particle label

  !  physics data
  integer :: ni, ne       !  # ions, electrons
  integer :: nep, nip     ! # particles/electrons/ions per PE
  real :: xl, yl, zl      ! box size
  real :: vte, vti       ! electron, ion thermal velocities
  real :: Te_keV, Ti_keV ! electron, ion emperatures in keV
  real :: T_scale = 1       ! factor for rescaling Te after restart 
  real :: force_const    ! force constant depending on unit system
  real :: mass_ratio     ! ion:electron mass ratio
  real :: qe, qi         ! electron, ion charge
  real :: mass_e, mass_i   ! electron, ion mass
  real :: r_sphere       ! initial radius of plasma sphere
  real :: x_plasma       ! initial plasma length (slab or disc targets)
  real :: y_plasma       ! initial plasma y-width (slab)
  real :: z_plasma       ! initial plasma z-width (slab)
  real :: plasma_centre(3) ! vector defining centre of plasma target
  real :: x_crit         ! critical surface
  real :: rho0           ! electron density (1)
  real :: Vplas          ! plasma volume
  real :: a_ii           ! mean ion spacing
  real :: eps            ! potential/force law cutoff
  real :: q_factor       ! Charge factor
  real*8 :: Ukine          ! Electron kinetic energy
  real*8 :: Ukini          ! Ion kinetic energy
  real*8 :: Umagnetic           ! Magnetic energy
  real*8 :: Ubeam          ! Beam energy

  ! tree stuff
  real :: theta       ! Clumping parameter
  integer :: mac = 0  ! MAC (default=BH)

  !  Variables needing 'copy' for tree routines
  integer :: npart_total  ! Total # particles (npart)
  integer :: np_local 
  integer :: nppm  ! Total # particles (npart)
  real :: np_mult=1.5

  !  Associated MPI stuff
  integer :: my_rank       ! Rank of current task
  integer :: n_cpu   ! # cpus used by program

  ! Control stuff
  integer :: system_config = 1  ! Switch for initial configuration (positions, velocities)
  integer :: target_geometry = 0  ! Geometry for plasma target
  integer :: idim=3  ! # dimensions (velocity and position updates)
  integer :: ispecial       ! Switch to select special electron configs 
  integer :: weighted
  integer :: debug_level =0 ! Debug level for printed O/P

   real :: dt             ! timestep
   real :: trun           ! total run time including restarts
   real :: convert_keV     ! conversion factor from norm energy to keV/particle
   integer :: nt, itime   ! # timesteps and current timestep
   integer :: itime_in    ! timestep to read mpi-io checkpoint from in case of ispecial==-1
   integer :: idump, idump_vtk, idump_checkpoint, idump_binary ! output frequency (timesteps): ascii, vtk and mpi-io-checkpoint
   integer :: db_level = 1  ! printed o/p debug level
   integer :: curve_type = 0 !< type of space-filling curve, 0=z-curve, 1=Hilbert-curve

   integer :: ifile_cpu    ! O/P stream

   ! constrain
   integer :: number_faces    ! # faces to constrain the particles in

  character*10 :: plasma_configs(0:2)= (/ &
       'no plasma ','plasma    ','special   ' /)

  character*10 :: geometries(0:10)= (/ &
       'slab      ', & 
       'sphere    ','disc      ','wire      ','ellipsoid ','wedge     ', &
       'hemisphere','hollow sph','hollow hsp','          ','special   ' /)

end module physvars





