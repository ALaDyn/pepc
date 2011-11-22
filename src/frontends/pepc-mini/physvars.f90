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
  real :: force_const    ! force constant depending on unit system
  real :: qe = -1.0
  real :: qi = 1.0         ! electron, ion charge
  real :: mass_e = 1.0, mass_i = 1856.0   ! electron, ion mass
  real :: r_sphere       ! initial radius of plasma sphere
  real :: eps            ! potential/force law cutoff
  real :: q_factor = 1.0       ! Charge factor

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
  integer :: idim=3  ! # dimensions (velocity and position updates)
  integer :: ispecial       ! Switch to select special electron configs 
  integer :: weighted
  integer :: debug_level =0 ! Debug level for printed O/P

   real :: dt             ! timestep
   real :: trun           ! total run time including restarts
   integer :: nt, itime   ! # timesteps and current timestep
   integer :: db_level = 1  ! printed o/p debug level
   integer :: curve_type = 0 !< type of space-filling curve, 0=z-curve, 1=Hilbert-curve

   integer :: ifile_cpu    ! O/P stream


end module physvars





