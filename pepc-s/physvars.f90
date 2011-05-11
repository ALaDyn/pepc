module physvars

  real, parameter :: pi=3.141592654

! particle arrays
  real*8, allocatable ::   work(:)          ! work load (interaction list length)
  integer, allocatable ::  pelabel(:)  ! particle label
    
  !  physics data
  integer :: ni, ne       !  # ions, electrons
  integer :: nep, nip     ! # particles/electrons/ions per PE
  real :: force_const    ! force constant depending on unit system
  real :: eps            ! potential/force law cutoff

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
  integer :: weighted
  integer :: debug_level =0 ! Debug level for printed O/P

   integer :: itime   ! # timesteps and current timestep
   integer :: db_level = 1  ! printed o/p debug level
   integer :: curve_type = 0 !< type of space-filling curve, 0=z-curve, 1=Hilbert-curve

   integer :: ifile_cpu    ! O/P stream


end module physvars





