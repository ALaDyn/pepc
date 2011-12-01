module physvars

  use treetypes

  type(t_particle), allocatable :: vortex_particles(:)

  real*8, parameter :: pi=3.141592654
    
  !  physics data
  integer  :: n   ! # vortices
  integer  :: np  ! # vortices per PE
  real     :: force_const    ! force constant depending on unit system
  real     :: eps            ! potential/force law cutoff
  real     :: h, m_h         ! initial particle distance and mesh width for remeshing
  integer  :: rem_freq       ! remeshing frequence
  real*8   :: kernel_c       ! mod. remeshing kernel parameter
  real*8   :: thresh         ! vorticity threshold: particles with lower vorticity mag. will be kicked out (mandatory to avoid zero abs_charge)
  real     :: nu             ! viscosity parameter
  real     :: rmax    ! radius of torus segment
  real     :: rl      ! temp radius of current circle
  real     :: r_torus ! radius of torus
  integer  :: nc      ! # circles per torus segment
  integer  :: nphi    ! # torus segments
  integer  :: ns      ! # particles per torus segment
  real     :: g       ! # vorticity amplifier
  real, dimension(3) :: torus_offset  ! shifts coords of both tori (one with +, one with -)


 ! tree stuff
  real :: theta       ! Clumping parameter
  integer :: mac      ! MAC (default=BH)

!  Variables needing 'copy' for tree routines
  integer :: npart_total  ! Total # particles (npart)
  integer :: np_local 
  integer :: nppm  ! Total # particles (npart)
  real    :: np_mult

!  Associated MPI stuff
  integer :: my_rank       ! Rank of current task
  integer :: n_cpu   ! # cpus used by program

! Control stuff
  integer :: ispecial       ! Switch to select special electron configs 
  integer :: weighted       ! load balancing switch
  real :: dt, ts, te        ! timestep, start-time, end-time
  integer :: nt             ! # timesteps and current timestep
  integer :: rk_stages      ! # Runge-Kutta stages
  integer :: db_level       ! printed o/p debug level
  integer :: curve_type     !< type of space-filling curve, 0=z-curve, 1=Hilbert-curve

  integer :: ifile_cpu    ! O/P stream


end module physvars





