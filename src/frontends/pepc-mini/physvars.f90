module physvars
  use module_pepc_types, only : t_particle
  use module_interaction_specific_types, only : t_particle_results
  implicit none

  real, parameter :: pi=3.141592654

  real*8, allocatable :: u(:,:), m(:)
  type(t_particle), allocatable :: particles(:)

  !  physics data
  integer :: ni, ne       !  # ions, electrons
  integer :: nep, nip     ! # particles/electrons/ions per PE
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

!  Associated MPI stuff

  integer :: my_rank       ! Rank of current task
  integer :: n_cpu   ! # cpus used by program

! Control stuff
  integer :: ispecial       ! Switch to select special electron configs 

  real :: dt             ! timestep
  real :: trun           ! total run time including restarts
  integer :: nt, itime   ! # timesteps and current timestep

  integer :: ifile_cpu    ! O/P stream


end module physvars





