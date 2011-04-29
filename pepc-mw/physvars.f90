module physvars

! particle arrays
  real*8, allocatable :: x(:),  y(:),  z(:), &     ! position
                      ux(:), uy(:), uz(:), &     ! velocity
                              q(:),  m(:), &     ! charge and mass
            Ex(:), Ey(:), Ez(:), &   ! E-field
            Ax(:), Ay(:), Az(:), &   ! E-field
            Bx(:), By(:), Bz(:), &   ! E-field
			pot(:), &	         ! scalar potential
			work(:)	         ! work load (interaction list length)

  integer, allocatable ::  pelabel(:)  ! particle label
    
  real, allocatable ::  rhoe_loc(:,:,:), rhoi_loc(:,:,:)  ! field arrays for time-averages
  real, allocatable ::  rhoi(:,:,:), rhoe(:,:,:)
  real, allocatable ::  ex_loc(:,:,:), ey_loc(:,:,:), ez_loc(:,:,:)  ! E-field 
  real, allocatable ::  bx_loc(:,:,:), by_loc(:,:,:), bz_loc(:,:,:)  ! B-field 
  real, allocatable ::  jxe_loc(:,:,:), jye_loc(:,:,:), jze_loc(:,:,:)  ! elec current

  real, allocatable :: xslice(:),  yslice(:),  zslice(:), &     ! Rezoning slice
                      uxslice(:), uyslice(:), uzslice(:), &     ! velocity
                              qslice(:),  mslice(:)    ! charge and mass

  !  physics data

  integer :: ni, ne       !  # ions, electrons
  integer :: nep, nip     ! # particles/electrons/ions per PE
  real*8 :: maxdt(4)       ! maximum allowed dt from different constraints
  real*8 :: xl, yl, zl      ! box size
  integer :: ngx, ngy, ngz  ! Plot grid dimensions
  real*8 :: vte, vti       ! electron, ion thermal velocities
  real*8 :: Te = 0., Ti = 0. ! electron, ion emperatures in program units
  real*8 :: Te_eV = 0., Ti_eV = 0. ! electron, ion emperatures in electron Volts
  real*8 :: Te_K = 0., Ti_K = 0. ! electron, ion emperatures in Kelvin
  real :: force_const    ! force constant depending on unit system
  real*8 :: rhoe_nm3 = 1., rhoi_nm3 = 0.       ! number of electrons and ions per nm^3
  real*8 :: qe, qi         ! electron, ion charge
  real*8 :: mass_e, mass_i   ! electron, ion mass
  real*8 :: wpl_e, wpl_i !< electron and ion plasma frequency
  real*8 :: lambdaD_e, lambdaD_i !< electron and ion Debye length
  real*8 :: r_sphere       ! initial radius of plasma sphere
  real*8 :: x_plasma       ! initial plasma length (slab or disc targets)
  real*8:: y_plasma       ! initial plasma y-width (slab)
  real*8 :: z_plasma       ! initial plasma z-width (slab)
  real*8 :: plasma_centre(3) ! vector defining centre of plasma target
  real*8 :: Vplas          ! plasma volume
  real*8 :: a_ii, a_ee           ! mean ion and electron spacing
  real*8 :: a_i            ! ion sphere radius
  real*8 :: physGamma      ! coupling parameter
  real :: eps            ! potential/force law cutoff
  integer :: Zion=1, Aion=1       ! ion charge and mass number

  real*8 :: Ukine          ! Electron kinetic energy
  real*8 :: Ukini          ! Ion kinetic energy

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
  integer :: MPI_COMM_PEPC ! MPI communicator to be used (will usually be a copy of MPI_COMM_WORLD)

! Control stuff
  integer :: idim=3  ! # dimensions (velocity and position updates)
  integer :: ispecial       ! Switch to select special electron configs 
  integer :: weighted
  integer :: debug_level = 0 ! Debug level for printed O/P

   logical, public :: restart = .false. !< Restart switch: config read from parts_all.in
   real*8 :: dt             ! timestep
   real*8 :: trun           ! total run time including restarts
   integer :: nt, itime   ! # timesteps and current timestep
   integer :: itime_in    ! timestep to read mpi-io checkpoint from in case of ispecial==-1
   integer :: idump, idump_vtk, idump_checkpoint, idump_binary ! output frequency (timesteps): ascii, vtk and mpi-io-checkpoint
   integer :: db_level = 1  ! printed o/p debug level

   integer :: ifile_cpu    ! O/P stream

end module physvars





