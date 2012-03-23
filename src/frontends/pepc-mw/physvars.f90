! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2012 Juelich Supercomputing Centre, 
!                         Forschungszentrum Juelich GmbH,
!                         Germany
! 
! PEPC is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! PEPC is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
! 
! You should have received a copy of the GNU Lesser General Public License
! along with PEPC.  If not, see <http://www.gnu.org/licenses/>.
!

module physvars
  use module_pepc_types

! particle arrays
  type(t_particle), allocatable :: particles(:)     ! position
  real*8, allocatable :: energy(:,:)         ! potential, kinetic, and total energy
    
  real, allocatable ::  rhoe_loc(:,:,:), rhoi_loc(:,:,:)  ! field arrays for time-averages
  real, allocatable ::  rhoi(:,:,:), rhoe(:,:,:)
  real, allocatable ::  ex_loc(:,:,:), ey_loc(:,:,:), ez_loc(:,:,:)  ! E-field 
  real, allocatable ::  bx_loc(:,:,:), by_loc(:,:,:), bz_loc(:,:,:)  ! B-field 
  real, allocatable ::  jxe_loc(:,:,:), jye_loc(:,:,:), jze_loc(:,:,:)  ! elec current

  !  physics data

  integer :: rngseed = 13
  integer :: ne = 200     !  # electrons
  integer :: ni = 0       !  # ions
  real*8 :: maxdt(4)       ! maximum allowed dt from different constraints
  real*8 :: xl = 1.      ! box size
  real*8 :: yl = 1.      ! box size
  real*8 :: zl = 1.      ! box size
  integer :: ngx = 25, ngy = 25, ngz = 25  ! Plot grid dimensions
  real*8 :: vte, vti       ! electron, ion thermal velocities
  real*8 :: Te = 0., Ti = 0. ! electron, ion emperatures in program units
  real*8 :: Te_eV = 50., Ti_eV = 10. ! electron, ion emperatures in electron Volts
  real*8 :: Te_K = 0., Ti_K = 0. ! electron, ion emperatures in Kelvin
  real*8 :: force_const = 1.   ! force constant depending on unit system
  real*8 :: rhoe_nm3 = 1., rhoi_nm3 = 0.       ! number of electrons and ions per nm^3
  real*8 :: qe, qi         ! electron, ion charge
  real*8 :: mass_e, mass_i   ! electron, ion mass
  real*8 :: wpl_e, wpl_i !< electron and ion plasma frequency
  real*8 :: lambdaD_e, lambdaD_i !< electron and ion Debye length
  real*8 :: r_sphere = 4.      ! initial radius of plasma sphere
  real*8 :: x_plasma = 1.      ! initial plasma length (slab or disc targets)
  real*8 :: y_plasma = 1.      ! initial plasma y-width (slab)
  real*8 :: z_plasma = 1.      ! initial plasma z-width (slab)
  real*8 :: particle_shift = 0.025
  real*8 :: plasma_centre(3) ! vector defining centre of plasma target
  real*8 :: Vplas          ! plasma volume
  real*8 :: a_ii, a_ee           ! mean ion and electron spacing
  real*8 :: a_i            ! ion sphere radius
  real*8 :: physGamma      ! coupling parameter
  real*8 :: V0_eV = 0.       ! desired potential at distance r=0 from an ion --> eps is adjusted to match this value
  real*8 :: eps = 1.           ! potential/force law cutoff
  integer :: Zion=1, Aion=1       ! ion charge and mass number
  integer :: setup_type = 0 !< for computing volume, interparticle distance, etc: 0-cubic, 1-spherical
  integer :: momentum_acf_from_timestep = 0

  real*8 :: Ukine          ! Electron kinetic energy
  real*8 :: Ukini          ! Ion kinetic energy

!  Variables needing 'copy' for tree routines
  integer :: npart_total  ! Total # particles (npart)
  integer :: np_local 

!  Associated MPI stuff

  integer :: my_rank       ! Rank of current task
  integer :: n_cpu   ! # cpus used by program
  integer :: MPI_COMM_PEPC ! MPI communicator to be used (will usually be a copy of MPI_COMM_WORLD)

! Control stuff
  integer :: idim=3  ! # dimensions (velocity and position updates)
  integer :: ispecial = 1      ! Switch to select special electron configs
  integer :: debug_level = 0 ! Debug level for printed O/P
  logical :: treediags = .false.

   logical, public :: restart = .false. !< Restart switch: config read from parts_all.in
   real*8 :: dt = 0.01            ! timestep
   real*8 :: trun = 0.          ! total run time including restarts
   integer :: nt = 100
   integer :: itime = 0   ! # timesteps and current timestep
   integer :: itime_in    ! timestep to read mpi-io checkpoint from in case of ispecial==-1
   integer :: idump = 0, idump_vtk = 0, idump_checkpoint = 0, idump_binary = 0 ! output frequency (timesteps): ascii, vtk and mpi-io-checkpoint

   integer :: ifile_cpu    ! O/P stream

   logical :: directforce = .false. !< if set to true, the frontend only performs a direct force computation instead of utilizing the treecode

end module physvars





