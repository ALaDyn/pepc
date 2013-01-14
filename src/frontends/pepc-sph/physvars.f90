! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2013 Juelich Supercomputing Centre, 
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
  use module_pepc_types, only: &
       t_particle

  real, parameter :: pi=3.141592654

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

  real*8 :: thermal_constant = 1.           !< TODO: is this the boltzmann constant?
  real*8 :: kappa            = 1.4          !< kappa = gamma: Heat capacity ratio for gas
  real*8 :: art_vis_alpha    = 1.           !< artificial viscosity parameter
  real*8 :: art_vis_beta     = 2.           !< artificial viscosity parameter


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

   logical :: initialized_v_minus_half = .false.  ! set this to .true. then v_minus_half is initialized in particle_setup or when data are read. otherwhise v_minus_half is computed from v befor entering the loop over timesteps

   ! I/O stuff
   integer :: dump_time, cp_time   ! When to dump, when to do a checkpoint (read-in)
   integer :: input_itime          ! Which step shpuld we read in?
   character(50) :: mpifile        ! MPI-IO-file
   












end module physvars





