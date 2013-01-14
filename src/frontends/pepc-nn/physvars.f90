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
  use module_pepc_types

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
  integer :: debug_level =0 ! Debug level for printed O/P

   real :: dt             ! timestep
   real :: trun           ! total run time including restarts
   integer :: nt, itime   ! # timesteps and current timestep

   integer :: ifile_cpu    ! O/P stream


end module physvars





