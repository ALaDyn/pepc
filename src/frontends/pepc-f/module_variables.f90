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

!!!!!!!!!!!!!!!!!!!!
!! helper module
!!!!!!!!!!!!!!!!!!!!

module variables
  use module_pepc_types
  implicit none

  ! output variables
  integer,parameter :: err_output_id=666

  ! MPI variables
  integer :: my_rank, n_ranks
  logical :: root

  ! time variables
  real*8 :: dt
  integer :: step

  ! control variables
  integer :: nt               ! number of timesteps
  integer :: tnp              ! total number of particles (wall +plasma)
  integer :: np               ! local number of particles (wall +plasma)

  ! type of source
  integer :: quelltyp

  
  real*8  :: delx      ! length of plasma (multiples of lambda_debye)
  real*8  :: dely       ! width of plasma (multiples of lamor radius)
  real*8  :: delz       ! width of plasma (multiples of lamor radius)
  real*8  :: dx,dy,dz           ! lenght/width in m
  
  ! physics constants
  real*8, parameter  ::  e=1.602176565e-19
  real*8, parameter  ::  me=9.1093829e-31
  real*8, parameter  ::  mp=1.67262177e-27
  real*8, parameter  ::  kb=1.3806488e-23
  real*8, parameter  ::  md=3.34358320e-2
  real*8, parameter  ::  eps0=8.854187e-12
  real*8             ::  pi
  real*8             ::  fc                  !force constant: SI: 1/4/Pi/eps0

  ! Parameters
  real*8 :: ti_ev  
  real*8 :: te_ev 
  real*8 :: Bx               !für räumlich
  real*8 :: By               !konstantes
  real*8 :: Bz               !Magnetfeld
  real*8 :: B
  real*8 :: ni  
  real*8 :: ne 
  real*8 :: l_debye, omega_p, r_lamor
  real*8 :: fsup             !superparticle factor

  ! particle data (position, velocity, mass, charge)
  type(t_particle), allocatable :: particles(:)
  integer                       :: next_label

  !wall particles
  type(t_particle), allocatable :: wall_particles(:)
  integer                       :: tnwpy                          !total number of wall particles in y dir
  integer                       :: tnwpz                          !total number of wall particles in z dir
  integer                       :: nwp                            !local number of wall particles
  integer                       :: tnwp
  real*8,allocatable            :: wall_pos(:,:)

  !plasma particles
  type(t_particle), allocatable :: plasma_particles(:)
  integer :: tnpp              ! total number of plasma particles 
  integer :: npp               ! local number of plasma particles 

end module