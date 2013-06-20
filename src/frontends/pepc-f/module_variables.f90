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

!!!!!!!!!!!!!!!!!!!!
!! variables module
!!!!!!!!!!!!!!!!!!!!

module variables
  use module_pepc_types
  use module_geometry_types
  use module_species_types
  implicit none

  integer,parameter :: err_output_id=666

  ! MPI variables
  integer :: my_rank, n_ranks
  logical :: root
  integer :: ierr

  ! filehandles
  integer :: out=199
  integer :: timing_out=200
  integer :: recycling_out=201

  ! time variables
  real*8 :: dt
  integer :: step,startstep

  ! control variables
  integer :: nt               ! number of timesteps
  integer :: tnp              ! total number of particles (all species)
  integer*8 :: npart          ! total number of particles, needed as int8 for checkpoints
  integer :: np               ! local number of particles (all species)
  integer :: diag_interval=0  ! interval for writing probe data
  integer :: vtk_interval=0   ! interval for writing vtk output
  integer :: checkp_interval=0! interval for setting checkpoints
  logical :: diags
  logical :: interaction_partner_diags

  ! type of rng (0=standard fortran,1=par_rand from module_zufall)
  integer :: rng
  ! treat electrons in guiding centre approximation
  logical :: guiding_centre_electrons


  real*8  :: dx,dy,dz           ! lenght/width in m
  real*8  :: xmin,ymin,zmin     ! used if the domain does not start at (0,0,0)
  real*8  :: xmax,ymax,zmax     ! 
  
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
  type(t_particle), allocatable :: mirror_particles(:)
  type(t_particle), allocatable :: all_particles(:)
  integer(kind_particle)        :: next_label

  !variables for reflux in every 2nd timestep
  integer                        :: new_e_r_last_ts=0
  integer                        :: new_i_r_last_ts=0
  integer                        :: last_reflux_step=0
  logical                        :: need_to_reflux=.false.


  integer,allocatable  :: tnpps(:)  !total number of particles per species
  integer,allocatable  :: npps(:)   !local number of particles per species


  !aux strings
  character(255) :: filename
  !aux ints
  integer :: cmd_args

  !other
  integer :: chunk_size_default

  !test
  type(t_boundary), allocatable :: boundaries(:)
  integer :: nb          ! number of boundaries
  type(t_species), allocatable :: species(:)
  integer :: nspecies    ! number of species

  integer :: spiegelladung=0

  real(KIND=8),allocatable :: probe_start_x(:), probe_start_y(:), probe_start_z(:)
  real(KIND=8),allocatable :: probe_end_x(:), probe_end_y(:), probe_end_z(:)

  namelist /probe_positions/ probe_start_x, probe_start_y, probe_start_z,probe_end_x, probe_end_y, probe_end_z
  namelist /pepcf/ fsup,guiding_centre_electrons, nt, dt, Bx, By, Bz, dx ,dy, dz,diag_interval, checkp_interval, vtk_interval,spiegelladung
  namelist /walk_para_smpss/ chunk_size_default


end module
