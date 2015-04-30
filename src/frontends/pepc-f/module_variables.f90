! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
!
! Copyright (C) 2002-2014 Juelich Supercomputing Centre,
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
  use module_pepc_kinds
  use module_pepc_types
  use module_geometry_types
  use module_species_types
  implicit none

  integer,parameter :: err_output_id=666

  ! physics constants
  real*8, parameter  ::  e=1.602176565e-19
  real*8, parameter  ::  me=9.1093829e-31
  real*8, parameter  ::  mp=1.67262177e-27
  real*8, parameter  ::  kb=1.3806488e-23
  real*8, parameter  ::  md=3.34358320e-2
  real*8, parameter  ::  eps0=8.854187e-12
  real*8             ::  pi
  real*8             ::  fc                  !force constant: SI: 1/4/Pi/eps0


  ! MPI variables
  integer :: my_rank, n_ranks
  logical :: root
  integer :: ierr


  ! filehandles
  integer :: out=199
  integer :: detailed_timing_out=200
  integer :: ph_timing_out=201
  integer :: recycling_out=202


  ! time variables
  real*8 :: dt
  integer :: step,startstep
  integer :: nt               ! number of timesteps


  ! output intervals
  integer :: diag_interval    ! interval for writing probe data
  integer :: vtk_interval     ! interval for writing vtk output
  integer :: checkp_interval  ! interval for setting checkpoints
  integer :: npy_interval     ! interval for dumping particles in npy format
  logical :: checkpoint_now
  logical :: vtk_now
  logical :: diag_now
  logical :: npy_now


  ! Physcial System
  real(KIND=8) :: dx,dy,dz                         ! lenght/width in m
  real(KIND=8) :: xmin,ymin,zmin
  real(KIND=8) :: xmax,ymax,zmax
  real(KIND=8) :: Bx                               ! für räumlich
  real(KIND=8) :: By                               ! konstantes
  real(KIND=8) :: Bz                               ! Magnetfeld
  real(KIND=8) :: B                                ! norm((Bx,By,Bz))
  real(KIND=8) :: fsup                             ! superparticle factor
  integer :: nb                                    ! number of boundaries
  type(t_boundary), allocatable :: boundaries(:)
  integer :: nspecies                              ! number of species
  type(t_species), allocatable :: species(:)


  ! variables for output of averaged physical quantities
  integer :: diag_bins_x
  integer :: diag_bins_y
  integer :: diag_bins_z
  logical :: bool_diag_bins_cylinder
  logical :: bool_avg_btwn_diag_steps
  real(KIND=8), allocatable :: data_bins(:,:,:,:,:) !(nspecies,38,diag_bins_x,diag_bins_y,diag_bins_z)
  integer, allocatable      :: n_bins(:,:,:,:)      !(nspecies,diag_bins_x,diag_bins_y,diag_bins_z)

  ! variables for output of averaged physical quantities in velocity space
  integer :: diag_bins_vx
  integer :: diag_bins_v2
  real(KIND=8) :: v_grid_max
  real(KIND=8), allocatable :: data_bins_v(:,:,:,:) !(nspecies,38,diag_bins_vx,diag_bins_v2)
  integer, allocatable      :: n_bins_v(:,:,:)      !(nspecies,diag_bins_vx,diag_bins_v2)


  ! variables for detailed boundary hit statistics
  logical :: bool_energy_resolved_hits
  logical :: bool_angle_resolved_hits
  logical :: bool_space_resolved_hits
  logical :: bool_age_resolved_hits
  integer :: last_diag_output
  integer :: nbins_energy_resolved_hits
  real(KIND=8) :: ehit_max_in_T, agehit_max_in_t_trav_ion
  real(KIND=8), allocatable :: ehit_max(:)
  real(KIND=8), allocatable :: agehit_max(:)
  integer, allocatable :: energy_resolved_hits(:,:,:)
  integer :: nbins_angle_resolved_hits
  integer, allocatable :: angle_resolved_hits(:,:,:)
  integer :: nbins_e1_space_resolved_hits
  integer :: nbins_e2_space_resolved_hits
  integer, allocatable :: space_resolved_hits(:,:,:,:)
  integer :: nbins_age_resolved_hits
  integer, allocatable :: age_resolved_hits(:,:,:)


  ! variables for delayed refluxing (Benjamin apparently did this every second step)
  integer :: last_reflux_step
  integer :: reflux_interval


  ! particle arrays
  type(t_particle), allocatable :: particles(:)                                    !particles
  type(t_particle), allocatable :: all_particles(:)                                !all particles if spiegelladung is set (includes mirror_particles)
  integer,allocatable  :: tnpps(:)                                                 !total number of particles per species
  integer,allocatable  :: npps(:)                                                  !local number of particles per species
  integer(kind_particle)        :: next_label
  real(KIND=8),allocatable :: probe_start_x(:), probe_start_y(:), probe_start_z(:) !positions of probe particles
  real(KIND=8),allocatable :: probe_end_x(:), probe_end_y(:), probe_end_z(:)       !positions of probe particles


  ! aux variable
  real(KIND=8),allocatable :: maxw_flux_table_F(:,:) ! needed for drifting maxwellian flux output
  real(KIND=8),allocatable :: maxw_flux_table_v(:,:) ! needed for drifting maxwellian flux output
  character(255) :: filename
  integer :: cmd_args                                ! number of command line args


  ! control variables (some are only temporary/for testing)
  integer :: spiegelladung = 0  !temp
  integer :: retherm = 0        !temp
  logical :: diags !temp
  logical :: guiding_centre_electrons  ! treat electrons in guiding centre approximation (scheme by Benjamin, not sure if correct)
  integer :: rng   !type of rng (0=standard fortran,1=par_rand from module_zufall)
  integer :: rngseed
  logical :: bool_hockney_diag


  ! other
  integer :: chunk_size_default


  ! namelists
  namelist /probe_positions/ probe_start_x, probe_start_y, probe_start_z,probe_end_x, probe_end_y, probe_end_z
  namelist /pepcf/ fsup,guiding_centre_electrons, nt, dt, Bx, By, Bz, xmin, rngseed ,&
                   xmax, ymin, ymax, zmin, zmax, diag_interval, checkp_interval,npy_interval,&
                   vtk_interval,spiegelladung, diag_bins_x,diag_bins_y,diag_bins_z,retherm,&
                   bool_angle_resolved_hits, bool_energy_resolved_hits, bool_space_resolved_hits, &
                   bool_age_resolved_hits, nbins_age_resolved_hits, &
                   nbins_angle_resolved_hits, nbins_energy_resolved_hits, nbins_e1_space_resolved_hits, &
                   nbins_e2_space_resolved_hits, bool_diag_bins_cylinder, bool_avg_btwn_diag_steps, reflux_interval
  namelist /walk_para_smpss/ chunk_size_default

end module
