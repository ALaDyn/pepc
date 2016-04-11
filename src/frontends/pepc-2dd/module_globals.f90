! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
!
! Copyright (C) 2002-2016 Juelich Supercomputing Centre,
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

module module_globals
  use module_pepc_kinds
  use module_pepc_types
  implicit none
  save

  !integer(kind_dim), parameter :: dim = 2

  ! MPI variables
  integer :: my_rank, n_ranks,communicator
  logical :: root

  ! time variables
  real(kind_particle)    :: dt
  integer                :: step

  ! control variables
  integer :: nt               ! number of timesteps
  integer(kind_particle) :: tnp ! total number of particles
  integer(kind_particle) :: np  ! local number of particles
  real(kind_particle)    :: eps2!,theta2
  integer(kind_dim)      :: ixdim
  integer(kind_dim)      :: ivdim
  integer(kind_particle) :: initial_setup,normal
  logical                :: particle_output  ! turn vtk output on/off
  logical                :: domain_output    ! turn vtk outpuit on/off
  logical                :: particle_test    ! turn direct summation on/off
  logical                :: periodicity_particles    ! particle periodicity on/off
  logical                :: flag_classic             ! classic/relativistic time integrator 
  integer                :: diag_interval    ! number of timesteps between all diagnostics and IO
  integer                :: restart_step     ! number of timesteps between all diagnostics and IO
  integer                :: nsp              ! number of species
  integer(kind_particle) :: adv              ! numeical scheme
  real(kind_particle)    :: Volume           ! Volume
  real(kind_particle)    :: wpe              ! Plasma density
  real(kind_particle)    :: we               ! Ratio of the number of real particles and simulated particles 
  real(kind_particle)    :: vmax             ! Max velocity - for velocity normalized by vth

  real(kind_particle)    :: norm_factor
            !!!!!  CGS System  !!!!!
  real(kind_particle), parameter  :: c  = 2.998e+10     ! speed of light
  real(kind_particle), parameter  :: qe = 4.803e-10     ! proton charge
  real(kind_particle), parameter  :: me = 9.109e-28     ! electron mass
  real(kind_particle), parameter  :: mi = 1.673e-24     ! proton mass
  real(kind_particle), parameter  :: G  = 6.673e-8      ! Gravitational constant
  real(kind_particle), parameter  :: k  = 1.381e-16     ! Boltzmann constant
  
            !!!!!  Conversion Parameters From CGS To SI Check NRL Plasma Formulary  !!!!!
  real(kind_particle), parameter  :: alpha   = 1.0e+2                           ! [cm/m]
  real(kind_particle), parameter  :: beta    = 1.0e+7                           ! [erg/J]
  real(kind_particle), parameter  :: epsilon0= 8.8542e-12                       ! [F/m]
  real(kind_particle), parameter  :: mu0     = 8.0_8*acos(0.0_8)*1.0e-7         ! [H/m]
  real(kind_particle), parameter  :: hcut    = 1.0546e-34                       ! [Js]

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer(kind_physics)               :: woutput
  integer(kind_physics), dimension(3) :: n       = [ 0 , 0  , 0  ]
  integer(kind_physics), dimension(3) :: nppd    = [ 0 , 0  , 0  ] ! number of particle per direction
  real(kind_physics)   , dimension(3) :: offset  = [ 0., 0. , 0. ]
  real(kind_physics)   , dimension(3) :: extent  = [ 0., 0. , 0. ]
  real(kind_physics)   , dimension(3) :: veth    = [ 0., 0. , 0. ]
  real(kind_physics)   , dimension(3) :: vedrift = [ 0., 0. , 0. ]
  real(kind_physics)   , dimension(3) :: vith    = [ 0., 0. , 0. ]
  real(kind_physics)   , dimension(3) :: vidrift = [ 0., 0. , 0. ]

  real(kind_physics)                  :: xtilde,vtilde,ttilde,qtilde,mtilde,&
                                         etilde,btilde,jtilde,rhotilde,atilde,&
                                         phitilde,lorentz_tilde

  character(*), parameter :: FRONTEND_NAME = "pepc-2dd"
  character(255)          :: ischeme,restart_file,folder


end module
