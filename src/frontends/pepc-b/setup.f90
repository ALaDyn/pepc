! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2017 Juelich Supercomputing Centre, 
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

!  ================================
!
!         SETUP
!
!   $Revision: 1347 $
!
!     Initialise physics constants and 
!      simulation variables
!
!  ================================


subroutine setup
  use module_pepc
  use module_mirror_boxes
  use module_fmm_framework
  use module_physvars
  use mpi

  implicit none
  integer :: ne_rest, ni_rest, ierr

  character(len=255) :: parameterfile
  logical :: read_param_file

  integer*4 IARGC

!  include 'namelist.h' 
  !  Default input set

  ! switches

  plasma_config = 1  ! plasma target
  target_geometry = 1         ! random sphere


  ! particles
  nep = 0 ! # plasma electrons per PE
  nip = 0
  ne = 100  ! Total # plasma electrons
  ni = 100  ! total # plasma ions
  mc_steps = 10

  xl = 2
  yl = 2
  zl = 2
  
  ! physics stuff
  force_const = 1./3.
  bond_const = 0.1
  rho0 = 1.0
  mac = 0        ! Multipole acceptance criterion (BH by default)
  theta = 0.5
  Te_keV = 1.
  Ti_keV = Te_keV/10.
  mass_ratio = 10.
  Zion = 1.
  uthresh = -1.

  r_sphere = 0.5
  x_plasma = 0.1    ! plasma disc thickness (2) or wire length (3)
  y_plasma = 1.     ! plasma width (slab target)
  z_plasma = 1.     ! plasma width (slab target)
  eps = 0.1
  fnn = 5  ! Neighbour search radius multiplier (x a_ii)
  delta_mc = r_sphere/5.
  displace(1:3) = (/0.,0.,0./)
  ! beam

  !  beam_config = 1  ! fixed beam, initialised at start
  ! beam_config = 2  ! user-controlled, real-time particle source
  beam_config_in = 0 ! beam off
  ! beam_config = 4  ! laser ponderomotive force

  r_beam = 0.15
  u_beam = 0.2
  theta_beam = 0.
  phi_beam = 0.3
  x_beam = .04
  start_beam = 0.4
  rho_beam = -5.
  mass_beam = 3.
  np_beam = 0 ! initial # beam particles/ dt

  ! laser stuff
  sigma = 1.
  tpulse = 10.
  vosc = 0.1
  omega = 0.5
  theta_inc = 0.
  x_offset = 0.
  z_offset = 0.
  tlaser = 0.    ! time since laser switched on
  lambda = 1.0
  rho_track = 1.5
  lolam = 0.1
  rho_min = 0.1

  ! control
  nt = 600
  dt = 0.2
  trun = 0.
  ivis = 1
  ivis_fields = 1
  ivis_domains = 1
  itime_start = 0
  itrack = 10
  domain_cut = zl

  ngx = 25   ! Grid size for plots
  ngy = 25
  ngz = 25
  ! constrain
!!$    len_tripod = .001
  constrain_proof = .001
  struct_step = 0

! Vectors for periodic bcs
  t_lattice_1(1:3)=(/1.,0.,0./)
  t_lattice_2(1:3)=(/0.,1.,0./)
  t_lattice_3(1:3)=(/0.,0.,1./)
  periodicity(1:3) = (/.false., .false., .false./)  !< boolean switches for determining periodicity directions
  fmm_extrinsic_correction = FMM_EXTRINSIC_CORRECTION_REDLACK

 ! read in first command line argument
  call pepc_get_para_file(read_param_file, parameterfile, my_rank)

  if (read_param_file) then
     if(my_rank .eq. 0) write(*,*) "reading parameter file, section pepcb: ", parameterfile
     open(10,file=parameterfile)
     read(10,NML=pepcb)
     close(10)
  else
     if(my_rank .eq. 0) write(*,*) "##### using default parameter #####"
  end if

 ! Scale periodic lattice vectors by box lengths
  t_lattice_1 = t_lattice_1*xl
  t_lattice_2 = t_lattice_2*yl
  t_lattice_3 = t_lattice_3*zl

 ! # particles in primary target component
  npart_total = ni+ne

! Adjust local numbers if total non-multiple of # PEs
! Nov 29 2009 PG: Changed from previous assignment, which placed all
! remaining particles on root cpu
! Now spread evenly over 1st ne_rest cpus (each gets one extra)

  nep = ne/n_cpu  ! local # electrons and ions - may be adjusted later
  nip = ni/n_cpu 
  ne_rest = mod(ne,n_cpu)
  ni_rest = mod(ni,n_cpu)

! Assign remainder of each particle species 
  if (my_rank.lt.ne_rest) nep=nep+1  
  if (my_rank.lt.ni_rest) nip=nip+1 

  np_local = nep+nip ! initial estimate for total local # particles
  new_label = npart_total  ! Rezone label


  !if (target_dup .or. scheme==5) np_mult = np_mult*2  ! double up particle array size if multi-target or ion config mode 

  if (.not. restart .and. plasma_config .ge. 10) npart_total = npart_total + 2*SUM(n_layer) ! Include particles from remaining layers for array setup.

  np_alloc = 2*np_local  ! Allow for load balancing in memory allocation

  beam_config=mod(beam_config_in,10)  ! derived config from s,p variations
  lolam = lolam*2.*pi/omega  ! normalise scale-length

end subroutine setup




