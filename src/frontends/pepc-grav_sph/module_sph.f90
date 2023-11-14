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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> A module to apply the sph-method to pepc particles
!>

!> \author Andreas Breslau
!> \date 2011.12
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module module_sph

   use module_pepc_kinds
   use module_interaction_specific_types, only: &
      num_neighbour_particles

   implicit none
   private

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!  public type declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer :: idim = 3

   real*8 :: thermal_constant
   real*8 :: kappa
   real*8 :: art_vis_alpha
   real*8 :: art_vis_beta
   logical :: use_artificial_viscosity = .true.

   real*8 :: pi = acos(-1.0)

   logical :: sph_debug = .false.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!  subroutine-implementation  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   public sph
   public sph_kernel_tests
!  public sph_initialize
!  public sph_kernel
!  public sph_density
!  public sph_grad_kernel
!  public sph_sum_force
!  public update_particle_props

contains

   subroutine sph_initialize(idim_tmp, thermal_constant_tmp, kappa_tmp, art_vis_alpha_tmp, art_vis_beta_tmp)

      use mpi
      implicit none

      integer, intent(in) :: idim_tmp
      real*8, intent(in) :: thermal_constant_tmp
      real*8, intent(in) :: kappa_tmp
      real*8, intent(in) :: art_vis_alpha_tmp
      real*8, intent(in) :: art_vis_beta_tmp

      integer :: ierr

      ! set module variables
      idim = idim_tmp
      thermal_constant = thermal_constant_tmp
      kappa = kappa_tmp
      art_vis_alpha = art_vis_alpha_tmp
      art_vis_beta = art_vis_beta_tmp

      use_artificial_viscosity = .true.

   end subroutine sph_initialize

   ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !> \brief output the values of the sph kernel and grad_kernel for 100 steps from -2h to 2h into a file
   !>

   !> \author Andreas Breslau
   !> \date 2011.12

   !   param[in,out]  Name      Description
   !> \param[in]      idim_     the dimension, changes the normalisation factor in the kernel
   subroutine sph_kernel_tests(idim_)

      implicit none

      integer, intent(in) :: idim_

      real*8 :: r
      real*8 :: h
      real*8 :: kernel
      real*8 :: grad_kernel
      integer :: i, steps

      idim = idim_

      steps = 100
      h = 1.

      open (77, FILE='sph_kernel_test.dat', STATUS='NEW')

      do i = 0, steps

         r = -2.+4./steps * i
         call sph_kernel(abs(r), h, kernel)
         call sph_grad_kernel(abs(r), h, grad_kernel)

         grad_kernel = grad_kernel * r

         write (77, *) i, r, kernel, grad_kernel

      end do

      close (77)

   end subroutine sph_kernel_tests

   ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !> \brief compute the kernel for sph
   !>
   !>

   !> \author Andreas Breslau
   !> \date 2011.12.01

   !   param[in,out]  Name      Description
   !> \param[in]      distance
   !> \param[in]      h
   !> \param[out]     kernel

   subroutine sph_kernel(r, h, kernel)

      use mpi
      implicit none

      real*8, intent(in) :: r                                 !< scalar distance between two particles, should allways be >= 0
      real*8, intent(in) :: h                                 !< sph smoothing length h
      real*8, intent(out) :: kernel                           !< scalar sph kernel, W

      real*8 :: q                                             !< r/h
      real*8 :: sph_kernel_factor

      integer :: ierr

      q = r / h

      if (idim .eq. 3) then
         sph_kernel_factor = 1._8 / pi / h**3                                     ! see monaghan 1992, S. 554, 555
      else if (idim .eq. 2) then
         sph_kernel_factor = 10._8 / 7._8 / pi / h**2
      else if (idim .eq. 1) then
         sph_kernel_factor = 2._8 / 3._8 / h
      else
         write (*, *) "idim not one of 1, 2, 3. Terminating:", idim

         call MPI_ABORT(MPI_COMM_WORLD, 666, ierr)
      end if

      ! see monaghan 1992, S. 554, 555
      if (q .ge. 0._8 .and. q .le. 1._8) then
         kernel = sph_kernel_factor * (1._8 - 3._8 / 2._8 * (r / h)**2 + 3._8 / 4._8 * (r / h)**3)
      else if (q .le. 2._8) then
         kernel = sph_kernel_factor * (1._8 / 4._8 * (2._8 - (r / h))**3)
      else ! q > 2, or negative
         write (*, *) "SPH kernel: q not in [0,2]. Should never happen:", q

         call MPI_ABORT(MPI_COMM_WORLD, 666, ierr)
      end if

   end subroutine sph_kernel

   ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !> \brief compute sph density
   !>
   !>

   !> \author Andreas Breslau
   !> \date 2011.12.01

   !   param[in,out]  Name      Description
   !> \param[in]
   !> \param[in]
   !> \param[in]
   ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine sph_density(np_local, particles, itime, num_neighbour_boxes, neighbour_boxes)

      use module_pepc, only: &
         t => global_tree

      use module_pepc_types, only: &
         t_particle, &
         t_tree_node

      use mpi
      implicit none

      integer, intent(in) :: np_local    !< # particles on this CPU
      type(t_particle), intent(inout) :: particles(:)
      integer, intent(in) :: itime  ! timestep
      integer, intent(in) :: num_neighbour_boxes !< number of shift vectors in neighbours list (must be at least 1 since [0, 0, 0] has to be inside the list)
      integer, intent(in) :: neighbour_boxes(3, num_neighbour_boxes) ! list with shift vectors to neighbour boxes that shall be included in interaction calculation, at least [0, 0, 0] should be inside this list

      integer :: ierr

      integer :: local_particle_index
      integer :: actual_neighbour
      real*8 :: h
      real*8 :: kernel
      type(t_tree_node), pointer :: actual_node

!     IF( me == 0) then

!        write(87, *) me, idim

!        DO actual_particle = 1, num_part
!           WRITE(87, *) part_indices( actual_particle ), r_nn( part_indices( actual_particle ) ), pelabel( part_indices( actual_particle ) )
!           DO actual_neighbour = 1, n_nn
!              WRITE(87, *) charge( nodelist( actual_neighbour, actual_particle ) ), dist2_list( actual_neighbour, actual_particle )
!           END DO
!        END DO

!     end if

      ! compute smoothing length (h) for all local particles
      ! TODO: move this out of sph_density
      ! TODO: parallelize with OpenMP

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(local_particle_index)
      do local_particle_index = 1, np_local
         particles(local_particle_index)%results%h = sqrt(particles(local_particle_index)%results%maxdist2) / 2._8
      end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(local_particle_index, h, kernel, actual_neighbour, actual_node)
      do local_particle_index = 1, np_local

         h = particles(local_particle_index)%results%h

         call sph_kernel(0._8, h, kernel) ! particle self
         particles(local_particle_index)%results%rho = particles(local_particle_index)%data%q * kernel

         do actual_neighbour = 1, num_neighbour_particles
            actual_node => t%nodes(particles(local_particle_index)%results%neighbour_nodes(actual_neighbour))

            call sph_kernel(sqrt(particles(local_particle_index)%results%dist2(actual_neighbour)), h, kernel)

            particles(local_particle_index)%results%rho = particles(local_particle_index)%results%rho + &
                                                          actual_node%interaction_data%charge * kernel

         end do

      end do
!$OMP END PARALLEL DO

   end subroutine sph_density

   ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !> \brief compute nabla_kernel for sph
   !>
   !>

   !> \author Andreas Breslau
   !> \date 2011.12.01

   !   param[in,out]  Name      Description
   !> \param[in]      r
   !> \param[in]      h
   !> \param[out]     grad_kernel

   subroutine sph_grad_kernel(r, h, grad_kernel)

      use mpi
      implicit none

      ! \bug ab: using idim for dimensions, be carefull with idim .ne. 3 because of pepc
      real*8, intent(in) :: r                                 !< scalar distance between two particles, should be >= 0
      real*8, intent(in) :: h                                 !< sph smoothing length h
      real*8, intent(out) :: grad_kernel                      !< scalar part of gradient of sph kernel, Nabla W,
      ! has to be multiplied by the distance vector by the calling function

      real*8 :: q                                             !< r/h
      real*8 :: sph_kernel_factor

      integer :: ierr

      q = r / h

      if (idim .eq. 3) then
         sph_kernel_factor = 1._8 / pi / h**3                                     ! see monaghan 1992, S. 554, 555
      else if (idim .eq. 2) then
         sph_kernel_factor = 10._8 / 7._8 / pi / h**2
      else if (idim .eq. 1) then
         sph_kernel_factor = 2._8 / 3._8 / h
      else
         write (*, *) "idim not one of 1, 2, 3. Terminating."

         call MPI_ABORT(MPI_COMM_WORLD, 666, ierr)
      end if

      if (q .ge. 0._8 .and. q .le. 1._8) then
         grad_kernel = sph_kernel_factor * (9._8 * r - 12._8 * h) / (4._8 * h**3)
      else if (q .gt. 1._8 .and. q .le. 2._8) then
         grad_kernel = sph_kernel_factor * (-3._8) * (2._8 - r / h)**2 / (4._8 * h * r)
      else ! q > 2, or negative
         write (*, *) "SPH grad kernel: q not in [0,2]. Should never happen:", q

         ! ab: contribution of last particle is always 0, because h is defined as half of the distance to this particle, so we are really using n-1 neighbours

         call MPI_ABORT(MPI_COMM_WORLD, 666, ierr)
      end if

   end subroutine sph_grad_kernel

   ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Subroutine to sum the sph acceleration, not the force. The tree_walk returns the field,
   ! which, in case of gravity, equals the acceleration.
   !
   ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine sph_sum_force(np_local, particles, itime, num_neighbour_boxes, neighbour_boxes)

      use module_pepc, only: &
         t => global_tree

      use module_pepc_types, only: &
         t_particle, &
         t_tree_node

      use physvars, only: &
         my_rank

      use mpi
      implicit none

      integer, intent(in) :: np_local    !< # particles on this CPU
      type(t_particle), intent(inout) :: particles(:)
      integer, intent(in) :: itime  ! timestep
      integer, intent(in) :: num_neighbour_boxes !< number of shift vectors in neighbours list (must be at least 1 since [0, 0, 0] has to be inside the list)
      integer, intent(in) :: neighbour_boxes(3, num_neighbour_boxes) ! list with shift vectors to neighbour boxes that shall be included in interaction calculation, at least [0, 0, 0] should be inside this list

      integer :: ierr

      integer :: local_particle_index
      integer :: actual_neighbour
      real*8 :: h1, h2                          !< smoothing length of particle and current interaction partner
      real*8 :: grad_kernel_1, grad_kernel_2    !< kernel gradient for particle and current interaction partner with h1 and h2
      type(t_tree_node), pointer :: actual_node

      real*8 :: const
      real*8, dimension(3) :: dist              !< distance vector from particle to current interaction partner
      real*8 :: distance                        !< scalar distance from particle to current interaction partner
      real*8 :: artificial_viscosity
      real*8 :: vr, rr, mu, eta
      real*8 :: sound_speed
      real*8 :: scalar_force
      real*8 :: thermal_energy_sum
      real*8 :: thermal_energy_factor
      real*8, dimension(3) :: dv                !< velocity difference between particle and current interaction partner
      real*8 :: energy_factor
      integer :: dim

      CHARACTER(100) :: forcefile

      do local_particle_index = 1, np_local
         particles(local_particle_index)%results%sph_force = [0._8, 0._8, 0._8]
      end do

      const = thermal_constant  ! m**2 s**-2 K
      ! this is the specific gas constant R_s = k_B / m_M (with m_M the mass of a molecule)

      ! p = rho *const*T , quick n dirty

      energy_factor = 1._8 * (kappa - 1._8) / const / 2._8

      if (sph_debug) then
         write (forcefile, '(a,i6.6,a,i6.6,a)') "sph_debug_", itime, "_", my_rank, ".list"
         open (50, FILE=forcefile, STATUS='NEW')
      end if

      ! thermal energy change:
      ! du_i/dt = 1/2 SUM (m_j ( P_b/rho_b^2 + P_a/rho_a^2 ) (v_i -v_j ) grad_kernel ) (Monaghan 1992 S. 549)
      ! dT/dt = (gamma -1) /k_B du/dt
      ! gamma = kappa (adiabatic exponent)
      ! \todo ab: correction factor f_i (siehe Springel 2010, SEREN 2011) ??

      ! acceleration:
      ! dv_i/dt = -SUM (m_j ( P_b/rho_b^2 + P_a/rho_a^2 ) grad_kernel ) (Monaghan 1992 S. 546)

      ! for artificial viscosity (art_vis) see Monaghan 1992, S. 550

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(local_particle_index, h1, h2, thermal_energy_sum, actual_neighbour, actual_node, grad_kernel_1, grad_kernel_2, &
!$OMP& rr, vr, dim, dist, dv, thermal_energy_factor, eta, sound_speed, mu, artificial_viscosity, scalar_force, distance)
      do local_particle_index = 1, np_local

         h1 = particles(local_particle_index)%results%h

         thermal_energy_sum = 0._8

         do actual_neighbour = 1, num_neighbour_particles

            actual_node => t%nodes(particles(local_particle_index)%results%neighbour_nodes(actual_neighbour))

            ! scalar distance
            distance = sqrt(particles(local_particle_index)%results%dist2(actual_neighbour))

            h2 = actual_node%interaction_data%h

            call sph_grad_kernel(distance, h1, grad_kernel_1)
            !          call sph_grad_kernel( distance , h2, grad_kernel_2 )

            ! PEPCs sum_force returns the electric field components ex, ey, ez. In velocities the accelerations are then computed as
            ! a = charge / mass * e(x|y|z). For gravity charge/mass = 1, so the fields and accelerations are equal.
            ! Because of this the accelerations caused by sph_force can be added to the fields.

            ! distance vector
            ! pointing from actual_particle to actual_neighbour
            ! the dist_vector stored by the walk points form actual_neighbour to actual_particle, therefore the "-"
            dist = -particles(local_particle_index)%results%dist_vector(:, actual_neighbour)

            if (sph_debug) then
               if (abs(sqrt(dist(1)**2 + dist(2)**2 + dist(3)**2) - distance) .gt. 1E-7) then
                  write (50, *) '|dist| - distance > 1E-7:', my_rank, local_particle_index, actual_neighbour, distance, dist
               end if
            end if

            ! rr: art_vis, distance squared
            rr = particles(local_particle_index)%results%dist2(actual_neighbour)

            vr = 0._8

            ! for all dimension
            do dim = 1, idim
               dv(dim) = actual_node%interaction_data%v(dim) - particles(local_particle_index)%data%v(dim)

               ! vr: art_vis, scalar product of velocity difference and distance
               vr = vr + dist(dim) * dv(dim)
            end do

            thermal_energy_factor = vr

            ! TODO: make eta parameter???
            eta = 0.1_8 * h1                                                ! art_vis

            ! TODO: make this a bit faster
            sound_speed = (sqrt(const * particles(local_particle_index)%data%temperature) &
                           + sqrt(const * actual_node%interaction_data%temperature)) / 2. ! mean sound_speed

            if (use_artificial_viscosity .and. (vr .lt. 0._8)) then

               mu = (h1 * vr) / (rr + eta * eta)                        ! art_vis

               artificial_viscosity = (-art_vis_alpha * sound_speed * mu + art_vis_beta * mu * mu) / &
                                      ((actual_node%interaction_data%rho + particles(local_particle_index)%results%rho) / 2.)

            else

               artificial_viscosity = 0._8
            end if

            ! compute scalar part of the force: mass * ( p1/rho1^2 + p2/rho2^2 + art_vis ) * grad_kernel
            scalar_force = actual_node%interaction_data%charge * &
                           ( &
                           const * actual_node%interaction_data%temperature / actual_node%interaction_data%rho + &
                           const * particles(local_particle_index)%data%temperature / particles(local_particle_index)%results%rho + &
                           artificial_viscosity &
                           ) * grad_kernel_1
            ! + grad_kernel_2 ) / 2._8

            ! see diploma thesis andreas breslau, eq. 3.10 (note: there is a missing '-' in the equation before: a = - grad(p)/rho
            ! Here is no "-" because in eq. 3.10 the distance vector (r_a - r_b) is used, here (r_b - r_a)
            particles(local_particle_index)%results%sph_force = particles(local_particle_index)%results%sph_force + scalar_force * dist

            if (sph_debug) then
               write (50, *) my_rank, local_particle_index, particles(local_particle_index)%x(1), h1, &
                  particles(local_particle_index)%data%temperature, particles(local_particle_index)%results%rho, &
                  actual_neighbour, actual_node%interaction_data%coc(1), distance, dist(1), grad_kernel_1, actual_node%interaction_data%charge, &
                  actual_node%interaction_data%temperature, actual_node%interaction_data%h, actual_node%interaction_data%rho, &
                  scalar_force * dist(1), particles(local_particle_index)%results%sph_force(1)
            end if

            ! No "-" is correct. see explanation for force
            thermal_energy_sum = thermal_energy_sum + vr * scalar_force

         end do ! end of loop over neighbours

         particles(local_particle_index)%results%temperature_change = energy_factor * thermal_energy_sum
         ! TODO: add art_visc term to temp_change

         if (sph_debug) then
            write (50, *) particles(local_particle_index)%results%temperature_change
         end if

      end do
!$OMP END PARALLEL DO

      if (sph_debug) then
         close (50)
      end if

   end subroutine sph_sum_force

   subroutine sph(np_local, particles, itime, num_neighbour_boxes, neighbour_boxes, idim)

      use module_pepc_types, only: &
         t_particle

      use physvars, only: &
         thermal_constant

      use mpi
      implicit none

      integer, intent(in) :: np_local    !< # particles on this CPU
      type(t_particle), intent(inout) :: particles(:)
      integer, intent(in) :: itime  ! timestep
      integer, intent(in) :: num_neighbour_boxes !< number of shift vectors in neighbours list (must be at least 1 since [0, 0, 0] has to be inside the list)
      integer, intent(in) :: neighbour_boxes(3, num_neighbour_boxes) ! list with shift vectors to neighbour boxes that shall be included in interaction calculation, at least [0, 0, 0] should be inside this list
      integer, intent(in) :: idim    ! dimension

      integer :: local_particle_index

      call sph_initialize(idim, thermal_constant, 1.4_8, 2._8, 1._8)

      call sph_density(np_local, particles, itime, num_neighbour_boxes, neighbour_boxes)

      call update_particle_props(np_local, particles)

!    call tree_build_upwards(particles(1:np_local)%key, np_local)

      particles(1:np_local)%results%maxidx = 1

      call sph_sum_force(np_local, particles, itime, num_neighbour_boxes, neighbour_boxes)

   end subroutine sph

   subroutine update_particle_props(np_local, particles)

      use module_pepc, only: &
         t => global_tree

      use module_pepc_types, only: &
         t_particle, &
         t_tree_node

      ! only for sort test
      use module_sort

      use physvars, only: &
         my_rank, &
         n_cpu

      use module_tree_node, only: &
         tree_node_is_leaf

      use module_tree, only: &
         tree_traverse_to_key

      use module_debug

      use mpi
      implicit none

      ! Data structure for shipping updated sph properties
      type t_property_update
         sequence
         integer*8 :: key                                                  !< key
         integer   :: owner                                                !< owner
         real*8    :: smoothing_length                                     !< \bug ab: comments needed
         real*8    :: rho                                                  !<
         real*8    :: v(1:3)                                               !< velocity
         REAL*8    :: temperature                                          !< SPH temperature
      end type t_property_update

      integer, parameter :: nprops_property_update = 6

      integer, intent(in) :: np_local    !< # particles on this CPU
      type(t_particle), intent(inout) :: particles(:)

      integer(kind_node) :: nleaf_non_local
      integer*8, allocatable :: non_local_node_keys(:)
      integer*8, allocatable :: non_local_nodes(:)
      integer(kind_pe), allocatable :: non_local_node_owner(:)
      integer, allocatable :: node_arr_cp(:)
      integer, allocatable :: int_arr(:)
      integer :: num_request
      integer :: i, j
      integer :: ierr
      integer, allocatable :: requests_per_process(:)
      integer, allocatable :: requests_from_process(:)
      integer*8, allocatable :: requested_keys(:)
      integer, allocatable :: sdispls(:)
      integer, allocatable :: rdispls(:)
      integer :: total_num_requests_from_others
      integer :: disp
      integer :: actual_address
      integer*8 :: actual_key
      type(t_tree_node), pointer :: actual_node
      integer(kind_node) :: iactual_node

      type(t_property_update), allocatable :: packed_updates(:)
      type(t_property_update), allocatable :: received_updates(:)

      type(t_property_update)  :: dummy_property_update
      integer :: mpi_type_property_update
      integer, parameter :: max_props = nprops_property_update
      ! address calculation
      integer, dimension(1:max_props) :: blocklengths, displacements, types
      integer(KIND=MPI_ADDRESS_KIND), dimension(0:max_props) :: address
      logical :: found

      nleaf_non_local = t%nleaf - t%nleaf_me ! bigger than necessary, TODO: find a better estimation for this

      allocate (non_local_node_keys(nleaf_non_local), non_local_node_owner(nleaf_non_local), requests_per_process(n_cpu), &
                non_local_nodes(nleaf_non_local), node_arr_cp(nleaf_non_local), int_arr(nleaf_non_local), requests_from_process(n_cpu), &
                sdispls(n_cpu), rdispls(n_cpu), STAT=ierr)
      ! TODO: remove int_arr after sort test
      ! TODO: test STAT

      num_request = 0

      ! get leafs from hashtabel with owner .ne. my_rank
      do iactual_node = 1, t%nodes_nentries
         actual_node => t%nodes(iactual_node)
         if ((actual_node%owner .ne. my_rank) .and. tree_node_is_leaf(actual_node)) then
            if (actual_node%owner .gt. n_cpu - 1) write (*, *) 'strange owner:', my_rank, actual_node%owner, actual_node%key

            num_request = num_request + 1
            non_local_nodes(num_request) = iactual_node
            non_local_node_keys(num_request) = actual_node%key
            non_local_node_owner(num_request) = actual_node%owner
         end if
      end do

      ! a test for debugging
      if (sph_debug) then
         do i = 1, num_request
            if (non_local_node_owner(i) .gt. n_cpu - 1) write (*, *) 'owner > n_cpu:', i, non_local_node_owner(i), non_local_node_keys(i)
         end do
      end if

      ! sort keys accorting to owner
      call sort(non_local_node_owner(1:num_request), int_arr(1:num_request))

      ! a test for debugging
      if (sph_debug) then
         do i = 1, num_request
            if (non_local_node_owner(i) .gt. n_cpu - 1) write (*, *) 'after sort, owner > n_cpu:', i, non_local_node_owner(i), non_local_node_keys(i)
         end do
      end if

      ! sort keys according to owners
      non_local_node_keys(1:num_request) = non_local_node_keys(int_arr(1:num_request))
      non_local_nodes(1:num_request) = non_local_nodes(int_arr(1:num_request))

      requests_per_process = 0

      do i = 1, num_request
         ! use owner + 1, because owner is from 0 and array index from 1
         requests_per_process(non_local_node_owner(i) + 1) = requests_per_process(non_local_node_owner(i) + 1) + 1
      end do

      if (requests_per_process(my_rank + 1) .ne. 0) write (*, *) 'on rank', my_rank, 'requests for self is non-zero:', &
         requests_per_process(my_rank + 1)

!    write(*,*) 'requests 1:', my_rank, requests_per_process

      call MPI_ALLTOALL(requests_per_process, 1, MPI_INTEGER, requests_from_process, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

      total_num_requests_from_others = sum(requests_from_process)

!    write(*,*) 'requests 2:', my_rank, requests_from_process

      allocate (requested_keys(total_num_requests_from_others))

      disp = 0
      do i = 1, n_cpu
         sdispls(i) = disp
         disp = disp + requests_per_process(i)
      end do

      disp = 0
      do i = 1, n_cpu
         rdispls(i) = disp
         disp = disp + requests_from_process(i)
      end do

!    write (*,*) 'sdispls:', my_rank, sdispls
!    write (*,*) 'rdispls:', my_rank, rdispls

      call MPI_ALLTOALLV(non_local_node_keys, requests_per_process, sdispls, MPI_INTEGER8, &
                         requested_keys, requests_from_process, rdispls, MPI_INTEGER8, MPI_COMM_WORLD, ierr)

      !&<
      ! register propertyupdate data type
      blocklengths(1:nprops_property_update)  = [1, 1, 1, 1, 3, 1]
      types(1:nprops_property_update)         = [MPI_INTEGER8, MPI_INTEGER, MPI_REAL8, MPI_REAL8, MPI_REAL8, MPI_REAL8]
      call MPI_GET_ADDRESS(dummy_property_update,                  address(0), ierr)
      call MPI_GET_ADDRESS(dummy_property_update%key,              address(1), ierr)
      call MPI_GET_ADDRESS(dummy_property_update%owner,            address(2), ierr)
      call MPI_GET_ADDRESS(dummy_property_update%smoothing_length, address(3), ierr)
      call MPI_GET_ADDRESS(dummy_property_update%rho,              address(4), ierr)
      call MPI_GET_ADDRESS(dummy_property_update%v,                address(5), ierr)
      call MPI_GET_ADDRESS(dummy_property_update%temperature,      address(6), ierr)
      displacements(1:nprops_property_update) = int(address(1:nprops_property_update) - address(0))
      call MPI_TYPE_STRUCT(nprops_property_update, blocklengths, displacements, types, mpi_type_property_update, ierr)
      call MPI_TYPE_COMMIT(mpi_type_property_update, ierr)
      !&>

      allocate (packed_updates(total_num_requests_from_others), received_updates(num_request))

      ! test whether requested keys are locally known and pack everything together
      do i = 1, total_num_requests_from_others
         !TODO: it would be more elegant to include the node_index inside the owner`s t%nodes(:) array into the t_tree_node data
         !      structure for faster lookup. This is essentially done in the tree_communicator by keeping %first_child as a pointer into the
         !      owner`s nodes-array. Here, this is not really possible until now as first_child is occupied with real data here.
         if (.not. tree_traverse_to_key(t, requested_keys(i), iactual_node)) then
            DEBUG_ERROR(*, 'Did not find requested key in my local tree')
         end if

         actual_node => t%nodes(iactual_node)

         found = .false.
         do j = 1, np_local
            if (particles(j)%node_leaf .eq. iactual_node) then
               packed_updates(i) = t_property_update(actual_node%key, actual_node%owner, particles(j)%results%h, &
                                                     particles(j)%results%rho, actual_node%interaction_data%v, actual_node%interaction_data%temperature)

               ! test whether requested key is parent of particle key
               ! TODO: write a test function for this?
               ! write(*,'(i3,a,O30,O30)') my_rank, 'packing:', htable(actual_address)%key, particles(actual_node)%key
               found = .true.
               exit
            end if
         end do

         if (.not. found) then
            DEBUG_ERROR(*, 'Did not find node in my local keylist')
         end if

      end do

      disp = 0
      do i = 1, n_cpu
         sdispls(i) = disp
         disp = disp + requests_from_process(i)
      end do

      disp = 0
      do i = 1, n_cpu
         rdispls(i) = disp
         disp = disp + requests_per_process(i)
      end do

      call MPI_ALLTOALLV(packed_updates, requests_from_process, sdispls, mpi_type_property_update, &
                         received_updates, requests_per_process, rdispls, mpi_type_property_update, MPI_COMM_WORLD, ierr)

      do i = 1, num_request
         if (received_updates(i)%key .ne. non_local_node_keys(i)) write (*, *) 'Error in update on', my_rank, 'key mismatch', received_updates(i)%key, non_local_node_keys(i)

         actual_node => t%nodes(non_local_nodes(i))

         !&<
         actual_node%interaction_data%rho         = received_updates(i)%rho
         actual_node%interaction_data%temperature = received_updates(i)%temperature
         actual_node%interaction_data%v           = received_updates(i)%v
         actual_node%interaction_data%h           = received_updates(i)%smoothing_length
         !&>

      end do

      do i = 1, np_local
         !&<
         actual_node => t%nodes(particles(i)%node_leaf)
         actual_node%interaction_data%rho         = particles(i)%results%rho
         actual_node%interaction_data%temperature = particles(i)%data%temperature
         actual_node%interaction_data%h           = particles(i)%results%h
         !&>
      end do

      ! TODO: update tree_nodes%h for all parents

      deallocate (non_local_nodes, non_local_node_keys, non_local_node_owner, requests_per_process, int_arr, node_arr_cp, STAT=ierr)

      deallocate (requested_keys, STAT=ierr)

      deallocate (packed_updates, received_updates, STAT=ierr)

      call MPI_TYPE_FREE(mpi_type_property_update, ierr)

   end subroutine update_particle_props

end module module_sph

