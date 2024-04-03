! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
!
! Copyright (C) 2002-2024 Juelich Supercomputing Centre,
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

!>
!> helper module
!>
module helper
   use module_pepc_kinds
   use module_pepc_types
   use module_timings
   use module_random
   implicit none

   ! timing variables
   integer, parameter :: t_user_total       = t_userdefined_first     !&
   integer, parameter :: t_user_init        = t_userdefined_first + 1 !&
   integer, parameter :: t_user_step        = t_userdefined_first + 2 !&
   integer, parameter :: t_user_directsum   = t_userdefined_first + 3 !&
   integer, parameter :: t_user_particleio  = t_userdefined_first + 4 !&

   ! MPI variables
   integer(kind_pe) :: my_rank, n_ranks
   logical :: root

   ! time variables
   real*8 :: dt
   integer :: step

   ! control variables
   integer :: nt                   ! number of timesteps
   integer(kind_particle) :: tnp   ! total number of particles
   integer(kind_particle) :: np    ! local number of particles
   character(len=10) :: setup      ! string to define testing/benchmarking setup
   logical :: particle_output      ! turn vtk output on/off
   logical :: domain_output        ! turn vtk output on/off
   logical :: diag_test            ! check diagnostics for correct physics
   integer :: check_step           ! timestep to check histogram
   logical :: particle_test        ! check tree code results against direct summation
   logical :: reflecting_walls     ! reflect particles at walls
   integer :: diag_interval        ! number of timesteps between diagnostics
   integer :: io_interval          ! number of timesteps between IO
   real(kind_physics) :: plasma_dimensions(3) ! size of the simulation box

   integer :: particle_direct = -1 ! number of particle for direct summation

   ! particle data (position, velocity, mass, charge)
   type(t_particle), allocatable   :: particles(:)
   real(kind_physics), allocatable :: direct_L2(:)
   ! summaries of some particle data
   real(kind_physics)              :: e_kin, e_pot, v_max, r_max
   integer, parameter              :: n_bins = 50 ! number of bins for the histogram
   integer, dimension(0:n_bins)    :: velocity_histo, global_velocity_histo

   ! PRNG state, PER MPI RANK, NOT THREAD!
   type(random_state_t) :: rng_state

contains

   subroutine set_parameter()

      use module_pepc
      use module_interaction_specific, only: theta2, eps2, force_law, include_far_field_if_periodic
      implicit none

      integer, parameter :: fid = 12
      character(255)     :: para_file
      logical            :: read_para_file

      namelist /pepcbenchmark/ tnp, nt, dt, setup, check_step, diag_test, particle_test, particle_output, domain_output, reflecting_walls, particle_direct, diag_interval, io_interval, plasma_dimensions

      ! set default parameter values
      tnp               = 10000                   !&
      nt                = 25                      !&
      dt                = 1e-2                    !&
      setup             = 'benchmark'             !&
      check_step        = 7000                    !&
      diag_test         = .false.                 !&
      particle_test     = .false.                 !&
      particle_output   = .false.                 !&
      domain_output     = .false.                 !&
      reflecting_walls  = .false.                 !&
      diag_interval     = 1                       !&
      io_interval       = 1                       !&
      plasma_dimensions = (/1.0_8, 1.0_8, 1.0_8/) !&

      ! read in namelist file
      call pepc_read_parameters_from_first_argument(read_para_file, para_file)

      if (read_para_file) then
         if (root) write (*, '(a)') " == reading parameter file, section pepcbenchmark: ", para_file
         open (fid, file=para_file)
         read (fid, NML=pepcbenchmark)
         close (fid)
      else
         if (root) write (*, '(a)') " == no param file, using default parameters and writing template to 'params.template' "
         open (fid, file='params.template', status='replace')
         write (fid, *) "!"
         write (fid, *) "!==============================================================================="
         write (fid, *) "! FRONTEND PARAMETERS FOR BENCHMARK"
         write (fid, *) "!"
         write (fid, *) "! tnp                  : total number of particles [10000] <module_helper>"
         write (fid, *) "! nt                   : number of timesteps [25] <module_helper>"
         write (fid, *) "! dt                   : time step size [0.02] <module_helper>"
         write (fid, *) "! setup                : string to define testing/benchmarking setup ['benchmark'] <module_helper>"
         write (fid, *) "! check_step           : timestep to check histogram [7000] <module_helper>"
         write (fid, *) "!                        this is highly dependant on number of ranks and setup, so needs to evaluated"
         write (fid, *) "! diag_test            : check diagnostics for correct physics [.false.] <module_helper>"
         write (fid, *) "!                        this is highly dependant on number of ranks and setup, so needs to evaluated"
         write (fid, *) "! particle_test        : check tree code results against direct summation [.false.] <module_helper>"
         write (fid, *) "! particle_output      : turn vtk output on/off [.false.] <module_helper>"
         write (fid, *) "! domain_output        : turn vtk output on/off [.false.] <module_helper>"
         write (fid, *) "! reflecting_walls     : reflect particles at walls [.false.] <module_helper>"
         write (fid, *) "!                        switch between Coulomb explosion and more deterministic benchmark"
         write (fid, *) "! particle_direct      : number of particle for direct summation [-1] <module_helper>"
         write (fid, *) "!                        -1 to test all particles"
         write (fid, *) "! diag_interval        : number of timesteps between diagnostics [1] <module_helper>"
         write (fid, *) "! io_interval          : number of timesteps between I/O [1] <module_helper>"
         write (fid, *) "! plasma_dimensions(3) : size of the simulation box [1.0, 1.0, 1.0] <module_helper>"
         write (fid, *) "!"
         write (fid, NML=pepcbenchmark)
         call pepc_write_parameters(fid)
         close (fid)

      end if

      if (root) then
         write (*, '(a,i12)')       " == total number of particles : ", tnp               !&
         write (*, '(a,i12)')       " == number of time steps      : ", nt                !&
         write (*, '(a,es12.4)')    " == time step                 : ", dt                !&
         write (*, '(a,a12)')       " == setup                     : ", trim(setup)       !&
         write (*, '(a,i12)')       " == diag interval             : ", diag_interval     !&
         write (*, '(a,i12)')       " == IO interval               : ", io_interval       !&
         write (*, '(a,l12)')       " == diag test                 : ", diag_test         !&
         write (*, '(a,l12)')       " == particle test             : ", particle_test     !&
         write (*, '(a,l12)')       " == particle output           : ", particle_output   !&
         write (*, '(a,l12)')       " == domain output             : ", domain_output     !&
         write (*, '(a,l12)')       " == reflecting walls          : ", reflecting_walls  !&
         write (*, '(a,3(es12.4))') " == plasma dimensions         : ", plasma_dimensions !&
      end if

      call pepc_prepare(3_kind_dim)
   end subroutine set_parameter

   subroutine init_particles(p)
      implicit none

      type(t_particle), allocatable, intent(inout) :: p(:)
      integer(kind_particle) :: ip, np_pad
      integer :: rc, prior_ranks
      real*8 :: dummy
      real(kind_physics) :: pi, phi, theta, rnd(1:2), s

      if (root) write (*, '(a)') " == [init] init particles "

      ! set initial number of local particles
      np = tnp / n_ranks
      ! fill up to tnp on last rank - this does not load-balance, but makes sure labels are ok
      if ((my_rank + 1) .eq. n_ranks) then
         np_pad = tnp - n_ranks * np
      else
         np_pad = 0
      end if

      allocate (particles(np + np_pad), stat=rc)
      if (rc .ne. 0) write (*, *) " === particle allocation error!"

      allocate (direct_L2(np + np_pad), stat=rc)
      if (rc .ne. 0) write (*, *) " === direct_L2 allocation error!"
      direct_L2 = -1.0_8

      ! set random seed to MPI rank
      rng_state%idum = -(my_rank + 1)
      call random(dummy, rng_state)

      pi = 2.0_kind_physics * acos(0.0_kind_physics)

      select case (trim(setup))
      case ('test')
         ! setup for Coulomb explosion, physics test case
         ! have fixed seed for more reproducible results
         rng_state%idum = -1234
         call random(dummy, rng_state)
         ! get enough random numbers on each rank to have matching streams, no matter how many ranks we have
         ! this does mean higher ranks draw a lot of random numbers for /dev/null
         do prior_ranks = 0, my_rank - 1
            do ip = 1, np
               call random(rnd, rng_state)
            end do
         end do
         ! now start 'locally'
         do ip = 1, np + np_pad
            p(ip)%label = my_rank * (np) + ip - 1
            p(ip)%data%q = -1.0_8 * 2.0_8 * &
                           plasma_dimensions(1) * plasma_dimensions(2) * &
                           plasma_dimensions(3) / tnp
            p(ip)%data%m = 1.d-3
            if (p(ip)%data%q .gt. 0.0) p(ip)%data%m = p(ip)%data%m * 100.0_8

            ! obtain 2 random numbers
            call random(rnd, rng_state)

            ! get angles to point into sphere
            phi = rnd(1) * 2 * pi
            theta = acos(rnd(2) * 2 - 1)

            ! compute a radial coordinate from distribution Eq.(9) in Kaplan`s paper, \mu = 1
            s = (1./(1 - real(p(ip)%label) / tnp) - 1)**(1./3.)

            ! fill cartesian coordinates
            p(ip)%x = s*[cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta)] * plasma_dimensions

            p(ip)%data%v = 0.0_8

            p(ip)%work = 1.0_8
         end do
         np = np + np_pad
      case ('benchmark')
         ! setup for random qubic particle cloud, benchmarking case
         do ip = 1, np
            p(ip)%label = my_rank * np + ip - 1
            p(ip)%data%q = (-1.0_8 + 2.0_8 * MOD(p(ip)%label, 2_kind_particle)) * 2.0_8 * &
                           plasma_dimensions(1) * plasma_dimensions(2) * &
                           plasma_dimensions(3) / tnp
            p(ip)%data%m = 1.0_8
            if (p(ip)%data%q .gt. 0.0) p(ip)%data%m = p(ip)%data%m * 100.0_8

            call random(p(ip)%x, rng_state)
            p(ip)%x = p(ip)%x * plasma_dimensions

            call random_gauss(p(ip)%data%v, rng_state)
            p(ip)%data%v = p(ip)%data%v / sqrt(p(ip)%data%m)

            p(ip)%work = 1.0_8
         end do
      case default
         stop 'wrong/no setup chosen'
      end select
   end subroutine init_particles

   subroutine push_particles(p)
      use module_mirror_boxes
      implicit none

      type(t_particle), allocatable, intent(inout) :: p(:)
      integer(kind_particle) :: ip

      if (root) write (*, '(a)') " == [pusher] push particles "

      do ip = 1, np
         p(ip)%data%v = p(ip)%data%v + dt * p(ip)%data%q / p(ip)%data%m * p(ip)%results%e
         p(ip)%x = p(ip)%x + dt * p(ip)%data%v
      end do
   end subroutine push_particles

   subroutine check_energies_local(p)
      use module_mirror_boxes
      implicit none

      type(t_particle), allocatable, intent(in) :: p(:)

      integer(kind_particle) :: ip
      real(kind_physics)     :: v_(size(p(1)%data%v)), v, r, rel_gamma

      if (root) write (*, '(a)') " == [energies] compute local net energies "

      e_kin = 0._kind_physics
      e_pot = 0._kind_physics
      v_max = -1._kind_physics
      r_max = -1._kind_physics

      do ip = 1, np
         ! velocity at -1/2*dt to sync w/ potential
         v_ = p(ip)%data%v - dt * p(ip)%data%q / p(ip)%data%m * p(ip)%results%e * 0.5_kind_physics
         rel_gamma = sqrt(1 + dot_product(v_, v_))

         ! compute kinetic and potential energies (a.u.)
         e_kin = e_kin + (rel_gamma - 1) * p(ip)%data%m
         e_pot = e_pot + p(ip)%data%q * p(ip)%results%pot * 0.5_kind_physics

         ! keep maximum velocity and distance to origin
#ifdef __PGI
         v = sqrt(dot_product(p(ip)%data%v, p(ip)%data%v))
         if (v .gt. v_max) v_max = v
         r = sqrt(dot_product(p(ip)%x, p(ip)%x))
         if (r .gt. r_max) r_max = r
#else
         v = norm2(p(ip)%data%v)
         if (v .gt. v_max) v_max = v
         r = norm2(p(ip)%x)
         if (r .gt. r_max) r_max = r
#endif
      end do
   end subroutine check_energies_local

   subroutine check_energies()
      use mpi
      implicit none

      integer            :: info
      real(kind_physics) :: local(2), global(2)

      if (root) write (*, '(a)') " == [energies] gather global energies "

      local = [e_kin, e_pot]
      global = 0._kind_physics
      call MPI_ALLREDUCE(local, global, size(local), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
      e_kin = global(1)
      e_pot = global(2)
      local = [v_max, r_max]
      global = 0._kind_physics
      call MPI_ALLREDUCE(local, global, size(local), MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, info)
      v_max = global(1)
      r_max = global(2)

      if (root) then
         write (*, '(a, es12.4)') " == [energies]          kinetic energy: ", e_kin
         write (*, '(a, es12.4)') " == [energies]        potential energy: ", e_pot
         write (*, '(a, es16.8)') " == [energies]                  energy: ", e_pot + e_kin
         write (*, '(a, es12.4)') " == [energies]        maximum velocity: ", v_max
         write (*, '(a, es12.4)') " == [energies] maximum radial distance: ", r_max
      end if

   end subroutine check_energies

   subroutine histogram_local(p)
      implicit none

      type(t_particle), allocatable, intent(in) :: p(:)

      integer                :: v_pos
      integer(kind_particle) :: ip
      real(kind_physics)     :: r_bin, v_bin, r_bin_min, r_bin_max

      if (root) write (*, '(a)') " == [histogram] compute local histogram "

      ! init histogram bin sizes
      r_bin = r_max / n_bins
      v_bin = v_max / n_bins
      r_bin_min = 7.d0  ! this region for test particles has been id'd experimentally
      r_bin_max = 8.5d0

      ! clear histogram
      velocity_histo = 0

      do ip = 1, np
         ! check radial position
#ifdef __PGI
         if (sqrt(dot_product(p(ip)%x, p(ip)%x)) .ge. r_bin_min .and. sqrt(dot_product(p(ip)%x, p(ip)%x)) .le. r_bin_max) then
            v_pos = int(sqrt(dot_product(p(ip)%data%v, p(ip)%data%v)) / v_bin)
            velocity_histo(v_pos) = velocity_histo(v_pos) + 1
         end if
#else
         if (norm2(p(ip)%x) .ge. r_bin_min .and. norm2(p(ip)%x) .le. r_bin_max) then
            v_pos = int(norm2(p(ip)%data%v) / v_bin)
            velocity_histo(v_pos) = velocity_histo(v_pos) + 1
         end if
#endif
      end do
   end subroutine histogram_local

   subroutine test_histogram(check)
      use mpi
      implicit none

      logical, intent(in)           :: check
      logical                       :: peak
      integer                       :: v_pos, c_peaks, info, hist_unit

      peak = .false.
      c_peaks = 0

      ! gather histogram data from other ranks
      call MPI_REDUCE(velocity_histo, global_velocity_histo, size(velocity_histo), MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, info)

      ! estimate number of peaks in histogram
      if (root) then
         !! save data to file
         !open(newunit=hist_unit, file='histograms.log', position='append')
         !write(hist_unit, *) global_velocity_histo
         !close(hist_unit)
         ! loop over all bins
         do v_pos = 0, n_bins
            if (global_velocity_histo(v_pos) .gt. 0) then
               ! we have particles in this bin
               if (.not. peak) then
                  ! we did not have particles in bin before
                  peak = .true.
                  c_peaks = c_peaks + 1
               end if
            else
               ! we have no particles in this bin
               if (peak) then
                  ! we previously had particles
                  peak = .false.
               end if
            end if
         end do

         write (*, '(a, i0)') " == [histogram] number of peaks found: ", c_peaks
         if (check) then
            if (c_peaks .eq. 2) then
               write (*, '(a)') " == [histogram] check                : passed"
            else
               write (*, '(a)') " == [histogram] check                : failed"
            end if
         end if
      end if
   end subroutine test_histogram

   subroutine filter_particles(p)
      use mpi
      implicit none

      type(t_particle), allocatable, intent(inout) :: p(:)
      integer(kind_particle) :: ip
      integer :: id, ncoll, ncoll_total, ierr

      ncoll = 0

      do ip = 1, np
         do id = 1, 3
            if (p(ip)%x(id) .lt. 0.0_8) then
               p(ip)%x(id) = -p(ip)%x(id)
               p(ip)%data%v(id) = -p(ip)%data%v(id)
               ncoll = ncoll + 1
            else if (p(ip)%x(id) .gt. plasma_dimensions(id)) then
               p(ip)%x(id) = 2 * plasma_dimensions(id) - p(ip)%x(id)
               p(ip)%data%v(id) = -p(ip)%data%v(id)
               ncoll = ncoll + 1
            end if
         end do
      end do

      call mpi_reduce(ncoll, ncoll_total, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      if (root) write (*, '(a,i12)') " == [filter] total number of wall collisions      : ", ncoll_total
   end subroutine filter_particles

   subroutine test_particles()
      use module_pepc_types
      use module_directsum
      use mpi
      implicit none

      integer(kind_particle), allocatable   :: tindx(:)
      real*8, allocatable                   :: trnd(:)
      type(t_particle_results), allocatable :: trslt(:)
      integer(kind_particle)                :: tn, tn_global, ti
      integer                               :: rc
      real(kind_physics)                    :: L2sum_local, L2sum_global, L2

      call timer_start(t_user_directsum)

      if (allocated(direct_L2)) then
         deallocate (direct_L2)
      end if
      allocate (direct_L2(np))
      direct_L2 = -1.0_8

      if (particle_direct .eq. -1) then
         tn = np
      else
         tn = particle_direct / n_ranks
         if (my_rank .eq. (n_ranks - 1)) tn = tn + MOD(particle_direct, n_ranks)
      end if

      allocate (tindx(tn), trnd(tn), trslt(tn))

      if (particle_direct .eq. -1) then
         do ti = 1, tn
            tindx(ti) = ti
         end do
      else
         call random(trnd(1:tn), rng_state)
         tindx(1:tn) = int(trnd(1:tn) * (np - 1)) + 1
      end if

      call directforce(particles, tindx, tn, trslt, MPI_COMM_WORLD)

      L2sum_local = 0.0
      L2sum_global = 0.0
      do ti = 1, tn
         L2 = &
            (particles(tindx(ti))%results%e(1) - trslt(ti)%e(1))**2 + &
            (particles(tindx(ti))%results%e(2) - trslt(ti)%e(2))**2 + &
            (particles(tindx(ti))%results%e(3) - trslt(ti)%e(3))**2
         L2sum_local = L2sum_local + L2
         direct_L2(tindx(ti)) = L2
      end do

      call MPI_ALLREDUCE(tn, tn_global, 1, MPI_KIND_PARTICLE, MPI_SUM, MPI_COMM_WORLD, rc)
      call MPI_ALLREDUCE(L2sum_local, L2sum_global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)

      L2sum_global = sqrt(L2sum_global) / tn_global

      call timer_stop(t_user_directsum)
      if (root) then
         write (*, '(a,i12)') " == [direct test] number tested particles         : ", tn_global
         write (*, '(a,es12.4)') " == [direct test] L2 error in probed particles    : ", L2sum_global
         write (*, '(a,es12.4)') " == [direct test] time in test [s]                : ", timer_read(t_user_directsum)
      end if

      deallocate (tindx)
      deallocate (trnd)
      deallocate (trslt)
   end subroutine test_particles

   integer function vtk_step_of_step(step) result(vtk_step)
      use module_vtk
      implicit none

      integer, intent(in) :: step

      if (step .eq. 0) then
         vtk_step = VTK_STEP_FIRST
      else if (nt - 1 - step .lt. io_interval) then
         vtk_step = VTK_STEP_LAST
      else
         vtk_step = VTK_STEP_NORMAL
      end if
   end function vtk_step_of_step

   subroutine write_particles(p)
      use module_vtk_helpers
      use mpi
      implicit none

      type(t_particle), intent(in) :: p(:)

      integer :: vtk_step

      call timer_start(t_user_particleio)
      vtk_step = vtk_step_of_step(step)
      call vtk_write_particles("particles", MPI_COMM_WORLD, step, dt * step, vtk_step, p, coulomb_and_l2)
      call timer_stop(t_user_particleio)
      if (root) write (*, '(a,es12.4)') " == [write particles] time in vtk output [s]      : ", timer_read(t_user_particleio)

   contains

      subroutine coulomb_and_l2(d, r, vtkf)
         use module_vtk
         use module_interaction_specific_types
         implicit none

         type(t_particle_data), intent(in) :: d(:)
         type(t_particle_results), intent(in) :: r(:)
         type(vtkfile_unstructured_grid), intent(inout) :: vtkf

         call vtk_write_particle_data_results(d, r, vtkf)
         if (particle_test) call vtkf%write_data_array("L2 error", direct_L2(:))
      end subroutine
   end subroutine write_particles

   subroutine write_domain(p)
      use module_vtk
      use module_vtk_helpers
      use module_pepc, only: global_tree
      implicit none

      type(t_particle), allocatable, intent(in) :: p(:)

      integer :: vtk_step

      ! output of tree diagnostics
      vtk_step = vtk_step_of_step(step)
      call vtk_write_branches(step, dt * step, vtk_step, global_tree)
      call vtk_write_leaves(step, dt * step, vtk_step, global_tree)
      call vtk_write_spacecurve(step, dt * step, vtk_step, p)
   end subroutine write_domain

   subroutine random_gauss(list, state)
      implicit none

      real(kind_physics), intent(inout) :: list(:)
      type(random_state_t), intent(inout) :: state

      real(kind_physics) :: v(2), pi, r, p
      integer :: n, i

      pi = 2.0_kind_physics * acos(0.0_kind_physics)
      n = size(list)

      do i = 1, n, 2

         call random(v, state)

         r = sqrt(-2.0_8 * log(v(1)))
         p = 2.0_8 * pi * v(2)

         list(i) = r * sin(p)
         if ((i + 1) .le. n) list(i + 1) = r * cos(p)

      end do
   end subroutine

end module
