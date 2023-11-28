! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
!
! Copyright (C) 2002-2023 Juelich Supercomputing Centre,
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

module manipulate_particles

   use module_pepc_kinds
   implicit none

contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !>
   !>  Remeshing with a Gaussian diffusion process using local grids and parallel sorting
   !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine remeshing()

      use physvars
      use module_user_timings
      use omp_lib
      use treevars, only: num_threads
      use mpi
      implicit none

      integer :: ierr, m_dim, omp_num_threads
      integer(kind_particle) :: i1, i2, i3
      integer(kind_particle) :: i, k, l
      integer(kind_particle) :: proximity(3), n_remesh_points, n_max_remesh_points
      integer(kind_particle), allocatable :: mapping_indices(:), index_map(:, :, :)
      real(kind_physics)     :: x, y, z, deno
      real(kind_physics)     :: dist2, alpha_sum
      real(kind_physics)     :: frac, pos(3), alpha(3), work
      real(kind_physics)     :: local_extent_min(3), local_extent_max(3)
      real(kind_physics), dimension(3) :: total_vort, total_vort_full_pre, total_vort_full_mid, total_vort_full_after
      real(kind_physics)               :: total_vortmod, total_vortmod_full_pre, total_vortmod_full_mid, total_vortmod_full_after
      type(t_particle_short), allocatable :: m_part(:)
      real(kind_physics), allocatable :: m_part_reduction(:, :)
      logical(1), allocatable :: grid_mask(:, :, :)

      if (vort_check) then
         total_vort = 0.
         total_vortmod = 0.
         do i = 1, np
            total_vort(1:3) = total_vort(1:3) + vortex_particles(i)%data%alpha(1:3)
            total_vortmod = total_vortmod + dot_product(vortex_particles(i)%data%alpha, vortex_particles(i)%data%alpha)
         end do
         total_vortmod = sqrt(total_vortmod)

         call MPI_ALLREDUCE(total_vort, total_vort_full_pre, 3, MPI_KIND_PHYSICS, MPI_SUM, MPI_COMM_WORLD, ierr)
         call MPI_ALLREDUCE(total_vortmod, total_vortmod_full_pre, 1, MPI_KIND_PHYSICS, MPI_SUM, MPI_COMM_WORLD, ierr)
      end if

      deno = (pi * kernel_c)**1.5

      !DANILO: PARADIGM SHIFT. WITH CLASSIC POPULATION CONTROL THE ALLOCATION OF TEMP PARTICLES MAY INDUCE STRONG PAGINATION.
      !        WE USE SOMETHING SIMILAR TO DVH 2D: DEFINE A MESH CONTAINING WHOLE PARTICLES AND CYLCING OVER MESH POINTS.

      call timer_start(t_remesh_interpol)

      local_extent_min(1) = minval(vortex_particles(1:np)%x(1))
      local_extent_max(1) = maxval(vortex_particles(1:np)%x(1))
      local_extent_min(2) = minval(vortex_particles(1:np)%x(2))
      local_extent_max(2) = maxval(vortex_particles(1:np)%x(2))
      local_extent_min(3) = minval(vortex_particles(1:np)%x(3))
      local_extent_max(3) = maxval(vortex_particles(1:np)%x(3))

      ! Safety margin - put buffer region around particles
      local_extent_max = local_extent_max + (nDeltar + 1) * m_h
      local_extent_min = local_extent_min - (nDeltar + 1) * m_h

      allocate(grid_mask(nint(local_extent_min(1)/m_h):nint(local_extent_max(1)/m_h), & !&
                         nint(local_extent_min(2)/m_h):nint(local_extent_max(2)/m_h), & !&
                         nint(local_extent_min(3)/m_h):nint(local_extent_max(3)/m_h) )) !&

      if (my_rank .eq. 0) write (*, *) 'Storage size (Mbytes) of grid_mask', real(storage_size(grid_mask)) * size(grid_mask) / (8 * 1024 * 1024)
      grid_mask = .false.

      omp_num_threads = 1
      ! Set number of OpenMP threads to the same number as pthreads used in the walk +1
      ! We use `+1` here since the walk and communication thread are not active here, so can use
      ! all available threads.
!$    call omp_set_num_threads(num_threads + 1)

      ! Inform the user that OpenMP is used, and with how many threads
!$    omp_num_threads = num_threads + 1
!$    if (my_rank .eq. 0) write (*, *) 'Using OpenMP with', omp_num_threads, 'threads.'

      ! Identify grid points that support the remeshing with new particles
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
!$OMP SHARED(vortex_particles, m_h, np, nDeltar) &
!$OMP PRIVATE(i1, i2, i3, k, proximity) &
!$OMP REDUCTION(.or.:grid_mask)
      do k = 1, np

         ! Find indices of nearest global grid point
         proximity = nint(vortex_particles(k)%x / m_h)

         do i3 = proximity(3) - nDeltar, proximity(3) + nDeltar
            do i2 = proximity(2) - nDeltar, proximity(2) + nDeltar
               do i1 = proximity(1) - nDeltar, proximity(1) + nDeltar

                  grid_mask(i1, i2, i3) = .true.

               end do
            end do
         end do

      end do
!$OMP END PARALLEL DO

      ! $OMP PARALLEL WORKSHARE DEFAULT(NONE)
      ! Would `n_remesh_points` be in a reduction clause? Which one?
      n_remesh_points = count(grid_mask)
      ! $OMP END PARALLEL WORKSHARE

      allocate (mapping_indices(n_remesh_points))
      allocate(index_map(lbound(grid_mask, 1):ubound(grid_mask, 1), & !&
                         lbound(grid_mask, 2):ubound(grid_mask, 2), & !&
                         lbound(grid_mask, 3):ubound(grid_mask, 3) )) !&
      if (my_rank .eq. 0) write (*, *) 'Storage size (Mbytes) of index_map', real(storage_size(index_map)) * size(index_map) / (8 * 1024 * 1024)

      index_map = 0

      ! $OMP PARALLEL WORKSHARE DEFAULT(NONE)
      ! Would rev_lg be in a reduction clause?
      mapping_indices = (/(i, i=1, n_remesh_points)/)
      index_map = unpack(mapping_indices, grid_mask, index_map)
      ! $OMP END PARALLEL WORKSHARE

      deallocate (mapping_indices)

      ! Find new (guessed) maximum number of particles after load-balancing
      ! The guess assumes the local number will NOT be higher than the current max number of local particles
      ! across all ranks, while allowing for a 5% increase.
      call MPI_ALLREDUCE(n_remesh_points, n_max_remesh_points, 1, MPI_KIND_PARTICLE, MPI_MAX, MPI_COMM_WORLD, ierr)
      n_max_remesh_points = ceiling(1.05 * n_max_remesh_points)
      allocate (m_part(n_max_remesh_points))
      if (my_rank .eq. 0) write (*, *) 'Storage size (Mbytes) of m_part', real(storage_size(m_part)) * size(m_part) / (8 * 1024 * 1024)
      m_part(:)%data = t_particle_data_short([0.d0, 0.d0, 0.d0], 0.d0)
      m_part(:)%work = 0.d0
      allocate (m_part_reduction(1:4, 1:n_max_remesh_points))
      if (my_rank .eq. 0) write (*, *) 'Storage size (Mbytes) of m_part_reduction', real(storage_size(m_part_reduction)) * size(m_part_reduction) / (8 * 1024 * 1024)
      m_part_reduction = 0.d0

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP PRIVATE(pos, x, y, z, i1, i2, i3, l, alpha ,alpha_sum, work, k, frac, dist2, proximity) &
#ifdef USER_REDUCTION
!$OMP SHARED(vortex_particles, np, m_h, kernel_c, deno, index_map, nDeltar) &
!$OMP REDUCTION(+: m_part)
#else
!$OMP SHARED(vortex_particles, np, m_h, kernel_c, deno, index_map, m_part, nDeltar) &
!$OMP REDUCTION(+: m_part_reduction)
#endif

!$OMP DO SCHEDULE(STATIC)
      do k = 1, np
         pos = vortex_particles(k)%x

         ! Find indexes of nearest global grid point
         proximity = nint(pos / m_h)

         alpha = vortex_particles(k)%data%alpha
         work = vortex_particles(k)%work
         alpha_sum = 0.d0

         ! Compute kernel sum for circulation renormalization component-wise
         do i3 = proximity(3) - nDeltar, proximity(3) + nDeltar
            z = m_h * i3
            do i2 = proximity(2) - nDeltar, proximity(2) + nDeltar
               y = m_h * i2
               do i1 = proximity(1) - nDeltar, proximity(1) + nDeltar
                  x = m_h * i1

                  dist2 = (pos(1) - x)**2 + (pos(2) - y)**2 + (pos(3) - z)**2
                  frac = ip_kernel(dist2, kernel_c, deno)
                  alpha_sum = alpha_sum + frac

               end do
            end do
         end do

         ! Redistribute vortex circulation on new mesh points.
         ! Total circulation is conserved after remeshing.
         do i3 = proximity(3) - nDeltar, proximity(3) + nDeltar
            z = m_h * i3
            do i2 = proximity(2) - nDeltar, proximity(2) + nDeltar
               y = m_h * i2
               do i1 = proximity(1) - nDeltar, proximity(1) + nDeltar
                  x = m_h * i1

                  l = index_map(i1, i2, i3)

                  dist2 = (pos(1) - x)**2 + (pos(2) - y)**2 + (pos(3) - z)**2
                  frac = ip_kernel(dist2, kernel_c, deno)

                  m_part(l)%x = [x, y, z]
#ifdef USER_REDUCTION
                  m_part(l)%data%alpha = m_part(l)%data%alpha + frac * alpha / alpha_sum
                  m_part(l)%work       = m_part(l)%work + frac * work          !&
#else
                  m_part_reduction(1, l) = m_part_reduction(1, l) + frac * alpha(1) / alpha_sum
                  m_part_reduction(2, l) = m_part_reduction(2, l) + frac * alpha(2) / alpha_sum
                  m_part_reduction(3, l) = m_part_reduction(3, l) + frac * alpha(3) / alpha_sum
                  m_part_reduction(4, l) = m_part_reduction(4, l) + frac * work
#endif
               end do
            end do
         end do

      end do
!$OMP END DO
!$OMP END PARALLEL

#ifndef USER_REDUCTION
      m_part(:)%data%alpha(1) = m_part_reduction(1, :)      !&
      m_part(:)%data%alpha(2) = m_part_reduction(2, :)      !&
      m_part(:)%data%alpha(3) = m_part_reduction(3, :)      !&
      m_part(:)%work          = m_part_reduction(4, :)      !&
#endif
      m_part(:)%data%vol      = m_h**3                      !&

      deallocate (grid_mask, index_map, vortex_particles, m_part_reduction)

      if (vort_check) then
         ! Sum vorticity of remeshed points
         total_vort = 0.
         total_vortmod = 0.
         do i = 1, n_remesh_points
            total_vort(1:3) = total_vort(1:3) + m_part(i)%data%alpha(1:3)
            total_vortmod = total_vortmod + dot_product(m_part(i)%data%alpha, m_part(i)%data%alpha)
         end do
         total_vortmod = sqrt(total_vortmod)

         ! Sum net vorticity of remeshed points
         call MPI_ALLREDUCE(total_vort, total_vort_full_mid, 3, MPI_KIND_PHYSICS, MPI_SUM, MPI_COMM_WORLD, ierr)
         call MPI_ALLREDUCE(total_vortmod, total_vortmod_full_mid, 1, MPI_KIND_PHYSICS, MPI_SUM, MPI_COMM_WORLD, ierr)
      end if

      ! Reset the number of OpenMP threads to num_threads, the number of WORK threads.
!$    call omp_set_num_threads(num_threads)

      call timer_stop(t_remesh_interpol)
      if (my_rank .eq. 0) write (*, '("PEPC-DVH | ", a,f12.8,a)') 'Finished interpolation for remeshing after ', timer_read(t_remesh_interpol), ' seconds'

      call timer_start(t_remesh_sort)
      call sort_remesh(m_part, local_extent_min, local_extent_max, &
                       n_remesh_points, n_max_remesh_points) ! may change n_remesh_points to updated, balanced number
      call timer_stop(t_remesh_sort)

      ! Update module variable that stores number of vortices
      np = n_remesh_points
      if (my_rank .eq. 0) write (*, *) 'New total number of vortices ', np

      ! Copy remeshed points back onto module variable, free working space
      allocate (vortex_particles(np))
      vortex_particles(1:np) = m_part(1:n_remesh_points)
      deallocate (m_part)

      if (vort_check) then
         ! Sum local vorticity after balancing and filtering
         total_vort = 0.
         total_vortmod = 0.
         do i = 1, np
            total_vort(1:3) = total_vort(1:3) + vortex_particles(i)%data%alpha(1:3)
            total_vortmod = total_vortmod + dot_product(vortex_particles(i)%data%alpha, vortex_particles(i)%data%alpha)
         end do
         total_vortmod = sqrt(total_vortmod)

         ! Sum net vorticity after balancing and filtering
         call MPI_ALLREDUCE(total_vort, total_vort_full_after, 3, MPI_KIND_PHYSICS, MPI_SUM, MPI_COMM_WORLD, ierr)
         call MPI_ALLREDUCE(total_vortmod, total_vortmod_full_after, 1, MPI_KIND_PHYSICS, MPI_SUM, MPI_COMM_WORLD, ierr)

         ! Output stats to check vorticity
         if (my_rank .eq. 0) then
            write (*, *) '   Gamma before remeshing (norm,sum_mod):', sqrt(dot_product(total_vort_full_pre, total_vort_full_pre)), total_vortmod_full_pre
            write (*, *) ' Gamma after pop. control (norm,sum_mod):', sqrt(dot_product(total_vort_full_mid, total_vort_full_mid)), total_vortmod_full_mid
            write (*, *) '    Gamma after remeshing (norm,sum_mod):', sqrt(dot_product(total_vort_full_after, total_vort_full_after)), total_vortmod_full_after
         end if
      end if

      ! Check load-balance
      call MPI_ALLREDUCE(np, n, 1, MPI_KIND_PARTICLE, MPI_SUM, MPI_COMM_WORLD, ierr)
      if (1.25 * n / n_cpu .lt. np) then
         write (*, *) 'warning, rank', my_rank, ' appears to be heavily imbalanced:', 1.0 * np / (1.0 * n / n_cpu)
      end if

      call reset_labels()

   end subroutine remeshing

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !>
   !>   Interpolation function for remeshing
   !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elemental function ip_kernel(dist2, c, d)
      use physvars, only: pi, Rd

      implicit none

      real(kind_physics), intent(in) :: dist2, c, d
      real(kind_physics) :: ip_kernel

      ip_kernel = 0.D0
      if (dist2 .le. Rd * Rd) then
         ip_kernel = exp(-dist2 / c) / d
      end if

   end function ip_kernel

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !>
   !>   Sorting function for remeshing
   !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine sort_remesh(particles, local_extent_min, local_extent_max, m_np, m_nppm)

      use physvars
      use treevars, only: maxlevel, idim
      use module_sort, only: sort
      use mpi
      implicit none

      type(t_particle_short), intent(inout) :: particles(*)
      integer(kind_particle), intent(inout) :: m_np
      integer(kind_particle), intent(in) :: m_nppm

      integer :: ierr
      integer(kind_particle) :: i, j, k
      type(t_particle_short) :: bound_parts_loc(2), bound_parts(0:2 * n_cpu - 1), ship_parts(m_np), get_parts(m_nppm)
      integer :: prev, next, nbits
!TODO should this be kind_default or really kind_particle?!?!?
      integer :: irnkl2(m_nppm), indxl(m_np), irnkl(m_nppm)
      real(kind_physics) :: local_extent_min(3), local_extent_max(3), s, local_work(m_nppm)
      real(kind_physics) :: global_extent_max(3), global_extent_min(3), boxsize(3), boxsize_max
      real(kind_physics) :: omega, maxome, minome
      integer(kind_particle) :: ix(m_np), iy(m_np), iz(m_np)
      integer(kind_key) :: sorted_keys(m_nppm), local_keys(m_nppm)
      integer :: fposts(n_cpu + 1), gposts(n_cpu + 1), islen(n_cpu), irlen(n_cpu)
      integer(kind_particle) :: npnew, npold

      integer(kind_key) :: iplace
      real(kind_physics), dimension(3) :: total_vort, total_vort_full_after

      interface
         subroutine slsort_keys(nin, nmax, keys, workload, balance_weight, max_imbalance, nout, indxl, &
                                irnkl, scounts, rcounts, sdispls, rdispls, keys2, irnkl2, size, rank, comm)
            use module_pepc_kinds
            integer(kind_particle), intent(in) :: nin
            integer(kind_particle), intent(in) :: nmax
            integer(kind_key), intent(inout) :: keys(*)
            real*8, intent(inout) :: workload(*)
            integer(kind_default), intent(in) :: balance_weight
            real*8, intent(in) :: max_imbalance
            integer(kind_particle), intent(out) :: nout
            integer(kind_default), intent(out) :: indxl(*), irnkl(*), scounts(*), rcounts(*), sdispls(*), rdispls(*)
            integer(kind_key), intent(out) :: keys2(*)
            integer(kind_default), intent(out) :: irnkl2(*)
            integer(kind_pe), intent(in) :: size, rank
            integer(kind_default), intent(in) :: comm
         end subroutine slsort_keys
      end interface

      iplace = 2_8**(idim * maxlevel)

      indxl = 0
      irnkl = 0
      islen = 0
      irlen = 0

      ! Define wraps for ring network  0 -> 1 -> 2 -> ... ... -> num_pe-1 -> 0 ...
      if (my_rank .eq. 0) then
         prev = n_cpu - 1
      else
         prev = my_rank - 1
      end if

      if (my_rank .eq. n_cpu - 1) then
         next = 0
      else
         next = my_rank + 1
      end if

      ! Find global limits
      call MPI_ALLREDUCE(local_extent_min, global_extent_min, 3, MPI_KIND_PHYSICS, MPI_MIN, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(local_extent_max, global_extent_max, 3, MPI_KIND_PHYSICS, MPI_MAX, MPI_COMM_WORLD, ierr)

      boxsize = global_extent_max - global_extent_min

      boxsize_max = maxval(global_extent_max - global_extent_min)

      s = boxsize_max / 2**maxlevel       ! refinement length

        !! Start key generation
      ! TODO: Use module_spacefilling here, problem: this module uses treevars variables, but we do not (e.g. npp vs. m_np)

      ! ( global_extent_min ) is the translation vector from the tree box to the simulation region (in 1st octant)
      ix(1:m_np) = (particles(1:m_np)%x(1) - global_extent_min(1)) / s           ! partial keys
      iy(1:m_np) = (particles(1:m_np)%x(2) - global_extent_min(2)) / s           !
      iz(1:m_np) = (particles(1:m_np)%x(3) - global_extent_min(3)) / s

      ! construct keys by interleaving coord bits and add placeholder bit
      ! - note use of 64-bit constants to ensure correct arithmetic
      nbits = maxlevel + 1
      do j = 1, m_np
         local_keys(j) = iplace
         do i = 0, nbits - 1
            local_keys(j) = local_keys(j) &
                            + 8_8**i * (4_8 * ibits(iz(j), i, 1) + 2_8 * ibits(iy(j), i, 1) + 1_8 * ibits(ix(j), i, 1))
         end do
      end do

      ! Sort particles locally according to local_keys
      particles(1:m_np)%key = local_keys(1:m_np)
      call sort(local_keys(1:m_np), indxl(1:m_np))
      ship_parts(1:m_np) = particles(indxl(1:m_np))
      particles(1:m_np) = ship_parts(1:m_np)

      ! Use Parallel Sort by Parallel Search (SLSORT by M. Hofmann, Chemnitz)
      indxl = 0
      npold = m_np
      npnew = m_np
      local_keys(1:npold) = particles(1:npold)%key
      local_work(1:npold) = particles(1:npold)%work

      call slsort_keys(npold, m_nppm, local_keys, local_work, 0, 0.05D0, npnew, indxl, irnkl, islen, irlen, fposts, gposts, &
                       sorted_keys, irnkl2, n_cpu, my_rank, MPI_COMM_WORLD)

      ! Permute particles according to arrays from slsort
      m_np = npnew
      ship_parts(1:npold) = particles(indxl(1:npold))
      call MPI_ALLTOALLV(ship_parts, islen, fposts, MPI_TYPE_PARTICLE_SHORT_sca, &
                         get_parts, irlen, gposts, MPI_TYPE_PARTICLE_SHORT_sca, MPI_COMM_WORLD, ierr)
      particles(irnkl(1:m_np)) = get_parts(1:m_np)
      particles(1:m_np)%key = sorted_keys(1:m_np)

      ! TODO: MAYBE THIS ROUTINE IS USELESS. WE ARE NOW WORKING ON A SINGLE GLOBAL UNIFORM MESH
      ! Check if sort finished and find inner doublets
      k = 0
      do i = 2, m_np
         if (particles(i)%key .lt. particles(i - 1)%key) then
            write (*, '(a,i5,2z20)') 'Sorting globally failed: i,particles(i)%key,particles(i-1)%key=', i, particles(i)%key, particles(i - 1)%key
            call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
         else if (particles(i)%key .eq. particles(i - 1)%key) then
            particles(i)%data%alpha(1) = particles(i)%data%alpha(1) + particles(i - 1)%data%alpha(1)
            particles(i)%data%alpha(2) = particles(i)%data%alpha(2) + particles(i - 1)%data%alpha(2)
            particles(i)%data%alpha(3) = particles(i)%data%alpha(3) + particles(i - 1)%data%alpha(3)
            particles(i)%work = particles(i)%work + particles(i - 1)%work
         else
            k = k + 1
            particles(k) = particles(i - 1)
         end if
      end do
      ! Last original particle will always be transfered to next/final active particle
      k = k + 1
      particles(k) = particles(m_np)

      if (k .eq. 0) then
         write (*, *) 'something is wrong here: no particles left after compression', my_rank
         call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
      end if
      if (k .eq. 1) then
         write (*, *) 'Warning: Only one particle left on this process, candidate for m_np=0', my_rank
      end if

      bound_parts_loc(1) = particles(1)
      bound_parts_loc(2) = particles(k)
      call MPI_ALLGATHER(bound_parts_loc, 2, MPI_TYPE_PARTICLE_SHORT_sca, &
                         bound_parts, 2, MPI_TYPE_PARTICLE_SHORT_sca, MPI_COMM_WORLD, ierr)

      ! Eliminate right boundary, if doublet with at least on neighbor
      if (bound_parts(2 * my_rank + 1)%key .eq. bound_parts(2 * next)%key) then
         write (*, *) my_rank, 'is looking towards', next, 'for elimination'
         k = k - 1
      end if

      ! Accumulate boundary particles, looking left for doublets
      i = prev
      do while (bound_parts(2 * my_rank)%key .eq. bound_parts(2 * i + 1)%key)
         write (*, *) my_rank, 'is looking towards', i, 'for accumulation'
         particles(1)%data%alpha(1) = particles(1)%data%alpha(1) + bound_parts(2 * i + 1)%data%alpha(1)
         particles(1)%data%alpha(2) = particles(1)%data%alpha(2) + bound_parts(2 * i + 1)%data%alpha(2)
         particles(1)%data%alpha(3) = particles(1)%data%alpha(3) + bound_parts(2 * i + 1)%data%alpha(3)
         particles(1)%work = particles(1)%work + bound_parts(2 * i + 1)%work
         i = i - 1
         if (i .lt. 0) exit
      end do

      m_np = k

      ! KICK OUT CRITERION MUST BE ENFORCED ON vorticity (omega), not on circulation (alpha)
      ! Kick out particles (cannot use subroutinee here, since we work on a temp-array)
      k = 0
      maxome = -1.d2
      minome = 1.d2
      do i = 1, m_np
         omega = norm2(particles(i)%data%alpha) / particles(i)%data%vol
         if (omega .gt. thresh) then
            k = k + 1
            particles(k) = particles(i)
         else
            maxome = max(maxome, omega)
            minome = min(minome, omega)
         end if
      end do

      ! How many particles are kicked out? Write for check on threshold limit
      if(my_rank.eq.0) write(*,*) 'Particles kicked-out',m_np - k

      m_np = k

      ! Rebalance: Use Parallel Sort by Parallel Search (SLSORT by M. Hofmann, Chemnitz)
      indxl = 0
      npold = m_np
      npnew = m_np
      local_keys(1:npold) = particles(1:npold)%key
      local_work(1:npold) = particles(1:npold)%work

      call slsort_keys(npold, m_nppm, local_keys, local_work, 0, 0.05D0, npnew, indxl, irnkl, islen, irlen, fposts, gposts, &
                       sorted_keys, irnkl2, n_cpu, my_rank, MPI_COMM_WORLD)

      ! Permute particles according to arrays from slsort
      m_np = npnew
      ship_parts(1:npold) = particles(indxl(1:npold))
      call MPI_ALLTOALLV(ship_parts, islen, fposts, MPI_TYPE_PARTICLE_SHORT_sca, &
                         get_parts, irlen, gposts, MPI_TYPE_PARTICLE_SHORT_sca, MPI_COMM_WORLD, ierr)
      particles(irnkl(1:m_np)) = get_parts(1:m_np)
      particles(1:m_np)%key = sorted_keys(1:m_np)
      particles(1:m_np)%work = 1. !TODO: is this elegant?

   end subroutine sort_remesh

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !>
   !>   Kick out particles according to threshold (compare with |alpha|)
   !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine kick_out_particles()

      use physvars
      use mpi
      implicit none

      integer :: ierr
      integer(kind_particle) :: i, k
      real(kind_physics) :: omega
      type(t_particle), allocatable :: temp_particles(:)

      ! kick out particles with vorticity magnitude below threshold
      k = 0
      do i = 1, np
         omega = norm2(vortex_particles(i)%data%alpha) / vortex_particles(i)%data%vol
         if (omega .gt. thresh) then
            k = k + 1
            vortex_particles(k) = vortex_particles(i)
         end if
      end do
      np = k
      call MPI_ALLREDUCE(np, n, 1, MPI_KIND_PARTICLE, MPI_SUM, MPI_COMM_WORLD, ierr)

      ! shrink vortex_particles to new size
      call move_alloc(vortex_particles, temp_particles)
      vortex_particles = temp_particles(1:np)

      if (1.25 * n / n_cpu .lt. np) then
         write (*, *) 'warning, rank', my_rank, ' appears to be heavily imbalanced:', 1.0 * np / (1.0 * n / n_cpu)
      end if

   end subroutine kick_out_particles

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !>
   !>   Initialize/reset labels after manipulation, esp. after remeshing
   !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine reset_labels()

      use physvars
      use mpi
      implicit none

      integer(kind_particle) :: i, nscan
      integer :: ierr

      ! Define valid labels for all local particles
      nscan = 0
      call MPI_SCAN(np, nscan, 1, MPI_KIND_PARTICLE, MPI_SUM, MPI_COMM_WORLD, ierr)
      nscan = nscan - np
      do i = 1, np
         vortex_particles(i)%label = nscan + i
      end do

   end subroutine reset_labels

   subroutine direct_sum(np_local, particles, results, my_rank, n_cpu)

      use module_pepc_types
      use module_interaction_specific_types, only: t_particle_results
      use module_directsum
      use mpi
      implicit none

      integer(kind_particle), intent(in) :: np_local
      integer, intent(in) :: my_rank, n_cpu
      type(t_particle), intent(in) :: particles(1:np_local)
      type(t_particle_results), intent(out) :: results(1:np_local)

      integer(kind_particle) :: i
      type(t_particle_results), allocatable :: directresults(:)
      integer(kind_particle) :: indices(1:np_local)

      if (my_rank .eq. 0) write (*, '("PEPC-DVH | ", a)') 'Starting direct summation ...'

      do i = 1, np_local
         indices(i) = i
      end do

      allocate (directresults(int(np_local, kind=kind_particle)))
      call directforce(particles, indices, int(np_local, kind=kind_particle), directresults, MPI_COMM_WORLD)
      results(1:np_local) = directresults(1:np_local)

      deallocate (directresults)

      if (my_rank .eq. 0) write (*, '("PEPC-DVH | ", a)') '                          ... done!'

   end subroutine direct_sum

end module manipulate_particles
