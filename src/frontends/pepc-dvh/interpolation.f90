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

module interpolation_on_grid

   use module_pepc_kinds
   implicit none

contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !>
   !>  Interpolation of data on a Cartesian grid using Shepard renormalization
   !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine interpolation()

      use physvars
      use omp_lib
      use manipulate_particles
      use treevars, only: num_threads
      use mpi

      integer                             :: ierr, omp_num_threads
      integer(kind_particle)              :: i, k, l, i1, i2, i3
      integer(kind_particle), parameter   :: supp_int = 2     ! support for interpolation kernel
      integer(kind_particle)              :: proximity(3), n_interp_points, n_max_interp_points
      integer(kind_particle), allocatable :: mapping_indices(:), index_map(:, :, :), not_zeros(:)
      type(t_particle_short), allocatable :: m_part(:)
      real(kind_physics), parameter       :: c2 = 2.5d-1        ! interpolation factor
      real(kind_physics)                  :: maxximo, minnimo, ratio_ma, ratio_mi
      real(kind_physics)                  :: x, y, z, vol
      real(kind_physics)                  :: pos(3), omega(3)
      real(kind_physics)                  :: local_extent_min(3), local_extent_max(3)
      real(kind_physics)                  :: kernel, tentx, tenty, tentz
      real(kind_physics), allocatable     :: omega_reduction(:, :), positions(:,:), den_int(:), ome_app(:)
      real(kind_physics), dimension(3)    :: total_vort, total_vort_full_pre, total_vort_full_after
      real(kind_physics)                  :: total_vortmod, total_vortmod_full_pre, total_vortmod_full_after

      logical(1), allocatable :: grid_mask(:, :, :)
      character(len=2)        :: rank_label ! format descriptor

      local_extent_min(1) = minval(vortex_particles(1:np)%x(1))
      local_extent_max(1) = maxval(vortex_particles(1:np)%x(1))

      local_extent_min(2) = minval(vortex_particles(1:np)%x(2))
      local_extent_max(2) = maxval(vortex_particles(1:np)%x(2))

      local_extent_min(3) = minval(vortex_particles(1:np)%x(3))
      local_extent_max(3) = maxval(vortex_particles(1:np)%x(3))

      ! Safety margin - put buffer region around particles
      local_extent_max = local_extent_max + supp_int * m_h
      local_extent_min = local_extent_min - supp_int * m_h

      omp_num_threads = 1

      if (vort_check) then
         total_vort = 0.
         total_vortmod = 0.
         do i = 1, np
            total_vort(1:3) = total_vort(1:3) + vortex_particles(i)%data%alpha(1:3)
         end do

         call MPI_ALLREDUCE(total_vort, total_vort_full_pre, 3, MPI_KIND_PHYSICS, MPI_SUM, MPI_COMM_WORLD, ierr)
      end if

!$    call omp_set_num_threads(num_threads + 1)
!$    omp_num_threads = num_threads + 1
!$    if (my_rank .eq. 0) write (*, *) 'Using OpenMP with', omp_num_threads, 'threads.'

      allocate(grid_mask(nint(local_extent_min(1)/m_h):nint(local_extent_max(1)/m_h), &   !&
                         nint(local_extent_min(2)/m_h):nint(local_extent_max(2)/m_h), &   !&
                         nint(local_extent_min(3)/m_h):nint(local_extent_max(3)/m_h)   )) !&

      grid_mask = .false.

      ! Identify grid points that support the remeshing with new particles
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
!$OMP SHARED(vortex_particles, np, m_h) &
!$OMP PRIVATE(i1, i2, i3, k, proximity) &
!$OMP REDUCTION(.or.:grid_mask)
      do k = 1, np

         ! Find indices of nearest global grid point
         proximity = nint(vortex_particles(k)%x / m_h)

         do i3 = proximity(3) - supp_int, proximity(3) + supp_int
            do i2 = proximity(2) - supp_int, proximity(2) + supp_int
               do i1 = proximity(1) - supp_int, proximity(1) + supp_int

                  grid_mask(i1, i2, i3) = .true.

               end do
            end do
         end do

      end do
!$OMP END PARALLEL DO

      ! $OMP PARALLEL WORKSHARE DEFAULT(NONE)
      n_interp_points = count(grid_mask)
      ! $OMP END PARALLEL WORKSHARE

      call MPI_ALLREDUCE(n_interp_points, n_max_interp_points, 1, MPI_KIND_PARTICLE, MPI_MAX, MPI_COMM_WORLD, ierr)

      allocate(mapping_indices(n_interp_points))

      allocate(index_map(lbound(grid_mask,1):ubound(grid_mask,1), &
                         lbound(grid_mask,2):ubound(grid_mask,2), &
                         lbound(grid_mask,3):ubound(grid_mask,3) ))

      index_map = 0

      ! $OMP PARALLEL WORKSHARE DEFAULT(NONE)
      mapping_indices = (/(i, i=1, n_interp_points)/)
      index_map = unpack(mapping_indices, grid_mask, index_map)
      ! $OMP END PARALLEL WORKSHARE

      deallocate (mapping_indices)

      allocate (      positions(1:3,n_interp_points),  &
                omega_reduction(1:3,n_interp_points),  &
                        den_int(    n_interp_points) )

      positions       = 0.d0
      den_int         = 0.d0
      omega_reduction = 0.d0

!$OMP PARALLEL DEFAULT(NONE)                                      &
!$OMP PRIVATE(pos,x,y,z,i1,i2,i3,k,l,omega,vol,                   &
!$OMP         kernel,proximity,tentx,tenty,tentz)                 &
!$OMP SHARED(vortex_particles,np,m_h,index_map,positions)  &
!$OMP REDUCTION(+: omega_reduction, den_int)

!$OMP DO SCHEDULE(STATIC)
      do k = 1, np
         pos = vortex_particles(k)%x

         ! Find indexes of nearest global grid point
         proximity = nint(pos / m_h)

           vol = vortex_particles(k)%data%vol
         omega = vortex_particles(k)%data%alpha / vol

         do i3 = proximity(3) - supp_int, proximity(3) + supp_int
            z = m_h * i3
            tentz = tent(dabs(z - pos(3))/m_h,c2)

            do i2 = proximity(2) - supp_int, proximity(2) + supp_int
               y = m_h * i2
               tenty = tent(dabs(y - pos(2))/m_h,c2)

               do i1 = proximity(1) - supp_int, proximity(1) + supp_int
                  x = m_h * i1
                  tentx = tent(dabs(x - pos(1))/m_h,c2)

                  kernel = tentx * tenty * tentz

                  l = index_map(i1, i2, i3)

                  positions(:, l) = [x, y, z]

                  omega_reduction(:, l) = omega_reduction(:, l) + kernel * omega(:)

                  den_int(l) = den_int(l) + kernel

               end do
            end do
         end do

      end do
!$OMP END DO
!$OMP END PARALLEL

      deallocate (grid_mask, index_map, vortex_particles)

      not_zeros = pack([(i,i=1,n_interp_points)],den_int/=0.)

      allocate (m_part(n_max_interp_points))

      m_part(:)%data = t_particle_data_short([0.d0, 0.d0, 0.d0], 0.d0)

      n_interp_points = size(not_zeros)

! IN %WORK WE PUT THE DENOMINATOR OF INTERPOLATION: this is a smart way to have summation over two different quantities
      vol = m_h**3

!$OMP PARALLEL DEFAULT(NONE)                       &
!$OMP PRIVATE(i,l)                                 &
!$OMP SHARED(m_part, n_interp_points, vol,         &
!$OMP        positions, omega_reduction, den_int,  &
!$OMP        not_zeros)
!$OMP DO SCHEDULE(STATIC)
      do i = 1, n_interp_points
         l = not_zeros(i)
         m_part(i)%x(:)          = positions(:,l)
         m_part(i)%data%alpha(:) = omega_reduction(:,l) * vol
         m_part(i)%data%vol      = vol
         m_part(i)%work          = den_int(l)
!         m_part(i)%work          = 1.d0       ! uncomment this line do disable Shepard renormalization
      end do
!$OMP END DO
!$OMP END PARALLEL

      deallocate(not_zeros, positions, omega_reduction, den_int)

      call sort_remesh(m_part, local_extent_min, local_extent_max, n_interp_points, n_max_interp_points)

      np = n_interp_points

      if (vort_check) then
         total_vort = 0.
         do i = 1, np
            total_vort(1:3) = total_vort(1:3) + m_part(i)%data%alpha(1:3)
         end do

         ! Sum net vorticity after balancing and filtering
         call MPI_ALLREDUCE(total_vort, total_vort_full_after, 3, MPI_KIND_PHYSICS, MPI_SUM, MPI_COMM_WORLD, ierr)

         ! Output stats to check vorticity
         if (my_rank .eq. 0) then
            write (*, *) 'Gamma before interpolation :', norm2(total_vort_full_pre)   !, total_vortmod_full_pre
            write (*, *) 'Gamma after interpolation (scattering) :', norm2(total_vort_full_after) !, total_vortmod_full_after
         end if
      end if

      allocate(vortex_particles(np))
      vortex_particles(1:np) = m_part(1:np)

      ! Reset the number of OpenMP threads to num_threads, the number of WORK threads.
!$    call omp_set_num_threads(num_threads)

      deallocate (m_part)

      ! Check load-balance
      call MPI_ALLREDUCE(np, n, 1, MPI_KIND_PARTICLE, MPI_SUM, MPI_COMM_WORLD, ierr)
      if (1.25 * n / n_cpu .lt. np) then
         write (*, *) 'warning, rank', my_rank, ' appears to be heavily imbalanced:', 1.0 * np / (1.0 * n / n_cpu)
      end if

      if (my_rank .eq. 0) write (*, *) 'New total number of vortices ', np

      call reset_labels() ! works on vortex_particles

   end subroutine interpolation

   elemental function tent(adis,c2)

      real(kind_physics), intent(in) :: adis, c2
      real(kind_physics)             :: adis2, adis3
      real(kind_physics)             :: tent

      adis2 = adis*adis
      adis3 = adis*adis2

      tent = 0.d0

!     if(adis.lt.1.d0)                  tent = (1.d0 - adis)*(1.d0 + adis)*(2.d0 - adis) / 2.d0
!     if(adis.gt.1.d0.and.adis.lt.2.d0) tent = (1.d0 - adis)*(2.d0 - adis)*(3.d0 - adis) / 6.d0

      if(adis.lt.1.d0)       tent =       1.d0 - 2.5d0*adis2 + 1.5d0*adis3  &
                                  - c2 * (2.d0 - 9.0d0*adis2 + 6.0d0*adis3)

      if(adis.gt.1.d0 .and.                                                                     &
         adis.lt.2.d0      ) tent =      (2.d0 - adis)*(2.d0 - adis)*(1.d0 -      adis) / 2.d0  &
                                  - c2 * (2.d0 - adis)*(2.d0 - adis)*(1.d0 - 2.d0*adis)

   end function tent

end module interpolation_on_grid
