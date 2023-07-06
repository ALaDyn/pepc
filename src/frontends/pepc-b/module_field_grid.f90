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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>  Encapsulates gridded fields
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_field_grid

   implicit none
   save
   private

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   public densities        !< Densities
   public densities_2d     !< 2D densities
   public fields_2d        !< 2D densities
   public sum_fields       !< 3D fields
   public sum_fieldave     !< time-ave fields
   public sum_radial       !< 1D, radially symmetric fields
   public field_lineout    !< 1D lineouts

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

   !  =================================
   !
   !    3D Density gather for ion density
   !
   !  =================================

   subroutine densities

      use module_physvars
      use module_particle_props
      use mpi

      implicit none

      real :: rdx, rdy, rdz, dx, dy, dz, cweight
      real :: fx1, fx2, fy1, fy2, fz1, fz2, xa, ya, za
      integer :: i, ng, i1, i2, j1, j2, k1, k2
      integer :: ierr

      dx = xl / ngx
      dy = yl / ngy
      dz = zl / ngz
      rdx = 1./dx
      rdy = 1./dy
      rdz = 1./dz

      !  field box limits: (0-xl, 0-yl, 0-zl)
      !  Any particle outside gets put in ghost cells 0, ngx+1

      !      write(15,'(//a,3f12.3)') 'cw,dx,dy',cweight,dx,dy

      rhoi_loc(0:ngx + 1, 0:ngy + 1, 0:ngz + 1) = 0.

      do i = 1, np_local

         xa = x(i) * rdx
         ya = y(i) * rdy
         za = z(i) * rdz

         !  indices
         i1 = xa + 1
         i2 = i1 + 1
         j1 = ya + 1
         j2 = j1 + 1
         k1 = za + 1
         k2 = k1 + 1

         i1 = min(max(0, i1), ngx + 1)
         i2 = min(max(0, i2), ngx + 1)
         j1 = min(max(0, j1), ngy + 1)
         j2 = min(max(0, j2), ngy + 1)
         k1 = min(max(0, k1), ngz + 1)
         k2 = min(max(0, k2), ngz + 1)

         !  linear weighting
         fx2 = min(max(i1 - xa, 0.), 1.)  ! Prevent overflow/negative weighting for particles outside box
         fx1 = 1.-fx2
         fy2 = min(max(j1 - ya, 0.), 1.)
         fy1 = 1.-fy2
         fz2 = min(max(k1 - za, 0.), 1.)
         fz1 = 1.-fz2

         !  gather charge at nearest grid points
         if (q(i) .gt. 0) then
            cweight = q(i) * rdx * rdy * rdz       ! charge weighting factor

            rhoi_loc(i1, j1, k1) = rhoi_loc(i1, j1, k1) + cweight * fx1 * fy1 * fz1
            rhoi_loc(i2, j1, k1) = rhoi_loc(i2, j1, k1) + cweight * fx2 * fy1 * fz1
            rhoi_loc(i1, j2, k1) = rhoi_loc(i1, j2, k1) + cweight * fx1 * fy2 * fz1
            rhoi_loc(i2, j2, k1) = rhoi_loc(i2, j2, k1) + cweight * fx2 * fy2 * fz1
            rhoi_loc(i1, j1, k2) = rhoi_loc(i1, j1, k2) + cweight * fx1 * fy1 * fz2
            rhoi_loc(i2, j1, k2) = rhoi_loc(i2, j1, k2) + cweight * fx2 * fy1 * fz2
            rhoi_loc(i1, j2, k2) = rhoi_loc(i1, j2, k2) + cweight * fx1 * fy2 * fz2
            rhoi_loc(i2, j2, k2) = rhoi_loc(i2, j2, k2) + cweight * fx2 * fy2 * fz2
         else

         end if
      end do

      ng = (ngx + 2) * (ngy + 2) * (ngz + 2)                         ! total # gridpoints
      ! gather on root
      call MPI_REDUCE(rhoi_loc, rhoi, ng, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      !  call MPI_REDUCE(rhoe_loc, rhoe, ng, MPI_REAL, MPI_SUM, 0,  MPI_COMM_WORLD, ierr)

   end subroutine densities

   !  =================================
   !
   !    2D Density gather for ion & electron density
   !
   !  =================================

   subroutine densities_2d

      use module_physvars
      use module_particle_props
      use mpi

      implicit none

      real :: rdx, rdy, dx, dy, cweight
      real :: fx1, fx2, fy1, fy2, xa, ya
      integer :: i, ng, i1, i2, j1, j2

      dx = xl / ngx
      dy = yl / ngy

      rdx = 1./dx
      rdy = 1./dy

      !  field box limits: (0-xl, 0-yl)
      !  Any particle outside gets put in ghost cells 0, ngx+1

      !      write(15,'(//a,3f12.3)') 'cw,dx,dy',cweight,dx,dy

      rhoi2d_loc(0:ngx + 1, 0:ngy + 1) = 0.
      rhoe2d_loc(0:ngx + 1, 0:ngy + 1) = 0.

      do i = 1, np_local

         xa = x(i) * rdx
         ya = y(i) * rdy

         !  indices
         i1 = xa + 1
         i2 = i1 + 1
         j1 = ya + 1
         j2 = j1 + 1

         i1 = min(max(0, i1), ngx + 1)
         i2 = min(max(0, i2), ngx + 1)
         j1 = min(max(0, j1), ngy + 1)
         j2 = min(max(0, j2), ngy + 1)

         !  linear weighting
         fx2 = min(max(i1 - xa, 0.), 1.)  ! Prevent overflow/negative weighting for particles outside box
         fx1 = 1.-fx2
         fy2 = min(max(j1 - ya, 0.), 1.)
         fy1 = 1.-fy2
         cweight = q(i) * rdx * rdy      ! charge weighting factor
         !  gather charge at nearest grid points
         if (q(i) .gt. 0) then
            rhoi2d_loc(i1, j1) = rhoi2d_loc(i1, j1) + cweight * fx1 * fy1
            rhoi2d_loc(i2, j1) = rhoi2d_loc(i2, j1) + cweight * fx2 * fy1
            rhoi2d_loc(i1, j2) = rhoi2d_loc(i1, j2) + cweight * fx1 * fy2
            rhoi2d_loc(i2, j2) = rhoi2d_loc(i2, j2) + cweight * fx2 * fy2
         else
            rhoe2d_loc(i1, j1) = rhoe2d_loc(i1, j1) + cweight * fx1 * fy1
            rhoe2d_loc(i2, j1) = rhoe2d_loc(i2, j1) + cweight * fx2 * fy1
            rhoe2d_loc(i1, j2) = rhoe2d_loc(i1, j2) + cweight * fx1 * fy2
            rhoe2d_loc(i2, j2) = rhoe2d_loc(i2, j2) + cweight * fx2 * fy2
         end if
      end do

      ng = (ngx + 2) * (ngy + 2)                         ! total # gridpoints

   end subroutine densities_2d

   !  =================================
   !
   !    2D fields and potential
   !
   !  =================================

   subroutine fields_2d

      use module_physvars
      use module_particle_props
      use module_pepc_wrappers
      use module_interaction_specific, only: eps2
      use module_pepc_types
      use mpi
      implicit none

      integer :: ng_total ! total # grid points
      integer :: ngp ! local # gp
      integer, allocatable :: ngps(:), igap(:) ! gps and strides
      integer :: ng_rest
      integer :: i, ierr, j, k, grid_ind
      real :: dx, dy
      logical :: field_debug = .false.
      real*8, allocatable :: p_x(:), p_y(:), p_z(:), p_ex(:), p_ey(:), p_ez(:), p_pot(:)
      integer, allocatable :: p_label(:)

      ! Create grid of marker particles to extract fields and potential from tree

      ng_total = ngx * ngy
      ngp = ng_total / n_cpu
      ng_rest = mod(ng_total, n_cpu)

      dx = xl / ngx
      dy = yl / ngy

      if (my_rank .lt. ng_rest) ngp = ngp + 1

      allocate (ngps(n_cpu + 3), igap(n_cpu + 3))
      allocate (p_x(ngp), p_y(ngp), p_z(ngp), p_ex(ngp), p_ey(ngp), p_ez(ngp), p_pot(ngp), p_label(ngp))

      ! Create array of local gps for strides
      call mpi_allgather(ngp, 1, MPI_INTEGER, ngps, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

      ! work out stride lengths so that partial arrays placed sequentially in global array

      igap(1) = 0
      igap(2) = ngps(1)
      do i = 3, n_cpu
         igap(i) = SUM(ngps(1:i - 1))
      end do

      ! Compute my set of coordinates
      do k = 1, ngps(my_rank + 1)
         grid_ind = k + igap(my_rank + 1)
         i = mod(grid_ind - 1, ngx) + 1  ! x-index
         j = (grid_ind - 1) / ngx   ! y-index
         p_x(k) = i * dx - dx / 2.
         p_y(k) = j * dy + dy / 2.
         p_z(k) = plasma_centre(3)
         p_label(k) = grid_ind
      end do

      call pepc_grid_fields_coulomb_wrapper(ngp, p_x, p_y, p_z, p_label, &
                                            p_Ex, p_Ey, p_Ez, p_pot, &
                                            num_neighbour_boxes, neighbour_boxes, force_const)

      if (field_debug) then
         write (*, '(a7,a50/2i5,3f15.2,i2)') 'PEPC | ', 'Params: itime, mac, theta, eps:', &
            itime, mac, theta, sqrt(eps2)
         write (*, '(a7,a20/(i10,5(1pe14.5)))') 'PEPC | ', 'fields: ', (p_label(i), p_x(i), p_y(i), p_ex(i), p_ey(i), p_pot(i), i=1, ngp)
      end if

      ! TODO insert particle field results into 2d arrays on root for dump_fields_2d
      call mpi_allgatherv(p_Ex, ngp, MPI_DOUBLE_PRECISION, ex2d, ngps, igap, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
      call mpi_allgatherv(p_Ey, ngp, MPI_DOUBLE_PRECISION, ey2d, ngps, igap, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
      call mpi_allgatherv(p_pot, ngp, MPI_DOUBLE_PRECISION, pot2d, ngps, igap, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)

      !  if (field_debug) then
      !    do j=1,ngy
      !        do i=1,ngx
      !        write (*,*) j,i,ex2d(i,j)
      !        end do
      !    end do
      !  endif
   end subroutine fields_2d

   !  =================================
   !
   !    1D gather for time-averaged electric fields
   !
   !  =================================

   subroutine sum_fieldave

      use module_physvars
      use module_particle_props
      use mpi
      implicit none

      real*8, dimension(0:ngav + 1) :: ex_g
      real :: dx, dy
      integer :: ndum, i, p
      integer, dimension(ngav + 1) :: pshortl

      if (my_rank .eq. 0) then
         write (*, '(//a/a,f10.2)') '1D field average '
!          ' fields on grid 0-',xl
      end if

      !  field box limits

      dx = (xgav_end - xgav_start) / ngav

      ! Axial electric field
      ! --------------------
      ! - set up dummy particles at grid points & find forces directly from tree
      ! Dummies set up at end of particle arrays on root to ensure unique labelling

      if (my_rank .eq. 0) then
         do i = 1, ngav + 1

            p = np_local + i   !index
            pshortl(i) = p   !index
            x(p) = (i - 1) * dx + xgav_start   ! axial position in box
            y(p) = plasma_centre(2)  ! target centre
            z(p) = plasma_centre(3)
         end do

         ! Get interaction lists
         !     write (*,*) 'Doing lists for dummy particles'
         !     write (*,'((i8,f12.3))') (pshortl(i),x(pshortl(i)),i=1,ngav+1)
         ndum = ngav + 1
      else
         ndum = 0  ! Remainder of CPUs just have to provide multipole info
      end if

      ! all CPUs must call walk
      !subroutine tree_walk(pshort,npshort, !pass,theta,eps,itime,mac,twalk,tfetch,ex_nps,ey_nps,ez_nps,np_local)

      ! mac set to BH, theta to 0.5

      !   call tree_walk(pshortl(1:ndum),ndum,1,0.5,eps,itime,0,ttrav,tfetch)

      if (my_rank .eq. 0) then
         ! Fields
         do i = 1, ngav + 1
            p = pshortl(i)
            !        call sum_force(p, nterm(i), nodelist( 1:nterm(i),i), eps, &
            !             ex_g(i-1), ey_g(i-1), ez_g(i-1), phi_g(i-1), w_g(i-1))
         end do

         ex_ave = ex_ave + force_const * Ex_g / navcycle    ! Accumulate axial field
      end if

      ! Radial electric field
      ! ----------------------
      ! - set up dummy particles at grid points & find forces directly from tree
      ! Dummies set up at end of particle arrays on root to ensure unique labelling

      ! TODO needs new kernel routine !

      dy = yl / ngav
      if (my_rank .eq. 0) then
         do i = 1, ngav + 1

            p = np_local + i   !index
            pshortl(i) = p   !index
            x(p) = xgav_pos(1)   ! 1st axial position in box
            y(p) = dy * (i - 1)  ! y position
            z(p) = plasma_centre(3) ! midpoint in z
         end do
         ! Get interaction lists
         ndum = ngav + 1
      else
         ndum = 0  ! Remainder of CPUs just have to provide multipole info
      end if

      ! all CPUs must call walk

      !   call tree_walk(pshortl(1:ndum),ndum,1,0.5,eps,itime,0,ttrav,tfetch)

      if (my_rank .eq. 0) then
         ! Fields
         do i = 1, ngav + 1
            p = pshortl(i)
            !        call sum_force(p, nterm(i), nodelist( 1:nterm(i),i), eps, &
            !             ex_g(i-1), ey_g(i-1), ez_g(i-1), phi_g(i-1), w_g(i-1))
         end do

         ! Accumulate radial field
         ey_ave(1:ngav, 1) = ey_ave(1:ngav, 1) + force_const * Ex_g(1:ngav) / navcycle
      end if

   end subroutine sum_fieldave

   !  =================================
   !
   !    3D gather for DC fields
   !  =================================

   subroutine sum_fields

      use module_physvars
      use module_particle_props

      implicit none

      real :: rdx, rdy, rdz, dx, dy, dz, cweight, jxweight, jyweight, jzweight
      real :: tweight

      real :: fx1, fx2, fy1, fy2, fz1, fz2, xa, ya, za, gamma, fr1, fr2
      integer :: i, j, k, i1, i2, j1, j2, k1, k2, nelecs, nions
      real, dimension(0:ngx + 1, 0:ngy + 1, 0:ngz + 1) :: ex_w, ey_w, ez_w
      real, dimension(0:ngx + 1, 0:ngy + 1, 0:ngz + 1) :: bx_w, by_w, bz_w
      real :: gmin = 1.e-3

      dx = xl / ngx
      dy = yl / ngy
      dz = zl / ngz
      rdx = 1./dx
      rdy = 1./dy
      rdz = 1./dz

      do k = 1, ngz
         do j = 1, ngy
            do i = 1, ngx
               ex_w(i, j, k) = 0.
               ey_w(i, j, k) = 0.
               ez_w(i, j, k) = 0.
               bx_w(i, j, k) = 0.
               by_w(i, j, k) = 0.
               bz_w(i, j, k) = 0.
               g_ion(i, j, k) = 0.
               g_ele(i, j, k) = 0.
               rhoi_loc(i, j, k) = 0.
               rhoe_loc(i, j, k) = 0.
               jxe_loc(i, j, k) = 0.
               jye_loc(i, j, k) = 0.
               jze_loc(i, j, k) = 0.
               te_loc(i, j, k) = 0.
               ti_loc(i, j, k) = 0.
            end do
         end do
      end do

      !  field box limits: (0-xl, 0-yl, 0-zl)
      !  Any particle outside gets put in ghost cells 0, ngx+1

      !      write(15,'(//a,3f12.3)') 'cw,dx,dy',cweight,dx,dy

      do i = 1, np_local

         xa = x(i) * rdx
         ya = y(i) * rdy
         za = z(i) * rdz

         !  indices
         i1 = xa + 1
         i2 = i1 + 1
         j1 = ya + 1
         j2 = j1 + 1
         k1 = za + 1
         k2 = k1 + 1

         i1 = min(max(0, i1), ngx + 1)
         i2 = min(max(0, i2), ngx + 1)
         j1 = min(max(0, j1), ngy + 1)
         j2 = min(max(0, j2), ngy + 1)
         k1 = min(max(0, k1), ngz + 1)
         k2 = min(max(0, k2), ngz + 1)

         !  linear weighting
         fx2 = min(max(i1 - xa, 0.), 1.)  ! Prevent overflow/negative weighting for particles outside box
         fx1 = 1.-fx2
         fy2 = min(max(j1 - ya, 0.), 1.)
         fy1 = 1.-fy2
         fz2 = min(max(k1 - za, 0.), 1.)
         fz1 = 1.-fz2
         !  gather charge at nearest grid points
         gamma = sqrt(1.0 + ux(i)**2 + uy(i)**2 + uz(i)**2)
         cweight = abs(q(i)) * rdx * rdy * rdz       ! charge weighting factor

         jxweight = cweight * ux(i) / gamma
         jyweight = cweight * uy(i) / gamma
         jzweight = cweight * uz(i) / gamma
         tweight = (gamma - 1.)  ! K.E. of particle in keV
         fr1 = sqrt(fx1**2 + fy1**2 + fz1**2 + eps**2)
         fr2 = sqrt(fx2**2 + fy2**2 + fz2**2 + eps**2)

         if (q(i) .lt. 0) then
            g_ele(i1, j1, k1) = g_ele(i1, j1, k1) + fx1 * fy1 * fz1  ! weighted # electrons
            g_ele(i2, j1, k1) = g_ele(i2, j1, k1) + fx2 * fy1 * fz1
            g_ele(i1, j2, k1) = g_ele(i1, j2, k1) + fx1 * fy2 * fz1
            g_ele(i2, j2, k1) = g_ele(i2, j2, k1) + fx2 * fy2 * fz1
            g_ele(i1, j1, k2) = g_ele(i1, j1, k2) + fx1 * fy1 * fz2
            g_ele(i2, j1, k2) = g_ele(i2, j1, k2) + fx2 * fy1 * fz2
            g_ele(i1, j2, k2) = g_ele(i1, j2, k2) + fx1 * fy2 * fz2
            g_ele(i2, j2, k2) = g_ele(i2, j2, k2) + fx2 * fy2 * fz2

            rhoe_loc(i1, j1, k1) = rhoe_loc(i1, j1, k1) + cweight * fx1 * fy1 * fz1
            rhoe_loc(i2, j1, k1) = rhoe_loc(i2, j1, k1) + cweight * fx2 * fy1 * fz1
            rhoe_loc(i1, j2, k1) = rhoe_loc(i1, j2, k1) + cweight * fx1 * fy2 * fz1
            rhoe_loc(i2, j2, k1) = rhoe_loc(i2, j2, k1) + cweight * fx2 * fy2 * fz1
            rhoe_loc(i1, j1, k2) = rhoe_loc(i1, j1, k2) + cweight * fx1 * fy1 * fz2
            rhoe_loc(i2, j1, k2) = rhoe_loc(i2, j1, k2) + cweight * fx2 * fy1 * fz2
            rhoe_loc(i1, j2, k2) = rhoe_loc(i1, j2, k2) + cweight * fx1 * fy2 * fz2
            rhoe_loc(i2, j2, k2) = rhoe_loc(i2, j2, k2) + cweight * fx2 * fy2 * fz2

            jxe_loc(i1, j1, k1) = jxe_loc(i1, j1, k1) + jxweight * fx1 * fy1 * fz1
            jxe_loc(i2, j1, k1) = jxe_loc(i2, j1, k1) + jxweight * fx2 * fy1 * fz1
            jxe_loc(i1, j2, k1) = jxe_loc(i1, j2, k1) + jxweight * fx1 * fy2 * fz1
            jxe_loc(i2, j2, k1) = jxe_loc(i2, j2, k1) + jxweight * fx2 * fy2 * fz1
            jxe_loc(i1, j1, k2) = jxe_loc(i1, j1, k2) + jxweight * fx1 * fy1 * fz2
            jxe_loc(i2, j1, k2) = jxe_loc(i2, j1, k2) + jxweight * fx2 * fy1 * fz2
            jxe_loc(i1, j2, k2) = jxe_loc(i1, j2, k2) + jxweight * fx1 * fy2 * fz2
            jxe_loc(i2, j2, k2) = jxe_loc(i2, j2, k2) + jxweight * fx2 * fy2 * fz2

            jye_loc(i1, j1, k1) = jye_loc(i1, j1, k1) + jyweight * fx1 * fy1 * fz1
            jye_loc(i2, j1, k1) = jye_loc(i2, j1, k1) + jyweight * fx2 * fy1 * fz1
            jye_loc(i1, j2, k1) = jye_loc(i1, j2, k1) + jyweight * fx1 * fy2 * fz1
            jye_loc(i2, j2, k1) = jye_loc(i2, j2, k1) + jyweight * fx2 * fy2 * fz1
            jye_loc(i1, j1, k2) = jye_loc(i1, j1, k2) + jyweight * fx1 * fy1 * fz2
            jye_loc(i2, j1, k2) = jye_loc(i2, j1, k2) + jyweight * fx2 * fy1 * fz2
            jye_loc(i1, j2, k2) = jye_loc(i1, j2, k2) + jyweight * fx1 * fy2 * fz2
            jye_loc(i2, j2, k2) = jye_loc(i2, j2, k2) + jyweight * fx2 * fy2 * fz2

            jze_loc(i1, j1, k1) = jze_loc(i1, j1, k1) + jzweight * fx1 * fy1 * fz1
            jze_loc(i2, j1, k1) = jze_loc(i2, j1, k1) + jzweight * fx2 * fy1 * fz1
            jze_loc(i1, j2, k1) = jze_loc(i1, j2, k1) + jzweight * fx1 * fy2 * fz1
            jze_loc(i2, j2, k1) = jze_loc(i2, j2, k1) + jzweight * fx2 * fy2 * fz1
            jze_loc(i1, j1, k2) = jze_loc(i1, j1, k2) + jzweight * fx1 * fy1 * fz2
            jze_loc(i2, j1, k2) = jze_loc(i2, j1, k2) + jzweight * fx2 * fy1 * fz2
            jze_loc(i1, j2, k2) = jze_loc(i1, j2, k2) + jzweight * fx1 * fy2 * fz2
            jze_loc(i2, j2, k2) = jze_loc(i2, j2, k2) + jzweight * fx2 * fy2 * fz2

            Te_loc(i1, j1, k1) = Te_loc(i1, j1, k1) + tweight * fx1 * fy1 * fz1
            Te_loc(i2, j1, k1) = Te_loc(i2, j1, k1) + tweight * fx2 * fy1 * fz1
            Te_loc(i1, j2, k1) = Te_loc(i1, j2, k1) + tweight * fx1 * fy2 * fz1
            Te_loc(i2, j2, k1) = Te_loc(i2, j2, k1) + tweight * fx2 * fy2 * fz1
            Te_loc(i1, j1, k2) = Te_loc(i1, j1, k2) + tweight * fx1 * fy1 * fz2
            Te_loc(i2, j1, k2) = Te_loc(i2, j1, k2) + tweight * fx2 * fy1 * fz2
            Te_loc(i1, j2, k2) = Te_loc(i1, j2, k2) + tweight * fx1 * fy2 * fz2
            Te_loc(i2, j2, k2) = Te_loc(i2, j2, k2) + tweight * fx2 * fy2 * fz2

         else
            g_ion(i1, j1, k1) = g_ion(i1, j1, k1) + fx1 * fy1 * fz1  ! weighted # ions
            g_ion(i2, j1, k1) = g_ion(i2, j1, k1) + fx2 * fy1 * fz1
            g_ion(i1, j2, k1) = g_ion(i1, j2, k1) + fx1 * fy2 * fz1
            g_ion(i2, j2, k1) = g_ion(i2, j2, k1) + fx2 * fy2 * fz1
            g_ion(i1, j1, k2) = g_ion(i1, j1, k2) + fx1 * fy1 * fz2
            g_ion(i2, j1, k2) = g_ion(i2, j1, k2) + fx2 * fy1 * fz2
            g_ion(i1, j2, k2) = g_ion(i1, j2, k2) + fx1 * fy2 * fz2
            g_ion(i2, j2, k2) = g_ion(i2, j2, k2) + fx2 * fy2 * fz2

            rhoi_loc(i1, j1, k1) = rhoi_loc(i1, j1, k1) + cweight * fx1 * fy1 * fz1
            rhoi_loc(i2, j1, k1) = rhoi_loc(i2, j1, k1) + cweight * fx2 * fy1 * fz1
            rhoi_loc(i1, j2, k1) = rhoi_loc(i1, j2, k1) + cweight * fx1 * fy2 * fz1
            rhoi_loc(i2, j2, k1) = rhoi_loc(i2, j2, k1) + cweight * fx2 * fy2 * fz1
            rhoi_loc(i1, j1, k2) = rhoi_loc(i1, j1, k2) + cweight * fx1 * fy1 * fz2
            rhoi_loc(i2, j1, k2) = rhoi_loc(i2, j1, k2) + cweight * fx2 * fy1 * fz2
            rhoi_loc(i1, j2, k2) = rhoi_loc(i1, j2, k2) + cweight * fx1 * fy2 * fz2
            rhoi_loc(i2, j2, k2) = rhoi_loc(i2, j2, k2) + cweight * fx2 * fy2 * fz2
            ! ion temp
            Ti_loc(i1, j1, k1) = Ti_loc(i1, j1, k1) + tweight * fx1 * fy1 * fz1
            Ti_loc(i2, j1, k1) = Ti_loc(i2, j1, k1) + tweight * fx2 * fy1 * fz1
            Ti_loc(i1, j2, k1) = Ti_loc(i1, j2, k1) + tweight * fx1 * fy2 * fz1
            Ti_loc(i2, j2, k1) = Ti_loc(i2, j2, k1) + tweight * fx2 * fy2 * fz1
            Ti_loc(i1, j1, k2) = Ti_loc(i1, j1, k2) + tweight * fx1 * fy1 * fz2
            Ti_loc(i2, j1, k2) = Ti_loc(i2, j1, k2) + tweight * fx2 * fy1 * fz2
            Ti_loc(i1, j2, k2) = Ti_loc(i1, j2, k2) + tweight * fx1 * fy2 * fz2
            Ti_loc(i2, j2, k2) = Ti_loc(i2, j2, k2) + tweight * fx2 * fy2 * fz2

         end if
         ! Electric field - use ngp, softened 1/r^2 weight
         ex_w(i1, j1, k1) = ex_w(i1, j1, k1) + Ex(i) / fr1**2
         ey_w(i1, j1, k1) = ey_w(i1, j1, k1) + Ey(i) / fr1**2
         ez_w(i1, j1, k1) = ez_w(i1, j1, k1) + Ez(i) / fr1**2
         bz_w(i1, j1, k1) = bz_w(i1, j1, k1) + Bz(i) / fr1**2 ! B-field

         !        phig(i1,j1,k1)=phig(i1,j1,k1) + phi(i)/fr1

      end do

      ! normalise averaged quantities
      nelecs = SUM(g_ele(1:ngx, 1:ngy, 1:ngz))
      nions = SUM(g_ion(1:ngx, 1:ngy, 1:ngz))
      if (debug_level .gt. 1) write (ipefile, *) 'density integrals: ', nelecs, nions
      g_ele = max(gmin, g_ele)
      g_ion = max(gmin, g_ion)
      Te_loc = .511 * Te_loc   ! Temperature in MeV (KE per particle)
      Ti_loc = .511 * mass_ratio * Ti_loc   ! K.E. in MeV

      Ex_loc = Ex_loc + ex_w / (g_ion + g_ele) / navcycle    ! Accumulate normalised fields
      Ey_loc = Ey_loc + ey_w / (g_ion + g_ele) / navcycle
      Ez_loc = Ez_loc + ez_w / (g_ion + g_ele) / navcycle

      Bz_loc = bz_w / (g_ion + g_ele)  ! Instantaneous B

      ! laser fields

   end subroutine sum_fields

   !  =================================
   !
   !    1D gather for field lineout along 0-xl
   !
   !  =================================

   subroutine field_lineout(timestamp)

      use module_physvars
      use module_particle_props
      use module_utils
      use mpi
      implicit none

      integer, intent(in) :: timestamp
      real :: dr
      integer :: ndum, i, ngr, icall, p
      character(30) :: cfile
      character(6) :: cdump
      real*8, dimension(0:ngx + 1) :: ex_g, phi_g
      integer, dimension(ngx + 1) :: pshortl

      icall = timestamp / ivis_fields

      ! get filename suffix from dump counter
      do i = 0, 4
         cdump(6 - i:6 - i) = achar(mod(timestamp / 10**i, 10) + 48)
      end do
      cdump(1:1) = achar(timestamp / 10**5 + 48)

      if (my_rank .eq. 0) then
         write (*, '(//a/a,f10.2)') 'Field lineout:', &
            'writing out Ex and Phi on grid 0-', xl
      end if

      ngr = ngx  ! use box resolution in x for radial distn
      dr = xl / ngr

      ! Field and potential lineouts - set up dummy particles & find forces directly from tree
      ! Dummies set up at end of particle arrays on root to ensure unique labelling
      if (my_rank .eq. 0) then
         do i = 1, ngr + 1

            p = np_local + i   !index
            pshortl(i) = p   !index
            x(p) = (i - 1) * dr   ! start at x=0
            y(p) = plasma_centre(2)
            z(p) = plasma_centre(3)
         end do
         ! Get interaction lists
         !     write (*,*) 'Doing lists for dummy particles'
         !     write (*,'((i8,f12.3))') (pshortl(i),x(pshortl(i)),i=1,ngr+1)
         ndum = ngr + 1
      else
         ndum = 0
      end if

      ! TODO need replacement for walk, sum_force for ghost particles

      ! all PEs must call walk
      ! mac set to BH, theta to 0.5
      !  call tree_walk(pshortl(1:ndum),ndum,1,0.5,eps,itime,0,ttrav,tfetch)

      if (my_rank .eq. 0) then
         ! Fields
         do i = 1, ngr + 1
            p = pshortl(i)
            !       call sum_force(p, nterm(i), nodelist( 1:nterm(i),i), eps, &
            !            ex_g(i-1), ey_g(i-1), ez_g(i-1), phi_g(i-1), w_g(i-1))
         end do

      end if

      ! Write out to file

      if (my_rank .eq. 0) then
         call create_directory("fields")
         cfile = "fields/lineout."//cdump
         open (60, file=cfile)
         write (60, '(3(a12))') '!   x   ', 'Ex', 'Phi'
         write (60, '((3(1pe12.4)))') (i * dr, ex_g(i), phi_g(i), i=0, ngr)
         close (60)

      end if
      icall = icall + 1

   end subroutine field_lineout

   !  =================================
   !
   !    1D gather for spherically symmetric fields
   !
   !  =================================

   subroutine sum_radial(timestamp)

      use module_physvars
      use module_particle_props
      use module_utils
      use mpi
      implicit none

      integer, intent(in) :: timestamp
      real :: rdr, dr, vweight, cweight
      real :: tweight, erweight
      real :: fr1, fr2, ra, gamma, xt, yt, zt, rt

      integer :: ndum, i, i1, i2, nelecs, nions, ngr, icall, ierr, p
      character(30) :: cfile
      character(6) :: cdump
      real, dimension(0:ngx + 1) :: ve_loc, vi_loc, Er_loc, ni_loc, ne_loc, ge_loc, gi_loc
      real, dimension(0:ngx + 1) :: ve_glob, vi_glob, Er_glob, ni_glob, ne_glob, ge_glob, gi_glob
      real, dimension(0:ngx + 1) :: volume
      real*8, dimension(0:ngx + 1) :: ex_g, phi_g
      integer, dimension(ngx + 1) :: pshortl
      real :: gmin = 1.e-3

      icall = timestamp / ivis_fields

      ! get filename suffix from dump counter
      do i = 0, 4
         cdump(6 - i:6 - i) = achar(mod(timestamp / 10**i, 10) + 48)
      end do
      cdump(1:1) = achar(timestamp / 10**5 + 48)

      if (my_rank .eq. 0) then
         write (*, '(//a/a,f10.2)') 'Radial field dump:', &
            'writing out densities, fields on grid 0-', xl
      end if

      ngr = ngx  ! use box resolution in x for radial distn
      dr = xl / ngr

      rdr = 1./dr

      ve_loc = 0.
      vi_loc = 0.
      ge_loc = 0.
      gi_loc = 0.
      ne_loc = 0.
      ni_loc = 0.
      Er_loc = 0.

      !  field box limits assumed to be: (0-xl)
      !  Any particle outside gets put in ghost cells ngr+1

      !      write(15,'(//a,3f12.3)') 'cw,dx,dy',cweight,dx,dy

      !  Accumulate local radial moments

      do i = 1, np_local

         ! particle position relative to plasma centre
         xt = x(i) - plasma_centre(1)
         yt = y(i) - plasma_centre(2)
         zt = z(i) - plasma_centre(3)

         ! radius
         rt = sqrt(xt**2 + yt**2 + zt**2)
         ra = rt * rdr

         !  indices - include zero for r=0
         i1 = ra + .5
         i2 = i1 + 1

         i1 = min(i1, ngr + 1)
         i2 = min(i2, ngr + 1)

         !  linear weighting
         !     fr2=ra-i1  ! Prevent overflow/negative weighting for particles outside box
         !     fr1=1.-fr1
         ! NGP
         fr2 = 0.
         fr1 = 1.
         !  gather charge at nearest grid points
         gamma = sqrt(1.0 + ux(i)**2 + uy(i)**2 + uz(i)**2)
         cweight = abs(q(i))    !  charge weighting factor - densities computed later
         vweight = sqrt(ux(i)**2 + uy(i)**2 + uz(i)**2)
         erweight = sqrt(ex(i)**2 + ey(i)**2 + ez(i)**2)

         tweight = (gamma - 1.)  ! K.E. of particle in keV

         if (q(i) .lt. 0) then
            ve_loc(i1) = ve_loc(i1) + fr1 * vweight
            ve_loc(i2) = ve_loc(i2) + fr2 * vweight

            ne_loc(i1) = ne_loc(i1) + fr1 * cweight
            ne_loc(i2) = ne_loc(i2) + fr2 * cweight

            ge_loc(i1) = ge_loc(i1) + fr1  ! weighted # electrons at r
            ge_loc(i2) = ge_loc(i2) + fr2

         else
            vi_loc(i1) = vi_loc(i1) + fr1 * vweight  ! TODO: This wont give multi-valued vel
            vi_loc(i2) = vi_loc(i2) + fr2 * vweight  ! beyond shock front - need phase space!

            ni_loc(i1) = ni_loc(i1) + fr1 * cweight
            ni_loc(i2) = ni_loc(i2) + fr2 * cweight

            gi_loc(i1) = gi_loc(i1) + fr1  ! weighted # ions at r
            gi_loc(i2) = gi_loc(i2) + fr2

         end if

      end do

      volume(1:ngr + 1) = (/(4./3.*pi * ((i * dr + .5 * dr)**3 - (i * dr - .5 * dr)**3), i=1, ngr + 1)/)
      !  ni_loc(1) = ni_loc(1) + ni_loc(0)  ! Fold charge at r=0 onto r=dr
      !  ne_loc(1) = ne_loc(1) + ne_loc(0)
      volume(0) = 4./3.*pi * (dr / 2.)**3  ! 1/2-vol weight for r=0

      !  Renormalise charge densities with volume-weighting
      ni_loc = ni_loc / volume
      ne_loc = ne_loc / volume

      ! Gather partial sums together to get global moments - arrays run from (0:ngr)

      call MPI_ALLREDUCE(ge_loc, ge_glob, ngr + 1, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(gi_loc, gi_glob, ngr + 1, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(ve_loc, ve_glob, ngr + 1, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(vi_loc, vi_glob, ngr + 1, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(ne_loc, ne_glob, ngr + 1, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(ni_loc, ni_glob, ngr + 1, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)

      ge_glob = max(gmin, ge_glob)
      gi_glob = max(gmin, gi_glob)

      ! normalise averaged quantities
      nelecs = SUM(ge_glob(1:ngr))
      nions = SUM(gi_glob(1:ngr))
      if (my_rank .eq. 0) then
         write (*, *) "Charge check:", nelecs, nions
      end if
      ve_glob = ve_glob / ge_glob
      vi_glob = vi_glob / gi_glob

      ! Radial fields - set up dummy particles & find forces directly from tree
      ! Dummies set up at end of particle arrays on root to ensure unique labelling
      if (my_rank .eq. 0) then
      do i = 1, ngr + 1

         p = np_local + i   !index
         pshortl(i) = p   !index
         x(p) = (i - 1) * dr + plasma_centre(1)   ! radius - include r=0
         y(p) = plasma_centre(2)
         z(p) = plasma_centre(3)
      end do
      ! Get interaction lists
      !     write (*,*) 'Doing lists for dummy particles'
      !     write (*,'((i8,f12.3))') (pshortl(i),x(pshortl(i)),i=1,ngr+1)
      ndum = ngr + 1
      else
      ndum = 0
      end if

      ! TODO need to replace these with new kernel routine avoiding tree build

      ! all PEs must call walk
      ! mac set to BH, theta to 0.5
      !   call tree_walk(pshortl(1:ndum),ndum,1,0.5,eps,itime,0,ttrav,tfetch)

      if (my_rank .eq. 0) then
         ! Fields
         do i = 1, ngr + 1
            p = pshortl(i)
            !        call sum_force(p, nterm(i), nodelist( 1:nterm(i),i), eps, &
            !            ex_g(i-1), ey_g(i-1), ez_g(i-1), phi_g(i-1), w_g(i-1))
         end do

         Er_glob(0:ngr) = Ex_g(0:ngr)    !  Radial field (leave off norm constant)
      end if

      ! Write out to file

      if (my_rank .eq. 0) then
         write (*, *) 'number density integrals: ', nelecs, nions
         call create_directory("fields")
         cfile = "fields/radial."//cdump
         open (60, file=trim(cfile))
         write (60, '(9(a12))') '!   r      ', ' r/r0   ', 'ne    ', 'ni   ', 'rhoe   ', 'rhoi   ', 've   ', 'vi   ', 'er'
         write (60, '((10(1pe12.4)))') &
            (i * dr, i * dr, ge_glob(i) / max(1, ne), gi_glob(i) / max(1, ni), max(ne_glob(i), 1.e-10), max(ni_glob(i), 1.e-10), &
             ve_glob(i), vi_glob(i), max(er_glob(i), 1.e-10), phi_g(i), i=0, ngr)
         close (60)

      end if
      icall = icall + 1

   end subroutine sum_radial

   !  =================================
   !
   !    3D gather for potential and E-field
   !
   !  =================================

end module module_field_grid
