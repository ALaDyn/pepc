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

module init_particles

   use module_pepc_kinds
   use manipulate_particles
   implicit none

contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !>
   !>   Initialize particles with different setup (choose via ispecial)
   !>
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine special_start()

      use physvars
      use files
      use mpi
      implicit none

      integer :: ierr
      integer(kind_particle) :: j, k, ind, ind0, i, m, l
      real :: par_rand_res
      real(kind_physics) :: part_2d, rc, xi1, xi2, xi, part_3d, eta1, eta2, eta, v(3), xt, yt, zt
      real(kind_physics), dimension(3, 3) :: D1, D2, D3, D4   !< rotation matrices for ring setup
      real(kind_physics), allocatable :: xp(:), yp(:), zp(:), volp(:), wxp(:), wyp(:), wzp(:)  !< helper arrays for ring setups

      ! weird helper variables for sphere setup (ask Holger Dachsel)
      real(kind_physics) ::  a, b, c, cth, sth, cphi, sphi, s, rr, rr1, rr2, expo, stheta
      real(kind_physics), parameter :: zero = 0.d0, one = 1.d0, mone = -one, two = 2.d0, five = 5.d0, nine = 9.d0
      real(kind_physics), parameter :: kappa = 2.24182**2 / 4.
      real(kind_physics), parameter :: r_core = 0.35d0
      real(kind_physics) :: sumvol, rrmax, rrmin

      !     interface
      !         subroutine par_rand(res, iseed)
      !            real, intent(inout) :: res
      !            integer, intent(in), optional :: iseed
      !         end subroutine
      !     end interface

      ! Set up particle data
      config: select case (ispecial)
      case (1)                               ! Vortex ring setup, side-by-side

         allocate (xp(ns), yp(ns), zp(ns), volp(ns), wxp(ns), wyp(ns), wzp(ns))

         j = 0

         do k = 1, nc
            part_2d = 2 * pi / (8 * k)
            rc = (1 + 12 * k**2) / (6 * k) * rl

            do l = 1, 8 * k
               j = j + 1
               xi1 = part_2d * (l - 1)
               xi2 = part_2d * l
               xi = (xi2 - xi1) / 2 + xi1
               xp(j) = rc * cos(xi)
               yp(j) = rc * sin(xi)
               zp(j) = 0
               volp(j) = (2 * pi**2 * (r_torus + (2 * k + 1) * rl) * ((2 * k + 1) * rl)**2 - 2 * pi**2 * (r_torus + (2 * k - 1) * rl) * ((2 * k - 1) * rl)**2) / (8 * k * Nphi)
               wxp(j) = 0.
               wyp(j) = 0.
               wzp(j) = g * exp(-(rc / rmax)**2)
            end do
         end do

         xp(ns) = 0.
         yp(ns) = 0.
         zp(ns) = 0.
         wxp(ns) = 0.
         wyp(ns) = 0.
         wzp(ns) = g
         volp(ns) = 2 * pi**2 * (r_torus + rl) * rl**2 / Nphi

         j = 0
         ind0 = 0
         ind = 0
         part_3d = 2 * pi / Nphi
         do m = 1, Nphi
            eta1 = part_3d * (m - 1)
            eta2 = part_3d * m
            eta = (eta2 - eta1) / 2 + eta1
            v(1) = cos(eta2)
            v(2) = sin(eta2)
            v(3) = 0
            D1 = reshape((/-1.0 + 2 * v(1)**2, 2 * v(2) * v(1), 0.0D0, 2 * v(1) * v(2), -1.0 + 2 * v(2)**2, 0.0D0, 0.0D0, 0.0D0, -1.0D0/), (/3, 3/))
            v(1) = cos(eta)
            v(2) = sin(eta)
            v(3) = 0
            D2 = reshape((/v(1)**2, v(2) * v(1), -v(2), v(1) * v(2), v(2)**2, v(1), v(2), -v(1), 0.0D0/), (/3, 3/))
            D3 = reshape((/-1.0 + 2 * v(1)**2, 2 * v(2) * v(1), 0.0D0, 2 * v(1) * v(2), -1.0 + 2 * v(2)**2, 0.0D0, 0.0D0, 0.0D0, -1.0D0/), (/3, 3/))
            D4 = matmul(D1, D3)
            do i = 1, Ns
               if (m .eq. 1) then
                  v(1) = xp(i) + (r_torus + rmax) * cos(eta)
                  v(2) = yp(i) + (r_torus + rmax) * sin(eta)
                  v(3) = zp(i)
                  xp(i) = dot_product(v, D2(1:3, 1))
                  yp(i) = dot_product(v, D2(1:3, 2))
                  zp(i) = dot_product(v, D2(1:3, 3))
                  v(1) = wxp(i)
                  v(2) = wyp(i)
                  v(3) = wzp(i)
                  wxp(i) = dot_product(v, D2(1:3, 1))
                  wyp(i) = dot_product(v, D2(1:3, 2))
                  wzp(i) = dot_product(v, D2(1:3, 3))
               else
                  v(1) = xp(i)
                  v(2) = yp(i)
                  v(3) = zp(i)
                  xp(i) = dot_product(v, D4(1:3, 1))
                  yp(i) = dot_product(v, D4(1:3, 2))
                  zp(i) = dot_product(v, D4(1:3, 3))
                  v(1) = wxp(i)
                  v(2) = wyp(i)
                  v(3) = wzp(i)
                  wxp(i) = dot_product(v, D4(1:3, 1))
                  wyp(i) = dot_product(v, D4(1:3, 2))
                  wzp(i) = dot_product(v, D4(1:3, 3))
               end if
               ind0 = ind0 + 1
               if (mod(ind0 - 1, n_cpu) .eq. my_rank) then
                  ind = ind + 1
                  if (ind .gt. np) then
                     write (*, *) 'something is wrong here: to many particles in init', my_rank, ind, np, n
                     call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
                  end if
                  vortex_particles(ind)%x(1) = xp(i) - (torus_offset(1) - (rmax - (1 + 12 * nc**2) / (6 * nc) * rl)) / 2.0
                  vortex_particles(ind)%x(2) = yp(i)
                  vortex_particles(ind)%x(3) = zp(i)
                  vortex_particles(ind)%data%alpha(1) = wxp(i) * volp(i)
                  vortex_particles(ind)%data%alpha(2) = wyp(i) * volp(i)
                  vortex_particles(ind)%data%alpha(3) = wzp(i) * volp(i)
                  ind = ind + 1
                  vortex_particles(ind)%x(1) = xp(i) + (torus_offset(1) - (rmax - (1 + 12 * nc**2) / (6 * nc) * rl)) / 2.0
                  vortex_particles(ind)%x(2) = yp(i)
                  vortex_particles(ind)%x(3) = zp(i)
                  vortex_particles(ind)%data%alpha(1) = wxp(i) * volp(i)
                  vortex_particles(ind)%data%alpha(2) = wyp(i) * volp(i)
                  vortex_particles(ind)%data%alpha(3) = wzp(i) * volp(i)
                  vortex_particles(ind)%data%vol = volp(i)
               end if
            end do
         end do
         np = ind
         deallocate (xp, yp, zp, volp, wxp, wyp, wzp)

      case (2) ! Vortex ring setup, offset collision

         allocate (xp(ns), yp(ns), zp(ns), volp(ns), wxp(ns), wyp(ns), wzp(ns))

         j = 0

         do k = 1, nc
            part_2d = 2 * pi / (8 * k)
            rc = (1 + 12 * k**2) / (6 * k) * rl

            do l = 1, 8 * k
               j = j + 1
               xi1 = part_2d * (l - 1)
               xi2 = part_2d * l
               xi = (xi2 - xi1) / 2 + xi1
               xp(j) = rc * cos(xi)
               yp(j) = rc * sin(xi)
               zp(j) = 0
               volp(j) = (2 * pi**2 * (r_torus + (2 * k + 1) * rl) * ((2 * k + 1) * rl)**2 - 2 * pi**2 * (r_torus + (2 * k - 1) * rl) * ((2 * k - 1) * rl)**2) / (8 * k * Nphi)
               wxp(j) = 0.
               wyp(j) = 0.
               wzp(j) = g * exp(-(rc / rmax)**2)
            end do
         end do

         xp(ns) = 0.
         yp(ns) = 0.
         zp(ns) = 0.
         wxp(ns) = 0.
         wyp(ns) = 0.
         wzp(ns) = g
         volp(ns) = 2 * pi**2 * (r_torus + rl) * rl**2 / Nphi

         j = 0
         ind0 = 0
         ind = 0
         part_3d = 2 * pi / Nphi
         do m = 1, Nphi
            eta1 = part_3d * (m - 1)
            eta2 = part_3d * m
            eta = (eta2 - eta1) / 2 + eta1
            v(1) = cos(eta2)
            v(2) = sin(eta2)
            v(3) = 0
            D1 = reshape((/-1.0 + 2 * v(1)**2, 2 * v(2) * v(1), 0.0D0, 2 * v(1) * v(2), -1.0 + 2 * v(2)**2, 0.0D0, 0.0D0, 0.0D0, -1.0D0/), (/3, 3/))
            v(1) = cos(eta)
            v(2) = sin(eta)
            v(3) = 0
            D2 = reshape((/v(1)**2, v(2) * v(1), -v(2), v(1) * v(2), v(2)**2, v(1), v(2), -v(1), 0.0D0/), (/3, 3/))
            D3 = reshape((/-1.0 + 2 * v(1)**2, 2 * v(2) * v(1), 0.0D0, 2 * v(1) * v(2), -1.0 + 2 * v(2)**2, 0.0D0, 0.0D0, 0.0D0, -1.0D0/), (/3, 3/))
            D4 = matmul(D1, D3)
            do i = 1, Ns
               if (m .eq. 1) then
                  v(1) = xp(i) + (r_torus + rmax) * cos(eta)
                  v(2) = yp(i) + (r_torus + rmax) * sin(eta)
                  v(3) = zp(i)
                  xp(i) = dot_product(v, D2(1:3, 1))
                  yp(i) = dot_product(v, D2(1:3, 2))
                  zp(i) = dot_product(v, D2(1:3, 3))
                  v(1) = wxp(i)
                  v(2) = wyp(i)
                  v(3) = wzp(i)
                  wxp(i) = dot_product(v, D2(1:3, 1))
                  wyp(i) = dot_product(v, D2(1:3, 2))
                  wzp(i) = dot_product(v, D2(1:3, 3))
               else
                  v(1) = xp(i)
                  v(2) = yp(i)
                  v(3) = zp(i)
                  xp(i) = dot_product(v, D4(1:3, 1))
                  yp(i) = dot_product(v, D4(1:3, 2))
                  zp(i) = dot_product(v, D4(1:3, 3))
                  v(1) = wxp(i)
                  v(2) = wyp(i)
                  v(3) = wzp(i)
                  wxp(i) = dot_product(v, D4(1:3, 1))
                  wyp(i) = dot_product(v, D4(1:3, 2))
                  wzp(i) = dot_product(v, D4(1:3, 3))
               end if
               ind0 = ind0 + 1
               if (mod(ind0 - 1, n_cpu) .eq. my_rank) then
                  ind = ind + 1
                  if (ind .gt. np) then
                     write (*, *) 'something is wrong here: to many particles in init of first ring', my_rank, ind, np, n
                     call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
                  end if
                  vortex_particles(ind)%x(1) = xp(i) - (torus_offset(1) - (rmax - (1 + 12 * nc**2) / (6 * nc) * rl)) / 2.0
                  vortex_particles(ind)%x(2) = yp(i) - (torus_offset(2) - (rmax - (1 + 12 * nc**2) / (6 * nc) * rl)) / 2.0
                  vortex_particles(ind)%x(3) = zp(i) - (torus_offset(3) - (rmax - (1 + 12 * nc**2) / (6 * nc) * rl)) / 2.0
                  vortex_particles(ind)%data%alpha(1) = wxp(i) * volp(i)
                  vortex_particles(ind)%data%alpha(2) = wyp(i) * volp(i)
                  vortex_particles(ind)%data%alpha(3) = wzp(i) * volp(i)
                  vortex_particles(ind)%data%vol = volp(i)
               end if
            end do
         end do

         j = 0

         do k = 1, nc
            part_2d = 2 * pi / (8 * k)
            rc = (1 + 12 * k**2) / (6 * k) * rl

            do l = 1, 8 * k
               j = j + 1
               xi1 = part_2d * (l - 1)
               xi2 = part_2d * l
               xi = (xi2 - xi1) / 2 + xi1
               xp(j) = rc * cos(xi)
               yp(j) = rc * sin(xi)
               zp(j) = 0
               volp(j) = (2 * pi**2 * (r_torus + (2 * k + 1) * rl) * ((2 * k + 1) * rl)**2 - 2 * pi**2 * (r_torus + (2 * k - 1) * rl) * ((2 * k - 1) * rl)**2) / (8 * k * Nphi)
               wxp(j) = 0.
               wyp(j) = 0.
               wzp(j) = G * exp(-(rc / rmax)**2)
            end do
         end do

         xp(ns) = 0.
         yp(ns) = 0.
         zp(ns) = 0.
         wxp(ns) = 0.
         wyp(ns) = 0.
         wzp(ns) = g
         volp(Ns) = 2 * pi**2 * (r_torus + rl) * rl**2 / Nphi

         j = 0
         ind0 = 0
         part_3d = 2 * pi / Nphi
         do m = 1, Nphi
            eta1 = part_3d * (m - 1)
            eta2 = part_3d * m
            eta = (eta2 - eta1) / 2 + eta1
            v(1) = cos(eta2)
            v(2) = sin(eta2)
            v(3) = 0
            D1 = reshape((/-1.0 + 2 * v(1)**2, 2 * v(2) * v(1), 0.0D0, 2 * v(1) * v(2), -1.0 + 2 * v(2)**2, 0.0D0, 0.0D0, 0.0D0, -1.0D0/), (/3, 3/))
            v(1) = cos(eta)
            v(2) = sin(eta)
            v(3) = 0
            D2 = reshape((/v(1)**2, v(2) * v(1), -v(2), v(1) * v(2), v(2)**2, v(1), v(2), -v(1), 0.0D0/), (/3, 3/))
            D3 = reshape((/-1.0 + 2 * v(1)**2, 2 * v(2) * v(1), 0.0D0, 2 * v(1) * v(2), -1.0 + 2 * v(2)**2, 0.0D0, 0.0D0, 0.0D0, -1.0D0/), (/3, 3/))
            D4 = matmul(D1, D3)
            do i = 1, Ns
               if (m .eq. 1) then
                  v(1) = xp(i) + (r_torus + rmax) * cos(eta)
                  v(2) = yp(i) + (r_torus + rmax) * sin(eta)
                  v(3) = zp(i)
                  xp(i) = dot_product(v, D2(1:3, 1))
                  yp(i) = dot_product(v, D2(1:3, 2))
                  zp(i) = dot_product(v, D2(1:3, 3))
                  v(1) = wxp(i)
                  v(2) = wyp(i)
                  v(3) = -wzp(i)
                  wxp(i) = dot_product(v, D2(1:3, 1))
                  wyp(i) = dot_product(v, D2(1:3, 2))
                  wzp(i) = dot_product(v, D2(1:3, 3))
               else
                  v(1) = xp(i)
                  v(2) = yp(i)
                  v(3) = zp(i)
                  xp(i) = dot_product(v, D4(1:3, 1))
                  yp(i) = dot_product(v, D4(1:3, 2))
                  zp(i) = dot_product(v, D4(1:3, 3))
                  v(1) = wxp(i)
                  v(2) = wyp(i)
                  v(3) = -wzp(i)
                  wxp(i) = dot_product(v, D4(1:3, 1))
                  wyp(i) = dot_product(v, D4(1:3, 2))
                  wzp(i) = dot_product(v, D4(1:3, 3))
               end if
               ind0 = ind0 + 1
               if (mod(ind0 - 1, n_cpu) .eq. my_rank) then
                  ind = ind + 1
                  if (ind .gt. np) then
                     write (*, *) 'something is wrong here: to many particles in init of second ring', my_rank, ind, np, n
                     call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
                  end if
                  vortex_particles(ind)%x(1) = xp(i) + (torus_offset(1) - (rmax - (1 + 12 * nc**2) / (6 * nc) * rl)) / 2.0
                  vortex_particles(ind)%x(2) = yp(i) + (torus_offset(2) - (rmax - (1 + 12 * nc**2) / (6 * nc) * rl)) / 2.0
                  vortex_particles(ind)%x(3) = zp(i) + (torus_offset(3) - (rmax - (1 + 12 * nc**2) / (6 * nc) * rl)) / 2.0
                  vortex_particles(ind)%data%alpha(1) = wxp(i) * volp(i)
                  vortex_particles(ind)%data%alpha(2) = wyp(i) * volp(i)
                  vortex_particles(ind)%data%alpha(3) = wzp(i) * volp(i)
                  vortex_particles(ind)%data%vol = volp(i)
               end if
            end do
         end do
         np = ind

         deallocate (xp, yp, zp, volp, wxp, wyp, wzp)

      case (3)  ! Sphere setup

         a = nine / five * sqrt(dble(n - 1) / dble(n))
         j = 0
         do i = 1, n

            if (i .eq. 1) then
               cth = mone
               sth = zero
               cphi = one
               sphi = zero
               b = cphi
               c = sphi
            elseif (i .eq. n) then
               cth = one
               sth = zero
               cphi = one
               sphi = zero
            else
               cth = dble(2 * i - n - 1) / dble(n - 1)
               sth = two * (sqrt(dble(i - 1) / dble(n - 1)) * sqrt(dble(n - i) / dble(n - 1)))
               s = a * sqrt(dble(n - 1) / (dble(i - 1) * dble(n - i)))
               cphi = b * cos(s) - c * sin(s)
               sphi = c * cos(s) + b * sin(s)
               b = cphi
               c = sphi
            end if
            if (mod(i + my_rank, n_cpu) .eq. 0) then
               if (j .gt. np - 1) then
                  write (*, *) 'something is wrong here: to many particles in init', my_rank, j, np, i
                  call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
               end if
               j = j + 1
               vortex_particles(j)%x(1) = sth * cphi
               vortex_particles(j)%x(2) = sth * sphi
               vortex_particles(j)%x(3) = cth
               vortex_particles(j)%data%alpha(1) = 0.3D01 / (0.8D01 * pi) * sth * sphi * h**2
               vortex_particles(j)%data%alpha(2) = 0.3D01 / (0.8D01 * pi) * sth * (-cphi) * h**2
               vortex_particles(j)%data%alpha(3) = 0.
            end if
         end do
         np = j

      case (4) ! Vortex wakes

         allocate (xp(ns), yp(ns), zp(ns), volp(ns), wxp(ns), wyp(ns), wzp(ns))

         j = 0

         do k = 1, nc
            part_2d = 2 * pi / (8 * k)
            rc = (1 + 12 * k**2) / (6 * k) * rl

            do l = 1, 8 * k
               j = j + 1
               xi1 = part_2d * (l - 1)
               xi2 = part_2d * l
               xi = (xi2 - xi1) / 2 + xi1
               xp(j) = rc * cos(xi)
               yp(j) = rc * sin(xi)
               zp(j) = 0
               volp(j) = (2 * pi**2 * (r_torus + (2 * k + 1) * rl) * ((2 * k + 1) * rl)**2 - 2 * pi**2 * (r_torus + (2 * k - 1) * rl) * ((2 * k - 1) * rl)**2) / (8 * k * Nphi)
               wxp(j) = 0.
               wyp(j) = 0.
               wzp(j) = g * exp(-(rc / rmax)**2)
            end do
         end do

         xp(ns) = 0.
         yp(ns) = 0.
         zp(ns) = 0.
         wxp(ns) = 0.
         wyp(ns) = 0.
         wzp(ns) = g
         volp(ns) = 2 * pi**2 * (r_torus + rl) * rl**2 / Nphi

         part_3d = 4 * pi / nphi
         ind0 = 0
         ind = 0
         zp(1:ns) = zp(1:ns) - 2 * pi - part_3d
         do k = 1, nphi

            zp(1:ns) = zp(1:ns) + part_3d

            do i = 1, ns
               ind0 = ind0 + 1
               if (mod(ind0 - 1, n_cpu) .eq. my_rank) then
                  ind = ind + 1
                  if (ind .gt. np) then
                     write (*, *) 'something is wrong here: to many particles in init', my_rank, ind, np, n
                     call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
                  end if
                  vortex_particles(ind)%x(1) = xp(i)
                  vortex_particles(ind)%x(2) = yp(i) - torus_offset(2)
                  vortex_particles(ind)%x(3) = zp(i)
                  vortex_particles(ind)%data%alpha(1) = 0.
                  vortex_particles(ind)%data%alpha(2) = 0.
                  vortex_particles(ind)%data%alpha(3) = -(exp(-(zp(i) - pi)**2) + exp(-(zp(i) + pi)**2)) * wzp(i) * volp(i)
                  ind = ind + 1
                  vortex_particles(ind)%x(1) = xp(i)
                  vortex_particles(ind)%x(2) = yp(i) + torus_offset(2)
                  vortex_particles(ind)%x(3) = zp(i)
                  vortex_particles(ind)%data%alpha(1) = 0.
                  vortex_particles(ind)%data%alpha(2) = 0.
                  vortex_particles(ind)%data%alpha(3) = +(exp(-(zp(i) - pi)**2) + exp(-zp(i)**2) + exp(-(zp(i) + pi)**2)) * wzp(i) * volp(i)
               end if
            end do

         end do

         np = ind

      case (5)  ! Different wakes

         ind0 = 0
         ind = 0
         do i = 1, ceiling(1.0 / h)
            do j = 1, ceiling(2 * pi / h)
               do k = 1, nc
                  ind0 = ind0 + 1
                  if (mod(ind0 + my_rank, n_cpu) .eq. 0) then
                     ind = ind + 1
                     if (ind .gt. np - 1) then
                        write (*, *) 'something is wrong here: to many particles in init', my_rank, ind, np, n, n_cpu
                        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
                     end if
                     xt = (i - 1) * h !+ h/2
                     yt = -pi + (j - 1) * h !+ h/2
                     zt = -pi + (k - 1) * h

                     vortex_particles(ind)%x(1) = xt + torus_offset(1)
                     vortex_particles(ind)%x(2) = yt + torus_offset(2)
                     vortex_particles(ind)%x(3) = zt + torus_offset(3)
                     vortex_particles(ind)%data%alpha(1) = 0.
                     vortex_particles(ind)%data%alpha(2) = 0.
                     vortex_particles(ind)%data%alpha(3) = -g / 2 * (1 - tanh(yt)**2) * h**3 !* (exp(-zt**2/2)+exp(-(zt-pi/2)**2/2)+exp(-(zt+pi/2)**2/2))
                     ind = ind + 1
                     vortex_particles(ind)%x(1) = xt - torus_offset(1)
                     vortex_particles(ind)%x(2) = yt - torus_offset(2)
                     vortex_particles(ind)%x(3) = zt - torus_offset(3)
                     vortex_particles(ind)%data%alpha(1) = 0.
                     vortex_particles(ind)%data%alpha(2) = 0.
                     vortex_particles(ind)%data%alpha(3) = +g / 2 * (1 - tanh(yt)**2) * h**3 !* (exp(-zt**2/2)+exp(-(zt-pi/2)**2/2)+exp(-(zt+pi/2)**2/2))
                  end if
               end do
            end do
         end do
         np = ind

      case (6) ! Single Vortex ring setup

         allocate (xp(ns), yp(ns), zp(ns), volp(ns), wxp(ns), wyp(ns), wzp(ns))

         j = 0
         sumvol = 0.d0

         do k = 1, nc
            part_2d = 2 * pi / (8 * k)
            rc = real(1 + 12 * k**2, kind(rc)) / real(6 * k, kind(rc)) * rl
            rrmax = rl * (2 * k + 1)
            rrmin = rl * (2 * k - 1)
            do l = 1, 8 * k
               j = j + 1
               xi1 = part_2d * (l - 1)
               xi2 = part_2d * l
               xi = (xi1 + xi2) * 5.d-1
               xp(j) = rc * cos(xi)
               yp(j) = rc * sin(xi)
               zp(j) = 0

               rr1 = xp(j) + r_torus
               rr2 = yp(j)
               rr = sqrt(rr1 * rr1 + rr2 * rr2)

               !volp(j) = (xi2 - xi1) * 0.5d0 * (rrmax**2 - rrmin**2) * 2.d0 * pi * rr1 / real(Nphi, kind(volp))
               volp(j) = (xi2 - xi1) * (rrmax**2 - rrmin**2) * pi * rr1 / real(Nphi, kind(volp))
               wxp(j) = 0.
               wyp(j) = 0.
               stheta = rr1 / rr
               expo = kappa * (r_torus * r_torus + rr * rr - 2.d0 * r_torus * rr * stheta)
               wzp(j) = kappa / pi * G / r_core / r_core * exp(-expo / r_core / r_core)
               sumvol = sumvol + volp(j)
            end do
         end do

         xp(ns) = 0.
         yp(ns) = 0.
         zp(ns) = 0.
         wxp(ns) = 0.
         wyp(ns) = 0.
         wzp(ns) = kappa / pi * G / r_core / r_core
         volp(ns) = pi * rl**2 * 2.d0 * pi * r_torus / real(Nphi, kind(volp))
         sumvol = sumvol + volp(ns)

         if (my_rank .eq. 0) write (*, *) 'Total discretized volume ', sumvol * real(Nphi), ' - theoretical', 2.d0 * pi * r_torus * pi * rmax * rmax

         j = 0
         ind0 = 0
         ind = 0
         part_3d = 2 * pi / Nphi
         do m = 1, Nphi
            eta1 = part_3d * (m - 1)
            eta2 = part_3d * m
            eta = (eta1 + eta2) * 5.d-1
            v(1) = cos(eta2)
            v(2) = sin(eta2)
            v(3) = 0
            D1 = reshape((/-1.0 + 2 * v(1)**2, 2 * v(2) * v(1), 0.0D0, 2 * v(1) * v(2), -1.0 + 2 * v(2)**2, 0.0D0, 0.0D0, 0.0D0, -1.0D0/), (/3, 3/))
            v(1) = cos(eta)
            v(2) = sin(eta)
            v(3) = 0
            D2 = reshape((/v(1)**2, v(2) * v(1), -v(2), v(1) * v(2), v(2)**2, v(1), v(2), -v(1), 0.0D0/), (/3, 3/))
            D3 = reshape((/-1.0 + 2 * v(1)**2, 2 * v(2) * v(1), 0.0D0, 2 * v(1) * v(2), -1.0 + 2 * v(2)**2, 0.0D0, 0.0D0, 0.0D0, -1.0D0/), (/3, 3/))
            D4 = matmul(D1, D3)
            do i = 1, Ns
               if (m .eq. 1) then
!                       v(1) = xp(i) + (r_torus+rmax)*cos(eta)
!                       v(2) = yp(i) + (r_torus+rmax)*sin(eta)
                  v(1) = xp(i) + r_torus * cos(eta)
                  v(2) = yp(i) + r_torus * sin(eta)
                  v(3) = zp(i)
                  xp(i) = dot_product(v, D2(1:3, 1))
                  yp(i) = dot_product(v, D2(1:3, 2))
                  zp(i) = dot_product(v, D2(1:3, 3))
                  v(1) = wxp(i)
                  v(2) = wyp(i)
                  v(3) = wzp(i)
                  wxp(i) = dot_product(v, D2(1:3, 1))
                  wyp(i) = dot_product(v, D2(1:3, 2))
                  wzp(i) = dot_product(v, D2(1:3, 3))
               else
                  v(1) = xp(i)
                  v(2) = yp(i)
                  v(3) = zp(i)
                  xp(i) = dot_product(v, D4(1:3, 1))
                  yp(i) = dot_product(v, D4(1:3, 2))
                  zp(i) = dot_product(v, D4(1:3, 3))
                  v(1) = wxp(i)
                  v(2) = wyp(i)
                  v(3) = wzp(i)
                  wxp(i) = dot_product(v, D4(1:3, 1))
                  wyp(i) = dot_product(v, D4(1:3, 2))
                  wzp(i) = dot_product(v, D4(1:3, 3))
               end if
               ind0 = ind0 + 1
               if (mod(ind0 - 1, n_cpu) .eq. my_rank) then
                  ind = ind + 1
                  if (ind .gt. np) then
                     write (*, *) 'something is wrong here: to many particles in init of first ring', my_rank, ind, np, n
                     call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
                  end if
                  vortex_particles(ind)%x(1) = xp(i) + (rmax - (1 + 12 * nc**2) / (6 * nc) * rl) / 2.0
                  vortex_particles(ind)%x(2) = yp(i) + (rmax - (1 + 12 * nc**2) / (6 * nc) * rl) / 2.0
                  vortex_particles(ind)%x(3) = zp(i) + (rmax - (1 + 12 * nc**2) / (6 * nc) * rl) / 2.0
                  vortex_particles(ind)%data%alpha(1) = wxp(i) * volp(i)
                  vortex_particles(ind)%data%alpha(2) = wyp(i) * volp(i)
                  vortex_particles(ind)%data%alpha(3) = wzp(i) * volp(i)
                  vortex_particles(ind)%data%vol = volp(i)
               end if
            end do
         end do
         np = ind

      case (7)

         block

            real(kind_physics)     :: vecp(3), volp, wxp, wyp, wzp, unit_vec(3), mod_w, r_loc
            integer(kind_particle) :: ngrid(3)

            !&<
            ngrid(1) = nint((r_torus + rmax) / m_h) + 1
            ngrid(2) = nint((r_torus + rmax) / m_h) + 1
            ngrid(3) = nint(           rmax  / m_h) + 1
            !&>

            ind = 0
            ind0 = 0
            volp = m_h**3
            do i = -ngrid(1), ngrid(1)
               do j = -ngrid(2), ngrid(2)
                  do k = -ngrid(3), ngrid(3)

                     ind0 = ind0 + 1
                     vecp(1) = i * m_h
                     vecp(2) = j * m_h
                     vecp(3) = k * m_h
                     rr = sqrt(dot_product(vecp, vecp))
                     rr1 = sqrt(vecp(1)**2 + vecp(2)**2)

                     if (rr .gt. 0.d0) then
                        stheta = rr1 / rr
                     else
                        stheta = 0.d0
                     end if

                     r_loc = sqrt((rr1 - r_torus)**2 + vecp(3)**2)

                     if (rr .gt. (r_torus - rmax) .and. r_loc .lt. rmax) then
                        expo = kappa * (r_torus * r_torus + rr * rr - 2.d0 * r_torus * rr * stheta)
                        mod_w = kappa / pi * G / r_core / r_core * exp(-expo / r_core / r_core)
                     else
                        mod_w = 0.d0
                     end if
                     if (rr .gt. 0.d0) then
                        unit_vec = vecp / rr
                     else
                        unit_vec = vecp
                     end if

                     wxp = -unit_vec(2) * mod_w
                     wyp = unit_vec(1) * mod_w
                     wzp = 0.d0
                     if (mod(ind0 - 1, n_cpu) .eq. my_rank) then
                        ind = ind + 1
                        if (ind .gt. np) then
                           write (*, *) 'something is wrong here: to many particles in init of first ring', my_rank, ind, np, n
                           call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
                        end if
                        !&<
                        vortex_particles(ind)%x             = vecp
                        vortex_particles(ind)%data%alpha(1) = wxp * volp
                        vortex_particles(ind)%data%alpha(2) = wyp * volp
                        vortex_particles(ind)%data%alpha(3) = wzp * volp
                        vortex_particles(ind)%data%vol      = volp
                        !&>
                     end if
                  end do
               end do
            end do

            np = ind

            if (norm2(torus_offset) .gt. 0.d0) then
               do i = 1, np
                  !&<
                  vortex_particles(i + np)%x             = vortex_particles(i)%x + torus_offset
                  vortex_particles(i + np)%data%alpha(1) = -vortex_particles(i)%data%alpha(1)
                  vortex_particles(i + np)%data%alpha(2) = -vortex_particles(i)%data%alpha(2)
                  vortex_particles(i + np)%data%alpha(3) = vortex_particles(i)%data%alpha(3)
                  !&>
               end do

               np = 2 * np
            end if

         end block

      case (98) ! Random cubic setup (for testing purpose only)

         j = 0

         do i = 1, n

            xt = 0.
            yt = 0.
            zt = 0.

            call par_rand(par_rand_res)
            xt = par_rand_res
            call par_rand(par_rand_res)
            yt = par_rand_res
            call par_rand(par_rand_res)
            zt = par_rand_res

            if (mod(i + my_rank, n_cpu) .eq. 0) then
               j = j + 1
               if (j .gt. np) then
                  write (*, *) 'something is wrong here: to many particles in init', my_rank, j, n
                  call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
               end if
               vortex_particles(j)%x(1) = xt
               vortex_particles(j)%x(2) = yt
               vortex_particles(j)%x(3) = zt
               call par_rand(par_rand_res)
               vortex_particles(j)%data%alpha(1) = par_rand_res * h**3
               call par_rand(par_rand_res)
               vortex_particles(j)%data%alpha(2) = -par_rand_res * h**3
               call par_rand(par_rand_res)
               vortex_particles(j)%data%alpha(3) = par_rand_res * h**3
            end if
         end do
         np = j

      case (99) ! Read-in MPI checkpoints

         call read_in_checkpoint()

      end select config

      ! shrink vortex_particles to correct size
      ! reasoning: vortex_particles may have been allocated larger than required
      ! but may be queried for its size at a later stage
      if (np .ne. size(vortex_particles)) then
         block
            type(t_particle), allocatable :: temp_particles(:)
            call move_alloc(vortex_particles, temp_particles)
            vortex_particles = temp_particles(1:np)
         end block
      end if

      ! initial dump if we did not read in via MPI
      if (ispecial .ne. 99) then
         n_out = 0
         call kick_out_particles()
         call reset_labels() ! works on vortex_particles
         call dump(0, ts)
         vortex_particles(1:np)%work = 1.
         n_out = 1
      end if

   end subroutine special_start

   ! portable random number generator, see numerical recipes
   ! check for the random numbers:
   ! the first numbers should be 0.2853809, 0.2533582 and 0.2533582
   subroutine par_rand(res, iseed)

      !use physvars
      implicit none

      integer :: idum, idum2, iy, j, k
      integer :: iv(32)

      real, intent(inout) :: res
      integer, intent(in), optional :: iseed

      integer :: IM1, IM2, IMM1, IA1, IA2, IQ1, IQ2, IR1, IR2, NTAB, NDIV
      real    :: AM, RNMX

      save

      data idum, idum2/-1, 123456789/

      IM1 = 2147483563
      IM2 = 2147483399
      AM = 1.0 / IM1
      IMM1 = IM1 - 1
      IA1 = 40014
      IA2 = 40692
      IQ1 = 53668
      IQ2 = 52774
      IR1 = 12211
      IR2 = 3791
      NTAB = 32
      NDIV = 1 + IMM1 / NTAB
      RNMX = 1.0 - 1.2e-7

      if (idum .lt. 0) then

         if (present(iseed)) then
            idum = iseed
         else
            idum = 1
         end if

         idum2 = idum

         do j = NTAB + 7, 0, -1
            k = idum / IQ1
            idum = IA1 * (idum - k * IQ1) - k * IR1
            if (idum .lt. 0) idum = idum + IM1

            if (j .lt. NTAB) iv(j + 1) = idum

         end do
         iy = iv(1)

      end if

      k = idum / IQ1
      idum = IA1 * (idum - k * IQ1) - k * IR1
      if (idum .lt. 0) idum = idum + IM1

      k = idum2 / IQ2
      idum2 = IA2 * (idum2 - k * IQ2) - k * IR2
      if (idum2 .lt. 0) idum2 = idum2 + IM2

      j = iy / NDIV + 1
      iy = iv(j) - idum2
      iv(j) = idum

      if (iy .lt. 1) iy = iy + IMM1
      res = AM * iy
      if (res .gt. RNMX) res = RNMX

   end subroutine par_rand

end module init_particles
