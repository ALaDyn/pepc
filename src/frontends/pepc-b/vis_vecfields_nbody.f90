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

! ======================
!
!   VIS_FIELDS
!
!   Send field data to VISIT for visualisation of scalar volumetric data
!
!
! ======================

subroutine vis_vecfields_nbody(timestamp)

   use physvars
   use treevars
   use mpi
   implicit none

   integer, intent(in) :: timestamp
   real*4, dimension(3*ngx*ngy*ngz) :: field1
   real*4, dimension(ngx, ngy, ngz) :: jelec_x, jelec_y, jelec_z
   real :: s, simtime, dummy, xd, yd, zd, dx, dz, dy, epond_max
   real :: xt, yt, zt, rt
   integer, parameter :: ngmax = 100
   integer :: i, j, k, ioffset, ixd, iyd, izd, ilev, lcount, iskip, itlas
   integer :: lvisit_active = 0, ierr
   integer :: npx, npy, npz, ng, jfoc, kfoc, nave
   real :: norm, j0
   integer :: iskip_x, iskip_y, iskip_z
   integer :: fselect = 0 ! field selector
   integer :: incdf
   real*4 :: grid_pars(24)  ! origins and mesh size of vis fields
   character(30) :: cfile
   character(5) :: cme
   character(6) :: cdump, cvis

   simtime = dt * (itime + itime_start)

#ifdef VISIT_NBODY
   if (me .eq. 0) call flvisit_nbody2_check_connection(lvisit_active)
   call MPI_BCAST(lvisit_active, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#endif

   if (lvisit_active .eq. 0) then
      if (me .eq. 0) write (*, '(a)') 'VIS_VECF     | No connection to visualization'
   end if

   ! Connected to vis, so proceed with field select & gather

   fselect = 2  ! Manual select - B-field
#ifdef VISIT_NBODY
   ! Fetch user-selected config from vis (TODO)
   if (me .eq. 0 .and. lvisit_active .ne. 0) then
      !        call flvisit_nbody2_selectvecfields_recv(fselect)
      write (*, '(a,i8)') "VIS_VECF    | Selected vector field", fselect
   end if
#endif

   ! get filename suffix from dump counter
   do i = 0, 4
      cdump(6 - i:6 - i) = achar(mod(timestamp / 10**i, 10) + 48)
   end do
   cdump(1:1) = achar(timestamp / 10**5 + 48)

   ng = ngx * ngy * ngz                         ! total # gridpoints
   ! Merge sums for gridded fields

   ! Should combine vec components into single array
   call MPI_ALLREDUCE(jxe_loc(1:ngx, 1:ngy, 1:ngz), jelec_x, ng, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
   call MPI_ALLREDUCE(jye_loc(1:ngx, 1:ngy, 1:ngz), jelec_y, ng, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
   call MPI_ALLREDUCE(jze_loc(1:ngx, 1:ngy, 1:ngz), jelec_z, ng, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)

   if (me .eq. 0) then

      lcount = 1

      ! limit size of field data to 100^3
      !     iskip_x = ngx/ngmax + 1
      !     iskip_y = ngy/ngmax + 1
      !     iskip_z = ngz/ngmax + 1
      iskip_x = 1
      iskip_y = 1
      iskip_z = 1

      !     npx = ngx/iskip_x + mod(ngx,2)
      !     npy = ngy/iskip_y + mod(ngy,2)
      !     npz = ngz/iskip_z + mod(ngz,2)
      npx = ngx
      npy = ngy
      npz = ngz

      dx = xl / ngx
      dz = zl / ngz
      dy = yl / ngy

      j0 = 1.0

      do k = 1, ngz, iskip_z
         do j = 1, ngy, iskip_y
            do i = 1, ngx, iskip_x

               ! Vector field selection
               f1: select case (fselect)
               case (1)  ! electron current
                  field1(lcount) = jelec_x(i, j, k)
                  field1(lcount + 1) = jelec_y(i, j, k)
                  field1(lcount + 2) = jelec_z(i, j, k)
                  !             field1(lcount+1)=0.
                  !             field1(lcount+2)=0.

               case (0)
                  xt = dx * (i - 0.5)
                  yt = dy * (j - 0.5) - yl / 2.
                  zt = dz * (k - 0.5) - zl / 2.
                  rt = sqrt(yt**2 + zt**2)
                  !                norm = rt*exp(-2*rt/yl)*xt/xl
                  norm = 1.
                  field1(lcount) = j0 / 3.
                  !                field1(lcount) = 0.
                  field1(lcount + 1) = -norm * zt / yl
                  field1(lcount + 2) = norm * yt / zl

               case (2)

                  if (beam_config_in .eq. 7) then ! Constant B in z-direction - either charge
                     field1(lcount) = 0.
                     field1(lcount + 1) = 0.
                     field1(lcount + 2) = vosc

                  else if (beam_config_in .eq. 17) then !  Z-pinch: circular B in x,y
                     xt = dx * (i - 0.5) - plasma_centre(1)
                     yt = dy * (j - 0.5) - plasma_centre(2)
                     rt = sqrt(xt**2 + yt**2)
                     field1(lcount) = -vosc * yt / rt
                     field1(lcount + 1) = vosc * xt / rt
                     field1(lcount + 2) = 0.

                  else if (beam_config_in .eq. 37) then !  Mirror in Bz
                     xt = dx * (i - 0.5) - plasma_centre(1)
                     yt = dy * (j - 0.5) - plasma_centre(2)
                     rt = sqrt(xt**2 + yt**2)
                     zt = dz * (k - 0.5) - zl / 2.
                     if (zt .gt. -zl / 2. .and. zt .lt. zl / 2.) then
                        field1(lcount) = -vosc / 5.*xt / rt * (2 * zt / zl)
                        field1(lcount + 1) = -vosc / 5.*yt / rt * (2 * zt / zl)
                        field1(lcount + 2) = vosc * (2 * zt / zl)**2 + vosc / 5.
                     else
                        field1(lcount:lcount + 2) = 0.

                     end if
                  end if
               end select f1
               lcount = lcount + 3
            end do
         end do
      end do

      !        write (*,*) 'Max currents:',maxval(jelec_x),maxval(jelec_y),maxval(jelec_z)
      write (*, *) 'Min/Max current:', minval(field1), maxval(field1)

#ifdef VISIT_NBODY
      call flvisit_nbody2_check_connection(lvisit_active)

      ! Tell vis which fields are coming

      if (lvisit_active .ne. 0) then
         !     call flvisit_nbody2_selectedvecfields_send(fselect)

#ifdef NETCDFLIB
         ! Netcdf write
         !         if (netcdf) call ncnbody_putselfield( ncid, simtime, fselect1, fselect2, fselect3, fselect4, incdf )
#endif

         !  Set up vis field grid - just one field at present
         !   do i=0,3
         i = 0
         grid_pars(6 * i + 1:6 * i + 3) = 0.
         grid_pars(6 * i + 4) = dx
         grid_pars(6 * i + 5) = dy
         grid_pars(6 * i + 6) = dz
         !   end do

         call flvisit_nbody2_vecfielddesc_send(grid_pars, 1, 6)

#ifdef NETCDFLIB
         !   if (netcdf) call ncnbody_putvecfielddesc( ncid, simtime, grid_pars, incdf )
#endif

         !   write(*,*) 'Grids: ',grid_pars

         !      if (fselect>0) then
         write (*, '(a,i8,a10,2(1pe12.3))') "VIS_VECF   | Shipping vector field", fselect, " min/max =", &
            minval(field1), maxval(field1)

         call flvisit_nbody2_vecfield1_send(field1, npx, npy, npz, 3)

#ifdef NETCDFLIB
         !         if (netcdf) call ncnbody_putvecfield( ncid, simtime, 1, npx, npy, npz, field1, incdf )
#endif

         !      endif
      end if
#endif

   end if

end subroutine vis_vecfields_nbody

