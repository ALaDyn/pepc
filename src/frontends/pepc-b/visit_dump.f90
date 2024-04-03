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
!   VISIT_DUMP
!
!   Write out particle data for
!   VISIT postprocessing
!   Each PE writes to its own file - to be concatenated later
!
!
! ======================

subroutine visit_dump(timestamp)

   use treevars
   use module_spacefilling
   use module_utils
   implicit none

   integer, parameter :: npart_visit_max = 250000  ! Max 25k data points for VIS
   real, dimension(npart_visit_max) :: xvis, yvis, zvis, vx, vy, vz, axvis, ayvis, azvis, qvis, mvis
   integer, dimension(npart_visit_max) :: ppid, plabel

   integer, dimension(npart_visit_max) :: icolour
   real, dimension(ngx) :: work1, work2

   real :: dx, dz, dy, xd, yd, zd, dummy, simtime, epondx, epondy, epondz, phipond, epond_max, box_max
   real :: Qtot, Qbox

   character(30) :: cfile
   character(5) :: cme
   character(6) :: cdump, cvis
   integer, intent(in) :: timestamp
   integer :: i, j, k, ioffset, idummy = 0, ilev, ixd, iyd, izd, npx, npz, npy
   integer :: icall, lcount
   integer :: jfoc, kfoc

   icall = timestamp / ivis
   simtime = timestamp * dt

   ! Filename (directory) prefix
   cme = "pe"//achar(me / 100 + 48)//achar(mod(me / 10, 10) + 48)//achar(mod(me, 10) + 48)

   ! get filename suffix from dump counter
   do i = 0, 4
      cdump(6 - i:6 - i) = achar(mod(timestamp / 10**i, 10) + 48)
   end do
   cdump(1:1) = achar(timestamp / 10**5 + 48)

   ! get filename suffix from call counter
   do i = 0, 4
      cvis(6 - i:6 - i) = achar(mod(icall / 10**i, 10) + 48)
   end do
   cvis(1:1) = achar(icall / 10**5 + 48)

   call create_directory(cme)

   cfile = cme//"/wf_info."//cdump
   open (60, file=cfile)
   write (60, '(a9,i8/4(a9,f12.5/),4(a9,i8/),a9,f8.3)') &
      'npp=', npp, 'xl=', xl, 'yl=', yl, 'zl=', zl, 'boxsize=', zl, &
      'ne=', ne, 'ni=', ni, 'np=', np_beam, 'itime=', timestamp, 't =', simtime
   close (60)

   cfile = cme//"/wf_ele."//cdump
   open (60, file=cfile)
   cfile = cme//"/wf_ion."//cdump
   open (61, file=cfile)

   do i = 1, npp
      if (pelabel(i) .le. ne) then
         write (60, '(10(1pe12.3),i6,i9)') &
            x(i), y(i), z(i), ux(i), uy(i), uz(i), q(i) &
            , ax(i) * m(i) / q(i), ay(i) * m(i) / q(i), az(i) * m(i) / q(i) &  ! electric field = m.a/q
            , pepid(i), pelabel(i)

      else if (pelabel(i) .le. ne + ni) then
         write (61, '(10(1pe12.3),i6,i9)') &
            x(i), y(i), z(i), ux(i), uy(i), uz(i), q(i) &
            , ax(i) * m(i) / q(i), ay(i) * m(i) / q(i), az(i) * m(i) / q(i) &
            , pepid(i), pelabel(i)
      else
      end if

   end do
   close (60)
   close (61)

   if (me .eq. 0) then
      !  Fields: electron and ion densities within xl*yl*zl
      cfile = "wf_field."//cdump
      open (62, file=cfile)

      dx = xl / ngx
      dy = yl / ngy
      dz = zl / ngz
      Qbox = 0.
      do k = 1, ngz
         do j = 1, ngy
            do i = 1, ngx
               Qbox = Qbox + rhoe(i, j, k) * dx * dy * dz
               write (62, '(3f13.5,2e13.3)') i * dx, j * dy, k * dz, abs(rhoe(i, j, k)), rhoi(i, j, k)
            end do
         end do
      end do
      Qtot = SUM(rhoe) * dx * dy * dz  ! including ghost cells
      if (debug .gt. 0) write (ipefile, '(4(a,f14.5/))') &
         'Total charge on grid:', Qbox, &
         '         ghost cells:', Qtot - Qbox, &
         '                 sum:', Qtot, &
         'Initial charge Q_s*Ne = rho0*V = ', Vplas * rho0

      close (62)

      cfile = "wf_field_slice."//cdump
      open (62, file=cfile)

      jfoc = focus(2) / dy
      kfoc = focus(3) / dz

      ! density average line-out along laser axis: 5x5 average, converted to n/nc

      work1 = 0.
      work2 = 0.
      do k = kfoc - 2, kfoc + 2
         do j = jfoc - 2, jfoc + 2
            work1(1:ngx) = work1(1:ngx) + rhoi(1:ngx, j, k) / 25./omega**2  ! density slice along laser axis: 5x5 average
            work2(1:ngx) = work2(1:ngx) + rhoe(1:ngx, j, k) / 25./omega**2
         end do
      end do

      write (62, '(3f13.5)') (i * dx, work1(i), work2(i), i=1, ngx)
      close (62)

      ! ship pond force as 'domain boxes' on grid - could carve this up as well
      if (beam_config .eq. 4) then
         box_max = xl / 25.
         epond_max = sqrt(1.+vosc**2 / 2.) / 2.
         lcount = 0
         dx = 2.*xl / ngx
         dz = zl / ngz
         dy = yl / ngy
         do i = 1, ngx
            do j = 1, ngy
               do k = 1, ngz
                  lcount = lcount + 1
                  xvis(lcount) = -1.5 * xl + (i - 0.5) * dx
                  yvis(lcount) = (j - 0.5) * dy
                  zvis(lcount) = (k - 0.5) * dz
                  xd = xvis(lcount) - focus(1) ! position relative to laser focus
                  yd = yvis(lcount) - focus(2)
                  zd = zvis(lcount) - focus(3)
                  call fpond(simtime, tpulse, sigma, vosc, omega, -xd, yd, zd, epondx, epondy, epondz, phipond)  ! evaluate pond force
                  mvis(lcount) = max(0.02, abs(epondx) / epond_max * box_max)  ! scale box size for plot
                  icolour(lcount) = 100 * epondx / epond_max * box_max
                  qvis(lcount) = phipond * 10
               end do
            end do
         end do

      else

         ! ship branch nodes to show domains
         do j = 1, nbranch_sum
            ilev = level_from_key(branch_key(j))
            ixd = SUM((/(2**i * ibits(branch_key(j), 3 * i, 1), i=0, ilev - 1)/))
            iyd = SUM((/(2**i * ibits(branch_key(j), 3 * i + 1, 1), i=0, ilev - 1)/))
            izd = SUM((/(2**i * ibits(branch_key(j), 3 * i + 2, 1), i=0, ilev - 1)/))
            mvis(j) = boxsize / 2**(ilev)          !  box length
            xvis(j) = ixd * mvis(j) + xmin
            yvis(j) = iyd * mvis(j) + ymin
            zvis(j) = izd * mvis(j) + zmin
            icolour(j) = branch_owner(j)
            lcount = nbranch_sum
         end do
         !       call flvisit_spk_domains_send( simtime, xvis,yvis,zvis,size,branch_owner,xvis,nbranch_sum)

      end if

      cfile = "wf_domain."//cdump
      open (60, file=cfile)
      write (60, '((4(f12.3),i6,2f12.3))') &
         (xvis(i), yvis(i), zvis(i), mvis(i), icolour(i), mvis(i), qvis(i), i=1, lcount)
      close (60)
   end if

   call MPI_BARRIER(MPI_COMM_WORLD, ierr)  ! Wait for everyone to catch up
   icall = icall + 1

end subroutine visit_dump
