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

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>
  !> 1D data gathering and output for spherically symmetric geometries
  !> prints ion density, charge density, radial velocity, radial field,
  !> and potential in dependence of distance r from the center into
  !> radial_fields.xxxxxx, where xxxxxx is the current timestep
  !>
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  =================================
!
!    1D gather for spherically symmetric fields
!
!  =================================
   !'r      ','ni   ','rhoi   ','vi   ','er',' phi'
subroutine sum_radial(timestamp)

  use physvars
  use module_units
  use module_utils
  use module_pepc_kinds
  use module_pepc_types
  use mpi
  implicit none


  integer, intent(in) :: timestamp
  real*8 :: rdr, dr, cweight
  real*8 :: fr1, fr2, ra, xt, yt, zt, rt
  integer :: i, i1, i2, ngr, ierr
  character(30) :: cfile

  real*8, dimension(0:ngx+1) :: vi_loc_min, vi_loc_max, ni_loc, gi_loc
  real*8, dimension(0:ngx+1) :: vi_glob_min, vi_glob_max, ni_glob, gi_glob
  real*8, dimension(0:ngx+1) :: volume
  real*8 :: vradial
  real*8 :: gmin=1.e-3

  if (my_rank==0) write(*,'("Radial fields: writing out densities, fields on grid 0.0 - ",f4.2)') xl

  ngr=ngx
  dr = xl/ngr

  rdr = 1./dr

  vi_loc_min = 1e10
  vi_loc_max = 0.
  gi_loc = 0.
  ni_loc = 0.

  !  Accumulate local radial moments
  do i=1,np_local-ngr
     ! particle position relative to box centre
     xt = particles(i)%x(1)-xl/2.
     yt = particles(i)%x(2)-yl/2.
     zt = particles(i)%x(3)-zl/2.
     ! radial velocity
     vradial=dot_product([particles(i)%data%v(1), particles(i)%data%v(2), particles(i)%data%v(3)], [xt, yt, zt])/sqrt(dot_product([xt, yt, zt],[xt, yt, zt]))
     ! radius
     rt = sqrt(xt**2+yt**2+zt**2)
     ra = rt*rdr
     ! indices - include zero for r=0
     i1 = min(int(ra),ngr+1)
     i2 = min(i1+1,   ngr+1)
     !  linear weighting - gather charge at nearest grid points
     fr1 = fraction(ra)
     fr2 = 1.- fr1

     cweight = abs(particles(i)%data%q)    !  charge weighting factor - densities computed later 

        if (fr1>fr2) then
          vi_loc_min(i1) = min(vi_loc_min(i1), vradial)
          vi_loc_max(i1) = max(vi_loc_max(i1), vradial)
        else
          vi_loc_min(i2) = min(vi_loc_min(i2), vradial)
          vi_loc_max(i2) = max(vi_loc_max(i2), vradial)
        end if

        ni_loc(i1) = ni_loc(i1) + fr1*cweight
        ni_loc(i2) = ni_loc(i2) + fr2*cweight

        gi_loc(i1) = gi_loc(i1) + fr1  ! weighted # ions at r
        gi_loc(i2) = gi_loc(i2) + fr2

  end do


  volume(1:ngr+1) = (/ (4./3.*pi*((i*dr+.5*dr)**3-(i*dr-.5*dr)**3), i=1,ngr+1) /)
  !  ni_loc(1) = ni_loc(1) + ni_loc(0)  ! Fold charge at r=0 onto r=dr
  !  ne_loc(1) = ne_loc(1) + ne_loc(0)
  volume(0) = 4./3.*pi*(dr/2.)**3  ! 1/2-vol weight for r=0

  !  Renormalise charge densities with volume-weighting
  ni_loc = ni_loc/volume

  where(vi_loc_min == 1e10) vi_loc_min = 0.

  ! Gather partial sums together to get global moments - arrays run from (0:ngr)
  call MPI_ALLREDUCE(gi_loc, gi_glob, ngr+1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(vi_loc_min, vi_glob_min, ngr+1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(vi_loc_max, vi_glob_max, ngr+1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(ni_loc, ni_glob, ngr+1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

  gi_glob = max(gmin,gi_glob)

  ! Fields ex and pot have been computed in pepc_fields() with the rest of the particles
  ! by just introducing ghost particles

  ! Write out to file
  if (my_rank == 0) then
     call create_directory("radial_fields")
     write(cfile,'("radial_fields/radial_fields.",i6.6)') timestamp
     open (60,file=trim(cfile))
     write(60,'(7(a12))') '#   r      ','ni   ','rhoi   ','vi_min','vi_max','er',' phi'
     write(60,'((7(1pe12.4)))') &
          ((i-1)*dr, gi_glob(i)/max(1_kind_particle,ni), ni_glob(i), &
		 vi_glob_min(i), vi_glob_max(i), particles(np_local-ngr+i)%results%e(1), particles(np_local-ngr+i)%results%pot, i=1,ngr)
     close(60)

  endif


end subroutine sum_radial





