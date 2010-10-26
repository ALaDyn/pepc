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
   !r      ','ni   ','rhoi   ','vi   ','er',' phi'
subroutine sum_radial(timestamp)

  use physvars
  implicit none
  include 'mpif.h'


  integer, intent(in) :: timestamp
  real*8 :: rdr, dr, cweight
  real*8 :: fr1, fr2, ra, xt, yt, zt, rt
  integer :: i, i1, i2, ngr, ierr
  character(30) :: cfile

  character(6) :: cdump
  real*8, dimension(0:ngx+1) :: vi_loc_min, vi_loc_max, ni_loc, gi_loc
  real*8, dimension(0:ngx+1) :: vi_glob_min, vi_glob_max, ni_glob, gi_glob
  real*8, dimension(0:ngx+1) :: volume
  real*8 :: vradial
  real*8 :: gmin=1.e-3


  ! get filename suffix from dump counter
  do i=0,4
     cdump(6-i:6-i) =  achar(mod(timestamp/10**i,10) + 48)  
  end do
  cdump(1:1) = achar(timestamp/10**5 + 48)

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
     xt = x(i)-xl/2.
     yt = y(i)-yl/2.
     zt = z(i)-zl/2.
     ! radial velocity
     vradial=dot_product([ux(i), uy(i), uz(i)], [xt, yt, zt])/sqrt(dot_product([xt, yt, zt],[xt, yt, zt]))
     ! radius
     rt = sqrt(xt**2+yt**2+zt**2)
     ra = rt*rdr
     ! indices - include zero for r=0
     i1 = min(int(ra),ngr+1)
     i2 = min(i1+1,   ngr+1)
     !  linear weighting - gather charge at nearest grid points
     fr1 = fraction(ra)
     fr2 = 1.- fr1

     cweight = abs(q(i))    !  charge weighting factor - densities computed later 

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
     cfile = "radial_fields."//cdump
     open (60,file=cfile)
     write(60,'(7(a12))') '#   r      ','ni   ','rhoi   ','vi_min','vi_max','er',' phi'
     write(60,'((7(1pe12.4)))') &
          ((i-1)*dr, gi_glob(i)/max(1,ni), ni_glob(i), &
		 vi_glob_min(i), vi_glob_max(i), ex(np_local-ngr+i), pot(np_local-ngr+i), i=1,ngr)
     close(60)

  endif


end subroutine sum_radial





