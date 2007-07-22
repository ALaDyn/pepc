!  =================================
!
!    1D gather for spherically symmetric fields
!
!  =================================

subroutine sum_radial(timestamp)

  use physvars
  use treevars
  implicit none
  include 'mpif.h'


  integer, intent(in) :: timestamp
  real :: rdr, dr, vweight, jweight, cweight
  real :: tweight, erweight, rhalf
  real :: xpdum, ypdum, zpdum
  real :: fr1, fr2, ra, gamma, xt, yt, zt, rt
  real :: ttrav, tfetch
  integer :: ndum, i, j, k, ng, i1, i2, nelecs, nions, ngr, icall, ierr, p
  character(30) :: cfile
  character(5) :: cme
  character(6) :: cdump, cvis
  real, dimension(0:ngx+1) :: ve_loc, vi_loc, Er_loc, ni_loc, ne_loc, ge_loc, gi_loc
  real, dimension(0:ngx+1) :: ve_glob, vi_glob, Er_glob, ni_glob, ne_glob, ge_glob, gi_glob
  real, dimension(0:ngx+1) :: volume
  real*8, dimension(0:ngx+1) :: ex_g, ey_g, ez_g, phi_g, w_g
  integer, dimension(ngx+1) :: pshortl
  real :: gmin=1.e-3

  icall = timestamp/ivis_fields


  ! get filename suffix from dump counter
  do i=0,4
     cdump(6-i:6-i) =  achar(mod(timestamp/10**i,10) + 48)  
  end do
  cdump(1:1) = achar(timestamp/10**5 + 48)

  if (my_rank==0) then
     write(*,'(//a/a,f10.2)') 'Radial field dump:', &
          'writing out densities, fields on grid 0-',xl
  endif

  ngr = ngx  ! use box resolution in x for radial distn
  dr = xl/ngr

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

  do i=1,npp

! particle position relative to plasma center
     xt = x(i)-plasma_center(1)
     yt = y(i)-plasma_center(2)
     zt = z(i)-plasma_center(3)

! radius
     rt=sqrt(xt**2+yt**2+zt**2)
     ra = rt*rdr

!  indices - include zero for r=0
     i1=ra+.5
     i2=i1+1

     i1 = min(i1,ngr+1)
     i2 = min(i2,ngr+1)

     !  linear weighting
!     fr2=ra-i1  ! Prevent overflow/negative weighting for particles outside box
!     fr1=1.-fr1
    ! NGP
     fr2=0.
     fr1=1.
     !  gather charge at nearest grid points
     gamma = sqrt(1.0+ux(i)**2+uy(i)**2+uz(i)**2)
     cweight = abs(q(i))    !  charge weighting factor - densities computed later 
     vweight = sqrt(ux(i)**2+uy(i)**2+uz(i)**2)
     erweight =  sqrt(ex(i)**2+ey(i)**2+ez(i)**2)

     tweight = (gamma-1.)  ! K.E. of particle in keV

     if (q(i)<0) then
        ve_loc(i1) = ve_loc(i1) + fr1*vweight    
        ve_loc(i2) = ve_loc(i2) + fr2*vweight  

        ne_loc(i1) = ne_loc(i1) + fr1*cweight
        ne_loc(i2) = ne_loc(i2) + fr2*cweight

        ge_loc(i1) = ge_loc(i1) + fr1  ! weighted # electrons at r
        ge_loc(i2) = ge_loc(i2) + fr2

     else
        vi_loc(i1) = vi_loc(i1) + fr1*vweight  ! TODO: This won't give multi-valued vel
        vi_loc(i2) = vi_loc(i2) + fr2*vweight  ! beyond shock front - need phase space!

        ni_loc(i1) = ni_loc(i1) + fr1*cweight
        ni_loc(i2) = ni_loc(i2) + fr2*cweight

        gi_loc(i1) = gi_loc(i1) + fr1  ! weighted # ions at r
        gi_loc(i2) = gi_loc(i2) + fr2

     endif

  end do


  volume(1:ngr+1) = (/ (4./3.*pi*((i*dr+.5*dr)**3-(i*dr-.5*dr)**3), i=1,ngr+1) /)
!  ni_loc(1) = ni_loc(1) + ni_loc(0)  ! Fold charge at r=0 onto r=dr
!  ne_loc(1) = ne_loc(1) + ne_loc(0)
  volume(0) = 4./3.*pi*(dr/2.)**3  ! 1/2-vol weight for r=0

!  Renormalise charge densities with volume-weighting
  ni_loc = ni_loc/volume
  ne_loc = ne_loc/volume


! Gather partial sums together to get global moments - arrays run from (0:ngr)
 
  call MPI_ALLREDUCE(ge_loc, ge_glob, ngr+1, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(gi_loc, gi_glob, ngr+1, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(ve_loc, ve_glob, ngr+1, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(vi_loc, vi_glob, ngr+1, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(ne_loc, ne_glob, ngr+1, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(ni_loc, ni_glob, ngr+1, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)

  ge_glob = max(gmin,ge_glob)
  gi_glob = max(gmin,gi_glob)

  ! normalise averaged quantities
  nelecs = SUM(ge_glob(1:ngr))
  nions = SUM(gi_glob(1:ngr))
 if (my_rank==0) then
	write(*,*) "Charge check:",nelecs,nions
 endif
  ve_glob = ve_glob/ge_glob
  vi_glob = vi_glob/gi_glob

! Radial fields - set up dummy particles & find forces directly from tree
! Dummies set up at end of particle arrays on root to ensure unique labelling
    if (my_rank==0) then
    do i=1,ngr+1
        
        p = npp+i   !index
        pshortl(i) = p   !index
        x(p) = (i-1)*dr + plasma_center(1)   ! radius - include r=0
        y(p) = plasma_center(2)
        z(p) = plasma_center(3)
     end do
! Get interaction lists
!     write (*,*) 'Doing lists for dummy particles'
!     write (*,'((i8,f12.3))') (pshortl(i),x(pshortl(i)),i=1,ngr+1)
     ndum = ngr+1
   else
     ndum=0      
   endif  

! all PEs must call walk

   call tree_walk(pshortl(1:ndum),ndum,1,theta,eps,current_step,mac,ttrav,tfetch)

   if (my_rank==0) then
! Fields
     do i=1,ngr+1
        p=pshortl(i)
        call sum_force(p, nterm(i), nodelist( 1:nterm(i),i), eps, &
             ex_g(i-1), ey_g(i-1), ez_g(i-1), phi_g(i-1), w_g(i-1))
     end do
 
     Er_glob(0:ngr) = Ex_g(0:ngr)    !  Radial field (leave off norm constant)
  endif

! Write out to file

  if (my_rank == 0) then
     write(*,*) 'number density integrals: ',nelecs, nions
     cfile = "fields/radial."//cdump
     open (60,file=cfile)
     write(60,'(9(a12))') '!   r      ',' r/r0   ','ne    ','ni   ','rhoe   ','rhoi   ','ve   ','vi   ','er'
     write(60,'((10(1pe12.4)))') &
          (i*dr, i*dr, ge_glob(i)/max(1,ne), gi_glob(i)/max(1,ni), max(ne_glob(i),1.e-10), max(ni_glob(i),1.e-10), ve_glob(i), vi_glob(i), max(er_glob(i),1.e-10), phi_g(i), i=0,ngr)
     close(60)

  endif
  icall = icall + 1

end subroutine sum_radial





