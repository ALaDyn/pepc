!  =================================
!
!    Tracking n_c for laser
!
!   $Revision 1.5$
!
!  =================================

subroutine track_nc

  use physvars

  implicit none
  include 'mpif.h'

  real, dimension(0:ngx+1) :: rho1d

  integer :: i, j, k, ng,jfoc, kfoc, icrit
  logical :: found
  real :: xc1, dx, dy, dz

  integer :: icm, icp, ixc, nover
  integer :: ierr

  call densities  ! Compute local ion density rhoi_loc on 3D grid

  ng = (ngx+2)*(ngy+2)*(ngz+2)                         ! total # gridpoints

  call MPI_ALLREDUCE(rhoi_loc, rhoi, ng, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)

  dy = yl/ngy
  dz = zl/ngz
  dx = xl/ngx

! density average line-out along laser axis: 5x5 average, converted to n/nc

  jfoc = max(2.0,focus(2)/dy) ! ensure index>=0
  kfoc = max(2.0,focus(3)/dz)
  rho1d(0:ngx+1) = 0.
  do k=kfoc-2,kfoc+2
     do j=jfoc-2,jfoc+2
        rho1d(1:ngx) = rho1d(1:ngx)+rhoi(1:ngx,j,k)/25./omega**2  
     end do
  end do


  ! Determine position of critical density along laser axis
 
! initial leading edge of plasma
  if (target_geometry == 1 .or. target_geometry ==3) then
     xc1 = plasma_centre(1)-r_sphere 
  else
     xc1 = plasma_centre(1)-x_plasma - x_layer(1)  
  endif

! abort tracking if particle number too low
  if (npart_total.le.100) then
	write(*,*) 'fixing xcrit at ',xc1
     x_crit=xc1
     rho_upper=1.0
     return
  endif

 !  start indice
  ixc = xc1/dx
  icm = min(max(ixc,1),ngx-1)
  icp = icm+1
  icrit = icm
  found = .false.
  rho_track = 1.5   ! tracking density

!  sweep down
  i = icm
  do while ( i > 0 .and. .not.found)
     if ( rho1d(i)<= rho_track .and. rho1d(i+1)> rho_track ) then
        found=.true.
        x_crit = i*dx + dx*(1.-rho1d(i))/(rho1d(i+1)-rho1d(i))
        icrit=i
     endif
     i = i-1
  end do

!  sweep up
  i = icm
  do while ( i < ngx .and. .not.found)
     if (rho1d(i)<= rho_track .and. rho1d(i+1)> rho_track) then
        found=.true.
        x_crit = i*dx + dx*(1.-rho1d(i))/(rho1d(i+1)-rho1d(i))
        icrit = i
     endif
     i = i+1
  end do

  if (.not.found .and. itime>0) then 
	beam_config=0 ! switch off laser
	if (my_rank==0) write(6,*) 'Target burnt through - switching off laser'
  endif

  rho_upper=0.
  nover = 5./dx

  do i=icrit+1,min(icrit+nover,ngx)
     rho_upper = rho_upper + rho1d(i)
  end do

  rho_upper = rho_upper/nover  ! ave. upper shelf density

!  if (.not.found) x_crit=xc1  ! original plasma edge 
!  if (.not.found) x_crit=x_crit + dt  ! vacuum propagation
 
if (my_rank==0) then
   write(*,'(/a15,f10.3,a15,f10.3,a15,f10.3)') &
        'plasma edge: ',xc1, ' x_crit: ',x_crit,' n_upper: ',rho_upper
   write(24,*) 'plasma edge: ',xc1, ' x_crit: ',x_crit

endif


end subroutine track_nc