!  =================================
!
!    Tracking n_c for laser
!
!  =================================

subroutine track_nc

  use treevars

  implicit none


  real, dimension(0:ngx+1) :: rho1d

  integer :: i, j, k, ng, i1, i2, j1, j2, k1, k2, jfoc, kfoc
  character(30) :: cfile
  character(5) :: cme
  character(6) :: cdump, cvis
  logical :: found
  real :: xc1, rho_track, dx, dy, dz

  integer :: icm, icp

  dy = yl/ngy
  dz = zl/ngz
  dx = xl/ngx

! density average line-out along laser axis: 5x5 average, converted to n/nc

  jfoc = focus(2)/dy
  kfoc = focus(3)/dz
  rho1d(0:ngx+1) = 0.
  do k=kfoc-2,kfoc+2
     do j=jfoc-2,jfoc+2
        rho1d(1:ngx) = rho1d(1:ngx)+rhoi(1:ngx,j,k)/25./omega**2  
     end do
  end do


  ! Determine position of critical density along laser axis
 
! initial leading edge of plasma
  if (initial_config == 1 .or.initial_config ==3) then
     xc1 = plasma_centre(1)-r_sphere 
  else
     xc1 = plasma_centre(1)-x_plasma/2.  
  endif

 !  start indices
  icm = xc1/dx
  icp = icm+1
  found = .false.
  rho_track = 1.5   ! tracking density

!  sweep down
  i = icm
  do while ( i > 0 .and. .not.found)
     if ( rho1d(i)<= rho_track .and. rho1d(i+1)> rho_track ) then
        found=.true.
        x_crit = i*dx + dx*(1.-rho1d(i))/(rho1d(i+1)-rho1d(i))
     endif
     i = i-1
  end do

!  sweep up
  i = icm
  do while ( i < ngx .and. .not.found)
     if (rho1d(i)<= rho_track .and. rho1d(i+1)> rho_track) then
        found=.true.
        x_crit = i*dx + dx*(1.-rho1d(i))/(rho1d(i+1)-rho1d(i))
     endif
     i = i+1
  end do

!  if (.not.found .and. itime>0) beam_config=0 ! switch off laser
  if (.not.found) x_crit=xc1  ! original plasma edge 
!  if (.not.found) x_crit=x_crit + dt  ! vacuum propagation
 
if (me==0) then
   write(*,*) 'plasma edge: ',xc1, ' x_crit: ',x_crit
   write(15,*) 'plasma edge: ',xc1, ' x_crit: ',x_crit

endif


end subroutine track_nc
