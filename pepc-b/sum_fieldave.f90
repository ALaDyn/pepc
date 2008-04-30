!  =================================
!
!    1D gather for time-averaged electric fields
!
!  =================================

subroutine sum_fieldave

  use physvars
  use treevars
  implicit none
  include 'mpif.h'


  real :: dx,dy
  real :: fr1, fr2, ra, gamma, xt, yt, zt, rt
  real :: ttrav, tfetch
  integer :: ndum, i, j, k, ng, i1, i2, nelecs, nions, ngr, icall, ierr, p
  character(30) :: cfile
  character(5) :: cme
  character(6) :: cdump, cvis

  real*8, dimension(0:ngav+1) :: ex_g, ey_g, ez_g, phi_g, w_g
  integer, dimension(ngav+1) :: pshortl


  if (my_rank==0) then
     write(*,'(//a/a,f10.2)') '1D field average '
!          ' fields on grid 0-',xl
  endif

  !  field box limits 

  dx = (xgav_end-xgav_start)/ngav



! Axial electric field 
! --------------------
! - set up dummy particles at grid points & find forces directly from tree
! Dummies set up at end of particle arrays on root to ensure unique labelling

  if (my_rank==0) then
      do i=1,ngav+1
        
        p = npp+i   !index
        pshortl(i) = p   !index
        x(p) = (i-1)*dx + xgav_start   ! axial position in box
        y(p) = plasma_centre(2)  ! target centre
        z(p) = plasma_centre(3)
      end do

! Get interaction lists
!     write (*,*) 'Doing lists for dummy particles'
!     write (*,'((i8,f12.3))') (pshortl(i),x(pshortl(i)),i=1,ngav+1)
     ndum = ngav+1
   else
     ndum=0  ! Remainder of CPUs just have to provide multipole info      
   endif  

! all CPUs must call walk
!subroutine tree_walk(pshort,npshort, !pass,theta,eps,itime,mac,twalk,tfetch,ex_nps,ey_nps,ez_nps,np_local)

! mac set to BH, theta to 0.5

   call tree_walk(pshortl(1:ndum),ndum,1,0.5,eps,itime,0,ttrav,tfetch)

   if (my_rank==0) then
! Fields
     do i=1,ngav+1
        p=pshortl(i)
        call sum_force(p, nterm(i), nodelist( 1:nterm(i),i), eps, &
             ex_g(i-1), ey_g(i-1), ez_g(i-1), phi_g(i-1), w_g(i-1))
     end do
 
     ex_ave = ex_ave + force_const*Ex_g/navcycle    ! Accumulate axial field
  endif


! Radial electric field 
! ----------------------
! - set up dummy particles at grid points & find forces directly from tree
! Dummies set up at end of particle arrays on root to ensure unique labelling

  dy = yl/ngav
  if (my_rank==0) then
      do i=1,ngav+1
        
        p = npp+i   !index
        pshortl(i) = p   !index
        x(p) = xgav_pos(1)   ! 1st axial position in box
        y(p) = dy*(i-1)  ! y position 
        z(p) = plasma_centre(3) ! midpoint in z
      end do
! Get interaction lists
     ndum = ngav+1
   else
     ndum=0  ! Remainder of CPUs just have to provide multipole info      
   endif  

! all CPUs must call walk

   call tree_walk(pshortl(1:ndum),ndum,1,0.5,eps,itime,0,ttrav,tfetch)

   if (my_rank==0) then
! Fields
     do i=1,ngav+1
        p=pshortl(i)
        call sum_force(p, nterm(i), nodelist( 1:nterm(i),i), eps, &
             ex_g(i-1), ey_g(i-1), ez_g(i-1), phi_g(i-1), w_g(i-1))
     end do
 
! Accumulate radial field
     ey_ave(1:ngav,1) = ey_ave(1:ngav,1) + force_const*Ex_g(1:ngav)/navcycle 
  endif



end subroutine sum_fieldave


