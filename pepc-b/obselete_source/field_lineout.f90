!  =================================
!
!    1D gather for field lineout along 0-xl
!
!  =================================

subroutine field_lineout(timestamp)

  use physvars
  use treevars
  implicit none
  include 'mpif.h'


  integer, intent(in) :: timestamp
  real :: dr

  real :: ttrav, tfetch
  integer :: ndum, i, ngr, icall, p
  character(30) :: cfile
  character(6) :: cdump
  real*8, dimension(0:ngx+1) :: ex_g, ey_g, ez_g, phi_g, w_g
  integer, dimension(ngx+1) :: pshortl


  icall = timestamp/ivis_fields


  ! get filename suffix from dump counter
  do i=0,4
     cdump(6-i:6-i) =  achar(mod(timestamp/10**i,10) + 48)  
  end do
  cdump(1:1) = achar(timestamp/10**5 + 48)

  if (my_rank==0) then
     write(*,'(//a/a,f10.2)') 'Field lineout:', &
          'writing out Ex and Phi on grid 0-',xl
  endif

  ngr = ngx  ! use box resolution in x for radial distn
  dr = xl/ngr

! Field and potential lineouts - set up dummy particles & find forces directly from tree
! Dummies set up at end of particle arrays on root to ensure unique labelling
    if (my_rank==0) then
    do i=1,ngr+1
        
        p = npp+i   !index
        pshortl(i) = p   !index
        x(p) = (i-1)*dr   ! start at x=0
        y(p) = plasma_centre(2)
        z(p) = plasma_centre(3)
     end do
! Get interaction lists
!     write (*,*) 'Doing lists for dummy particles'
!     write (*,'((i8,f12.3))') (pshortl(i),x(pshortl(i)),i=1,ngr+1)
     ndum = ngr+1
   else
     ndum=0      
   endif  

! all PEs must call walk
! mac set to BH, theta to 0.5
   call tree_walk(pshortl(1:ndum),ndum,1,0.5,eps,itime,0,ttrav,tfetch)

   if (my_rank==0) then
! Fields
     do i=1,ngr+1
        p=pshortl(i)
        call sum_force(p, nterm(i), nodelist( 1:nterm(i),i), eps, &
             ex_g(i-1), ey_g(i-1), ez_g(i-1), phi_g(i-1), w_g(i-1))
     end do
 
  endif

! Write out to file

  if (my_rank == 0) then

     cfile = "fields/lineout."//cdump
     open (60,file=cfile)
     write(60,'(3(a12))') '!   x   ','Ex','Phi'
     write(60,'((3(1pe12.4)))') &
          (i*dr, ex_g(i), phi_g(i), i=0,ngr)
     close(60)

  endif
  icall = icall + 1

end subroutine field_lineout





