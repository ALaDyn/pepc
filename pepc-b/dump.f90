! ======================
!
!   DUMP
!
!   Gather and write out particle data for restart 
!
!
! ======================

subroutine dump(timestamp)


  use physvars
  use treevars
  implicit none   
  include 'mpif.h'

  character(30) :: cfile
  character(6) :: cdump
  character(5) :: cme
  integer, intent(in) :: timestamp
  integer :: i, j, ioffset, idummy=0, ierr
  integer, save :: icall=0
  real :: simtime


  simtime = dt*timestamp


  ! get filename suffix from dump counter
  do i=0,4
     cdump(6-i:6-i) =  achar(mod(timestamp/10**i,10) + 48)  
  end do
  cdump(1:1) = achar(timestamp/10**5 + 48)

  cfile="data/pe"//csubme//"/parts_info."//cdump(1:6)


  open (60,file=cfile)    
  write(60,'(7(a9,i8/),10(a9,f12.5/),9(a9,1pe12.5/),2(a9,3f12.5/))')  &    ! info block
       'itime=',timestamp, 'npp=',npp, &
       'ne=',ne, 'ni=',ni, 'npbeam=',np_beam, 'geometry=', target_geometry, &
       'scheme=',scheme, &
       'xl=',xl, 'yl=',yl, 'zl=',zl, 'boxsize=',zl, &
       'eps=', eps, 'theta=',theta,' tlaser= ',tlaser,' trun= ',trun, &
       'omega=',omega,'lambda=',lambda,'  qe=',qe,'  qi=',qi, &
       'mass_e=',mass_e,'mass_i=',mass_i,'Zion=',Zion,'a_ii=',a_ii, &
       'Vplas=',Vplas,'Aplas=',Aplas,'Qplas=',Qplas, &
       'centre=',plasma_centre(1:3),'focus=',focus(1:3)
  

  if (me.eq.0) then

    cfile="parts_info.in"     ! copy to default restart block
    open (61,file=cfile)    
    write(61,'(7(a9,i8/),10(a9,f12.5/),9(a9,1pe12.5/),2(a9,3f12.5/))')  &    ! info block
       'itime=',timestamp, 'npp=',npp, &
       'ne=',ne, 'ni=',ni, 'npbeam=',np_beam, 'geometry=', target_geometry, &
       'scheme=',scheme, &
       'xl=',xl, 'yl=',yl, 'zl=',zl, 'boxsize=',zl, &
       'eps=', eps, 'theta=',theta,' tlaser= ',tlaser,' trun= ',trun, &
       'omega=',omega,'lambda=',lambda,'  qe=',qe,'  qi=',qi, &
       'mass_e=',mass_e,'mass_i=',mass_i,'Zion=',Zion,'a_ii=',a_ii, &
       'Vplas=',Vplas,'Aplas=',Aplas,'Qplas=',Qplas, &
       'centre=',plasma_centre(1:3),'focus=',focus(1:3)
  
    close (61)
    open (62,file="runstamp")  ! time stamp 
    write(62,'(a)') cdump(1:6)
    close (62)

  write(6,'(//a/7(a9,i8/),10(a9,f12.5/),9(a9,1pe12.5/),2(a9,3f12.5/))') 'PARTICLE DUMP:', &    ! info block
       'itime=',timestamp, 'npp=',npp, &
       'ne=',ne, 'ni=',ni, 'npbeam=',np_beam, 'geometry=', target_geometry, &
       'scheme=',scheme, &
       'xl=',xl, 'yl=',yl, 'zl=',zl, 'boxsize=',zl, &
       'eps=', eps, 'theta=',theta,' tlaser= ',tlaser,' trun= ',trun, &
       'omega=',omega,'lambda=',lambda,'  qe=',qe,'  qi=',qi, &
       'mass_e=',mass_e,'mass_i=',mass_i,'Zion=',Zion,'a_ii=',a_ii, &
       'Vplas=',Vplas,'Aplas=',Aplas,'Qplas=',Qplas, &
       'centre=',plasma_centre(1:3),'focus=',focus(1:3)

  endif
  close(60)


  
  cfile="data/pe"//csubme//"/parts_dump."//cdump(1:6)
  open (60,file=cfile) 
  write(60,'((12(1pe14.5),2i9))')  &
       (x(i), y(i), z(i), ux(i), uy(i), uz(i), q(i), m(i), &
        Ex(i), Ey(i), Ez(i), &  ! electric field
        pot(i), &  ! potential
        pepid(i), pelabel(i),i=1,npp)
  close(60)

  icall = icall + 1

end subroutine dump








