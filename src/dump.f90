! ======================
!
!   DUMP
!
!   Gather and write out particle data for restart 
!
!
! ======================

subroutine dump(timestamp)


  use treevars
  implicit none   


  character(30) :: cfile
  character(6) :: cdump
  character(5) :: cme
  integer, intent(in) :: timestamp
  integer :: i, j, ioffset, idummy=0
  integer, save :: icall=0
  real :: simtime


  simtime = dt*timestamp

  ! Filename (directory) prefix
  cme = "pe"// achar(me/100+48) // achar(mod(me/10,10)+48) // achar(mod(me,10)+48)  

  ! get filename suffix from dump counter
  do i=0,4
     cdump(6-i:6-i) =  achar(mod(timestamp/10**i,10) + 48)  
  end do
  cdump(1:1) = achar(timestamp/10**5 + 48)

  cfile=cme//"/parts_info."//cdump(1:6)


  open (60,file=cfile)    
  write(60,'(7(a9,i8/),7(a9,f12.5/))')  &    ! info block
       'itime=',timestamp, 'npp=',npp, &
       'ne=',ne, 'ni=',ni, 'npbeam=',np_beam, 'config=', initial_config, &
       'ensemble=',ensemble, &
       'xl=',xl, 'yl=',yl, 'zl=',zl, 'boxsize=',zl, &
       'eps=', eps, 'theta=',theta,'tlaser = ',tlaser   

  close(60)

  cfile=cme//"/parts_dump."//cdump(1:6)
  open (60,file=cfile) 
  write(60,'((12(1pe14.5),2i9))')  &
       (x(i), y(i), z(i), ux(i), uy(i), uz(i), q(i), m(i), &
        ax(i)*m(i)/q(i), ay(i)*m(i)/q(i), az(i)*m(i)/q(i), &  ! electric field = m.a/q
        pot(i), &  ! potential
        pepid(i), pelabel(i),i=1,npp)
  close(60)

  call MPI_BARRIER( MPI_COMM_WORLD, ierr)  ! Wait for everyone to catch up
  icall = icall + 1

end subroutine dump
