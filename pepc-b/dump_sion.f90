! ======================
!
!   Sion DUMP
!
!   Gather and write out particle data for restart
!   using sionlib
!
! ======================

subroutine dump(timestamp)

  use physvars
  use treevars
  implicit none
  include 'mpif.h'

  character(30) :: cfile, cinfofile
  character(6) :: cdump
  character(5) :: cme
  integer, intent(in) :: timestamp
  integer :: i, j, ioffset, idummy=0, ierr
  integer, save :: icall=0
  real :: simtime

  integer*8 chunksize, size_info, size1, size2, bwrote
  ! 2MB GPFS
  integer 	fsblksize
  integer 	sid

  real*8 realdummy

! Temp variable use for the info blockt
  TYPE, BIND(C) :: Infoblock
  integer, dimension(7) :: intarr
  real, dimension(18) :: realarr
  real, dimension(3) :: plasma_centre, focus
  END TYPE Infoblock

  TYPE(Infoblock) :: infoblk

  infoblk%intarr(1)=timestamp
  infoblk%intarr(2)=npp
  infoblk%intarr(3)=ne
  infoblk%intarr(4)=ni
  infoblk%intarr(5)=np_beam
  infoblk%intarr(6)=target_geometry
  infoblk%intarr(7)=scheme

  infoblk%realarr(1)=xl
  infoblk%realarr(2)=yl
  infoblk%realarr(3)=zl
  infoblk%realarr(4)=eps
  infoblk%realarr(5)=theta
  infoblk%realarr(6)=tlaser
  infoblk%realarr(7)=trun
  infoblk%realarr(8)=omega
  infoblk%realarr(9)=lambda
  infoblk%realarr(10)=qe
  infoblk%realarr(11)=qi
  infoblk%realarr(12)=mass_e
  infoblk%realarr(13)=mass_i
  infoblk%realarr(14)=Zion
  infoblk%realarr(15)=a_ii
  infoblk%realarr(16)=Vplas
  infoblk%realarr(17)=Aplas
  infoblk%realarr(18)=Qplas
  do i=1,3
  	infoblk%plasma_centre(i) = plasma_centre(i)
  	infoblk%focus(i)=focus(i)
  end do

  simtime = dt*timestamp


  ! get filename suffix from dump counter
  do i=0,4
     cdump(6-i:6-i) =  achar(mod(timestamp/10**i,10) + 48)
  end do
  cdump(1:1) = achar(timestamp/10**5 + 48)

! Write particles info file
  if(me == 0) then
	  cinfofile="parts_info.in"
	  open (60,file=cinfofile)
	  write(60,'(7(a9,i8/),10(a9,f12.5/),9(a9,1pe12.5/),2(a9,3(1pe12.5)/))')  &    ! info block
		   'itime=',timestamp, 'npp=',npp, &
		   'ne=',ne, 'ni=',ni, 'npbeam=',np_beam, 'geometry=', target_geometry, &
		   'scheme=',scheme, &
		   'xl=',xl, 'yl=',yl, 'zl=',zl, 'boxsize=',zl, &
		   'eps=', eps, 'theta=',theta,' tlaser= ',tlaser,' trun= ',trun, &
		   'omega=',omega,'lambda=',lambda,'  qe=',qe,'  qi=',qi, &
		   'mass_e=',mass_e,'mass_i=',mass_i,'Zion=',Zion,'a_ii=',a_ii, &
		   'Vplas=',Vplas,'Aplas=',Aplas,'Qplas=',Qplas, &
		   'centre=',plasma_centre(1:3),'focus=',focus(1:3)
	  close (60)
  endif

! Calc the size of the chunk
! Info block: 7*sizeof(integer)+24*sizeof(real) = size(infoblock)
! Particles: 12*npp*sizeof(real*8) + 2*npp*sizeof(integer)
! should be:

  size_info = sizeof(infoblk%intarr)+sizeof(infoblk%realarr)+sizeof(infoblk%plasma_centre)+sizeof(infoblk%focus)

  ! npp * sizepof(real*8)
  size1 = npp*sizeof(realdummy)
  ! npp*sizeof(integer)
  size2 = npp*sizeof(sid)
  chunksize = 12*size1 + 2*size2 + size_info
  fsblksize = 2*1024*1024

!  if(me == 0) then
!    write(6,*)'**************** npp=',npp,' **************'
!    write(6,*)'**************** size_info=',size_info,' **************'
!    write(6,*)'**************** size1    =',size1,' **************'
!    write(6,*)'**************** size2    =',size2,' **************'
!   	write(6,*)'**************** chunksize=',chunksize,' **************'
!  end if

! Particles dump file
  cfile="data/parts_all_dump."//cdump(1:6)

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call fsion_paropen_mpi(trim(cfile),'bw',MPI_COMM_WORLD,chunksize,fsblksize,me,sid)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  ! Write the inforblock
  call fsion_write(infoblk, 1, size_info, sid, bwrote)
  if (bwrote /= size_info) then
   write(*,*) 'Error writing Infoblock. Wrote: ', bwrote
  else
  	write(167,*)'Wrote infoblock: ',bwrote
  end if

  ! Write X
  call fsion_write(x, 1, size1, sid, bwrote)
  if (bwrote /= size1) write(167,*) 'Error writing x(). Wrote: ', bwrote

  ! Write Y
  call fsion_write(y, 1, size1, sid, bwrote)
  if (bwrote /= size1) write(167,*) 'Error writing y(). Wrote: ', bwrote

  ! Write Z
  call fsion_write(z, 1, size1, sid, bwrote)
  if (bwrote /= size1) write(167,*) 'Error writing z(). Wrote: ', bwrote

  ! Write UX
  call fsion_write(ux, 1, size1, sid, bwrote)
  if (bwrote /= size1) write(167,*) 'Error writing ux(). Wrote: ', bwrote

  ! Write UY
  call fsion_write(uy, 1, size1, sid, bwrote)
  if (bwrote /= size1) write(167,*) 'Error writing uy(). Wrote: ', bwrote

  ! Write UY
  call fsion_write(uz, 1, size1, sid, bwrote)
  if (bwrote /= size1) write(167,*) 'Error writing uz(). Wrote: ', bwrote

  ! Write Q
  call fsion_write(q, 1, size1, sid, bwrote)
  if (bwrote /= size1) write(167,*) 'Error writing q(). Wrote: ', bwrote

  ! Write M
  call fsion_write(m, 1, size1, sid, bwrote)
  if (bwrote /= size1) write(167,*) 'Error writing m(). Wrote: ', bwrote

  ! Write Ex
  call fsion_write(Ex, 1, size1, sid, bwrote)
  if (bwrote /= size1) write(167,*) 'Error writing Ex(). Wrote: ', bwrote

  ! Write Ey
  call fsion_write(Ey, 1, size1, sid, bwrote)
  if (bwrote /= size1) write(167,*) 'Error writing Ey(). Wrote: ', bwrote

  ! Write Ez
  call fsion_write(Ez, 1, size1, sid, bwrote)
  if (bwrote /= size1) write(167,*) 'Error writing Ez(). Wrote: ', bwrote

  ! Write pot
  call fsion_write(pot, 1, size1, sid, bwrote)
  if (bwrote /= size1) write(167,*) 'Error writing pot(). Wrote: ', bwrote

  ! Write pepid
  call fsion_write(pepid, 1, size2, sid, bwrote)
  if (bwrote /= size2) write(167,*) 'Error writing pepid(). Wrote: ', bwrote

  ! Write pelabel
  call fsion_write(pelabel, 1, size2, sid, bwrote)
  if (bwrote /= size2) write(167,*) 'Error writing pelabel(). Wrote: ', bwrote

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call fsion_parclose_mpi(sid,MPI_COMM_WORLD,ierr)

! add to runstamp the stamp and print to stdout
  if (me.eq.0) then
    open (62,file="runstamp")  ! time stamp
    write(62,'(a)') cdump(1:6)
    close (62)

!    write(166,'(2a)') 'Particle dump',cdump(1:6)
!  	write(166,'(//a/7(a9,i8/),10(a9,f12.5/),9(a9,1pe12.5/),2(a9,3(1pe12.5)/))') 'PARTICLE DUMP:', &    ! info block
!       'itime=',timestamp, 'npp=',npp, &
!       'ne=',ne, 'ni=',ni, 'npbeam=',np_beam, 'geometry=', target_geometry, &
!       'scheme=',scheme, &
!       'xl=',xl, 'yl=',yl, 'zl=',zl, 'boxsize=',zl, &
!       'eps=', eps, 'theta=',theta,' tlaser= ',tlaser,' trun= ',trun, &
!       'omega=',omega,'lambda=',lambda,'  qe=',qe,'  qi=',qi, &
!       'mass_e=',mass_e,'mass_i=',mass_i,'Zion=',Zion,'a_ii=',a_ii, &
!       'Vplas=',Vplas,'Aplas=',Aplas,'Qplas=',Qplas, &
!       'centre=',plasma_centre(1:3),'focus=',focus(1:3)

  endif
  icall = icall + 1

end subroutine dump
