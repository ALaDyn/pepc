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
  use physvars

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

  ! Filename (directory) prefix
  cme = "pe"// achar(me/100+48) // achar(mod(me/10,10)+48) // achar(mod(me,10)+48)  

  ! get filename suffix from dump counter
  do i=0,4
     cdump(6-i:6-i) =  achar(mod(timestamp/10**i,10) + 48)  
  end do
  cdump(1:1) = achar(timestamp/10**5 + 48)



  if (me.eq.0) then

     ! Write general header with global run parameters

     open (62,file="runstamp")  ! time stamp 
     write(62,'(a)') cdump(1:6)
     close (62)

     cfile="parts_info."//cdump(1:6)

     open (60,file=cfile)    
     write(60,'(i12,a20)') timestamp,'!itime'
     write(60,'(i12,a20)') npart,      '!npart  '
     write(60,'(i12,a20)') ne,       '!ndust   '
     write(60,'(i12,a20)') ni,       '!nstar   '
     write(60,'(i12,a20)') scheme,   '!scheme'
     write(60,'(f12.3,a20)') mdisc,  '!mdisc'
     write(60,'(f12.3,a20)') xl,     '!xl'
     write(60,'(f12.3,a20)') yl,     '!yl'
     write(60,'(f12.3,a20)') zl,     '!zl'
     write(60,'(f12.3,a20)') boxsize,'!boxsize'
     write(60,'(3f12.3,a20)') box_centre,'!shift vector'
     write(60,'(f12.3,a20)') eps,    '!eps'
     write(60,'(f12.3,a20)') theta,  '!theta'
     write(60,'(i12,a20)') nmerge,   '!nmerge'
     write(60,'(i12,a20)') nt,   '!nt'
     write(60,'(f12.3,a20)') dt,  '!timestep'
     write(60,'(i12,a20)') ivis,   '!particle vis'
     write(60,'(i12,a20)') ivis_fields,   '!field vis'
     write(60,'(i12,a20)') idump,   '!dump frequency'
     write(60,'(i12,a20)') iprot,   '!written o/p'


     write(6,*)'PARTICLE DUMP:'
     write(6,'(i8,a20)')timestamp,'itime'
     write(6,'(i8,a20)')npart,      'npart'
     write(6,'(i8,a20)')npp,      'npp'
     write(6,'(i8,a20)')ne,       'ndust'
     write(6,'(i8,a20)')ni,       'nstar'


  endif


  ! Write particle header local CPU

  cfile=cme//"/parts_info."//cdump(1:6)

  open (60,file=cfile)    
  write(60,'(i8,a20)')timestamp,'     !itime'
  write(60,'(i8,a18)')npp,      '     !npp  '
  write(60,'(i8,a17)')ne,       '     !ne   '
  write(60,'(i8,a17)')ni,       '     !ni   '
  close(60)

  !  Dump particles

  cfile=cme//"/parts_dump."//cdump(1:6)
  open (60,file=cfile) 
  write(60,'((12(1pe14.5),2i9))')  &
       (x(i), y(i), z(i), ux(i), uy(i), uz(i), q(i), m(i), &
       ax(i), ay(i), az(i), &  ! electric field = m.a/q
       pot(i), &  ! potential
       pepid(i), pelabel(i),i=1,npp)
  close(60)


  !  Dump stars

  if (me==0) then
     cfile="star_info."//cdump(1:6)
     open (81,file=cfile)
!     write(81,'(i12,a20)') timestamp,'!itime'
 
     do i=1,ni
        write(81,'(f15.5,10x,a20,i2)') m_star(i),  '  ! mass of star ',i
     end do
     do i=1,ni
        write(81,'(3f15.5,a20,i2)')x_star(i),y_star(i),z_star(i),'  ! position of star ',i
     end do
     do i=1,ni
        write(81,'(3f15.5,a20,i2)')ux_star(i),uy_star(i),uz_star(i),'  ! velocity of star ',i
     end do

  endif

  close(81)



  call MPI_BARRIER( MPI_COMM_WORLD, ierr)  ! Wait for everyone to catch up
  icall = icall + 1

end subroutine dump






