
! ======================
!
!   DRAW_TREE_2D 
!
!   Output tree structure for processing
!   by GLE 
!
! ======================

subroutine draw2d



  !  Tree written to ....  tree.box2dN
  !  Particles       ....  tree.partsN
  !  Centres         ....  tree.comasN


  use treevars
!VAMPINST include
      INCLUDE 'VTcommon.h'
      INTEGER VTIERR
!

  implicit none

  integer, dimension(nppm) :: ix, iy

  character(30) :: cfile
  character(1) :: csnap
  integer ic(5),jc(5),lc(5)
  real octx(9),octy(9)

  integer :: i, ip, j, ilev, isnap, ibt
  real :: s, xt, yt

  save isnap
  data isnap/1/
  !  square data
  data ic/0,1,1,0,0/
  data jc/0,0,1,1,0/
  data lc/1,0,0,0,0/

  !  coma shape
  data octx/-.333, -1., -1., -.333, .333, 1.,  1., .333, 0./
  data octy/-1., -.333, .333, 1., 1., .333, -.333, -1., 0./

!VAMPINST subroutine_start
       CALL VTENTER(IF_draw2d,VTNOSCL,VTIERR)
!      write(*,*) 'VT: draw2d S>',VTIERR,
!     *    IF_draw2d,ICLASSH
!
  csnap=achar(mod(isnap,10)+48)


  cfile="box2d_"//csnap//".gle"
  write (15,'(a)') cfile
  open(60,file=cfile)

  !  initialise graphics filter

  write (60,'(a,2(/a),2(/a,2f13.4))') &
       'size 15 15' &
       , 'set font rm' &
       , 'set lwidth 0.01' &
       , 'begin scale ',12./xl,12./xl &
       , 'begin translate ',.01*xl,.01*yl

  write (60,'(a/a,2f13.4)') 'amove 0 0','box ',xl,yl


  ! parent cells - in blue, no symbol
  write (60,'(a)') 'set color yellow  lwidth .005'

  do ilev = nlev,1,-1
     s = xl/2**(ilev)          !  box length
     ! recover box coordinates of parents
     ibt = nlev-ilev           ! bit shift factor (0=particle node, nlev-1 = root)
     do j = 1,npp
        ix(j) = SUM( (/ (2**i*ibits( ishft( pekey(j),-2*ibt ),2*i,1 ), i=0,nbits-1-ibt) /) )
        iy(j) = SUM( (/ (2**i*ibits( ishft( pekey(j),-2*ibt ),2*i+1,1 ), i=0,nbits-1-ibt) /) )
     end do

     do ip = 1,npp
        xt=ix(ip)*s
        yt=iy(ip)*s
        write (60,'(a,2f13.4)') 'amove ',xt,yt
        write (60,'(a,2f13.4)') 'box ',s,s
     end do

  end do

  write (60,'(a/a)') 'set color red','set lwidth .005'


  do i=1,npp
     write (60,'(a,2f13.4)') 'amove ',x(i),y(i)
     write (60,'(a)') 'circle .002 fill red'
  end do

  ! map of interleaved bits 
  write (60,'(a)') 'set color blue lwidth .002'

  ilev = nlev
  s = xl/2**(ilev)          !  box length
  ! recover box coordinates of parents
  ibt = nlev-ilev           ! bit shift factor (0=particle node, nlev-1 = root)
  do j = 1,npp
     ix(j) = SUM( (/ (2**i*ibits( ishft( pekey(j),-2*ibt ),2*i,1 ), i=0,nbits-1-ibt) /) )
     iy(j) = SUM( (/ (2**i*ibits( ishft( pekey(j),-2*ibt ),2*i+1,1 ), i=0,nbits-1-ibt) /) )
  end do

! first point
  xt=ix(1)*s+s/2
  yt=iy(1)*s+s/2
  write (60,'(a,2f13.4)') 'amove ',xt,yt

  do ip = 2,npp
     xt=ix(ip)*s+s/2
     yt=iy(ip)*s+s/2
     write (60,'(a,2f13.4)') 'aline ',xt,yt
  end do


  write (60,'(a/a)') 'end translate','end scale'
  close(60)


  isnap=isnap+1

!VAMPINST subroutine_end
       CALL VTLEAVE(ICLASSH,VTIERR)
!      write(*,*) 'VT: draw2d S<',VTIERR,ICLASSH
!
end subroutine draw2d
