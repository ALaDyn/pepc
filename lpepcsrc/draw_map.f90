
! ======================
!
!   DRAW_MAP
!
!   Output particles and space-filling curve for postprocessing
!   by GLE 
!
! ======================

subroutine draw_map


  use treevars

  implicit none

  integer*8 :: key_twig(ntwig), key_leaf(nleaf)

  integer, dimension(ntwig) :: level_twig, node_twig, owner_twig       ! twig-nodes
  integer, dimension(nleaf) :: level_leaf, plist_leaf, ind_leaf, owner_leaf       ! leaf-nodes

  character(30) :: cfile
  character(1) :: csnap
  character(3) :: cme
  character(7), parameter :: colors(0:9) = (/"orange ", "cyan   ", "magenta", "blue   ", "green  ", &
       "red    ","yellow ","grey20 ","violet ","brown  "/) 


  integer :: i, ip, j, ilev, isnap, ibt, ix, iy, nbits
  real :: s, xt, yt


  csnap=achar(mod(me,10)+48)


  !  Header file written out by root PE: does box and includes particle O/P from all PEs
  if ( me==0 ) then
     cfile="tree2d.gle"
     write (ipefile,'(a)') cfile
     open(60,file=cfile)

     !  initialise graphics filter

     write (60,'(a,3(/a),2(/a,2f13.4))') &
	  'size 18 18' &
	  , 'set font rm' &
	  , 'set lwidth 0.05 lstyle 1' &
	  , 'psize=0.04' &
	  , 'begin scale ',15./boxsize,15./boxsize &
	  , 'begin translate ',.01*boxsize,.01*yl

     write (60,'(a/a,2f13.4)') 'amove 0 0','box ',boxsize,yl

     do i=0, num_pe-1
	cme = achar(i/100+48) // achar(i/10+48) // achar(mod(i,10)+48)  ! Convert 3-digit PE number into character string
	cfile="tree_"//cme//".gle"
	write (60,'(2a)') 'include ',cfile
     end do
     write (60,'(a/a)') 'end translate','end scale'
     close(60)

  endif


  !  Now do particles and boxes belonging to each processor domain

  cme = achar(me/100+48) // achar(me/10+48) // achar(mod(me,10)+48)  ! Convert 3-digit PE number into character string
  cfile="tree_"//cme//".gle"
  open (60,file=cfile) 


  do i=1,npp
     write (60,'(a,a)') 'set color ',colors( mod(me,8) )

     write (60,'(a,2f13.4)') 'amove ',x(i),y(i)
     write (60,'(2a)') 'circle psize fill ',colors( mod(me,8))
!     write (60,'(2a)') 'circle psize fill ','black'
  end do


  ! map of interleaved bits 

  ilev = nlev
  s = boxsize/2**(ilev)          !  box length
  ! recover box coordinates of parents
  ibt = nlev-ilev           ! bit shift factor (0=particle node, nlev-1 = root)
  do j = 1,npp
     ix = SUM( (/ (2**i*ibits( ishft( pekey(j),-2*ibt ),2*i,1 ), i=0,nbits-2-ibt) /) )
     iy = SUM( (/ (2**i*ibits( ishft( pekey(j),-2*ibt ),2*i+1,1 ), i=0,nbits-2-ibt) /) )
  end do


  ! first point
  xt=x(1)
  yt=y(1)
  write (60,'(a,2f13.4)') 'amove ',xt,yt
  write (60,'(a,a)') 'set color ',colors( mod(me,8) )
  write (60,'(a,a)') 'set lwidth 0.03'

  do ip = 2,npp
     xt=x(ip)
     yt=y(ip)
     write (60,'(a,2f13.4)') 'aline ',xt,yt
  end do

  !  boundary links
  if (me /= lastpe) then
          write (60,'(a)') 'set lstyle 2'
          write (60,'(a,2f13.4)') 'aline ',x(npp+1),y(npp+1)
  endif



  close(60)

end subroutine draw_map
