
! ======================
!
!   DRAW_TREE2D
!
!   Output particles and tree structure for processing
!   by GLE 
!
! ======================

subroutine draw_tree2d(xl,yl)


  use treevars

  implicit none

  real, intent(in) :: xl,yl
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

  ! get keys of twig nodes from hash table
  key_twig = pack(htable%key,mask=htable%node<0)

  ! get levels of twigs
  level_twig = log( 1.*key_twig )/log(8.)
  node_twig = pack(htable%node,mask=htable%node<0)         ! twig index
  owner_twig = pack(htable%owner,mask=htable%node<0)         ! twig owner



  do j=2,ntwig
     write (60,'(a/a,a7)') 'set hei .02','set color ',colors( mod(owner_twig(j),10) )
     ! get coords from keys
     nbits = level_twig(j)
     ix = SUM( (/ (2**i*ibits( key_twig(j),idim*i,1 ), i=0,nbits-1) /) )
     iy = SUM( (/ (2**i*ibits( key_twig(j),idim*i+1,1 ), i=0,nbits-1) /) )
 
     s = boxsize/2**(level_twig(j))          !  box length
     xt=ix*s + xmin
     yt=iy*s + ymin
     write (60,'(a)') 'set color black'
     write (60,'(a,2f13.4)') 'amove ',xt,yt
     if (owner_twig(j) == me) then
        write (60,'(a,2f13.4)') 'box ',s,s
     else
        write (60,'(a,2f13.4,2a)') 'box ',s,s,' fill ',colors( mod(owner_twig(j),10) )
     endif

     !        write (60,'(a,2f13.4)') 'amove ',xcoc( node_twig(j) ),ycoc( node_twig(j) )        ! Centre of charge of twig node
     !         write (60,'(a,f10.5,a)') 'circle ',.005*sqrt(abs_charge( node_twig(j) )),' fill grey'
     !        write (60,'(a)') 'set hei 0.01'
     !
  !   write (60,'(a)') 'set color black'
 !    write (60,'(a,2f13.4)') 'amove ',xt+s/2,yt+s/2
!     write (60,'(a,i6)') 'text ',key_twig(j)         ! write out key

  end do



  ! get keys of leaf nodes from hash table
  key_leaf(1:nleaf) = pack(htable%key,mask=htable%node>0)
  ind_leaf(1:nleaf) = pack(htable%node,mask=htable%node>0)         ! particle index
  owner_leaf(1:nleaf) = pack(htable%owner,mask=htable%node>0)         ! leaf owner
  plist_leaf(1:nleaf) = pack(htable%childcode,mask=htable%node>0)   ! particle label

  ! get levels of leaves
  level_leaf(1:nleaf) = log(1.*key_leaf(1:nleaf))/log(8.)


  do j=1,nleaf
     ! get box coords from keys
     nbits = level_leaf(j)    ! # bits per ordinate
     ix = SUM( (/ (2**i*ibits( key_leaf(j),idim*i,1 ), i=0,nbits-1) /) )
     iy = SUM( (/ (2**i*ibits( key_leaf(j),idim*i+1,1 ), i=0,nbits-1) /) )

     s = boxsize/2**(level_leaf(j))          !  box length
     xt=ix*s + xmin
     yt=iy*s + ymin
!         write (60,'(a,a)') 'set color ',colors( mod(owner_leaf(j),10) )
     ! keys
     !     write (60,'(a)') 'set hei 0.01'
    write (60,'(a)') 'set color black'

     write (60,'(a,2f13.4)') 'amove ',xt,yt
     if (owner_leaf(j) == me) then
        write (60,'(a,2f13.4)') 'box ',s,s
     else
        write (60,'(a,2f13.4,2a)') 'box ',s,s,' fill ',colors( mod(owner_leaf(j),10) )
 !       write (60,'(a,2f13.4)') 'amove ',xt+s/2,yt+s/2
!        write (60,'(a,i6)') 'text ',key_leaf(j)
     endif
     !     write (60,'(a/a)') 'set color green','set lwidth .002'
     !     write (60,'(a,2f13.4)') 'amove ',xt,yt
     !     write (60,'(a,2f13.4)') 'circle .02'


  end do



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
!    write (60,'(a,2f13.4)') 'amove ',xt,yt
!     write (60,'(a,a)') 'set color ',colors( mod(me,8) )
 !    write (60,'(a,a)') 'set lwidth 0.03'

  do ip = 2,npp
     xt=x(ip)
     yt=y(ip)
 !    write (60,'(a,2f13.4)') 'aline ',xt,yt
  end do

  !  boundary links
  if (me /= num_pe-1) then
     !     write (60,'(a)') 'set color grey50'
     !     write (60,'(a,2f13.4)') 'aline ',x(npp+1),y(npp+1)
  endif



  close(60)

end subroutine draw_tree2d
