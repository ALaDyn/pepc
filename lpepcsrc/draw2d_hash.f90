
! ======================
!
!   DRAW_TREE_2D 
!
!   Output tree structure from hash table for processing
!   by GLE 
!
! ======================

subroutine draw2d_hash



  !  Tree written to ....  tree.box2dN
  !  Particles       ....  tree.partsN
  !  Centres         ....  tree.comasN


  use treevars

  implicit none

  integer*8 :: key_twig(ntwig), key_leaf(nleaf)

  integer, dimension(nppm) :: ix, iy

  integer, dimension(ntwig) :: level_twig, node_twig       ! twig-nodes
  integer, dimension(nleaf) :: level_leaf, plist_leaf, ind_leaf       ! leaf-nodes

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

  csnap=achar(mod(isnap,10)+48)


  cfile="box2d_"//csnap//".gle"
  write (ipefile,'(a)') cfile
  open(60,file=cfile)


  !  initialise graphics filter

  write (60,'(a,2(/a),2(/a,2f13.4))') &
       'size 15 15' &
       , 'set font psh hei .02' &
       , 'set lwidth 0.01' &
       , 'begin scale ',12./xl,12./xl &
       , 'begin translate ',.01*xl,.01*yl


  ! get keys of twig nodes from hash table
  key_twig = pack(htable%key,mask=htable%node<0)

  ! get levels of twigs
  level_twig = log( 1.*key_twig )/log(4.)
  node_twig = pack(htable%node,mask=htable%node<0)         ! twig index


  write (60,'(a/a)') 'set color blue','set lwidth .005'

  do j=1,ntwig
     ! get coords from keys
     nbits = level_twig(j)
     ix(j) = SUM( (/ (2**i*ibits( key_twig(j),3*i,1 ), i=0,nbits-1) /) )
     iy(j) = SUM( (/ (2**i*ibits( key_twig(j),3*i+1,1 ), i=0,nbits-1) /) )

     s = xl/2**(level_twig(j))          !  box length
     xt=ix(j)*s
     yt=iy(j)*s
     write (60,'(a,2f13.4)') 'amove ',xt,yt
     write (60,'(a,2f13.4)') 'box ',s,s
     write (60,'(a,2f13.4)') 'amove ',xcoc( node_twig(j) ),ycoc( node_twig(j) )        ! Centre of charge of twig node
     write (60,'(a,f10.5,a)') 'circle ',.005*sqrt(abs_charge( node_twig(j) )),' fill cyan'
     !     write (60,'(a)') 'set hei 0.03'
     !     write (60,'(a,b10)') 'text ',key_twig(j)         ! write out key



  end do


  ! get keys of twig nodes from hash table
  key_leaf(1:nleaf) = pack(htable%key,mask=htable%node>0)
  ind_leaf(1:nleaf) = pack(htable%node,mask=htable%node>0)         ! particle index
  plist_leaf(1:nleaf) = pack(htable%childcode,mask=htable%node>0)   ! particle label

  ! get levels of twigs
  level_leaf(1:nleaf) = log(1.*key_leaf(1:nleaf))/log(4.)



  write (60,'(a/a)') 'set color green','set lwidth .002'

  do j=1,nleaf
     ! get box coords from keys
     nbits = level_leaf(j)    ! # bits per ordinate
     ix(j) = SUM( (/ (2**i*ibits( key_leaf(j),3*i,1 ), i=0,nbits-1) /) )
     iy(j) = SUM( (/ (2**i*ibits( key_leaf(j),3*i+1,1 ), i=0,nbits-1) /) )

     s = xl/2**(level_leaf(j))          !  box length
     xt=ix(j)*s
     yt=iy(j)*s
     write (60,'(a/a)') 'set color black','set lwidth .002'
     write (60,'(a,2f13.4)') 'amove ',xt,yt
     write (60,'(a,2f13.4)') 'box ',s,s
     ! keys
     !     write (60,'(a,2f13.4)') 'amove ',xt+s/4,yt+s/2
     !     write (60,'(a)') 'set hei 0.02'
     !     write (60,'(a,b10)') 'text ',key_leaf(j)

     ! particle coords derived from leaf data

     !     write (60,'(a/a)') 'set color green','set lwidth .002'
     !     write (60,'(a,2f13.4)') 'amove ',xt,yt
     !     write (60,'(a,2f13.4)') 'circle .02'


  end do

  write (60,'(a/a)') 'set color red','set lwidth .002'


  do i=1,npp
     write (60,'(a,2f13.4)') 'amove ',x(i),y(i)
     if (q(i)<0) then
        write (60,'(a)') 'circle .005 fill red'
     else
        write (60,'(a)') 'circle .005 fill green'
     endif

  end do



  write (60,'(a/a)') 'end translate','end scale'
  close(60)


  isnap=isnap+1

end subroutine draw2d_hash
