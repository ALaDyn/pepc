
! ======================
!
!   DRAW_DOMAINS
!
!   Output particles and tree structure for processing
!   by GLE 
!
! ======================

subroutine draw_domains(timestamp)


  use treevars
!VAMPINST include
      INCLUDE 'VTcommon.h'
      INTEGER VTIERR
!

  implicit none

  integer*8 :: key_twig(ntwig), key_leaf(nleaf)
  integer, dimension(ntwig) :: level_twig, node_twig       ! twig-nodes
  integer, dimension(nleaf) :: level_leaf, plist_leaf, ind_leaf       ! leaf-nodes

  character(30) :: cfile
  character(1) :: csnap
  character(3) :: cme
  character(6) :: cvisit
  character(7), parameter :: colors(0:9) = (/"orange ", "cyan   ", "magenta", "blue   ", "green  ", &
       "red    ","yellow ","grey20 ","violet ","brown  "/) 


  integer :: i, ip, j, ilev, isnap, ibt, ix, iy, iz, nbits
  real :: s, xt, yt, zt
  integer :: key2addr        ! Mapping function to get hash table address from key
  integer :: icall
  integer, intent(in) :: timestamp

!VAMPINST subroutine_start
       CALL VTENTER(IF_draw_domains,VTNOSCL,VTIERR)
!      write(*,*) 'VT: draw_domains S>',VTIERR,
!     *    IF_draw_domains,ICLASSH
!
  icall = timestamp/ivis

  !  Header file written out by root PE: does box and includes particle O/P from all PEs
  if ( me==0 ) then
     cfile="domains.gle"
     write (ipefile,'(a)') cfile
     open(61,file=cfile)


     !  initialise graphics filter

     write (61,'(a,3(/a),2(/a,2f13.4))') &
	  'size 18 18' &
	  , 'set font rm' &
	  , 'set lwidth 0.001 lstyle 1' &
	  , 'psize=0.04' &
	  , 'begin scale ',15./boxsize,15./boxsize &
	  , 'begin translate ',.01*boxsize,.01*boxsize

!     write (61,'(a/a,2f13.4)') 'amove 0 0','box ',boxsize,boxsize

     do i=0, num_pe-1
	cme = achar(i/100+48) // achar(mod(i/10,10)+48) // achar(mod(i,10)+48)  ! Convert 3-digit PE number into character string
	cfile="domain_"//cme//".gle"
	write (61,'(2a)') 'include ',cfile
     end do
     write (61,'(a/a)') 'end translate','end scale'
     close(61)

  endif


  !  Now do particles and boxes belonging to each processor domain

  cme = achar(me/100+48) // achar(mod(me/10,10)+48) // achar(mod(me,10)+48)  ! Convert 3-digit PE number into character string
  cfile="domain_"//cme//".gle"
  open (60,file=cfile) 

  ! get keys of twig nodes from hash table
  key_twig = pack(htable%key,mask=htable%node<0)

  ! get levels of twigs
  level_twig = log( 1.*key_twig )/log(2.**idim)
  node_twig = pack(htable%node,mask=htable%node<0)         ! twig index


  !  write (60,'(a/a,a7)') 'set lwidth .001','set color ',colors( mod(me,10) )
  write (60,'(a/a)') 'set lwidth .001','set color black'

  do j=2,ntwig
     ! get coords from keys
     nbits = level_twig(j)
     ix = SUM( (/ (2**i*ibits( key_twig(j),idim*i,1 ), i=0,nbits-1) /) )
     iy = SUM( (/ (2**i*ibits( key_twig(j),idim*i+1,1 ), i=0,nbits-1) /) )

     s = boxsize/2**(level_twig(j))          !  box length
     xt=ix*s + xmin
     yt=iy*s + ymin
     !write (60,'(a,2f13.4)') 'box ',s,s
     !     write (60,'(a,2f13.4)') 'amove ',xcoc( node_twig(j) ),ycoc( node_twig(j) )        ! Centre of charge of twig node
     !     write (60,'(a,f10.5,a)') 'circle ',.005*sqrt(abs_charge( node_twig(j) )),' fill cyan'
     !     write (60,'(a)') 'set hei 0.01'
     !     write (60,'(a,b10)') 'text ',key_twig(j)         ! write out key

  end do


  ! Branch nodes

  write (60,'(a)') 'set lwidth .001 color white'

 ! Dump domain data in VISIT format

  if (me.eq.0) then
     do i=0,4
	cvisit(6-i:6-i) =  achar(mod(icall/10**i,10) + 48)
     end do
     cvisit(1:1) = achar(icall/10**5 + 48)

     cfile="wf_domain."//cvisit
     open (61,file=cfile) 
     write(61,'(f12.5)')  dt*timestamp

  endif 

!  branch_owner = (/ ( htable( key2addr( branch_key(j) ) )%owner, j=1,nbranch_sum) /)         ! Owner-PE of branch

  do j=1,nbranch_sum
     ilev = log( 1.*branch_key(j) )/log(2.**idim)
     ix = SUM( (/ (2**i*ibits( branch_key(j),idim*i,1 ), i=0,ilev-1) /) )
     iy = SUM( (/ (2**i*ibits( branch_key(j),idim*i+1,1 ), i=0,ilev-1) /) )
     iz = SUM( (/ (2**i*ibits( branch_key(j),idim*i+2,1 ), i=0,ilev-1) /) )

     s = boxsize/2**(ilev)          !  box length
     xt=ix*s + xmin
     yt=iy*s + ymin
     zt=iz*s + zmin
     write (60,'(a)') 'set color white'
     write (60,'(a,2f13.4)') 'amove ',xt,yt
     if (branch_owner(j) == me ) then
        write (60,'(a,2f13.4,a,a7)') 'box ',s,s,' fill ',colors( mod(me,10) )
     else
!        write (60,'(a,2f13.4,a)') 'box ',s,s,' fill grey5'
     endif

     if (me==0) write (61,'(4f15.5,i6,2f15.5)') xt, yt, zt, s, branch_owner(j), s, s

  !     write (60,'(a)') 'set hei 0.01 color black'
  !     write (60,'(a,i5)') 'text ',branch_key(j)
  end do

  ! get keys of leaf nodes from hash table
  key_leaf(1:nleaf) = pack(htable%key,mask=htable%node>0)
  ind_leaf(1:nleaf) = pack(htable%node,mask=htable%node>0)         ! particle index
  plist_leaf(1:nleaf) = pack(htable%childcode,mask=htable%node>0)   ! particle label

  ! get levels of leaves
  level_leaf(1:nleaf) = log(1.*key_leaf(1:nleaf))/log(2.**idim)

  write (60,'(a)') 'set lwidth .001 color black'

  do j=1,nleaf
     ! get box coords from keys
     nbits = level_leaf(j)    ! # bits per ordinate
     ix = SUM( (/ (2**i*ibits( key_leaf(j),idim*i,1 ), i=0,nbits-1) /) )
     iy = SUM( (/ (2**i*ibits( key_leaf(j),idim*i+1,1 ), i=0,nbits-1) /) )

     s = boxsize/2**(level_leaf(j))          !  box length
     xt=ix*s + xmin
     yt=iy*s + ymin
     !     write (60,'(a/a)') 'set color black','set lwidth .002'
     write (60,'(a,2f13.4)') 'amove ',xt,yt
     !     write (60,'(a,2f13.4)') 'box ',s,s
     ! keys
     !     write (60,'(a,2f13.4)') 'amove ',xt,yt
     !     write (60,'(a)') 'set hei 0.01'
     !     write (60,'(a,b14)') 'text ',key_leaf(j)



  end do



  do i=1,npp
     write (60,'(a,2f13.4)') 'amove ',x(i),y(i)
     !     write (60,'(2a)') 'circle .002 fill ',colors( mod(me,8) )
     write (60,'(2a)') 'circle psize fill black'
  end do


  ! map of interleaved bits 


  ! first point
  xt=x(1)
  yt=y(1)
  !  write (60,'(a,2f13.4)') 'amove ',xt,yt

  do ip = 2,npp
     xt=x(ip)
     yt=y(ip)
     !     write (60,'(a,2f13.4)') 'aline ',xt,yt
  end do

  !  boundary links
  if (me /= lastpe) then
     !     write (60,'(a)') 'set color grey50'
     !     write (60,'(a,2f13.4)') 'aline ',x(npp+1),y(npp+1)
  endif



  close(60)
  close(61)
  icall = icall + 1

!VAMPINST subroutine_end
       CALL VTLEAVE(ICLASSH,VTIERR)
!      write(*,*) 'VT: draw_domains S<',VTIERR,ICLASSH
!
end subroutine draw_domains






