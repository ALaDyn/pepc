!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!> Helper functions for gle postprocessing
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_gle
  implicit none
  private

    public draw_lists
    public draw_map
    public draw_domains
    public draw_tree2d
    public draw2d
    public draw2d_hash

  contains

	! ======================
	!
	!   DRAW_LISTS
	!
	!   Draw selected particle lists
	!
	! ======================

	subroutine draw_lists


	  use treevars
	  use module_htable
	  use module_pepc_wrappers
	  use module_spacefilling
	  use module_debug, only : ipefile

	  implicit none

	  integer*8 :: key_list(npp)

	  integer, dimension(npp) ::  addr_list, level_list, node_list, owner_list       ! list data

	  character(15) :: cfile
	  character(3) :: cme,cip
	  character(7), parameter :: colors(0:9) = (/"orange ", "cyan   ", "magenta", "blue   ", "green  ", &
	       "red    ","yellow ","grey20 ","violet ","brown  "/)


	  integer :: i, ip, j, ix, iy, nbits
	  real*8 :: xt, yt, scale, s
	  logical :: write_keys=.false.


	!  do ip=1,npp,npp/2
	  ip=1
	  !  Header file
	  cme = achar(me/100+48) // achar(me/10+48) // achar(mod(me,10)+48)  ! Convert 3-digit PE number into character string
	  cip = achar(ip/100+48) // achar(ip/10+48) // achar(mod(ip,10)+48)  ! Convert 3-digit particle number into character string
	  cfile="list_"//cme//cip//".gle"
	  write (ipefile,'(a)') cfile
	 open(60,file=cfile)

	  !  initialise graphics filter

	!  write (60,'(a,2(/a),2(/a,2f13.4))') &
	!       'size 18 18' &
	!       , 'set font rm' &
	!       , 'set lwidth 0.001 lstyle 1' &
	!       , 'begin scale ',15./boxsize,15./boxsize &
	!       , 'begin translate ',.01*boxsize,.01*yl

	!  write (60,'(a/a,2f13.4)') 'amove 0 0','box ',boxsize,yl

	  write (60,'(a)')  'psize=0.1'

	  scale = 15./boxsize
	  write (60,'(a/a,2f12.2)') 'set lwidth .01 hei .02', &
			     'begin scale ',scale,scale
	  write (60,'(a)') 'set color black hei .3'
	  write (60,'(a,2f13.4)') 'amove ',particles(ip)%x(1),particles(ip)%x(2)
	  write (60,'(a,f10.5,2a)') 'marker otimes '
	  write (60,'(a,2f13.4)') 'rmove ',.1,.1

	  write (60,'(a,i6)') 'text ',particles(ip)%label         ! write out label


	  !  Now do selected particle lists
	  !  TODO need optional reinstatement of interaction lists for diagnostic purposes

!	  nlist = nterm(ip)
	  nlist = 1
	  ! get keys of twig nodes from hash table
!	  key_list(1:nlist) = intlist(1:nlist,ip)
	  key_list(1) = 1  ! dummy
	  addr_list(1:nlist) = (/ ( key2addr( key_list(i),'DRAW_LISTS' ),i=1,nlist) /)   !  Table address

	  ! get levels of twigs
	  level_list(1:nlist) = level_from_key(key_list(1:nlist))
	  node_list(1:nlist) =   htable( addr_list(1:nlist))%node   !  nodes
	  owner_list(1:nlist) =   htable( addr_list(1:nlist))%owner      ! owner




	  do j=1,nlist
	    write (60,'(a,a7)') 'set color ',colors( mod(owner_list(j),10) )
	     ! get box coords from keys
	     nbits = level_list(j)
	     ix = int(SUM( (/ (2**i*ibits( key_list(j),3*i,1 ), i=0,nbits-1) /) ))
	     iy = int(SUM( (/ (2**i*ibits( key_list(j),3*i+1,1 ), i=0,nbits-1) /) ))

	     s = boxsize/2**(level_list(j))          !  box length
	     xt=ix*s + xmin
	     yt=iy*s + ymin
	     write (60,'(a,2f13.4)') 'amove ',xt,yt
	     write (60,'(a,2f13.4)') 'box ',s,s

	     write (60,'(a,2f13.4)') 'amove ',tree_nodes( node_list(j) )%coc(1),tree_nodes( node_list(j) )%coc(2)     ! Centre of charge of twig node
	     write (60,'(a,f10.5,2a)') 'circle ',(tree_nodes( node_list(j) )%abs_charge)**(.33), &
	     		'*psize fill ',colors( mod(owner_list(j),10) )
	    ! write (60,'(a,f10.5,2a)') 'circle ',psize,' fill ',colors( mod(owner_list(j),10) )
	     !
	     write (60,'(a)') 'set color black'
	     write (60,'(a,2f13.4)') 'amove ',xt+s/2,yt+s/2
	     if (write_keys) write (60,'(a,i6)') 'text ',key_list(j)  ! write out key

	  end do

	  write (60,'((/a))') 'end scale '

	  close(60)
	! end do
	end subroutine draw_lists


	! ======================
	!
	!   DRAW_MAP
	!
	!   Output particles and space-filling curve for postprocessing
	!   by GLE
	!
	! ======================

	subroutine draw_map(yl)

	  use treevars
	  use module_pepc_wrappers
      use module_debug, only : ipefile

	  implicit none

	  real*8, intent(in) :: yl

	  character(30) :: cfile
	  character(1) :: csnap
	  character(3) :: cme
	  character(7), parameter :: colors(0:9) = (/"orange ", "cyan   ", "magenta", "blue   ", "green  ", &
	       "red    ","yellow ","grey20 ","violet ","brown  "/)


	  integer :: i, ip, j, ilev, ibt, ix, iy, nbits
	  real*8 :: s, xt, yt


	  csnap=achar(mod(me,10)+48)
      nbits=nlev+1

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

	     write (60,'(a,2f13.4)') 'amove ',particles(i)%x(1),particles(i)%x(2)
	     write (60,'(2a)') 'circle psize fill ',colors( mod(me,8))
	!     write (60,'(2a)') 'circle psize fill ','black'
	  end do


	  ! map of interleaved bits

	  ilev = nlev
	  s = boxsize/2**(ilev)          !  box length
	  ! recover box coordinates of parents
	  ibt = nlev-ilev           ! bit shift factor (0=particle node, nlev-1 = root)
	  do j = 1,npp
	     ix = int(SUM( (/ (2**i*ibits( ishft( particles(j)%key,-2*ibt ),2*i  ,1 ), i=0,nbits-2-ibt) /) ))
	     iy = int(SUM( (/ (2**i*ibits( ishft( particles(j)%key,-2*ibt ),2*i+1,1 ), i=0,nbits-2-ibt) /) ))
	  end do


	  ! first point
	  xt=particles(1)%x(1)
	  yt=particles(1)%x(2)
	  write (60,'(a,2f13.4)') 'amove ',xt,yt
	  write (60,'(a,a)') 'set color ',colors( mod(me,8) )
	  write (60,'(a,a)') 'set lwidth 0.03'

	  do ip = 2,npp
	     xt=particles(ip)%x(1)
	     yt=particles(ip)%x(2)
	     write (60,'(a,2f13.4)') 'aline ',xt,yt
	  end do

	  !  boundary links
	  if (me /= num_pe-1) then
	          write (60,'(a)') 'set lstyle 2'
	          write (60,'(a,2f13.4)') 'aline ',particles(npp+1)%x(1),particles(npp+1)%x(2)
	  endif



	  close(60)

	end subroutine draw_map


	! ======================
	!
	!   DRAW_TREE2D
	!
	!   Output particles and tree structure for processing
	!   by GLE
	!
	! ======================

	subroutine draw_tree2d(yl)


	  use treevars
	  use module_htable
      use module_pepc_wrappers
      use module_spacefilling
      use module_debug, only : ipefile

	  implicit none

	  real, intent(in) :: yl
	  integer*8 :: key_twig(ntwig), key_leaf(nleaf)

	  integer, dimension(ntwig) :: level_twig, node_twig, owner_twig       ! twig-nodes
	  integer, dimension(nleaf) :: level_leaf, plist_leaf, ind_leaf, owner_leaf       ! leaf-nodes

	  character(30) :: cfile
	  character(1) :: csnap
	  character(3) :: cme
	  character(7), parameter :: colors(0:9) = (/"orange ", "cyan   ", "magenta", "blue   ", "green  ", &
	       "red    ","yellow ","grey20 ","violet ","brown  "/)


	  integer :: i, ip, j, ilev, ibt, ix, iy, nbits
	  real*8 :: s, xt, yt


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
	  level_twig = level_from_key(key_twig)
	  node_twig = pack(htable%node,mask=htable%node<0)         ! twig index
	  owner_twig = pack(htable%owner,mask=htable%node<0)         ! twig owner



	  do j=2,ntwig
	     write (60,'(a/a,a7)') 'set hei .02','set color ',colors( mod(owner_twig(j),10) )
	     ! get coords from keys
	     nbits = level_twig(j)
	     ix = int(SUM( (/ (2**i*ibits( key_twig(j),3*i,1 ), i=0,nbits-1) /) ))
	     iy = int(SUM( (/ (2**i*ibits( key_twig(j),3*i+1,1 ), i=0,nbits-1) /) ))

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

	     !        write (60,'(a,2f13.4)') 'amove ',tree_nodes( node_twig(j) )%coc(1),tree_nodes( node_twig(j) )%coc(2)    ! Centre of charge of twig node
	     !         write (60,'(a,f10.5,a)') 'circle ',.005*sqrt(tree_nodes( node_twig(j) )%abs_charge),' fill grey'
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
	  level_leaf(1:nleaf) = level_from_key(key_leaf(1:nleaf))


	  do j=1,nleaf
	     ! get box coords from keys
	     nbits = level_leaf(j)    ! # bits per ordinate
	     ix = int(SUM( (/ (2**i*ibits( key_leaf(j),3*i,1 ), i=0,nbits-1) /) ))
	     iy = int(SUM( (/ (2**i*ibits( key_leaf(j),3*i+1,1 ), i=0,nbits-1) /) ))

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

	     write (60,'(a,2f13.4)') 'amove ',particles(i)%x(1),particles(i)%x(2)
	     write (60,'(2a)') 'circle psize fill ',colors( mod(me,8))
	!     write (60,'(2a)') 'circle psize fill ','black'
	  end do


	  ! map of interleaved bits

	  ilev = nlev
	  s = boxsize/2**(ilev)          !  box length
	  ! recover box coordinates of parents
	  ibt = nlev-ilev           ! bit shift factor (0=particle node, nlev-1 = root)
	  do j = 1,npp
	     ix = int(SUM( (/ (2**i*ibits( ishft( particles(j)%key,-2*ibt ),2*i  ,1 ), i=0,nbits-2-ibt) /) ))
	     iy = int(SUM( (/ (2**i*ibits( ishft( particles(j)%key,-2*ibt ),2*i+1,1 ), i=0,nbits-2-ibt) /) ))
	  end do


	  ! first point
	  xt=particles(1)%x(1)
	  yt=particles(1)%x(2)
	!    write (60,'(a,2f13.4)') 'amove ',xt,yt
	!     write (60,'(a,a)') 'set color ',colors( mod(me,8) )
	 !    write (60,'(a,a)') 'set lwidth 0.03'

	  do ip = 2,npp
	     xt=particles(ip)%x(1)
	     yt=particles(ip)%x(2)
	 !    write (60,'(a,2f13.4)') 'aline ',xt,yt
	  end do

	  !  boundary links
	  if (me /= num_pe-1) then
	     !     write (60,'(a)') 'set color grey50'
	     !     write (60,'(a,2f13.4)') 'aline ',particles(npp+1)%x(1),particles(npp+1)%x(2)
	  endif



	  close(60)

	end subroutine draw_tree2d



	! ======================
	!
	!   DRAW_TREE_2D
	!
	!   Output tree structure from hash table for processing
	!   by GLE
	!
	! ======================

	subroutine draw2d_hash(xl, yl)
	  !  Tree written to ....  tree.box2dN
	  !  Particles       ....  tree.partsN
	  !  Centres         ....  tree.comasN

	  use treevars
	  use module_htable
      use module_pepc_wrappers
      use module_debug, only : ipefile

	  implicit none
	  real*8, intent(in) :: xl, yl

	  integer*8 :: key_twig(ntwig), key_leaf(nleaf)

	  integer, dimension(npp) :: ix, iy

	  integer, dimension(ntwig) :: level_twig, node_twig       ! twig-nodes
	  integer, dimension(nleaf) :: level_leaf, plist_leaf, ind_leaf       ! leaf-nodes

	  character(30) :: cfile
	  character(1) :: csnap
	  integer ic(5),jc(5),lc(5)

	  integer :: i, j, isnap, nbits
	  real*8 :: s, xt, yt

	  save isnap
	  data isnap/1/
	  !  square data
	  data ic/0,1,1,0,0/
	  data jc/0,0,1,1,0/
	  data lc/1,0,0,0,0/

	  csnap=achar(mod(isnap,10)+48)
	  nbits=nlev+1

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
	  level_twig = int(log( 1.*key_twig )/log(4.))
	  node_twig = pack(htable%node,mask=htable%node<0)         ! twig index


	  write (60,'(a/a)') 'set color blue','set lwidth .005'

	  do j=1,ntwig
	     ! get coords from keys
	     nbits = level_twig(j)
	     ix(j) = int(SUM( (/ (2**i*ibits( key_twig(j),3*i,1 ), i=0,nbits-1) /) ))
	     iy(j) = int(SUM( (/ (2**i*ibits( key_twig(j),3*i+1,1 ), i=0,nbits-1) /) ))

	     s = xl/2**(level_twig(j))          !  box length
	     xt=ix(j)*s
	     yt=iy(j)*s
	     write (60,'(a,2f13.4)') 'amove ',xt,yt
	     write (60,'(a,2f13.4)') 'box ',s,s
	     write (60,'(a,2f13.4)') 'amove ',tree_nodes( node_twig(j) )%coc(1),tree_nodes( node_twig(j) )%coc(2)      ! Centre of charge of twig node
	     write (60,'(a,f10.5,a)') 'circle ',.005*sqrt(tree_nodes( node_twig(j) )%abs_charge),' fill cyan'
	     !     write (60,'(a)') 'set hei 0.03'
	     !     write (60,'(a,b10)') 'text ',key_twig(j)         ! write out key



	  end do


	  ! get keys of twig nodes from hash table
	  key_leaf(1:nleaf) = pack(htable%key,mask=htable%node>0)
	  ind_leaf(1:nleaf) = pack(htable%node,mask=htable%node>0)         ! particle index
	  plist_leaf(1:nleaf) = pack(htable%childcode,mask=htable%node>0)   ! particle label

	  ! get levels of twigs
	  level_leaf(1:nleaf) = int(log(1.*key_leaf(1:nleaf))/log(4.))



	  write (60,'(a/a)') 'set color green','set lwidth .002'

	  do j=1,nleaf
	     ! get box coords from keys
	     nbits = level_leaf(j)    ! # bits per ordinate
	     ix(j) = int(SUM( (/ (2**i*ibits( key_leaf(j),3*i,1 ), i=0,nbits-1) /) ))
	     iy(j) = int(SUM( (/ (2**i*ibits( key_leaf(j),3*i+1,1 ), i=0,nbits-1) /) ))

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
	     write (60,'(a,2f13.4)') 'amove ',particles(i)%x(1),particles(i)%x(2)
	     if (particles(i)%data%q<0) then
	        write (60,'(a)') 'circle .005 fill red'
	     else
	        write (60,'(a)') 'circle .005 fill green'
	     endif

	  end do



	  write (60,'(a/a)') 'end translate','end scale'
	  close(60)


	  isnap=isnap+1

	end subroutine draw2d_hash

	! ======================
	!
	!   DRAW_TREE_2D
	!
	!   Output tree structure for processing
	!   by GLE
	!
	! ======================

	subroutine draw2d(xl, yl)

	  !  Tree written to ....  tree.box2dN
	  !  Particles       ....  tree.partsN
	  !  Centres         ....  tree.comasN

	  use treevars
      use module_pepc_wrappers

	  implicit none

	  real*8, intent(in) :: xl, yl
	  integer, dimension(npp) :: ix, iy

	  character(30) :: cfile
	  character(1) :: csnap
	  integer ic(5),jc(5),lc(5)
	  real octx(9),octy(9)

	  integer :: i, ip, j, ilev, isnap, ibt, nbits
	  real*8 :: s, xt, yt

	  save isnap
	  data isnap/1/
	  !  square data
	  data ic/0,1,1,0,0/
	  data jc/0,0,1,1,0/
	  data lc/1,0,0,0,0/

	  !  coma shape
	  data octx/-.333, -1., -1., -.333, .333, 1.,  1., .333, 0./
	  data octy/-1., -.333, .333, 1., 1., .333, -.333, -1., 0./

	  csnap=achar(mod(isnap,10)+48)
	  nbits=nlev+1

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
	        ix(j) = int(SUM( (/ (2**i*ibits( ishft( particles(j)%key,-2*ibt ),2*i  ,1 ), i=0,nbits-1-ibt) /) ))
	        iy(j) = int(SUM( (/ (2**i*ibits( ishft( particles(j)%key,-2*ibt ),2*i+1,1 ), i=0,nbits-1-ibt) /) ))
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
	     write (60,'(a,2f13.4)') 'amove ',particles(i)%x(1),particles(i)%x(2)
	     write (60,'(a)') 'circle .002 fill red'
	  end do

	  ! map of interleaved bits
	  write (60,'(a)') 'set color blue lwidth .002'

	  ilev = nlev
	  s = xl/2**(ilev)          !  box length
	  ! recover box coordinates of parents
	  ibt = nlev-ilev           ! bit shift factor (0=particle node, nlev-1 = root)
	  do j = 1,npp
	     ix(j) = int(SUM( (/ (2**i*ibits( ishft( particles(j)%key,-2*ibt ),2*i  ,1 ), i=0,nbits-1-ibt) /) ))
	     iy(j) = int(SUM( (/ (2**i*ibits( ishft( particles(j)%key,-2*ibt ),2*i+1,1 ), i=0,nbits-1-ibt) /) ))
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

	end subroutine draw2d


	! ======================
	!
	!   DRAW_DOMAINS
	!
	!   Output particles and tree structure for processing
	!   by GLE
	!
	! ======================

	subroutine draw_domains()


	  use treevars
	  use module_htable
      use module_pepc_wrappers
      use module_spacefilling
      use module_debug, only : ipefile

	  implicit none

	  integer*8 :: key_twig(ntwig), key_leaf(nleaf)
	  integer, dimension(ntwig) :: level_twig, node_twig       ! twig-nodes
	  integer, dimension(nleaf) :: level_leaf, plist_leaf, ind_leaf       ! leaf-nodes

	  character(30) :: cfile
	  character(3) :: cme
	  character(15), parameter :: colors(0:9) = (/"red           ", &
	                                              "bright_gold   ", &
	                                              "lime_green    ", &
	                                              "medium_blue   ", &
	                                              "magenta       ", &
	                                              "brass         ", &
	                                              "orange        ", &
	                                              "aqua          ", &
	                                              "violet_red    ", &
	                                              "yellow_green  "/)


	  integer :: i, ip, j, ilev, ix, iy, iz, nbits
	  real*8 :: s, xt, yt, zt

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
	      , 'begin translate ',.05*boxsize,.05*boxsize

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
	  level_twig = level_from_key(key_twig)
	  node_twig = pack(htable%node,mask=htable%node<0)         ! twig index


	  !  write (60,'(a/a,a15)') 'set lwidth .001','set color ',colors( mod(me,10) )
	  write (60,'(a/a)') 'set lwidth .001','set color black'

	  do j=2,ntwig
	     ! get coords from keys
	     nbits = level_twig(j)
	     ix = int(SUM( (/ (2**i*ibits( key_twig(j),3*i,1 ), i=0,nbits-1) /) ))
	     iy = int(SUM( (/ (2**i*ibits( key_twig(j),3*i+1,1 ), i=0,nbits-1) /) ))

	     s = boxsize/2**(level_twig(j))          !  box length
	     xt=ix*s + xmin
	     yt=iy*s + ymin
	     !write (60,'(a,2f13.4)') 'box ',s,s
	     !     write (60,'(a,2f13.4)') 'amove ',tree_nodes( node_twig(j) )%coc(1),tree_nodes( node_twig(j) )%coc(2)        ! Centre of charge of twig node
	     !     write (60,'(a,f10.5,a)') 'circle ',.005*sqrt(tree_nodes( node_twig(j) )%abs_charge),' fill cyan'
	     !     write (60,'(a)') 'set hei 0.01'
	     !     write (60,'(a,b10)') 'text ',key_twig(j)         ! write out key

	  end do


	  ! Branch nodes

	  write (60,'(a)') 'set lwidth .02 color white'


	!  branch_owner = (/ ( htable( key2addr( branch_key(j) ) )%owner, j=1,nbranch_sum) /)         ! Owner-PE of branch

	  do j=1,nbranch_sum
	     ilev = level_from_key(branch_key(j))
	     ix = int(SUM( (/ (2**i*ibits( branch_key(j),3*i,1 ), i=0,ilev-1) /) ))
	     iy = int(SUM( (/ (2**i*ibits( branch_key(j),3*i+1,1 ), i=0,ilev-1) /) ))
	     iz = int(SUM( (/ (2**i*ibits( branch_key(j),3*i+2,1 ), i=0,ilev-1) /) ))

	     s = boxsize/2**(ilev)          !  box length
	     xt=ix*s + xmin
	     yt=iy*s + ymin
	     zt=iz*s + zmin
	     write (60,'(a)') 'set color white'
	     write (60,'(a,2f13.4)') 'amove ',xt,yt
	     if (branch_owner(j) == me ) then
	        write (60,'(a,2f13.4,a,a15)') 'box ',s,s,' fill ',colors( mod(me,10) )
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
	  level_leaf(1:nleaf) = level_from_key(key_leaf(1:nleaf))

	  write (60,'(a)') 'set lwidth .001 color black'

	  do j=1,nleaf
	     ! get box coords from keys
	     nbits = level_leaf(j)    ! # bits per ordinate
	     ix = int(SUM( (/ (2**i*ibits( key_leaf(j),3*i,1 ), i=0,nbits-1) /) ))
	     iy = int(SUM( (/ (2**i*ibits( key_leaf(j),3*i+1,1 ), i=0,nbits-1) /) ))

	     s = boxsize/2**(level_leaf(j))          !  box length
	     xt=ix*s + xmin
	     yt=iy*s + ymin
	     !     write (60,'(a/a)') 'set color black','set lwidth .002'
	!     write (60,'(a,2f13.4)') 'amove ',xt,yt
	     !     write (60,'(a,2f13.4)') 'box ',s,s
	     ! keys
	     !     write (60,'(a,2f13.4)') 'amove ',xt,yt
	     !     write (60,'(a)') 'set hei 0.01'
	     !     write (60,'(a,b14)') 'text ',key_leaf(j)



	  end do



	  do i=1,npp
	!     write (60,'(a,2f13.4)') 'amove ',particles(i)%x(1),particles(i)%x(2)
	     !     write (60,'(2a)') 'circle .002 fill ',colors( mod(me,8) )
	 !    write (60,'(2a)') 'circle psize fill black'
	  end do


	  ! map of interleaved bits


	  ! first point
	  xt=particles(1)%x(1)
	  yt=particles(1)%x(1)
	  !  write (60,'(a,2f13.4)') 'amove ',xt,yt

	  do ip = 2,npp
	     xt=particles(ip)%x(1)
	     yt=particles(ip)%x(2)
	     !     write (60,'(a,2f13.4)') 'aline ',xt,yt
	  end do

	  !  boundary links
	  if (me /= num_pe-1) then
	     !     write (60,'(a)') 'set color grey50'
	     !     write (60,'(a,2f13.4)') 'aline ',particles(npp+1)%x(1),particles(npp+1)%x(2)
	  endif



	  close(60)
	  close(61)

	end subroutine draw_domains


end module module_gle
