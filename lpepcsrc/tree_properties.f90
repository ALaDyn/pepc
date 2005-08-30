! ===========================================
!
!           TREE_PROPERTIES
!
!   Fill in tree properties: charges, multipole moments etc
!   from bottom up
!
! ===========================================

subroutine tree_properties

  use treevars
  use tree_utils

  implicit none
  include 'mpif.h'

  integer*8, dimension(maxaddress) ::  search_key, resolve_key, parent_key, sum_key, key_twig
  integer*8, dimension(8) :: sub_key, key_child

  integer, dimension(maxaddress) ::  node_twig, child_twig    ! all twig-nodes, including branches and base nodes
  integer, dimension(maxaddress) ::  parent_node, parent_addr
  integer, dimension(maxaddress) :: local_node
  integer :: addr_twig

  integer, parameter :: n_moments = 23  ! # property arrays
  real, dimension(4*n_moments*nbranch_max/num_pe) :: local_moments      ! local branch properties    - size depends on # moments          
  real, dimension(n_moments*nbranch_max) :: branch_moments   ! global branch properties
  integer, dimension(num_pe) :: nbranchmoments ! array containing total # multipole terms*branch list length
  integer, dimension(num_pe) :: recv_strides, recv_counts
  integer, dimension(nbranch_max) :: bindex, branch_addr, branch_node, branch_level
  logical :: duplicate(maxaddress)

  integer :: addr_child
  integer, dimension(8) :: node_child  !  child nodes
  real, dimension(8) :: xs, ys, zs   ! multipole shift vector

  real :: xss, yss, zss, gamma, vx, vy, vz
  integer :: i, j, k, maxtwig, maxleaf, maxlevel, nchild, ncheck, ntwig_domain, nsearch, newsub, cchild 
  integer :: node_b, nuniq, nsub, nparent, ilevel, ibr, pno, bno
  integer ::  addr_leaf, p_leaf, node_leaf       ! local leaf-nodes
  integer :: ierr, nbuf
  integer :: key2addr        ! Mapping function to get hash table address from key


  if (tree_debug) write(ipefile,'(a)') 'TREE PROPERTIES'
  if (me==0 .and. tree_debug) write(*,'(a)') 'LPEPC | PROPERTIES'



  ! Zero multipole arrays
  do i=-ntwig_pw, nleaf_pw
     charge(i) = 0.
     abs_charge(i) = 0.
     xcoc(i) = 0.
     ycoc(i) = 0.
     zcoc(i) = 0.
     xdip(i) = 0.
     ydip(i) = 0.
     zdip(i) = 0.
     xxquad(i) = 0.
     yyquad(i) = 0.
     zzquad(i) = 0.
     xyquad(i) = 0.
     yzquad(i) = 0.
     zxquad(i) = 0.
     magmx(i) = 0.
     magmy(i) = 0.
     magmz(i) = 0.
     jx(i) = 0.
     jy(i) = 0.
     jz(i) = 0.
  end do

  !  Start with *local* leaf properties

  ! key_leaf(1:nleaf_me) =  pack(htable%key, mask = (htable%node > 0 .and. htable%owner == me) )  ! local leaf keys (could be predefined
  do i=1, nleaf_me
     addr_leaf = key2addr( leaf_key(i) )   !  Table address
     p_leaf = htable( addr_leaf )%node   !  Local particle index  - points to properties on PE
     node_leaf = p_leaf   !  Leaf node index is identical to particle index for *local* leaves 

     xcoc( node_leaf ) = x( p_leaf )         ! Centre of charge
     ycoc( node_leaf ) = y( p_leaf )
     zcoc( node_leaf ) = z( p_leaf )

     charge( node_leaf ) = q( p_leaf )       ! Charge
     abs_charge( node_leaf ) = abs( q(p_leaf) )   ! Absolute charge (needed for c.o.c)

     xdip( node_leaf ) = x(p_leaf) * q(p_leaf)   ! Dipole moment
     ydip( node_leaf ) = y(p_leaf) * q(p_leaf) 
     zdip( node_leaf ) = z(p_leaf) * q(p_leaf) 

     xxquad( node_leaf ) = q(p_leaf) * x(p_leaf)**2   ! Quadrupole moment
     yyquad( node_leaf ) = q(p_leaf) * y(p_leaf)**2   
     zzquad( node_leaf ) = q(p_leaf) * z(p_leaf)**2   
     xyquad( node_leaf ) = q(p_leaf) * x(p_leaf) * y(p_leaf) 
     yzquad( node_leaf ) = q(p_leaf) * y(p_leaf) * z(p_leaf)  
     zxquad( node_leaf ) = q(p_leaf) * z(p_leaf) * x(p_leaf)   

     !  get true velocities from momenta
     gamma = sqrt(1.0+ux(p_leaf)**2+uy(p_leaf)**2 + uz(p_leaf)**2)
     vx = ux(p_leaf)/gamma
     vy = uy(p_leaf)/gamma
     vz = uz(p_leaf)/gamma
     jx( node_leaf ) = vx * q(p_leaf)   ! Drift moment
     jy( node_leaf ) = vy * q(p_leaf) 
     jz( node_leaf ) = vz * q(p_leaf) 

     magmx( node_leaf ) = 0.5*q(p_leaf) * ( y(p_leaf)*vz - z(p_leaf)*vy)  ! Magnetic dipole moment
     magmy( node_leaf ) = 0.5*q(p_leaf) * ( z(p_leaf)*vx - x(p_leaf)*vz)  
     magmz( node_leaf ) = 0.5*q(p_leaf) * ( x(p_leaf)*vy - y(p_leaf)*vx)  


     !  Zero shift vector
     xshift( node_leaf ) = 0.
     yshift( node_leaf ) = 0.
     zshift( node_leaf ) = 0.
  end do


  !  Accumulate local twig properties

  ! Make list of twig nodes contained within local branch list (pebranch).  Recursive search adapted from make_branches.


  search_key(1:nbranch) = pebranch(1:nbranch)            ! start with branch list
  ncheck = 0       ! Checksum
  newsub = 0
  ntwig_domain = 0  ! # twigs contained in branch list
  nsearch = nbranch

  do while ( ncheck < nleaf_me )     ! Repeat until all local leaves found

     do i=1,nsearch
        if (  htable( key2addr( search_key(i) ) )%node > 0 ) then
           !  leaf,  so skip and increment checksum
           ncheck = ncheck +  1  

        else 
           ! twig: add to list and  subdivide
           ntwig_domain = ntwig_domain + 1
           key_twig(ntwig_domain) = search_key(i)

           cchild = htable( key2addr( search_key(i) ) )%childcode   !  Children byte-code
           nchild = SUM( (/ (ibits(cchild,j,1),j=0,7) /) ) ! # children = sum of bits in byte-code
           sub_key(1:nchild) = pack( bitarr, mask=(/ (btest(cchild,j),j=0,7) /) )  ! Extract child sub-keys from byte code

           resolve_key(newsub+1:newsub+nchild) = IOR( ishft( search_key(i),3 ), sub_key(1:nchild) ) ! Add keys of children to new search list
           newsub = newsub + nchild

        endif
     end do

     search_key(1:newsub) = resolve_key(1:newsub)        ! Put children into search list
     nsearch = newsub  ! # new nodes to search
     newsub = 0
  end do

  call sort(key_twig(1:ntwig_domain))

  if (tree_debug) then
!     write (ipefile,*) 'Twigs contained in local branch list: ',key_twig(1:ntwig_domain)
     write (ipefile,*) '# Twigs contained in local branch list: ',ntwig_domain
     write (ipefile,*) 'Found ',ncheck,' out of ',nleaf,' leaves'
  endif

  do i=1,ntwig_domain
     addr_twig = key2addr( key_twig(i) )   !  Table address
     node_twig(i) = htable( addr_twig )%node   !  Twig node index  
     child_twig(i) = htable( addr_twig )%childcode   !  Twig children byte-code 
  end do

  ! Go up through tree, starting at deepest level (largest key first)
  ! and accumulate multipole moments onto twig nodes

  do i = ntwig_domain,1,-1
     nchild = SUM( (/ (ibits(child_twig(i),j,1),j=0,7) /) )                 ! Get # children
     sub_key(1:nchild) = pack( bitarr, mask=(/ (btest(child_twig(i),j),j=0,7) /) )  ! Extract sub key from byte code

     do j=1,nchild
        key_child(j) = IOR( ishft( key_twig(i),3 ), sub_key(j) )      ! Construct keys of children
        addr_child = key2addr( key_child(j) )             ! Table address of children
        node_child(j) = htable( addr_child )%node                     ! Child node index  
     end do

     charge( node_twig(i) ) = SUM( charge( node_child(1:nchild) ) )                 ! Sum charge of child nodes and place
     ! result with parent node
     abs_charge( node_twig(i) ) = SUM( abs_charge( node_child(1:nchild) ) )                 ! Sum |q|

     ! Centres of charge
     xcoc( node_twig(i) ) = 0.
     ycoc( node_twig(i) ) = 0.
     zcoc( node_twig(i) ) = 0.

     do j=1,nchild
        xcoc( node_twig(i) ) = xcoc( node_twig(i) ) + ( xcoc(node_child(j)) * abs_charge(node_child(j))) &
             / abs_charge( node_twig(i) )
        ycoc( node_twig(i) ) = ycoc( node_twig(i) ) + ( ycoc(node_child(j)) * abs_charge(node_child(j))) &
             / abs_charge( node_twig(i) )
        zcoc( node_twig(i) ) = zcoc( node_twig(i) ) + ( zcoc(node_child(j)) * abs_charge(node_child(j))) &
             / abs_charge( node_twig(i) )
     end do

     ! Shifts and multipole moments
     xs(1:nchild) = xcoc( node_twig(i) ) - xshift( node_child(1:nchild) )     ! Shift vector for current node
     ys(1:nchild) = ycoc( node_twig(i) ) - yshift( node_child(1:nchild) )
     zs(1:nchild) = zcoc( node_twig(i) ) - zshift( node_child(1:nchild) )

     xshift( node_twig(i) ) = xcoc( node_twig(i) ) ! Shift variable for next level up
     yshift( node_twig(i) ) = ycoc( node_twig(i) ) 
     zshift( node_twig(i) ) = zcoc( node_twig(i) ) 



     do j = 1,nchild
        ! dipole moment

        xdip( node_twig(i) ) = xdip( node_twig(i) ) + xdip( node_child(j) ) - charge( node_child(j) )*xs(j) 
        ydip( node_twig(i) ) = ydip( node_twig(i) ) + ydip( node_child(j) ) - charge( node_child(j) )*ys(j) 
        zdip( node_twig(i) ) = zdip( node_twig(i) ) + zdip( node_child(j) ) - charge( node_child(j) )*zs(j) 

        ! quadrupole moment
        xxquad( node_twig(i) ) = xxquad( node_twig(i) ) + xxquad( node_child(j) ) - 2*xdip( node_child(j) )*xs(j) &
             + charge( node_child(j) )*xs(j)**2 
        yyquad( node_twig(i) ) = yyquad( node_twig(i) ) + yyquad( node_child(j) ) - 2*ydip( node_child(j) )*ys(j) &
             + charge( node_child(j) )*ys(j)**2 
        zzquad( node_twig(i) ) = zzquad( node_twig(i) ) + zzquad( node_child(j) ) - 2*zdip( node_child(j) )*zs(j) &
             + charge( node_child(j) )*zs(j)**2 

        xyquad( node_twig(i) ) = xyquad( node_twig(i) ) + xyquad( node_child(j) ) - xdip( node_child(j) )*ys(j) &
             - ydip( node_child(j) )*xs(j) + charge( node_child(j) )*xs(j)*ys(j) 
        yzquad( node_twig(i) ) = yzquad( node_twig(i) ) + yzquad( node_child(j) ) - ydip( node_child(j) )*zs(j) &
             - zdip( node_child(j) )*ys(j) + charge( node_child(j) )*ys(j)*zs(j) 
        zxquad( node_twig(i) ) = zxquad( node_twig(i) ) + zxquad( node_child(j) ) - zdip( node_child(j) )*xs(j) &
             - xdip( node_child(j) )*zs(j) + charge( node_child(j) )*zs(j)*xs(j) 

        ! magnetic dipole moment
        magmx( node_twig(i) ) = magmx( node_twig(i) ) + magmx( node_child(j) ) - 0.5*ys(j)*jz( node_child(j) )  - 0.5*zs(j)*jy( node_child(j) )
        magmy( node_twig(i) ) = magmy( node_twig(i) ) + magmy( node_child(j) ) - 0.5*zs(j)*jx( node_child(j) )  - 0.5*xs(j)*jz( node_child(j) )
        magmz( node_twig(i) ) = magmz( node_twig(i) ) + magmz( node_child(j) ) - 0.5*xs(j)*jy( node_child(j) )  - 0.5*ys(j)*jx( node_child(j) )

        jx( node_twig(i) ) = jx( node_twig(i) ) + jx( node_child(j) )       ! Sum currents of child nodes
        jy( node_twig(i) ) = jy( node_twig(i) ) + jy( node_child(j) )  
        jz( node_twig(i) ) = jz( node_twig(i) ) + jz( node_child(j) )  
     end do
  end do


  ! Should now have multipole information up to branch list level(s).
  ! By definition, this is complete: each branch node is self-contained.
  ! This information has to be broadcast to the other PEs so that the top levels can be filled in.

  local_node(1:nbranch) =  (/ ( htable( key2addr( pebranch(i) ) )%node, i=1,nbranch ) /)    ! Node #s of local branches

  ! Prepare properties for broadcast

  do i=1,nbranch
     ibr = (i-1)*n_moments
     local_moments(ibr+1) = charge( local_node(i) )
     local_moments(ibr+2) = abs_charge( local_node(i) )
     local_moments(ibr+3) = xcoc( local_node(i) )
     local_moments(ibr+4) = ycoc( local_node(i) )
     local_moments(ibr+5) = zcoc( local_node(i) )
     local_moments(ibr+6) = xshift( local_node(i) )
     local_moments(ibr+7) = yshift( local_node(i) )
     local_moments(ibr+8) = zshift( local_node(i) )
     local_moments(ibr+9) =  xdip( local_node(i) )
     local_moments(ibr+10) = ydip( local_node(i) )
     local_moments(ibr+11) = zdip( local_node(i) )
     local_moments(ibr+12) = xxquad( local_node(i) )
     local_moments(ibr+13) = yyquad( local_node(i) )
     local_moments(ibr+14) = zzquad( local_node(i) )
     local_moments(ibr+15) = xyquad( local_node(i) )
     local_moments(ibr+16) = yzquad( local_node(i) )
     local_moments(ibr+17) = zxquad( local_node(i) )
     local_moments(ibr+18) = magmx( local_node(i) )
     local_moments(ibr+19) = magmy( local_node(i) )
     local_moments(ibr+20) = magmz( local_node(i) )
     local_moments(ibr+21) = jx( local_node(i) )
     local_moments(ibr+22) = jy( local_node(i) )
     local_moments(ibr+23) = jz( local_node(i) )
  end do

  call MPI_BARRIER( MPI_COMM_WORLD, ierr)  ! Synchronize

  !  Collect multipole properties

  recv_counts(1:num_pe) = n_moments*nbranches(1:num_pe)   ! # terms in global array
  recv_strides(1:num_pe) = n_moments*igap(1:num_pe)       ! Stride lengths for global array
  nbuf = n_moments*nbranch


  call MPI_ALLGATHERV( local_moments, nbuf, MPI_REAL8, &
       branch_moments, recv_counts, recv_strides, MPI_REAL8, MPI_COMM_WORLD, ierr )


  !  write (ipefile,'(/a/a/(2i5,i10,f10.3))') 'Global branch properties:','   #     owner     key     abs(q) ', &
  !       (i,branch_owner(i), branch_key(i), branch_acharge(i),i=1,nbranch_sum)

  ! Copy info from *non-local* branches into property arrays of local tree

  do i=1,nbranch_sum
     ibr = (i-1)*n_moments
     if (branch_owner(i) /= me) then
        node_b = htable( key2addr( branch_key(i) ))%node   !  branch node index: already defined in MAKE_BRANCHES
        charge( node_b ) = branch_moments(ibr+1)      ! Copy collected remote branch props into main arrays
        abs_charge( node_b ) = branch_moments(ibr+2)  
        xcoc( node_b ) = branch_moments(ibr+3)  
        ycoc( node_b ) = branch_moments(ibr+4)  
        zcoc( node_b ) = branch_moments(ibr+5)  
        xshift( node_b ) = branch_moments(ibr+6)  
        yshift( node_b ) = branch_moments(ibr+7)  
        zshift( node_b ) = branch_moments(ibr+8)  
        xdip( node_b ) = branch_moments(ibr+9)  
        ydip( node_b ) = branch_moments(ibr+10)  
        zdip( node_b ) = branch_moments(ibr+11)  
        xxquad( node_b ) = branch_moments(ibr+12)  
        yyquad( node_b ) = branch_moments(ibr+13)  
        zzquad( node_b ) = branch_moments(ibr+14)  
        xyquad( node_b ) = branch_moments(ibr+15)  
        yzquad( node_b ) = branch_moments(ibr+16)  
        zxquad( node_b ) = branch_moments(ibr+17)  
        magmx( node_b ) = branch_moments(ibr+18)  
        magmy( node_b ) = branch_moments(ibr+19)  
        magmz( node_b ) = branch_moments(ibr+20)  
        jx( node_b ) = branch_moments(ibr+21)  
        jy( node_b ) = branch_moments(ibr+22)  
        jz( node_b ) = branch_moments(ibr+23)  
     endif
  end do

  ! Fill up top levels with tree properties: same algorithm as TREE_FILL

  maxlevel = 0
  do i=1,nbranch_sum
     branch_level(i) = log(1.*branch_key(i))/log(8.)       ! Get levels of branch nodes
     maxlevel = max( maxlevel,branch_level(i) )                                    ! Find maximum level
  end do

  nparent = 0

  do ilevel = maxlevel,1,-1                                           ! Start at deepest branch level and work up to root
     nsub = count( mask=branch_level(1:nbranch_sum) == ilevel )                      ! Count # branches at this level
     sum_key(1:nsub) =  pack(branch_key(1:nbranch_sum), mask = branch_level(1:nbranch_sum) == ilevel)       ! Pick out branches at current level

     sum_key(nsub+1:nsub+nparent) = parent_key(1:nparent)             ! Augment list with parent keys checked at previous level
     nsub = nsub + nparent
     call sort(sum_key(1:nsub))   ! Sort keys

     sum_key(nsub+1) = 0
     duplicate(1:nsub) = (/ (sum_key(i) /= sum_key(i+1), i=1,nsub) /)  ! Identify unique keys     
     nuniq= count(mask = duplicate(1:nsub))                            ! Count them
     sum_key(1:nuniq) = pack(sum_key(1:nsub), mask = duplicate(1:nsub))        ! Compress list

     do i=1,nuniq
        branch_addr(i) = key2addr( sum_key(i) )   !  branches' #table addresses
        branch_node(i) =  htable( branch_addr(i ) )%node

        parent_key(i) = ISHFT( sum_key(i),-3 )                ! parent keys
        parent_addr(i) = key2addr( parent_key(i) )   ! parents' #table addresses
        parent_node(i) =  htable( parent_addr(i) )%node          ! parents' node numbers
     end do

     ! Compute parent properties from children
     do i=nuniq,1,-1
        abs_charge( parent_node(i) ) = abs_charge( parent_node(i) ) + abs_charge( branch_node(i) )           ! Sum |q|
        charge( parent_node(i) ) = charge( parent_node(i) ) + charge( branch_node(i) )                       ! Sum q
     end do

     ! parent charges should be complete before computing coq's

     do i=nuniq,1,-1
        ! Centres of charge
        xcoc( parent_node(i) ) = xcoc(parent_node(i)) + (xcoc( branch_node(i) ) * abs_charge( branch_node(i) ) ) &
             / abs_charge( parent_node(i)) ! coq
        ycoc( parent_node(i) ) = ycoc(parent_node(i)) + (ycoc( branch_node(i) ) * abs_charge( branch_node(i) ) ) &
             / abs_charge( parent_node(i)) ! coq
        zcoc( parent_node(i) ) = zcoc(parent_node(i)) + (zcoc( branch_node(i) ) * abs_charge( branch_node(i) ) ) &
             / abs_charge( parent_node(i)) ! coq

     end do

     do i=nuniq,1,-1
        ! Shifts and multipole moments
        xss = xcoc( parent_node(i) ) - xshift( branch_node(i)  )     ! Shift vector for current child node
        yss = ycoc( parent_node(i) ) - yshift( branch_node(i) )
        zss = zcoc( parent_node(i) ) - zshift( branch_node(i) )

        xshift( parent_node(i) ) = xcoc( parent_node(i) ) ! Shift variable for next level up
        yshift( parent_node(i) ) = ycoc( parent_node(i) ) 
        zshift( parent_node(i) ) = zcoc( parent_node(i) ) 

        ! dipole moment
        xdip( parent_node(i) ) = xdip( parent_node(i) ) + xdip( branch_node(i) ) - charge( branch_node(i) )*xss 
        ydip( parent_node(i) ) = ydip( parent_node(i) ) + ydip( branch_node(i) ) - charge( branch_node(i) )*yss 
        zdip( parent_node(i) ) = zdip( parent_node(i) ) + zdip( branch_node(i) ) - charge( branch_node(i) )*zss 

        ! quadrupole moment
        xxquad( parent_node(i) ) = xxquad( parent_node(i) ) +  xxquad( branch_node(i) ) - 2*xdip( branch_node(i) )*xss &
             + charge( branch_node(i) )*xss**2 
        yyquad( parent_node(i) ) = yyquad( parent_node(i) ) +  yyquad( branch_node(i) ) - 2*ydip( branch_node(i) )*yss &
             + charge( branch_node(i) )*yss**2
        zzquad( parent_node(i) ) = zzquad( parent_node(i) ) +  zzquad( branch_node(i) ) - 2*zdip( branch_node(i) )*zss &
             + charge( branch_node(i) )*zss**2 
        xyquad( parent_node(i) ) = xyquad( parent_node(i) ) +  xyquad( branch_node(i) ) - xdip( branch_node(i) )*yss &
             - ydip( branch_node(i) )*xss + charge( branch_node(i) )*xss*yss
        yzquad( parent_node(i) ) = yzquad( parent_node(i) ) +  yzquad( branch_node(i) ) - ydip( branch_node(i) )*zss &
             - zdip( branch_node(i) )*yss + charge( branch_node(i) )*yss*zss
        zxquad( parent_node(i) ) = zxquad( parent_node(i) ) +  zxquad( branch_node(i) ) - zdip( branch_node(i) )*xss &
             - xdip( branch_node(i) )*zss + charge( branch_node(i) )*zss*xss

        magmx( parent_node(i) ) = magmx( parent_node(i) ) + magmx( branch_node(i) ) - 0.5*yss*jz( branch_node(i) )  - 0.5*zss*jy( branch_node(i) ) 
        magmy( parent_node(i) ) = magmy( parent_node(i) ) + magmy( branch_node(i) ) - 0.5*zss*jx( branch_node(i) )  - 0.5*xss*jz( branch_node(i) ) 
        magmz( parent_node(i) ) = magmz( parent_node(i) ) + magmz( branch_node(i) ) - 0.5*xss*jy( branch_node(i) )  - 0.5*yss*jx( branch_node(i) ) 

        jx( parent_node(i) ) = jx( parent_node(i) ) + jx( branch_node(i) )
        jy( parent_node(i) ) = jy( parent_node(i) ) + jy( branch_node(i) )
        jz( parent_node(i) ) = jz( parent_node(i) ) + jz( branch_node(i) )
     end do


     nparent = nuniq
  end do

  ! Rezero dipole and quadrupole sums of all local leaf nodes
  do i=1,nleaf_pw
     xdip(i) = 0.
     ydip(i) = 0.
     zdip(i) = 0.
     xxquad(i) = 0.
     yyquad(i) = 0.
     zzquad(i) = 0.
     xyquad(i) = 0.
     yzquad(i) = 0.
     zxquad(i) = 0.
     magmx(i) = 0.
     magmy(i) = 0.
     magmz(i) = 0.
  end do

  call MPI_BARRIER( MPI_COMM_WORLD, ierr)  ! Synchronize

end subroutine tree_properties
