
!  ================================
!
!         TREE_DOMAINS
!
!     Domain decomposition:
!     Share particle keys amoung PEs
!     - weighting according to load incurred on 
!       previous timestep
!
!  ================================


subroutine tree_domains(xl,yl,zl,indxl,irnkl,islen,irlen,fposts,gposts,npnew,npold,choose_sort,weighted)

  use treevars
  use tree_utils
  use utils
  use timings
  implicit none
  include 'mpif.h'

  real, intent(in) :: xl,yl,zl  ! initial box limits
  integer, intent(in) :: choose_sort, weighted
  integer, intent(out) :: indxl(nppm),irnkl(nppm)
  integer, intent(out) :: islen(num_pe),irlen(num_pe)
  integer, intent(out) :: fposts(num_pe+1),gposts(num_pe+1)
  integer :: npnew,npold

  integer*8, dimension(nppm) :: ix, iy, iz
  integer*8, dimension(nppm) :: ixd, iyd, izd
  integer*8, dimension(nppm) :: local_key

  real*8 :: ts1b=0., ts1e=0., ta1b=0., ta1e=0., ta2b=0., ta2e=0.

  integer ::  source_pe(nppm)
  integer :: i, j, ind_recv, inc, prev, next, handle(4)

  integer :: nbits

  real*8 :: s
  real*8 :: xmin_local, xmax_local, ymin_local, ymax_local, zmin_local, zmax_local
  logical :: boundary_debug=.false. 
  logical :: identical_keys=.false. 

  integer status(MPI_STATUS_SIZE), ierr, tag1

  ! arrays for parallel sort

  type (particle) :: ship_parts(nppm), get_parts(nppm)

  integer*8 :: xarray(nppm),keys(nppm),w1(nppm),wi2(nppm),wi3(nppm),compare(nppm)
  integer :: iteration, niterations, keycheck_pass, ipp
  integer :: errcount

  integer, dimension(nppm) ::  w2, w3 ! scratch arrays for integer*4 permute
  integer*8 :: tmp
  logical :: sort_debug
  real*8 :: xboxsize, yboxsize, zboxsize

  real*8 imba

  interface
     subroutine slsort_keys(nin,nmax,keys,workload,balance_weight,max_imbalance,nout,indxl,irnkl,scounts,rcounts,sdispls,rdispls,keys2,irnkl2,size,rank)
       integer,intent(in) :: nin,nmax,balance_weight,size,rank
       real*8,intent(in) :: max_imbalance
       integer,intent(out) :: nout,indxl(*),irnkl(*),scounts(*),rcounts(*),sdispls(*),rdispls(*),irnkl2(*)
       integer*8,intent(out) :: keys2(*)
       integer*8,intent(inout) :: keys(*)
       real*8,intent(inout) :: workload(*)
     end subroutine slsort_keys

     subroutine slsort_parts(n,nmax,keys,x,y,z,ux,uy,uz,q,m,work,ex,ey,ez,pelabel,balance_weight,max_imbalance,indxl,irnkl,scounts,rcounts,sdispls,rdispls,scratch0,scratch1,keys2,work2,irnkl2,size,rank)
       use treetypes
       integer,intent(inout) :: n
       integer,intent(in) :: nmax,balance_weight,size,rank
       integer*8,intent(inout) :: keys(*)
       real*8,intent(inout) :: x(*),y(*),z(*),ux(*),uy(*),uz(*),q(*),m(*),work(*),ex(*),ey(*),ez(*)
       integer,intent(inout) :: pelabel(*)
       real*8,intent(in) :: max_imbalance
       integer,intent(out) :: indxl(*),irnkl(*),scounts(*),rcounts(*),sdispls(*),rdispls(*)
       type (particle) :: scratch0(*), scratch1(*)
       integer*8,intent(out) :: keys2(*)
       real*8,intent(out) :: work2(*)
       integer,intent(out) :: irnkl2(*)
     end subroutine slsort_parts
  end interface

  integer*8 :: keys2(nppm)
  real*8 work2(nppm)
  integer irnkl2(nppm)

  integer local_count,count_stats(num_pe+1)
  real*8 local_work,work_stats(num_pe+1)
  real*8 d,minc,maxc,sumc,minw,maxw,sumw

  DOUBLE PRECISION :: ts, tpack, tunpack, talltoallv

  !POMP$ INST BEGIN(keys)

  ts1b = MPI_WTIME()
  ta1b = MPI_WTIME()

  sort_debug=domain_debug

  if (me==0 .and. tree_debug) write(*,'(a)') 'LPEPC | DOMAINS..'


  ! Find limits of local simulation region
  xmin_local = minval(x(1:npp))
  xmax_local = maxval(x(1:npp))
  ymin_local = minval(y(1:npp))
  ymax_local = maxval(y(1:npp))
  zmin_local = minval(z(1:npp))
  zmax_local = maxval(z(1:npp))

  ! Find global limits
  call MPI_ALLREDUCE(xmin_local, xmin, 1, MPI_REAL8, MPI_MIN,  MPI_COMM_WORLD, ierr )
  call MPI_ALLREDUCE(xmax_local, xmax, 1, MPI_REAL8, MPI_MAX,  MPI_COMM_WORLD, ierr )
  call MPI_ALLREDUCE(ymin_local, ymin, 1, MPI_REAL8, MPI_MIN,  MPI_COMM_WORLD, ierr )
  call MPI_ALLREDUCE(ymax_local, ymax, 1, MPI_REAL8, MPI_MAX,  MPI_COMM_WORLD, ierr )
  call MPI_ALLREDUCE(zmin_local, zmin, 1, MPI_REAL8, MPI_MIN,  MPI_COMM_WORLD, ierr )
  call MPI_ALLREDUCE(zmax_local, zmax, 1, MPI_REAL8, MPI_MAX,  MPI_COMM_WORLD, ierr )

  xboxsize = xmax-xmin
  yboxsize = ymax-ymin
  zboxsize = zmax-zmin

  ! Safety margin - put buffer region around particles
  xmax = xmax + xboxsize/10000.
  xmin = xmin - xboxsize/10000.
  ymax = ymax + yboxsize/10000.
  ymin = ymin - yboxsize/10000.
  zmax = zmax + zboxsize/10000.
  zmin = zmin - zboxsize/10000.

  boxsize = max(xmax-xmin, ymax-ymin, zmax-zmin)

  if (domain_debug) write (ipefile,'(4(a15,f12.4/))') &
       'xmin = ',xmin,'xmax = ',xmax, &
       'ymin = ',ymin,'ymax = ',ymax, &
       'zmin = ',zmin,'zmax = ',zmax, &
       'boxsize = ',boxsize

  s=boxsize/2**nlev       ! refinement length

  ! (xmin, ymin, zmin) is the translation vector from the tree box to the simulation region (in 1st octant)

  ix(1:npp) = ( x(1:npp) - xmin )/s           ! partial keys
  iy(1:npp) = ( y(1:npp) - ymin )/s           !
  iz(1:npp) = ( z(1:npp) - zmin )/s        

  ! construct keys by interleaving coord bits and add placeholder bit
  ! - note use of 64-bit constants to ensure correct arithmetic

  nbits = nlev+1
  do j = 1,npp
     !     local_key(j) = iplace + &
     !          SUM( (/ (8_8**i*(4_8*ibits( iz(j),i,1) + 2_8*ibits( iy(j),i,1 ) + 1_8*ibits( ix(j),i,1) ),i=0,nbits-1) /) )
     local_key(j) = iplace
     do i=0,nbits-1
        local_key(j) = local_key(j) &
             + 8_8**i*(4_8*ibits( iz(j),i,1) + 2_8*ibits( iy(j),i,1 ) + 1_8*ibits( ix(j),i,1) )
     end do
  end do

  !  if (.true.) then
  if (domain_debug) then
     write (ipefile,'(/a/a/(z21,i8,3f12.4,3i8,2f12.4))') 'Particle list before key sort:', &
          '  key,             label   coords     q ', &
          (local_key(i),pelabel(i),x(i),y(i),z(i),ix(i),iy(i),iz(i),q(i),work(i),i=1,npp) 
     !    write (ipefile,'(/a/a/(z21,i8,3f12.4,3i8))') '(last 10):', &
     !         '  key,                  label        coords              q ', &
     !         (local_key(i),pelabel(i),x(i),y(i),z(i),ix(i),iy(i),iz(i),q(i),work(i),i=max(1,npp-10),npp) 

     call blankn(ipefile)
  endif

  !POMP$ INST END(keys)

  !POMP$ INST BEGIN(sort)

  ! Use Parallel Sort by Regular Sampling (PSRS) 

  ta1e = MPI_WTIME()
  t_domains_keys = ta1e-ta1b  
  ta1b = MPI_WTIME()

  if (domain_debug .and. me==proc_debug) then
     write (*,*) 'MPI psrssort() commencing'
     write (*,*) 'iproc=',me
     write (*,*) 'num_pe=',num_pe
     write (*,*) 'npp=',npp
  endif

  ! Define wraps for ring network  0 -> 1 -> 2 -> ... ... -> num_pe-1 -> 0 ...
  if (me == 0) then
     prev = num_pe - 1
  else
     prev = me-1
  endif

  if (me == num_pe-1 ) then
     next = 0
  else
     next = me+1
  endif

  imba = 0.01

  npold = npp
  npnew = npp

  ts = MPI_WTIME()

  if (choose_sort == 2) then

     source_pe(1:npold) = pepid(1:npold)
     pekey(1:npp) = local_key(1:npp)

     ! sort particles with their data
     call slsort_parts(npp,nppm,&
          pekey,x,y,z,ux,uy,uz,q,m,work,ex,ey,ez,pelabel,&
          weighted,imba,&
          indxl,irnkl,islen,irlen,fposts,gposts,&
          ship_parts,get_parts,keys2,work2,irnkl2,&
          num_pe,me)

     ! set npnew to the new number of local particles
     npnew = npp

     t_domains_sort_pure = 0

  else

     iteration = 0
     niterations = 3  ! Max # iterations
     errcount = 1

     ! start permutation of local key list
     do i=1,npp
        xarray(i) = local_key(i)
     enddo

     do while (errcount /= 0 .or. iteration == niterations)
        iteration = iteration + 1

        if (iteration .lt. niterations) then
           do i=1,npold
              keys(i) = xarray(i)
           enddo
        else
           if (domain_debug .and. errcount .eq. 0) then
              write (ipefile,*) 'This time we are using already sorted data'
           endif
           npold = npnew
           do i=1,npold
              xarray(i) = w1(i)
              keys(i) = w1(i)
           enddo
        endif

        work2 = work

        ta2b = MPI_WTIME()

        ! perform index sort on keys
        if (choose_sort == 0) then
           call psrssort(nppm,npold,npnew,num_pe,me,keys, &
                indxl,irnkl,islen,irlen,fposts,gposts,w1)
!          call pswssort(nppm,npold,npnew,num_pe,me,keys, &
!               indxl,irnkl,islen,irlen,fposts,gposts,w1,work,key_box,load_balance,sort_debug)
        else if (choose_sort == 1) then
           !     call pbalsortr(nppm,npold,npnew,num_pe,me,keys, &
           !          indxl,irnkl,islen,irlen,fposts,gposts,pivots,w1,work,key_box,load_balance,sort_debug,work_local)
           call pbalsort(nppm,npold,npnew,num_pe,me,keys, &
                indxl,irnkl,islen,irlen,fposts,gposts,pivots,w1,work,nkeys_total,weighted,sort_debug,work_local)
        else ! Default: fastest
           call slsort_keys(npold,nppm,keys,work2,weighted,imba,npnew,indxl,irnkl,islen,irlen,fposts,gposts,w1,irnkl2,num_pe,me)
        end if

        ta2e = MPI_WTIME()
        t_domains_sort_pure = ta2e-ta2b

        !     write(*,*) me,npp,npold,npnew
        do i=1,npold
           w1(i) = xarray(i)
        enddo

        ! permute keys according to sorted indices
        call pll_permute(nppm,npold,npnew,num_pe,me,w1,wi2,wi3, &
             indxl,irnkl,islen,irlen,fposts,gposts)

        if (domain_debug) then
           write (ipefile,'(a/(i5,z20))') 'output array (1st 10): ',(i,w1(i),i=1,10)
           write (ipefile,'(a/(i5,z20))') 'output array (last 10): ',(i,w1(i),i=npnew-10,npnew)
        endif

        ! Check if sort finished
        errcount = 0
        do i=2,npnew
           if (w1(i) .lt. w1(i-1)) then
              errcount = errcount + 1
              if (domain_debug .and. errcount .lt. 10) then
                 write (ipefile,'(a,i5,2z20)') 'i,w1(i),w1(i-1)=',i,w1(i),w1(i-1)
              endif
           endif
        enddo
        !     write(*,*) me,errcount

        ! swap end items
        if (me .ne. 0) then
           call MPI_ISEND(w1, 1, MPI_INTEGER8, prev, 1, MPI_COMM_WORLD, handle(1), ierr)
           call MPI_REQUEST_FREE(handle(1),ierr)
        endif

        if (me .ne. num_pe-1) then
           call MPI_RECV(tmp, 1, MPI_INTEGER8, next, 1, MPI_COMM_WORLD, status, ierr)
        endif


        if (me .ne. num_pe-1) then
           !    if (me == 50 ) then
           if (domain_debug .and. tmp .lt. w1(npnew)) then          ! still something to sort
              write (*,'(a,i3,a1,2z20)') 'w1(npnew), w1(1) from',me+1, '=',w1(npnew),tmp
              errcount = errcount + 1  
           endif
        endif

        if (errcount .ne. 0 .and. domain_debug) then
           write (ipefile,*) 'errcount=',errcount
        endif

     enddo

     ta1e = MPI_WTIME()
     t_domains_sort = ta1e-ta1b
     ta1b = MPI_WTIME()

     npp = npnew
     pekey(1:npp) = w1(1:npp)

     ! Now permute remaining particle properties : x,y,z; vx,vy,vz; q,m, label, load

     source_pe(1:npold) = pepid(1:npold)   ! where particle came from

     ! Set up particle structure - keys and source_pe are dummies
     ! ( pekey is already sorted)

     tpack = MPI_WTIME()

     do i=1,npold
        ship_parts(i) = particle( x(indxl(i)), y(indxl(i)), z(indxl(i)), &
             ux(indxl(i)), uy(indxl(i)), uz(indxl(i)), &
             q(indxl(i)), m(indxl(i)), work(indxl(i)), &
             !          ax(indxl(i)), ay(indxl(i)), az(indxl(i)), &
             ex(indxl(i)), ey(indxl(i)), ez(indxl(i)), &
             keys(indxl(i)), pelabel(indxl(i)), source_pe(indxl(i))    )
     enddo

     tpack = MPI_WTIME() - tpack

     talltoallv = MPI_WTIME()

     ! perform permute
     call MPI_alltoallv(  ship_parts, islen, fposts, mpi_type_particle, &
          get_parts, irlen, gposts, mpi_type_particle, &
          MPI_COMM_WORLD,ierr)

     talltoallv = MPI_WTIME() - talltoallv

     tunpack = MPI_WTIME()

     do i=1,npp
        x(irnkl(i)) = get_parts(i)%x
        y(irnkl(i)) = get_parts(i)%y
        z(irnkl(i)) = get_parts(i)%z
        ux(irnkl(i)) = get_parts(i)%ux
        uy(irnkl(i)) = get_parts(i)%uy
        uz(irnkl(i)) = get_parts(i)%uz
        q(irnkl(i)) = get_parts(i)%q
        m(irnkl(i)) = get_parts(i)%m
        work(irnkl(i)) = get_parts(i)%work
        !     ax(irnkl(i)) = get_parts(i)%ax
        !     ay(irnkl(i)) = get_parts(i)%ay
        !     az(irnkl(i)) = get_parts(i)%az
        ex(irnkl(i)) = get_parts(i)%ax
        ey(irnkl(i)) = get_parts(i)%ay
        ez(irnkl(i)) = get_parts(i)%az
        pelabel(irnkl(i)) = get_parts(i)%label
     enddo

     tunpack = MPI_WTIME() - tunpack

  endif

  ts = MPI_WTIME() - ts

  ta1e = MPI_WTIME()
  t_domains_ship = ta1e - ta1b
  ta1b = MPI_WTIME()

  if (me .eq. -1) then
     write (*,*) 'sorting: ', ts, ' = ', t_domains_sort_pure, ' + ', t_domains_ship, '(', tpack, ' + ', talltoallv, ' + ', tunpack, ')'
  endif

  ! gather workload information
  if (0 .eq. -1) then

     ! new workloads
     local_count = npp;
     local_work = 0;
     do i=1,npp
!        write(*,*) i,work(i)
        local_work = local_work + work(i)
     enddo

     call MPI_Allgather(local_count, 1, MPI_INTEGER, count_stats, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr);
     call MPI_Allgather(local_work, 1, MPI_DOUBLE_PRECISION, work_stats, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr);

     work_stats(num_pe+1) = 0;
     count_stats(num_pe+1) = 0;
     do i=1,num_pe
        work_stats(num_pe+1) = work_stats(num_pe+1) + work_stats(i)
        count_stats(num_pe+1) = count_stats(num_pe+1) + count_stats(i)
     enddo

     if (me .eq. 0) then
        write(*,*) 'total_count:',count_stats(num_pe+1)
        sumc = 0;
        minc = count_stats(num_pe+1)/num_pe;
        maxc = 0;
        do i=1,num_pe
!           write(*,*) i,count_stats(i),real(count_stats(i))/count_stats(num_pe+1),'/',((real(count_stats(num_pe+1))/num_pe)-count_stats(i)),abs(1.0d0-(real(count_stats(i))*num_pe/count_stats(num_pe+1)))
           d = abs((count_stats(num_pe+1)/num_pe)-count_stats(i));
           minc = min(minc, d);
           maxc = max(maxc, d);
           sumc = sumc + d;
        enddo
        write(*,*) 'avg_count:',count_stats(num_pe+1)/num_pe
        write(*,*) 'avg_cdiff:',sumc/num_pe,sumc/count_stats(num_pe+1)
        write(*,*) 'min_cdiff:',minc,minc*num_pe/count_stats(num_pe+1)
        write(*,*) 'max_cdiff:',maxc,maxc*num_pe/count_stats(num_pe+1)
        write(*,*) 'total_weight:',work_stats(num_pe+1)
        sumw = 0;
        minw = work_stats(num_pe+1)/num_pe;
        maxw = 0;
        do i=1,num_pe
!           write(*,*) i,work_stats(i),real(work_stats(i))/work_stats(num_pe+1),'/',((real(work_stats(num_pe+1))/num_pe)-work_stats(i)),abs(1.0d0-(real(work_stats(i))*num_pe/work_stats(num_pe+1)))
           d = abs((work_stats(num_pe+1)/num_pe)-work_stats(i));
           minw = min(minw, d);
           maxw = max(maxw, d);
           sumw = sumw + d;
        enddo
        write(*,*) 'avg_weight:',work_stats(num_pe+1)/num_pe
        write(*,*) 'avg_wdiff:',sumw/num_pe,sumw/work_stats(num_pe+1)
        write(*,*) 'min_wdiff:',minw,minw*num_pe/work_stats(num_pe+1)
        write(*,*) 'max_wdiff:',maxw,maxw*num_pe/work_stats(num_pe+1)
     endif

  endif

  pepid(1:npp) = me  ! new owner

  if (domain_debug) then
     do j=1,npp
        ixd(j) = SUM( (/ (2_8**i*ibits( pekey(j)-iplace,3*i,1 ), i=0,nbits-2) /) )
        iyd(j) = SUM( (/ (2_8**i*ibits( pekey(j)-iplace,3*i+1,1 ), i=0,nbits-2) /) )
        izd(j) = SUM( (/ (2_8**i*ibits( pekey(j)-iplace,3*i+2,1 ), i=0,nbits-2) /) )
     end do
     write (ipefile,'(/a/a/(z21,2i6,a2,i8,6f12.4))') 'Particle list after key sort:', &
          '  key,                owner,    from PE  |  label  Fetched coords      derived from key', &
          (pekey(i),pepid(i),source_pe(i),'|', &
          pelabel(i),x(i),y(i),z(i),ixd(i)*s+xmin,iyd(i)*s+ymin,izd(i)*s+zmin,i=1,npp) 

     call blankn(ipefile)
  endif


  ! Each PE now has sorted segment of particles of length npp
  ! Note that now npp /= npart/num_pe, only approx depending on key distribution, or target shape.

  ! Check for identical keys
  keycheck_pass=1

  do i=2,npp-1
     if (pekey(i) == pekey(i-1)) then
        pekey(i) = pekey(i) + 1  ! Augment higher key
        pekey(i+1) = pekey(i+1) + 2  ! Augment next higher key to avoid 'triplet'
        !  Tweak momenta and positions to ensure particles drift apart
        ux(i-1) = ux(i-1)*.99999
        ux(i+1) = ux(i+1)*1.0001
        x(i-1) = x(i-1)*.99999
        x(i+1) = x(i+1)*1.0001

        ! TODO: need 'ripple' here up to next large gap in keys i+1->npp

        write(*,'(a15,i5,a8,i3,a30,2i6,3i10,a25,o25,a12,o25)') 'LPEPC | PE ',me,' pass ',keycheck_pass, &
             ' WARNING: identical keys found for particles  ',i,npp,pelabel(i-1),pelabel(i),pelabel(i+1), &
             ' - upper increased to: ',pekey(i),' next key: ',pekey(i+1)
        !        if (x(i) == x(i-1)) write(*,*) "HELP"
     endif
  end do


  ! Special case for last pair to avoid possibility of identical keys across processor boundary
  ! - work back down from last key

  ipp=npp-1
  identical_keys=.true.
  keycheck_pass=2

  do while (identical_keys .and. ipp.gt.2)
     identical_keys=.false.
     if (pekey(ipp+1) == pekey(ipp)) then
        pekey(ipp) = pekey(ipp)-1
        write(*,'(a15,i5,a8,i3,a30,2i15/a25,o30)') 'LPEPC | PE ',me,' pass ',keycheck_pass, &
             ' WARNING: identical keys found for particles  ',pelabel(ipp+1),pelabel(ipp), &
             'LPEPC | Lower key decreased to:  ',pekey(ipp)
        identical_keys=.true.
     endif
     ipp=ipp-1
  end do

  ! Copy boundary particles to adjacent PEs to ensure proper tree construction
  !  - if we don't do this, can get two particles on separate PEs 'sharing' a leaf


  !  ship_props = particle ( x(1), y(1), z(1), ux(1), uy(1), uz(1), q(1), m(1), work(1), &
  !       ax(1),ay(1),az(1), pekey(1), pelabel(1), pepid(1) )

  ship_props = particle ( x(1), y(1), z(1), ux(1), uy(1), uz(1), q(1), m(1), work(1), &
       ax(1),ay(1),az(1), pekey(1), pelabel(1), pepid(1) )

  !  write (*,'(9f12.3,z20,2i6)') ship_props

  ! Ship 1st particle data to end of list of LH neighbour PE

  if (me /= 0 ) then
     call MPI_ISEND( ship_props, 1, mpi_type_particle, prev, 1, MPI_COMM_WORLD, handle(1), ierr ) 
     call MPI_REQUEST_FREE(handle(1),ierr) 
  endif

  ! Place incoming data at end of array
  if ( me /= num_pe-1) then
     !     call MPI_IRECV( get_props, 1, mpi_type_particle, next, 1,  MPI_COMM_WORLD, handle(2), ierr )
     call MPI_RECV( get_props, 1, mpi_type_particle, next, 1,  MPI_COMM_WORLD, status, ierr )
     x(npp+1) = get_props%x
     y(npp+1) = get_props%y
     z(npp+1) = get_props%z
     ux(npp+1) = get_props%ux
     uy(npp+1) = get_props%uy
     uz(npp+1) = get_props%uz
     q(npp+1) = get_props%q
     m(npp+1) = get_props%m
     work(npp+1) = get_props%work
     ex(npp+1) = get_props%ax
     ey(npp+1) = get_props%ay
     ez(npp+1) = get_props%az
     !     ax(npp+1) = get_props%ax
     !     ay(npp+1) = get_props%ay
     !     az(npp+1) = get_props%az
     pekey(npp+1) = get_props%key
     pelabel(npp+1) = get_props%label
     pepid(npp+1) = get_props%pid
  endif

  ! Ship  end particle data to start of list of RH neighbour PE

  !  ship_props = particle ( x(npp), y(npp), z(npp), ux(npp), uy(npp), uz(npp), q(npp), m(npp), work(npp), &
  !       ax(npp),ay(npp),az(npp), pekey(npp), pelabel(npp), pepid(npp) )

  ship_props = particle ( x(npp), y(npp), z(npp), ux(npp), uy(npp), uz(npp), q(npp), m(npp), work(npp), &
       ex(npp),ey(npp),ez(npp), pekey(npp), pelabel(npp), pepid(npp) )

  if (me /= num_pe-1 ) then
     call MPI_ISEND( ship_props, 1, mpi_type_particle, next, 2, MPI_COMM_WORLD, handle(3), ierr )
     call MPI_REQUEST_FREE(handle(3),ierr) 
  endif

  ! Place incoming data at end of array

  if (me == num_pe-1) then
     ind_recv = npp+1   ! PEn array hasn't yet received boundary value
  else
     ind_recv = npp+2
  endif

  if ( me /= 0) then
     !     call MPI_IRECV( get_props, 1, mpi_type_particle, prev, 2,  MPI_COMM_WORLD, handle(4), ierr )
     call MPI_RECV( get_props, 1, mpi_type_particle, prev, 2,  MPI_COMM_WORLD, status, ierr )
     x(ind_recv) = get_props%x
     y(ind_recv) = get_props%y
     z(ind_recv) = get_props%z
     ux(ind_recv) = get_props%ux
     uy(ind_recv) = get_props%uy
     uz(ind_recv) = get_props%uz
     q(ind_recv) = get_props%q
     m(ind_recv) = get_props%m
     work(ind_recv) = get_props%work
     !     ax(ind_recv) = get_props%ax
     !     ay(ind_recv) = get_props%ay
     !     az(ind_recv) = get_props%az
     ex(ind_recv) = get_props%ax
     ey(ind_recv) = get_props%ay
     ez(ind_recv) = get_props%az
     pekey(ind_recv) = get_props%key
     pelabel(ind_recv) = get_props%label
     pepid(ind_recv) = get_props%pid
  endif


  !  if (me /= 0) call MPI_WAIT(handle(4),status,ierr)

  if (boundary_debug) then
     if (me /= 0 .and. me /= num_pe-1) then
        inc = 1
     else
        inc = 0
     endif

     write (ipefile,'(/a/a/(i5,z21,2i8,3f12.5,2f12.3))') 'Particle list after boundary swap:', &
          ' index   key,     label,   on PE,    x      y     q       m', &
          (i,pekey(i),pelabel(i),pepid(i),x(i),y(i),z(i),q(i),m(i),i=1,npp+1+inc) 

     call blankn(ipefile)
  endif

  ta1e = MPI_WTIME()
  t_domains_bound = ta1e - ta1b
  ts1e = MPI_WTIME()
  t_domains = ts1e-ts1b

  if (me==0 .and. tree_debug) write(*,'(a)') 'LPEPC | ..done'
  !POMP$ INST END(sort)

end subroutine tree_domains
