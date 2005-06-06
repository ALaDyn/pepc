!  ===================================================================
!
!                              FIELDS
!
!   $Revision$
!
!   Calculate fields and potential from coordinates x,y,z:
!
!
!   ** Returns fields Ex, Ey, Ez and potential pot excluding external terms **
!
!
!  ===================================================================


subroutine pepc_fields(np_local, p_x, p_y, p_z, p_vx, p_vy, p_vz, p_q, p_m, p_w, p_label, &
     Ex, Ey, Ez, pot, &
     mac, theta, eps, force_const, err_f, xl, yl, zl, itime, &
     t_domain,t_build,t_prefetch, t_walk, t_walkc, t_force)

  ! use physvars
  use treevars
  use utils
  implicit none
  include 'mpif.h'

  integer, intent(in) :: np_local  ! # particles on this CPU
  real, intent(in) :: theta       ! multipole opening angle
  real, intent(in) :: err_f       ! max tolerated force error (rms)
  real, intent(in) :: force_const       ! scaling factor for fields & potential
  real, intent(in) :: eps         ! potential softening distance
  real, intent(in) :: xl, yl, zl         ! box dimensions
  integer, intent(in) :: itime  ! timestep
  integer, intent(in) :: mac  ! choice of mac
  real, intent(in), dimension(np_local) :: p_x, p_y, p_z  ! coords and velocities: x1,x2,x3, y1,y2,y3, etc 
  real, intent(in),  dimension(np_local) :: p_vx, p_vy, p_vz  ! coords and velocities: x1,x2,x3, y1,y2,y3, etc 
  real, intent(in), dimension(np_local) :: p_q, p_m ! charges, masses
  real, dimension(np_local) :: p_w ! work loads
  integer, intent(in), dimension(np_local) :: p_label  ! particle label 
  real, intent(out), dimension(np_local) :: ex, ey, ez, pot  ! fields and potential to return



  integer, parameter :: npassm=100000 ! Max # passes - will need npp/nshortm

  integer :: p, i, j, npass, jpass, ip1, nps,  max_npass,nshort_list, ipe
  real :: t_domain, t_build, t_prefetch, t_walk, t_walkc, t_force, ttrav, tfetch, t1, t2, t3  ! timing integrals
  real :: tb1, tb2, td1, td2, tp1, tp2
  integer :: pshortlist(nshortm),nshort(npassm),pstart(npassm) ! work balance arrays
  integer :: hashaddr ! Key address 

  integer :: max_local,  timestamp
  integer :: ierr
  integer :: iprot = 50  ! frequency for load balance dump

  real :: fsx, fsy, fsz, phi, phi_coul, ex_coul, ey_coul, ez_coul
  real :: ax_ind, ay_ind, az_ind, bx_ind, by_ind, bz_ind
  real :: Epon_x, Epon_y, Epon_z, Phipon, ex_em, ey_em, ez_em, bx_em, by_em, bz_em
  real :: xd, yd, zd  ! positions relative to centre of laser spot
  real :: work_local, load_average, load_integral, total_work, average_work
  integer :: total_parts
  character(30) :: cfile, ccol1, ccol2
  character(4) :: cme
  integer :: key2addr        ! Mapping function to get hash table address from key

!  force_debug=.true.
!  tree_debug=.false.
!  build_debug=.false.
!  domain_debug = .false.
!  branch_debug=.false.
!  prefetch_debug=.false.
!  walk_debug=.false.
!  walk_summary=.true.
!  dump_tree=.true.
  npp = np_local  ! assumed lists matched for now
  
  if (force_debug) then
     write (*,'(a7,a40,2i5,4f15.2)') 'PEPC | ','Params itime, mac, theta, eps, force_const, err:',itime, mac, theta, eps, force_const, err_f
     write (*,'(a7,a20/(i16,5f15.3))') 'PEPC | ','Initial buffers: ',(p_label(i), p_x(i), p_y(i), p_z(i), p_vx(i), p_q(i),i=1,npp) 
  endif

 ! Copy particle buffers to tree arrays
  x(1:npp) = p_x(1:npp)
  y(1:npp) = p_y(1:npp)
  z(1:npp) = p_z(1:npp)
  ux(1:npp) = p_vx(1:npp)
  uy(1:npp) = p_vy(1:npp)
  uz(1:npp) = p_vz(1:npp)
  q(1:npp) = p_q(1:npp)
  m(1:npp) = p_m(1:npp)
  pelabel(1:npp) = p_label(1:npp)
  pepid(1:npp) = me
  work(1:npp) = p_w(1:npp)
  ax(1:npp) = 0.
  ay(1:npp) = 0.
  az(1:npp) = 0.
 
  call cputime(td1)
  call tree_domains(xl,yl,zl)    ! Domain decomposition: allocate particle keys to PEs

! particles now sorted according to keys assigned in tree_domains.
! Serial mode: incoming labels can serve to restore order - now kept in pelabel 
! Parallel mode: have to retain sort vectors from tree_domains

  call cputime(tb1)
  call tree_build      ! Build trees from local particle lists
  call tree_branches   ! Determine and concatenate branch nodes
  call tree_fill       ! Fill in remainder of local tree
  call tree_properties ! Compute multipole moments for local tree
  call cputime(tp1)
  if (num_pe>1) call tree_prefetch(itime)
  call cputime(tp2)

  t_domain = tb1-td1
  t_build = tp1-tb1
  if (num_pe>1) then
     t_prefetch = tp2-tp1
  else
     t_prefetch = 0.
  endif
  t_walk=0.
  t_walkc=0.
  t_force=0.
  max_local = 0   ! max length of interaction list
  max_list_length = 0
  work_local = 0  ! total workload
  maxtraverse=0   ! max # traversals
  maxships=0      ! max # multipole shipments/traversal
  sumships=0      ! total # multipole shipments/iteration

  nfetch_total=0     ! Zero key fetch counters if prefetch mode off
  nreqs_total=0

  !  # passes needed to process all particles
  nshort_list =nshortm/4 
  npass = max(1,npart/num_pe/nshort_list)   ! ave(npp)/nshort_list   - make nshort_list a power of 2
  load_average = SUM(work(1:npp))/npass   ! Ave. workload per pass: same for all PEs if already load balanced
  nshort(1:npass+1) = 0

  load_integral = 0.
  jpass = 1
  pstart(jpass) = 1
  do i=1,npp
     load_integral = load_integral + work(i)   ! integrate workload

     if (i-pstart(jpass) + 1 == nshortm) then ! Need to check that nshort < nshortm
        write(*,*) 'Warning from PE: ',me,' # parts on pass ',jpass,' in shortlist exceeds array limit ',nshortm
        write(*,*) 'Putting spill-over into following pass'
        nshort(jpass) = nshortm
        jpass = jpass + 1
        pstart(jpass) = i+1

     else if (load_integral >= load_average * jpass .or. i==npp ) then
        nshort(jpass) = i-pstart(jpass) + 1
        jpass = jpass + 1
        pstart(jpass) = i+1

     endif
  end do


  if (jpass-1 > npass ) then
     write(*,*) 'PE',me,' missed some:',nshort(npass+1)
     if (nshort(npass) + nshort(npass+1) <= nshortm) then
        nshort(npass) = nshort(npass) + nshort(npass+1)
     else
        npass = npass+1
     endif
  endif

  if (force_debug)   write (ipefile,*) 'Shortlists: ',(nshort(j),j=1,npass+1)

  max_npass = npass

  ip1 = 1


  do jpass = 1,max_npass
     !  make short-list
     nps = nshort(jpass)
     ip1 = pstart(jpass)
     pshortlist(1:nps) = (/ (ip1+i-1, i=1,nps) /)

     if (force_debug) then
       	write(*,*) 'pass ',jpass,' of ',max_npass,': # parts ',ip1,' to ',ip1+nps-1
        write(ipefile,*) 'pass ',jpass,' # parts ',ip1,' to ',ip1+nps-1
     endif

     !  build interaction list: 
     ! tree walk creates intlist(1:nps), nodelist(1:nps) for particles on short list

     call tree_walk(pshortlist,nps,jpass,theta,itime,mac,ttrav,tfetch)
     t_walk = t_walk + ttrav  ! traversal time (serial)
     t_walkc = t_walkc + tfetch  ! multipole swaps

     call cputime(t2)   ! timing
     do i = 1, nps

        p = pshortlist(i)    ! local particle index


        !  compute Coulomb fields and potential of particle p from its interaction list
        call sum_force(p, nterm(i), nodelist( 1:nterm(i),i), eps, ex_coul, ey_coul, ez_coul, phi_coul, work(p))

! restore initial particle order specified by calling routine
! - only works for serial mode at present 
        pot(pelabel(p)) = force_const * phi_coul
        Ex(pelabel(p)) = force_const * ex_coul
        Ey(pelabel(p)) = force_const * ey_coul
        Ez(pelabel(p)) = force_const * ez_coul
        p_w(pelabel(p)) = work(p)  ! send back work load for next iteration

        work_local = work_local+nterm(i)
     end do

     call cputime(t3)   ! timing
     t_force = t_force + t3-t2

     max_local = max( max_local,maxval(nterm(1:nps)) )  ! Max length of interaction list

     if (dump_tree) call diagnose_tree

  end do





  call MPI_ALLREDUCE(max_local, max_list_length, 1, MPI_INTEGER, MPI_MAX,  MPI_COMM_WORLD, ierr )
  call MPI_GATHER(work_local, 1, MPI_REAL8, work_loads, 1, MPI_REAL8, 0,  MPI_COMM_WORLD, ierr )  ! Gather work integrals
  call MPI_GATHER(npp, 1, MPI_INTEGER, npps, 1, MPI_INTEGER, 0,  MPI_COMM_WORLD, ierr )  ! Gather particle distn

!  timestamp = itime + itime_start
  timestamp = itime

  if (me ==0 .and. mod(itime,iprot)==0) then
     total_work = SUM(work_loads)
     average_work = total_work/num_pe
     cme = achar(timestamp/1000+48) // achar(mod(timestamp/100,10)+48) &
          // achar(mod(timestamp/10,10)+48) // achar(mod(timestamp,10)+48) 
     cfile="load_"//cme//".dat"
     total_parts=SUM(npps)
     open(60, file=cfile)
     write(60,'(a/a,i8,2(a,1pe15.6))')  '! Full balancing','Parts: ',total_parts,' Work: ',total_work, &
          ' Ave. work:',average_work        
     write(60,'(2i8,f12.3))')  (i-1,npps(i),work_loads(i)/average_work,i=1,num_pe)
     close(60)
  endif

  if (force_debug) then
     write (ipefile,101)
     write (*,101)
     write (ipefile,101) force_const

     do i=1,npp
        write (ipefile,102) pelabel(i), & 
             q(i), m(i), ux(i), pot(i), ex(i)
        write (*,102) pelabel(i), x(i), & 
             q(i), m(i), ux(i), pot(i)
     end do

101  format('Tree forces:'/'   p    q   m   ux   pot  ',f8.2)
102  format(1x,i7,5(1pe14.5))

  endif

end subroutine pepc_fields









