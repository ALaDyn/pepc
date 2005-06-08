!  ===================================================================
!
!                              FIELDS - internal parallel version
!
!   $Revision$
!
!   Calculate fields and potential from INTERNAL coordinates x,y,z, ux, uy, uz in treevars module:
!
!
!   ** Returns fields Ex, Ey, Ez and potential pot excluding external terms **
!
!
!  ===================================================================


subroutine pepc_fields_p(np_local, mac, theta, eps, force_const, bond_const, err_f, delta_t,  xl, yl, zl, itime, &
     coulomb, bfield_on, bonds, lenjones, &
     t_domain,t_build,t_prefetch, t_walk, t_walkc, t_force)

  use treevars
  use utils
  implicit none
  include 'mpif.h'

  integer, intent(in) :: np_local  ! # particles on this CPU
  real, intent(in) :: theta       ! multipole opening angle
  real, intent(in) :: err_f       ! max tolerated force error (rms)
  real, intent(in) :: delta_t       ! timestep 
  real, intent(in) :: force_const, bond_const       ! scaling factors for Coulomb/Len-Jones fields & potential
  real, intent(in) :: eps         ! potential softening distance
  logical, intent(in) :: coulomb  ! Compute Coulomb fields
  logical, intent(in) :: bfield_on  ! Switch for including B-fields
  logical, intent(in) :: bonds  ! Include bond forces
  logical, intent(in) :: lenjones ! Include Lennard-Jones potential
  real, intent(in) :: xl, yl, zl         ! box dimensions
  integer, intent(in) :: itime  ! timestep
  integer, intent(in) :: mac  ! choice of mac



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
     write (*,'(a7,a20/(i16,4f15.3))') 'PEPC | ','Initial buffers: ',(pelabel(i), x(i), y(i), z(i), q(i),i=1,npp) 
  endif


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

        ! zero field sums
        Ex(p) = 0.
        Ey(p) = 0.
        Ez(p) = 0.
        pot(p) = 0.
        Axo(p) = Ax(p)
        Ayo(p) = Ay(p)
        Azo(p) = Az(p)
        Ax(p) = 0.
        Ay(p) = 0.
        Az(p) = 0.
        Bx(p) = 0.
        By(p) = 0.
        Bz(p) = 0.

        if (coulomb) then
           !  compute Coulomb fields and potential of particle p from its interaction list
           call sum_force(p, nterm(i), nodelist( 1:nterm(i),i), eps, ex_coul, ey_coul, ez_coul, phi_coul, work(p))


           pot(p) = pot(p) + force_const * phi_coul
           Ex(p) = Ex(p) + force_const * ex_coul
           Ey(p) = Ey(p) + force_const * ey_coul
           Ez(p) = Ez(p) + force_const * ez_coul
        endif

        if (bfield_on) then
           !  compute magnetic fields and vector potential 
           call sum_bfield(p, nterm(i), nodelist( 1:nterm(i),i), eps, bx_ind, by_ind, bz_ind, ax_ind, ay_ind, az_ind )

           Ax(p) =  force_const * ax_ind
           Ay(p) =  force_const * ay_ind
           Az(p) =  force_const * az_ind
           Bx(p) =  force_const * bx_ind
           By(p) =  force_const * by_ind
           Bz(p) =  force_const * bz_ind
           ! Adjust E-field with inductive term
	   Ex(p) = Ex(p) - (Ax(p)-Axo(p))/delta_t
	   Ey(p) = Ey(p) - (Ay(p)-Ayo(p))/delta_t
   	   Ez(p) = Ez(p) - (Az(p)-Azo(p))/delta_t
        endif
        if (bonds) then
           !  compute short-range forces and potential of particle p from its interaction list
           !          call sum_bond(p, nterm(i), nodelist( 1:nterm(i),i ), fsx, fsy, fsz, phi )

           !          pot(p) = pot(p) + bond_const * phi
           !          fx(p) = fx(p) + bond_const * fsx
           !          fy(p) = fy(p) + bond_const * fsy
           !          fz(p) = fz(p) + bond_const * fsz
        endif

        if (lenjones) then
           !  compute short-range Lennard-Jones forces and potential of particle p from its interaction list
           call sum_lennardjones(p, nterm(i), nodelist( 1:nterm(i),i ), fsx, fsy, fsz, phi )

           pot(p) = pot(p) + bond_const * phi
           Ex(p) = Ex(p) + bond_const * fsx
           Ey(p) = Ey(p) + bond_const * fsy
           Ez(p) = Ez(p) + bond_const * fsz
        endif

        work(p) = nterm(i)        ! Should really compute this in sum_force to allow for leaf/twig terms
        work_local = work_local+nterm(i)
     end do

     call cputime(t3)   ! timing
     t_force = t_force + t3-t2

     max_local = max( max_local,maxval(nterm(1:nps)) )  ! Max length of interaction list


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

end subroutine pepc_fields_p









