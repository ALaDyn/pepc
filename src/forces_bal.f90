!  ===================================================================
!
!                              FORCES
!
!   $Revision$
!
!   Calculate forces from interaction list:
!     nterm(i) = terms
!   Pseudoparticles are given by: intlist(j,i), j=1,nterm(i);
!
!
!   ** Returns accelerations ax,ay,az including external field term **
!
!   fx,fy,fz   -    local forces
!
!  ===================================================================


subroutine forces_bal(p_start,p_finish,delta_t, t_walk, t_force)

  use treevars
  use utils
  implicit none

  real, intent(in) :: delta_t
  integer, intent(in) :: p_start,p_finish  ! min, max particle nos.

  integer, parameter :: npassm=100000 ! Max # passes - will need npp/nshortm

  integer :: p, i, j, npass, jpass, ip1, nps,  max_npass,nshort_list, ipe
  real :: t_walk, t_force, t1, t2, t3  ! timing integrals
  integer pshortlist(nshortm),nshort(npassm),pstart(npassm)
  real :: work_loads(num_pe)  ! Load balance array
  integer :: npps(num_pe)  ! Particle distrib amoung PEs
  integer :: max_local,  timestamp
  real, dimension(nppm) :: fx, fy, fz

  real :: fsx, fsy, fsz, phi, Ex, Ey, Ez, phi_coul, ex_coul, ey_coul, ez_coul
  real :: Epon_x, Epon_y, Epon_z, Phipon
  real :: xd, yd, zd  ! positions relative to centre of laser spot
  real :: work_local, load_average, load_integral, total_work, average_work
  integer :: total_parts
  character(30) :: cfile, ccol1, ccol2
  character(4) :: cme

  Ex = 0.
  Ey = 0.
  Ez = 0.
  fx=0.
  fy=0.
  fz=0.
  t_walk=0.
  t_force=0.
  max_local = 0   ! max length of interaction list
  work_local = 0  ! total workload

  !  # passes needed to process all particles
  nshort_list =nshortm/5 
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
     if (force_debug)     write(*,*) 'PE',me,' missed some:',nshort(npass+1)
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

     if (walk_debug) then
        !	write(*,*) 'pass ',jpass,' of ',max_npass,': # parts ',ip1,' to ',ip1+nps-1
        write(ipefile,*) 'pass ',jpass,' # parts ',ip1,' to ',ip1+nps-1
     endif

     !  build interaction list: 
     ! tree walk returns intlist(1:nps), nodelist(1:nps) for particles on short list

     call cputime(t1)
     call tree_walk(pshortlist,nps)
     call cputime(t2)
     t_walk = t_walk + t2-t1

     do i = 1, nps

        p = pshortlist(i)    ! local particle index
        fx(p) = 0.
        fy(p) = 0.
        fz(p) = 0.
        pot(p) = 0.

        if (coulomb) then
           !  compute Coulomb forces and potential of particle p from its interaction list
           call sum_force(p, nterm(i), nodelist( 1:nterm(i),i ), ex_coul, ey_coul, ez_coul, phi_coul )

           pot(p) = pot(p) + force_const * phi_coul
           fx(p) = fx(p) + force_const * q(p) * ex_coul
           fy(p) = fy(p) + force_const * q(p) * ey_coul
           fz(p) = fz(p) + force_const * q(p) * ez_coul
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
           fx(p) = fx(p) + bond_const * fsx
           fy(p) = fy(p) + bond_const * fsy
           fz(p) = fz(p) + bond_const * fsz
        endif

        work(p) = nterm(i)        ! Should really compute this in sum_force to allow for leaf/twig terms
        work_local = work_local+nterm(i)
     end do

     call cputime(t3)   ! timing
     t_force = t_force + t3-t2

     max_local = max( max_local,maxval(nterm(1:nps)) )  ! Max length of interaction list

  end do



  ! Include ponderomotive force from laser on electrons 
  ! - select either standing wave (overdense) or propagating (underdense) fpond


  if (beam_config ==4 .or. beam_config ==5) then
     do p = p_start, p_finish
        if (q(p)<0) then
           xd = x(p)-focus(1)   ! position relative to laser focus
           yd = y(p)-focus(2)
           zd = z(p)-focus(3)

           laser_model: select case(beam_config)

           case(4)  ! standing wave fpond
              call fpond( tlaser, tpulse,sigma,vosc,omega,xd,yd,zd,epon_x,epon_y,epon_z,phipon)
           case(5)  ! propagating fpond
              call laser_bullet( tlaser, focus(1), tpulse,sigma,vosc,omega, & 
                   xd,yd,zd,epon_x,epon_y,epon_z,phipon)
           case default
              Epon_x=0
              Epon_y=0
              Epon_z=0

           end select laser_model

           fx(p) = fx(p) + q(p) * Epon_x
           fy(p) = fy(p) + q(p) * Epon_y
           fz(p) = fz(p) + q(p) * Epon_z
        endif

     end do
  endif

  ! store accelerations


  do p = p_start, p_finish
     ax(p) = fx(p)/m(p)
     ay(p) = fy(p)/m(p)
     az(p) = fz(p)/m(p)
  end do

  call MPI_ALLREDUCE(max_local, max_list_length, 1, MPI_INTEGER, MPI_MAX,  MPI_COMM_WORLD, ierr )
  call MPI_GATHER(work_local, 1, MPI_REAL8, work_loads, 1, MPI_REAL8, 0,  MPI_COMM_WORLD, ierr )  ! Gather work integrals
  call MPI_GATHER(npp, 1, MPI_INTEGER, npps, 1, MPI_INTEGER, 0,  MPI_COMM_WORLD, ierr )  ! Gather particle distn

  timestamp = itime + itime_start

  if (me ==0 .and. mod(itime,idump)==0) then
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
     write (ipefile,101) force_const, delta_t

     do i=p_start,p_finish
        write (ipefile,102) pelabel(i), pepid(i), ax(i), ay(i), az(i)

        !        write (ipefile,102) i, fx(i), fy(i), fz(i)
     end do

101  format('Tree forces:'/'   p     owner    ax         ay      az  ',2f8.2)
102  format(1x,2i7,3(1pe14.5))

  endif



end subroutine forces_bal
