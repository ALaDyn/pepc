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
!   ** Returns fields Ex, Ey, Ez and potential pot excluding external terms **
!
!
!  ===================================================================


subroutine forces(p_start,p_finish,delta_t, t_walk, t_force)

  use physvars
  use treevars
  use utils
  implicit none

  real, intent(in) :: delta_t
  integer, intent(in) :: p_start,p_finish  ! min, max particle nos.

  integer, parameter :: npassm=100000 ! Max # passes - will need npp/nshortm

  integer :: p, i, j, npass, jpass, ip1, nps,  max_npass,nshort_list, ipe
  real :: t_walk, t_force, t1, t2, t3  ! timing integrals
  integer :: pshortlist(nshortm),nodlist(nintmax),nshort(npassm),pstart(npassm)
  integer :: hashaddr ! Key address 
  real :: work_loads(num_pe)  ! Load balance array
  integer :: npps(num_pe)  ! Particle distrib amoung PEs
  integer :: max_local,  timestamp
  integer :: ierr

  real :: fsx, fsy, fsz, phi, phi_coul, ex_coul, ey_coul, ez_coul
  real :: ax_ind, ay_ind, az_ind, bx_ind, by_ind, bz_ind
  real :: Epon_x, Epon_y, Epon_z, Phipon
  real :: xd, yd, zd  ! positions relative to centre of laser spot
  real :: work_local, load_average, load_integral, total_work, average_work
  integer :: total_parts
  character(30) :: cfile, ccol1, ccol2
  character(4) :: cme
  integer :: key2addr        ! Mapping function to get hash table address from key


  t_walk=0.
  t_force=0.
  max_local = 0   ! max length of interaction list
  work_local = 0  ! total workload

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

     if (walk_debug) then
       	write(*,*) 'pass ',jpass,' of ',max_npass,': # parts ',ip1,' to ',ip1+nps-1
        write(ipefile,*) 'pass ',jpass,' # parts ',ip1,' to ',ip1+nps-1
     endif

     !  build interaction list: 
     ! tree walk returns intlist(1:nps), nodelist(1:nps) for particles on short list

     call cputime(t1)
     call tree_walk(pshortlist,nps,jpass,theta,itime)
     call cputime(t2)
     t_walk = t_walk + t2-t1

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
           call sum_force(p, nterm(i), nodelist( 1:nterm(i),i), eps, ex_coul, ey_coul, ez_coul, phi_coul )

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


  if (itime>0 .and. beam_config==4) focus(1) = x_crit  ! laser tracks n_c

  ! Include force from laser on electrons - ES scheme
  if (beam_config >= 3 .and. beam_config <=5 ) then
     do p = p_start, p_finish
        if (q(p)<0) then
           xd = x(p)-focus(1)
           yd = y(p)-focus(2)
           zd = z(p)-focus(3)

           laser_model: select case(beam_config)

           case(3)  ! Uniform sinusoid in z (s-pol)
              Epon_x = 0.
              Epon_y = 0.
              Epon_z = vosc*omega*sin(omega*tlaser)

           case(4)  ! standing wave fpond
              call fpond( tlaser, tpulse,sigma,vosc,omega,rho_upper, &
                xd,yd,zd,epon_x,epon_y,epon_z,phipon)

           case(5)  ! propagating fpond
              call laser_bullet( tlaser, focus(1), tpulse,sigma,vosc,omega, & 
                   xd,yd,zd,epon_x,epon_y,epon_z,phipon)

           case default  ! no laser
              Epon_x=0
              Epon_y=0
              Epon_z=0

           end select laser_model


           Ex(p) = Ex(p) + Epon_x
           Ey(p) = Ey(p) + Epon_y
           Ez(p) = Ez(p) + Epon_z
        endif

     end do
  endif


  call MPI_ALLREDUCE(max_local, max_list_length, 1, MPI_INTEGER, MPI_MAX,  MPI_COMM_WORLD, ierr )
  call MPI_GATHER(work_local, 1, MPI_REAL8, work_loads, 1, MPI_REAL8, 0,  MPI_COMM_WORLD, ierr )  ! Gather work integrals
  call MPI_GATHER(npp, 1, MPI_INTEGER, npps, 1, MPI_INTEGER, 0,  MPI_COMM_WORLD, ierr )  ! Gather particle distn

  timestamp = itime + itime_start

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
     write (ipefile,101) force_const, delta_t

     do i=p_start,p_finish
        write (ipefile,102) pelabel(i), pepid(i), & 
             q(i), m(i), ux(i), Ex(i),  Ax(i), Axo(i), Bx(i)
     end do

101  format('Tree forces:'/'   p     owner    ux   Ex         Ax      Axo      Bx  ',2f8.2)
102  format(1x,2i7,7(1pe14.5))

  endif

end subroutine forces




