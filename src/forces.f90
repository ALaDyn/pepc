!  ===================================================================
!
!                              FORCES
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


subroutine forces(p_start,p_finish,delta_t, t_walk, t_force)


  use treevars
  use utils

  implicit none
  real, intent(in) :: delta_t
  integer, intent(in) :: p_start,p_finish  ! min, max particle nos.

  integer, parameter :: npassm=100000 ! Max # passes - will need npp/nshortm

  integer :: p, i, j, npass, jpass, ip1, nps, max_local,  max_npass, timestamp
  real :: t_walk, t_force, t1, t2, t3  ! timing integrals
  integer pshortlist(nshortm),nshort(0:npassm)
  real :: work_loads(num_pe)  ! Load balance array
  integer :: npps(num_pe)  ! Particle distrib amoung PEs
  real :: work_local, load_average, load_integral, total_work, average_work
  integer :: total_parts
  character(30) :: cfile
  character(4) :: cme


  real, dimension(nppm) :: fx, fy, fz

  real :: fsx, fsy, fsz,  Ex, Ey, Ez, phi
  real :: Epon_x, Epon_y, Epon_z, Phipon
  real :: xd, yd, zd  ! positions relative to centre of laser spot

  !  Sum forces:  uses interaction list from INTLIST
  !  ------------------------------------------------

  Ex = 0.
  Ey = 0.
  Ez = 0.
  fx=0.
  fy=0.
  fz=0.
  t_walk=0.
  t_force=0.
  max_local = 0   ! max length of interaction list
  work_local = 0.  ! total workload on this PE
 
  !  # passes needed to process all particles

  if  (mod(npp,nshortm).eq.0) then
     npass = npp/nshortm
     !  process nlim particles at a time
     do j=0,npass-1
        nshort(j)=nshortm
     end do
  else
     npass = npp/nshortm + 1
     !  process nlim particles at a time
     do j=0,npass-2
        nshort(j)=nshortm
     end do
     !  last bunch is 'remainder' of npp/nshortm
     nshort(npass-1) = mod(npp,nshortm)
  endif


! determine global max # passes (npass need not be the same across PEs if npp varies)
  call MPI_ALLREDUCE(npass, max_npass, one, MPI_INTEGER, MPI_MAX,  MPI_COMM_WORLD, ierr )

  if (npass < max_npass) then
     do j=npass, max_npass-1
	nshort(j) = 0  ! zero # particles for dummy passes
     end do
  endif

  do jpass = 0,max_npass-1
     nps = nshort(jpass)
     ip1 = jpass*nshortm + 1
     !  make short-list
     pshortlist(1:nps) = (/ (ip1+i-1, i=1,nps) /)

     if (force_debug) then
	write(*,*) 'pass ',jpass+1,' # parts ',ip1,' to ',ip1+nps-1
	write(ipefile,*) 'pass ',jpass+1,' # parts ',ip1,' to ',ip1+nps-1
     endif

     !  build interaction list: 
     ! tree walk returns intlist(1:nps), nodelist(1:nps) for particles on short list

     call cputime(t1)
     call tree_walk(pshortlist,nps)
     call cputime(t2)
     t_walk = t_walk + t2-t1

     do i = 1, nps
	fsx=0. ! partial forces for particle p
	fsy=0.
	fsz=0.
        p = pshortlist(i)    ! local particle index

        !  compute forces of particle p from its interaction list
        call sum_force(p, nterm(i), nodelist( 1:nterm(i),i ), fsx, fsy, fsz, phi )

        pot(p) = force_const * phi
        fx(p) = force_const * q(p) * fsx
        fy(p) = force_const * q(p) * fsy
        fz(p) = force_const * q(p) * fsz
        work(p) = nterm(i)        ! Should really compute this in sum_force to allow for leaf/twig terms
        work_local = work_local+nterm(i) 
     end do
     call cputime(t3)
     t_force = t_force + t3-t2
     max_local = max( max_local,maxval(nterm(1:nps)) )  ! Max length of interaction list

  end do



  ! Include ponderomotive force from laser on electrons
  if (beam_config ==4 ) then
     do p = p_start, p_finish
        if (q(p)<0) then
           xd = x(p)-focus(1)
           yd = y(p)-focus(2)
           zd = z(p)-focus(3)

           call fpond( tlaser, tpulse,sigma,vosc,omega,-xd,yd,zd,epon_x,epon_y,epon_z,phipon)

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

  if (me ==0 .and. mod(itime,idump)==0) then
     total_work = SUM(work_loads)
     average_work = total_work/num_pe
     timestamp = itime + itime_start
     cme = achar(timestamp/1000+48) // achar(mod(timestamp/100,10)+48) &
          // achar(mod(timestamp/10,10)+48) // achar(mod(timestamp,10)+48)
     cfile="load_"//cme//".dat"
     total_parts=SUM(npps)
     open(60, file=cfile)
     write(60,'(a,i8,2(a,1pe12.5))')  '! Parts: ',total_parts,' Work: ',total_work, &
          ' Ave. work:',average_work
     write(60,'(2i8,f15.3))')  (i-1,npps(i),work_loads(i)/average_work,i=1,num_pe)
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


end subroutine forces
