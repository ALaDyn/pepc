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
!VAMPINST include
      INCLUDE 'VTcommon.h'
      INTEGER VTIERR
!

  real, intent(in) :: delta_t
  integer, intent(in) :: p_start,p_finish  ! min, max particle nos.

  integer, parameter :: npassm=1000 ! Max # passes - will need npp/nshortm

  integer :: p, i, j, npass, jpass, ip1, nps, max_local, max_npass
  real :: t_walk, t_force, t1, t2  ! timing integrals
  integer pshortlist(nshortm),nshort(0:npassm)


  real, dimension(nppm) :: fx, fy, fz

  real :: fsx, fsy, fsz,  Ex, Ey, Ez
  real :: xd, yd, zd  ! positions relative to centre of laser spot
  logical :: force_debug=.true.
  real :: load_average, load_integral

!VAMPINST subroutine_start
       CALL VTENTER(IF_forces,VTNOSCL,VTIERR)
!      write(*,*) 'VT: forces S>',VTIERR,
!     *    IF_forces,ICLASSH
!
  Ex = 0.
  Ey = 0.
  Ez = 0.
  fx=0.
  fy=0.
  fz=0.
  t_walk=0.
  t_force=0.

  !  # passes needed to process all particles
  npass = npart/num_pe/nshortm   ! ave(npp)/nshortm   - make nshortm a power of 2
  load_average = SUM(work(1:npp))/npass   ! Ave. workload per pass: same for all PEs if already load balanced
  nshort(1:npass+1) = 0

  load_integral = 0.
  jpass = 1
  do i=1,npp
     load_integral = load_integral + work(i)   ! integrate workload
     if (load_integral >= load_average * jpass .or. i==npp) then
        nshort(jpass) = i
        jpass = jpass + 1
     endif
  end do

  if (jpass-1 > npass) then
     write(*,*) 'PE',me,' missed some:'
     nshort(npass) = nshort(npass+1)
  endif

  write (ipefile,*) 'Shortlists: ',(nshort(j),j=1,npass+1)
  max_npass = npass

! determine global max # passes (npass need not be the same across PEs if npp varies)
!  call MPI_ALLREDUCE(npass, max_npass, 1, MPI_INTEGER, MPI_MAX,  MPI_COMM_WORLD, ierr )

!  if (npass < max_npass) then
!     do j=npass, max_npass-1
!	nshort(j) = 0  ! zero # particles for dummy passes
!     end do
!  endif

  ip1 = 1
  do jpass = 1,max_npass
     !  make short-list
     nps = nshort(jpass)
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
	call sum_force(p, nterm(i), nodelist( 1:nterm(i),i ), fsx, fsy, fsz )

	fx(p) = force_const * q(p) * fsx
	fy(p) = force_const * q(p) * fsy
	fz(p) = force_const * q(p) * fsz
        work(p) = nterm(i)        ! Should really compute this in sum_force to allow for leaf/twig terms
     end do
     call cputime(t3)
     t_force = t_force + t3-t2
     ip1 = ip1 + nshort(jpass)    ! start number for next pass

  end do



  ! Include ponderomotive force from laser on electrons
  if (beam_config ==4 ) then
     do p = p_start, p_finish
        if (q(p)<0) then
           xd = x(p)-focus(1)
           yd = y(p)-focus(2)
           zd = z(p)-focus(3)

           call fpond( (itime+itime_start)*delta_t,tpulse,sigma,vosc,omega,-xd,yd,zd,epon_x,epon_y,epon_z,phipon)

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


!VAMPINST subroutine_end
       CALL VTLEAVE(ICLASSH,VTIERR)
!      write(*,*) 'VT: forces S<',VTIERR,ICLASSH
!
end subroutine forces
