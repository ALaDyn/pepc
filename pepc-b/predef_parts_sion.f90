! ==============================================
!
!                PREDEF_PARTS
!
!  Read particle distribution from a sion file
!
! ==============================================

subroutine predef_parts

  use physvars
  use treevars
  implicit none
  include 'mpif.h'

  integer :: i, j, idummy=0, ierr
  character(30) :: cinfile, cdump, cfile, cinfofile
  integer :: ner, nir, np_beamr, npartr, iconf, iens, timestamp
  real :: epsr, thetar, xlr, ylr, zlr, boxr
  real :: omegar, lambdar
  real :: axdum, aydum, azdum,phidum, bdum
  integer :: ioffset, i1, i2, npp_partial, npp_total, ipass, me_read, nrest, nadd
  integer :: nslice_e, nslice_i
  logical :: stopflag=.false.


  integer*8 chunksize, size_info, size1, size2, bread
  ! 2MB GPFS
  integer 	fsblksize, feof
  integer 	sid, globalrank
  real*8 realdummy

! Temp variable use for the info blockt
  TYPE, BIND(C) :: Infoblock
  integer, dimension(7) :: intarr
  real, dimension(18) :: realarr
  real, dimension(3) :: plasma_centre, focus
  END TYPE Infoblock

  TYPE(Infoblock) :: infoblk

  !get the last timestep
  if (me == 0) then
     ! input file for PE me
     cinfile="parts_info.in"
     open(80,file=cinfile)
     ! Root reads info block in run directory to check that run parameters correct
     read(80,'(9x,i8/)') itime_start  ! itime_start from the info block - skip variable names
     close(80)
  endif
  call MPI_BCAST( itime_start, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)

  ! get filename suffix from dump counter
  do i=0,4
     cdump(6-i:6-i) =  achar(mod(itime_start/10**i,10) + 48)
  end do
  cdump(1:1) = achar(itime_start/10**5 + 48)

  cinfile="dumps/parts_dump."//cdump(1:6)
  call fsion_paropen_mpi(trim(cinfile),"br",MPI_COMM_WORLD,chunksize,fsblksize,me,sid)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  size_info = sizeof(infoblk%intarr)+sizeof(infoblk%realarr)+sizeof(infoblk%plasma_centre)+sizeof(infoblk%focus)

  call fsion_read(infoblk, 1, size_info, sid, bread)
  if (bread /= size_info) then
   write(*,*) 'Error reading Infoblock. Read: ', bread
  else
  	write(167,*)'Read infoblock: ',bread
  end if

  itime_start = infoblk%intarr(1)
  npartr = infoblk%intarr(2)
  ner = infoblk%intarr(3)
  nir = infoblk%intarr(4)
  np_beamr = infoblk%intarr(5)
  iconf = infoblk%intarr(6)
  iens = infoblk%intarr(7)
  xlr = infoblk%realarr(1)
  ylr = infoblk%realarr(2)
  zlr = infoblk%realarr(3)
  boxr= infoblk%realarr(3)
  thetar = infoblk%realarr(5)
  epsr= infoblk%realarr(4)
  tlaser = infoblk%realarr(6)
  trun = infoblk%realarr(7)
  omegar = infoblk%realarr(8)
  lambdar = infoblk%realarr(9)
  qe = infoblk%realarr(10)
  qi = infoblk%realarr(11)
  mass_e = infoblk%realarr(12)
  mass_i = infoblk%realarr(13)
  Zion = infoblk%realarr(14)
  a_ii = infoblk%realarr(15)
  Vplas = infoblk%realarr(16)
  Aplas = infoblk%realarr(17)
  Qplas = infoblk%realarr(18)
  do i=1,3
  	plasma_centre(i) = infoblk%plasma_centre(i)
  	focus(i) = infoblk%focus(i)
  end do

  if (me == 0) then
     write(6,'(/a/a/a,i5)') 'RESTART:','Reading run data from info block: parts_info.in','Timestep:',itime_start
!     write(*,'(7(a12,i12/),9(a12,f12.5/),9(a12,1pe12.5/),2(a12,3(1pe12.5)/))')  &    ! info block - skip variable names
!          'Start time: ', itime_start, &
!          '# particles: ', npartr, &
!	  	  '# electrons: ', ner, &
!          '# ions: ', nir, &
!          '# beam: ', np_beamr, &
!          'Geom  : ', iconf, &
!          'Scheme: ', iens, &
!	  	  'Box_x: ',xlr,&
!          'Box_y: ', ylr, &
!          'Box_z: ', zlr, &
!          'eps: ', epsr,  &
!  	  	  'theta: ', thetar, &
!          'tlaser: ',tlaser, &
!          'trun: ',trun, &
!          'omega:',omegar, &
!          'lambda:', lambdar, &
!          'qe:',qe, &
!          'qi:',qi, &
!          'mass_e:',mass_e, &
!          'mass_i:',mass_i, &
!          'Zion:',Zion, &
!          'a_ii:', a_ii, &
!          'Vplas: ',Vplas, &
!          'Aplas: ',Aplas, &
!          'Qplas: ',Qplas, &
!          'centre: ',plasma_centre(1:3), &
!          'focus: ',focus(1:3)

     if (ner /= ne .or. nir /= ni) then
        write(*,*) '*** Particle nos. in input deck do not match those in restart file parts_info.in'
        stopflag=.true.
     endif

     if (epsr /= eps) then
        write(*,*) '*** Warning: potential cutoff eps changed - check inputs( eps=',eps,' epsr=',epsr, ')'
     endif

     if (thetar /= theta) then
        write(*,*) '*** Warning: MAC changed - check inputs( theta=',theta,' thetar=',thetar, ')'
     endif

     if (iconf /= target_geometry) then
        write(*,*) '*** Warning: Target geometry in restart file ',iconf, &
             ' does not match value ',target_geometry,' in run.h - check inputs'
     endif
  endif

  call MPI_BARRIER( MPI_COMM_WORLD, ierr)   ! Synchronize first

  if (stopflag) then
    if (me==0) write (*,*) "*** Stopping program"
    call MPI_FINALIZE(ierr)
    call closefiles
    stop
  endif

!  if (ncpu_merge < 0 ) then

!     ! num_pe > # data sets, so carve up restart data
!     if (ncpu_merge==num_pe) npp_total = ne+ni  ! Override total # if dataset to be split amoung all PEs

!     npp = npp_total/ncpu_merge
!     nrest = mod(npp_total,ncpu_merge)  ! Remainder
!
!     if (mod(me,ncpu_merge)==ncpu_merge-1) then
!        nadd = nrest  ! Put remainder on last PE
!     else
!        nadd = 0
!     endif
!
!     if(mod(me,1000)==0) write(*,*) 'PE ',me,': Reading ',npp+nadd,' particles out of ',npp_total,' from ',cinfile
! ! Record # particles read in openfiles.out
!     write(90,*) 'PE ',me,': Reading ',npp+nadd,' particles out of ',npp_total,' from ',cinfile

     !  Skip dummy blocks up to previous PEs
!     nslice_e=0
!     nslice_i = 0
!     nslice = 0
!     if (debug_level>0) write(ipefile,*) 'skip pass ',j

	npp = npartr

  ! npp * sizepof(real*8)
	size1 = npp*sizeof(realdummy)
  ! npp*sizeof(integer)
	size2 = npp*sizeof(sid)

!    if(me == 0) then
!      write(6,*)'**************** npp=',npp,' **************'
!  	  write(6,*)'**************** size_info=',size_info,' **************'
!      write(6,*)'**************** size1    =',size1,' **************'
!      write(6,*)'**************** size2    =',size2,' **************'
!   	  write(6,*)'**************** chunksize=',chunksize,' **************'
!    end if

  	! Read X
  	call fsion_read(x, 1, size1, sid, bread)
  	if (bread /= size1) write(167,*) 'Error reading x(). Read: ', bread

	! Read Y
  	call fsion_read(y, 1, size1, sid, bread)
  	if (bread /= size1) write(167,*) 'Error reading y(). Read: ', bread

	! Read Z
  	call fsion_read(z, 1, size1, sid, bread)
  	if (bread /= size1) write(167,*) 'Error reading z(). Read: ', bread

	! Read UX
  	call fsion_read(ux, 1, size1, sid, bread)
  	if (bread /= size1) write(167,*) 'Error reading ux(). Read: ', bread

	! Read UY
  	call fsion_read(uy, 1, size1, sid, bread)
  	if (bread /= size1) write(167,*) 'Error reading uy(). Read: ', bread

	! Read UY
  	call fsion_read(uz, 1, size1, sid, bread)
  	if (bread /= size1) write(167,*) 'Error reading uz(). Read: ', bread

  	! Read Q
  	call fsion_read(q, 1, size1, sid, bread)
  	if (bread /= size1) write(167,*) 'Error reading q(). Read: ', bread

  	! Read M
  	call fsion_read(m, 1, size1, sid, bread)
  	if (bread /= size1) write(167,*) 'Error reading m(). Read: ', bread

  	! NOT NEEDED Read Ex
  	call fsion_read(Ex, 1, size1, sid, bread)
  	if (bread /= size1) write(167,*) 'Error reading Ex(). Read: ', bread

  	! NOT NEEDED Read Ey
  	call fsion_read(Ey, 1, size1, sid, bread)
  	if (bread /= size1) write(167,*) 'Error reading Ey(). Read: ', bread

  	! NOT NEEDED Read Ez
  	call fsion_read(Ez, 1, size1, sid, bread)
  	if (bread /= size1) write(167,*) 'Error reading Ez(). Read: ', bread

  	! NOT NEEDED Read pot
  	call fsion_read(pot, 1, size1, sid, bread)
  	if (bread /= size1) write(167,*) 'Error reading pot(). Read: ', bread

  	! NOT NEEDED Read pepid
  	call fsion_read(pepid, 1, size2, sid, bread)
  	if (bread /= size2) write(167,*) 'Error reading pepid(). Read: ', bread

  	! Read pelabel
  	call fsion_read(pelabel, 1, size2, sid, bread)
  	if (bread /= size2) write(167,*) 'Error reading pelabel(). Read: ', bread

!    do i=1,npp
!       if (beam_config.eq.5 .and. q(i)>0 .and. x(i) < window_min+dt .and. x(i) > window_min) then
! 		! create rezoning slice for wakefield mode - first few blocks should be sufficient
!          nslice = nslice+1
!          xslice(nslice) = x(i)+x_plasma ! Include offset for new slice
!          yslice(nslice) = y(i)
!          zslice(nslice) = z(i)
!       endif
!    end do
!  else
!  	 ! num_pe < # data sets, so merge together
!     ioffset = 0
!  endif

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call fsion_parclose_mpi(sid,MPI_COMM_WORLD,ierr)

!  if (me==0) write(*,*) "Dumping particle data for checking"
!  cfile="data/pe"//csubme//"/parts_dump."//cdump(1:6)
!  open (60,file=cfile)
!  write(60,'((12(1pe14.5),2i9))')  &
!       (x(i), y(i), z(i), ux(i), uy(i), uz(i), q(i), m(i), &
!        Ex(i), Ey(i), Ez(i), &  ! electric field
!        pot(i), &  ! potential
!        pepid(i), pelabel(i), i=1,npp)
!  close(60)

!  call MPI_BARRIER( MPI_COMM_WORLD, ierr)   ! Synchronize first
!  call MPI_BCAST( nslice, 1, MPI_INTEGER, num_pe-1, MPI_COMM_WORLD,ierr)

  pepid(1:npp) = me                ! processor ID
  Ex(1:npp) = 0.       ! zero fields until first force comp.
  Ey(1:npp) = 0.
  Ez(1:npp) = 0.
  work(1:npp) = 1.

  ! Rescale velocities if different temperature required
  if (T_scale /= 1) then
     if (me==0) write(*,*) 'Rescaling temperature by ',T_scale,' to ',Te_keV
     do i=1,npp
        if (q(i)<0) then
           ux(i) = ux(i)*sqrt(T_scale)
           uy(i) = uy(i)*sqrt(T_scale)
           uz(i) = uz(i)*sqrt(T_scale)
        endif
     end do
  endif

  if (me==0) write(*,*) "Finished reading particle data"

end subroutine predef_parts

