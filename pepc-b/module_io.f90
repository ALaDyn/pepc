!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>  Encapsulates any file i/o, opening, closing etc.
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_io
      implicit none
      save
      private

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      character*4, public :: csubme   !< Character string of data subdirectory 'data/peXXXX'

      integer, parameter, public :: file_tree_out    = 15
      integer, parameter, public :: file_pepc_out    = 24
      integer, parameter, public :: file_domains_dat = 70
      integer, parameter, public :: file_laser_dat   = 71
      integer, parameter, public :: file_energy_dat  = 75
      integer, parameter, public :: file_parts_info_in = 80
      integer, parameter, public :: file_memory_dat  = 59

      integer, parameter, public :: file_ipefile     = 20

      integer, parameter, public :: file_stdout = 6

      integer, parameter, public :: file_tempfile_1  = 60
      integer, parameter, public :: file_tempfile_2  = 61
      integer, parameter, public :: file_tempfile_3  = 62

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      public openfiles
      public flushfiles
      public closefiles
      public stamp
      public dump		!< Particle dump
      public predef_parts	!< Particle read from checkpoint
      public dump_fields	!< Write out 3D fields and cuts
      public slices		!< Write out 1D cuts

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      contains

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!>
		!> Open files for result and debugging output
		!> - Rank 0:
		!>          - tree.out (\e 15): Tree stats
		!>          - pepc.out (\e 24): Physics log
		!>          - domains.dat (\e 70):
		!>          - laser.dat (\e 71): Laser parameters
        !>          - energy.dat (\e 75): energies
        !>          - memory.dat (\e 59): memory status information
		!>          .
		!> - all ranks (only if \code(debug_level>2 .and. idump>0)\endcode) :
		!>          - data/out.1234 (\e 20):
		!>
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		subroutine openfiles

		  use module_physvars
		  implicit none
		  character(30) :: cfile


		  if (my_rank == 0) then
		     !  master diagnostics output
		     open(file_tree_out, file='tree.out')  ! Tree stats
		     open(file_pepc_out, file='pepc.out') ! Physics log

		     open(file_domains_dat, file='domains.dat')
		     open(file_laser_dat,   file='laser.dat')       ! laser parameters
		     open(file_energy_dat,  file='energy.dat')      ! energies

		     write(*,*) 'debug level: ',debug_level,' idump',idump
		  endif

		  !  stdout for PE my_rank
		  if (debug_level > 0) then
            call system("mkdir -p " // "diag")
		    write(cfile,'("diag/diag_",i6.6,".dat")') my_rank
		    open(file_ipefile, file=cfile,STATUS='UNKNOWN', POSITION = 'APPEND')
		  endif

		  ipefile = file_ipefile ! copy file handle to core

		end subroutine openfiles


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine flushfiles
          use module_physvars
          implicit none

          if (my_rank == 0) then
             !  master diagnostics output
             flush(file_tree_out)
             flush(file_pepc_out)

             flush(file_domains_dat)
             flush(file_laser_dat)
             flush(file_energy_dat)

             flush(file_memory_dat)
          endif

        end subroutine flushfiles


		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!>
		!> Tidy up O/P files
		!>
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		subroutine closefiles
		  use module_physvars
		  implicit none

		  if (my_rank == 0) then
		     close(file_tree_out)
             close(file_pepc_out)
		     close(file_domains_dat)
		     close(file_laser_dat)
             close(file_energy_dat)
             close(file_memory_dat)
		     close(90)
		  endif
		  if (debug_level > 0) close(file_ipefile)
		  close(file_parts_info_in)  ! initial particle data


		end subroutine closefiles


		!     =========================
		!
		!     Time stamp
		!
		!     =========================

		subroutine stamp(istream,ibegin)
		  implicit none

		  character :: cdate*8, ctime*10, czone*5
		  integer :: ibegin
		  integer :: istream

		     !      call DATE_AND_TIME(cdate,ctime,czone,vals)
		     call DATE_AND_TIME(cdate,ctime,czone)

		     if (ibegin.eq.1) then

		        write(istream,'(//a20,a12/a20,a12/a20,a12//)') 'PEPC run on ' &
		             ,cdate(7:8)//'/'//cdate(5:6)//'/'//cdate(1:4) &
		             ,'Time: ',ctime(1:2)//':'//ctime(3:4)//':'//ctime(5:6),' GMT+',czone

		     else
		        write(istream,'(a,a9)') 'Finished run at time: ',ctime(1:2)//':'//ctime(3:4)//':'//ctime(5:6)
		     endif
		end subroutine stamp



#ifdef SION
! ======================
!
! Sion DUMP 
!
! Gather and write out particle data for restart
! using sionlib
!
! Only compiled if SION flag set - needs sionlib!
! ======================

subroutine dump(timestamp)

  use module_physvars
  use module_particle_props

  implicit none
  include 'mpif.h'

  character(30) :: cfile, cinfofile, cfile_new
  character(6) :: cdump
  character(5) :: cme
  integer, intent(in) :: timestamp
  integer :: i, j, ioffset, idummy=0, ierr
  integer, save :: icall=0
  integer :: numfiles=0
  real :: simtime

  integer*8 chunksize, size_info, size1, size2, bwrote
  ! 2MB GPFS
  integer fsblksize
  integer sid

  real*8 realdummy

! Temp variable use for the info blockt
  TYPE, BIND(C) :: Infoblock
  integer, dimension(7) :: intarr
  real, dimension(18) :: realarr
  real, dimension(3) :: plasma_centre, focus
  END TYPE Infoblock

  TYPE(Infoblock) :: infoblk

  infoblk%intarr(1)=timestamp
  infoblk%intarr(2)=np_local
  infoblk%intarr(3)=ne
  infoblk%intarr(4)=ni
  infoblk%intarr(5)=np_beam
  infoblk%intarr(6)=target_geometry
  infoblk%intarr(7)=scheme

  infoblk%realarr(1)=xl
  infoblk%realarr(2)=yl
  infoblk%realarr(3)=zl
  infoblk%realarr(4)=eps
  infoblk%realarr(5)=theta
  infoblk%realarr(6)=tlaser
  infoblk%realarr(7)=trun
  infoblk%realarr(8)=omega
  infoblk%realarr(9)=lambda
  infoblk%realarr(10)=qe
  infoblk%realarr(11)=qi
  infoblk%realarr(12)=mass_e
  infoblk%realarr(13)=mass_i
  infoblk%realarr(14)=Zion
  infoblk%realarr(15)=a_ii
  infoblk%realarr(16)=Vplas
  infoblk%realarr(17)=Aplas
  infoblk%realarr(18)=Qplas
  do i=1,3
   infoblk%plasma_centre(i) = plasma_centre(i)
   infoblk%focus(i)=focus(i)
  end do

  simtime = dt*timestamp


  ! get filename suffix from dump counter
  do i=0,4
     cdump(6-i:6-i) = achar(mod(timestamp/10**i,10) + 48)
  end do
  cdump(1:1) = achar(timestamp/10**5 + 48)

! Write particles info file
  if(my_rank == 0) then
   cinfofile="dumps/parts_info."//cdump(1:6)
   open (60,file=cinfofile)
   write(60,'(7(a9,i8/),9(a9,f12.5/),9(a9,1pe12.5/),2(a9,3(1pe12.5)/))') & ! info block
     'itime=',timestamp, 'np_local=',np_local, &
     'ne=',ne, 'ni=',ni, 'npbeam=',np_beam, 'geometry=', target_geometry, &
     'scheme=',scheme, &
     'xl=',xl, 'yl=',yl, 'zl=',zl,  &
     'eps=', eps, 'theta=',theta,' tlaser= ',tlaser,' trun= ',trun, &
     'omega=',omega,'lambda=',lambda,'  qe=',qe,'  qi=',qi, &
     'mass_e=',mass_e,'mass_i=',mass_i,'Zion=',Zion,'a_ii=',a_ii, &
     'Vplas=',Vplas,'Aplas=',Aplas,'Qplas=',Qplas, &
     'centre=',plasma_centre(1:3),'focus=',focus(1:3)
   close (60)
  endif

! Calc the size of the chunk
! Info block: 7*sizeof(integer)+24*sizeof(real) = size(infoblock)
! Particles: 12*np_local*sizeof(real*8) + 2*np_local*sizeof(integer)
! should be:

  size_info = sizeof(infoblk%intarr)+sizeof(infoblk%realarr)+sizeof(infoblk%plasma_centre)+sizeof(infoblk%focus)

  ! np_local * sizepof(real*8)
  size1 = np_local*sizeof(realdummy)
  ! np_local*sizeof(integer)
  size2 = np_local*sizeof(sid)
  chunksize = 12*size1 + 2*size2 + size_info
  fsblksize = 2*1024*1024

! if(my_rank == 0) then
! write(6,*)'**************** np_local=',np_local,' **************'
! write(6,*)'**************** size_info=',size_info,' **************'
! write(6,*)'**************** size1    =',size1,' **************'
! write(6,*)'**************** size2    =',size2,' **************'
! write(6,*)'**************** chunksize=',chunksize,' **************'
! end if

! Particles dump file
  cfile="dumps/parts_dump."//cdump(1:6)

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call fsion_paropen_mpi(trim(cfile),'bw',numfiles,MPI_COMM_WORLD,MPI_COMM_WORLD,chunksize,fsblksize,my_rank,cfile_new,sid)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  ! Write the inforblock
  call fsion_write(infoblk, 1, size_info, sid, bwrote)
  if (bwrote /= size_info) then
   write(*,*) 'Error writing Infoblock. Wrote: ', bwrote
  else
   write(167,*)'Wrote infoblock: ',bwrote
  end if

  ! Write X
  call fsion_write(x, 1, size1, sid, bwrote)
  if (bwrote /= size1) write(167,*) 'Error writing x(). Wrote: ', bwrote

  ! Write Y
  call fsion_write(y, 1, size1, sid, bwrote)
  if (bwrote /= size1) write(167,*) 'Error writing y(). Wrote: ', bwrote

  ! Write Z
  call fsion_write(z, 1, size1, sid, bwrote)
  if (bwrote /= size1) write(167,*) 'Error writing z(). Wrote: ', bwrote

  ! Write UX
  call fsion_write(ux, 1, size1, sid, bwrote)
  if (bwrote /= size1) write(167,*) 'Error writing ux(). Wrote: ', bwrote

  ! Write UY
  call fsion_write(uy, 1, size1, sid, bwrote)
  if (bwrote /= size1) write(167,*) 'Error writing uy(). Wrote: ', bwrote

  ! Write UY
  call fsion_write(uz, 1, size1, sid, bwrote)
  if (bwrote /= size1) write(167,*) 'Error writing uz(). Wrote: ', bwrote

  ! Write Q
  call fsion_write(q, 1, size1, sid, bwrote)
  if (bwrote /= size1) write(167,*) 'Error writing q(). Wrote: ', bwrote

  ! Write M
  call fsion_write(m, 1, size1, sid, bwrote)
  if (bwrote /= size1) write(167,*) 'Error writing m(). Wrote: ', bwrote

  ! Write Ex
  call fsion_write(Ex, 1, size1, sid, bwrote)
  if (bwrote /= size1) write(167,*) 'Error writing Ex(). Wrote: ', bwrote

  ! Write Ey
  call fsion_write(Ey, 1, size1, sid, bwrote)
  if (bwrote /= size1) write(167,*) 'Error writing Ey(). Wrote: ', bwrote

  ! Write Ez
  call fsion_write(Ez, 1, size1, sid, bwrote)
  if (bwrote /= size1) write(167,*) 'Error writing Ez(). Wrote: ', bwrote

  ! Write pot
  call fsion_write(pot, 1, size1, sid, bwrote)
  if (bwrote /= size1) write(167,*) 'Error writing pot(). Wrote: ', bwrote

  ! Write pepid
  call fsion_write(pepid, 1, size2, sid, bwrote)
  if (bwrote /= size2) write(167,*) 'Error writing pepid(). Wrote: ', bwrote

  ! Write pelabel
  call fsion_write(pelabel, 1, size2, sid, bwrote)
  if (bwrote /= size2) write(167,*) 'Error writing pelabel(). Wrote: ', bwrote

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call fsion_parclose_mpi(sid,MPI_COMM_WORLD,ierr)

! add to runstamp the stamp and print to stdout
  if (my_rank.eq.0) then
    open (62,file="runstamp") ! time stamp
    write(62,'(a)') cdump(1:6)
    close (62)

! write(166,'(2a)') 'Particle dump',cdump(1:6)
! write(166,'(//a/7(a9,i8/),9(a9,f12.5/),9(a9,1pe12.5/),2(a9,3(1pe12.5)/))') 'PARTICLE DUMP:', & ! info block
! 'itime=',timestamp, 'np_local=',np_local, &
! 'ne=',ne, 'ni=',ni, 'npbeam=',np_beam, 'geometry=', target_geometry, &
! 'scheme=',scheme, &
! 'xl=',xl, 'yl=',yl, 'zl=',zl, &
! 'eps=', eps, 'theta=',theta,' tlaser= ',tlaser,' trun= ',trun, &
! 'omega=',omega,'lambda=',lambda,'  qe=',qe,'  qi=',qi, &
! 'mass_e=',mass_e,'mass_i=',mass_i,'Zion=',Zion,'a_ii=',a_ii, &
! 'Vplas=',Vplas,'Aplas=',Aplas,'Qplas=',Qplas, &
! 'centre=',plasma_centre(1:3),'focus=',focus(1:3)

  endif
  icall = icall + 1

end subroutine dump

! ==============================================
!
!                PREDEF_PARTS
!
!  Read particle distribution from a sion file
!
! ==============================================

subroutine predef_parts

  use module_physvars
  use module_particle_props

  implicit none
  include 'mpif.h'

  integer :: i, j, idummy=0, ierr
  character(30) :: cinfile, cdump, cfile, cinfofile, cfile_new
  integer :: ner, nir, np_beamr, npartr, iconf, iens, timestamp
  real :: epsr, thetar, xlr, ylr, zlr, boxr
  real :: omegar, lambdar
  real :: axdum, aydum, azdum,phidum, bdum
  integer :: ioffset, i1, i2, npp_partial, npp_total, ipass, me_read, nrest, nadd
  integer :: nslice_e, nslice_i
  integer :: numfiles = 0
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
  if (my_rank == 0) then
     ! input file for PE my_rank
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
  call fsion_paropen_mpi(trim(cinfile),"br",numfiles,MPI_COMM_WORLD,MPI_COMM_WORLD,chunksize,fsblksize,my_rank,cfile_new,sid)
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

  if (my_rank == 0) then
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
    if (my_rank==0) write (*,*) "*** Stopping program"
    call MPI_FINALIZE(ierr)
    call closefiles
    stop
  endif

!  if (ncpu_merge < 0 ) then

!     ! n_cpu > # data sets, so carve up restart data
!     if (ncpu_merge==n_cpu) np_local_total = ne+ni  ! Override total # if dataset to be split amoung all PEs

!     npp = npp_total/ncpu_merge
!     nrest = mod(npp_total,ncpu_merge)  ! Remainder
!
!     if (mod(me,ncpu_merge)==ncpu_merge-1) then
!        nadd = nrest  ! Put remainder on last PE
!     else
!        nadd = 0
!     endif
!
!     if(mod(my_rank,1000)==0) write(*,*) 'PE ',my_rank,': Reading ',npp+nadd,' particles out of ',npp_total,' from ',cinfile
! ! Record # particles read in openfiles.out
!     write(90,*) 'PE ',my_rank,': Reading ',npp+nadd,' particles out of ',npp_total,' from ',cinfile

     !  Skip dummy blocks up to previous PEs
!     nslice_e=0
!     nslice_i = 0
!     nslice = 0
!     if (debug_level>0) write(ipefile,*) 'skip pass ',j

	np_local = npartr

  ! npp * sizepof(real*8)
	size1 = npp*sizeof(realdummy)
  ! npp*sizeof(integer)
	size2 = np_local*sizeof(sid)

!    if(my_rank == 0) then
!      write(6,*)'**************** np_local=',np_local,' **************'
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

!    do i=1,np_local
!       if (beam_config.eq.5 .and. q(i)>0 .and. x(i) < window_min+dt .and. x(i) > window_min) then
! 		! create rezoning slice for wakefield mode - first few blocks should be sufficient
!          nslice = nslice+1
!          xslice(nslice) = x(i)+x_plasma ! Include offset for new slice
!          yslice(nslice) = y(i)
!          zslice(nslice) = z(i)
!       endif
!    end do
!  else
!  	 ! n_cpu < # data sets, so merge together
!     ioffset = 0
!  endif

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  call fsion_parclose_mpi(sid,MPI_COMM_WORLD,ierr)

!  if (my_rank==0) write(*,*) "Dumping particle data for checking"
!  cfile="data/pe"//csubme//"/parts_dump."//cdump(1:6)
!  open (60,file=cfile)
!  write(60,'((12(1pe14.5),2i9))')  &
!       (x(i), y(i), z(i), ux(i), uy(i), uz(i), q(i), m(i), &
!        Ex(i), Ey(i), Ez(i), &  ! electric field
!        pot(i), &  ! potential
!        pepid(i), pelabel(i), i=1,np_local)
!  close(60)

!  call MPI_BARRIER( MPI_COMM_WORLD, ierr)   ! Synchronize first
!  call MPI_BCAST( nslice, 1, MPI_INTEGER, n_cpu-1, MPI_COMM_WORLD,ierr)

  pepid(1:np_local) = my_rank                ! processor ID
  Ex(1:np_local) = 0.       ! zero fields until first force comp.
  Ey(1:np_local) = 0.
  Ez(1:np_local) = 0.
  work(1:np_local) = 1.

  ! Rescale velocities if different temperature required
  if (T_scale /= 1) then
     if (my_rank==0) write(*,*) 'Rescaling temperature by ',T_scale,' to ',Te_keV
     do i=1,np_local
        if (q(i)<0) then
           ux(i) = ux(i)*sqrt(T_scale)
           uy(i) = uy(i)*sqrt(T_scale)
           uz(i) = uz(i)*sqrt(T_scale)
        endif
     end do
  endif

  if (my_rank==0) write(*,*) "Finished reading particle data"

end subroutine predef_parts

#else

!  ---- No SION installed ---------
! ======================
!
!   DUMP
!
!   Gather and write out particle data for restart 
!
!
! ======================

subroutine dump(timestamp)


  use module_physvars
  use module_particle_props

  implicit none   
  include 'mpif.h'

  character(30) :: cfile
  character(6) :: cdump
  integer, intent(in) :: timestamp
  integer :: i
  integer, save :: icall=0
  real :: simtime


  simtime = dt*timestamp


  ! get filename suffix from dump counter
  do i=0,4
     cdump(6-i:6-i) =  achar(mod(timestamp/10**i,10) + 48)  
  end do
  cdump(1:1) = achar(timestamp/10**5 + 48)

!  cfile="data/pe"//csubme//"/parts_info."//cdump(1:6)
  call system("mkdir -p " // "dumps")
  cfile="dumps/info_p"//csubme//"."//cdump(1:6)

  open (60,file=cfile)    
  write(60,'(7(a9,i8/),9(a9,f12.5/),9(a9,1pe12.5/),2(a9,3(1pe12.5)/))')  &    ! info block
       'itime=',timestamp, 'np_local=',np_local, &
       'ne=',ne, 'ni=',ni, 'npbeam=',np_beam, 'geometry=', target_geometry, &
       'scheme=',scheme, &
       'xl=',xl, 'yl=',yl, 'zl=',zl, &
       'eps=', eps, 'theta=',theta,' tlaser= ',tlaser,' trun= ',trun, &
       'omega=',omega,'lambda=',lambda,'  qe=',qe,'  qi=',qi, &
       'mass_e=',mass_e,'mass_i=',mass_i,'Zion=',Zion,'a_ii=',a_ii, &
       'Vplas=',Vplas,'Aplas=',Aplas,'Qplas=',Qplas, &
       'centre=',plasma_centre(1:3),'focus=',focus(1:3)
  

  if (my_rank.eq.0) then

    cfile="parts_info.in"     ! copy to default restart block
    open (61,file=cfile)    
    write(61,'(7(a9,i8/),9(a9,f12.5/),9(a9,1pe12.5/),2(a9,3(1pe12.5)/))')  &    ! info block
       'itime=',timestamp, 'np_local=',np_local, &
       'ne=',ne, 'ni=',ni, 'npbeam=',np_beam, 'geometry=', target_geometry, &
       'scheme=',scheme, &
       'xl=',xl, 'yl=',yl, 'zl=',zl, &
       'eps=', eps, 'theta=',theta,' tlaser= ',tlaser,' trun= ',trun, &
       'omega=',omega,'lambda=',lambda,'  qe=',qe,'  qi=',qi, &
       'mass_e=',mass_e,'mass_i=',mass_i,'Zion=',Zion,'a_ii=',a_ii, &
       'Vplas=',Vplas,'Aplas=',Aplas,'Qplas=',Qplas, &
       'centre=',plasma_centre(1:3),'focus=',focus(1:3)
  
    close (61)
    open (62,file="runstamp")  ! time stamp 
    write(62,'(a)') cdump(1:6)
    close (62)
    write(ipefile,'(2a)') 'Particle dump',cdump(1:6)

  write(6,'(//a/7(a9,i8/),9(a9,f12.5/),9(a9,1pe12.5/),2(a9,3(1pe12.5)/))') 'PARTICLE DUMP:', &    ! info block
       'itime=',timestamp, 'np_local=',np_local, &
       'ne=',ne, 'ni=',ni, 'npbeam=',np_beam, 'geometry=', target_geometry, &
       'scheme=',scheme, &
       'xl=',xl, 'yl=',yl, 'zl=',zl, &
       'eps=', eps, 'theta=',theta,' tlaser= ',tlaser,' trun= ',trun, &
       'omega=',omega,'lambda=',lambda,'  qe=',qe,'  qi=',qi, &
       'mass_e=',mass_e,'mass_i=',mass_i,'Zion=',Zion,'a_ii=',a_ii, &
       'Vplas=',Vplas,'Aplas=',Aplas,'Qplas=',Qplas, &
       'centre=',plasma_centre(1:3),'focus=',focus(1:3)

  endif
  close(60)


  
!  cfile="data/pe"//csubme//"/parts_dump."//cdump(1:6)
  cfile="dumps/parts_p"//csubme//"."//cdump(1:6)
  open (60,file=cfile) 
  write(60,'((12(1pe14.5),2i9))')  &
       (x(i), y(i), z(i), ux(i), uy(i), uz(i), q(i), m(i), &
        Ex(i), Ey(i), Ez(i), &  ! electric field
        pot(i), &  ! potential
        pepid(i), pelabel(i),i=1,np_local)
  close(60)

  icall = icall + 1

end subroutine dump



! ==============================================
!
!                PREDEF_PARTS
!
!  Read particle distribution from file
!
! ==============================================

subroutine predef_parts

  use module_physvars
  use module_particle_props
  implicit none
  include 'mpif.h'

  integer :: i, idummy, ierr

  character(30) :: cinfile, cdump, cfile
  character(5) :: cme
  integer :: ner, nir, np_beamr, npartr, iconf, iens, timestamp
  real :: epsr, thetar, xlr, ylr, zlr, boxr
  real :: omegar, lambdar
  real :: axdum, aydum, azdum,phidum
  integer :: ioffset, i1, i2, npp_partial, npp_total, ipass, me_read, j, nrest, nadd
  integer :: nslice_e, nslice_i
  integer :: lendir  ! length of directory string
  logical :: stopflag=.false.

  if (my_rank == 0) then 


     ! input file for PE my_rank

     cinfile="parts_info.in"
     open(80,file=cinfile)

     ! Root reads info block in run directory to check that run parameters correct

     read(80,'(7(9x,i8/),10(9x,f12.5/),9(9x,1pe12.5/),2(9x,3f12.5/))')  &    ! info block - skip variable names
          itime_start, npartr, &
	  ner, nir, np_beamr, iconf, iens, &
	  xlr, ylr, zlr, boxr, &
  	  epsr, thetar, &
          tlaser, trun, omegar, lambdar, &
          qe, qi, mass_e, mass_i, Zion, a_ii, Vplas, Aplas, Qplas, &
          plasma_centre(1:3), focus(1:3)
     close(80)
     write(6,'(/a/a/a,i5)') 'RESTART:','Reading run data from info block: parts_info.in','Timestep:',itime_start

     if (my_rank==0) write(*,'(7(a12,i12/),9(a12,f12.5/),9(a12,1pe12.5/),2(a12,3(1pe12.5)/))')  &    ! info block - skip variable names
          'Start time: ', itime_start, & 
          '# particles: ', npartr, &
	  '# electrons: ', ner, &
          '# ions: ', nir, & 
          '# beam: ', np_beamr, &
          'Geom  : ', iconf, & 
          'Scheme: ', iens, &
	  'Box_x: ',xlr,&
          'Box_y: ', ylr, &
          'Box_z: ', zlr, &
          'eps: ', eps,  &
  	  'theta: ', thetar, &
          'tlaser: ',tlaser, &
          'trun: ',trun, &
          'omega:',omegar, &
          'lambda:', lambdar, &
          'qe:',qe, &
          'qi:',qi, &
          'mass_e:',mass_e, &
          'mass_i:',mass_i, &
          'Zion:',Zion, &
          'a_ii:', a_ii, &
          'Vplas: ',Vplas, &
          'Aplas: ',Aplas, &
          'Qplas: ',Qplas, &
          'centre: ',plasma_centre(1:3), &
          'focus: ',focus(1:3)

     if (ner /= ne .or. nir /= ni) then
        write(*,*) '*** Particle nos. in input deck do not match those in restart file parts_info.in'
        stopflag=.true.
     endif

     if (epsr /= eps) then
        write(*,*) '*** Warning: potential cutoff eps changed - check inputs'
     endif

     if (thetar /= theta) then
        write(*,*) '*** Warning: MAC changed - check inputs'
     endif

     if (iconf /= target_geometry) then
        write(*,*) '*** Warning: Target geometry in restart file ',iconf, &
             ' does not match value ',target_geometry,' in run.h - check inputs'
     endif
  endif



  call MPI_BARRIER( MPI_COMM_WORLD, ierr)   ! Synchronize first

  if (stopflag) then
    if (my_rank==0) write (*,*) "*** Stopping program"
    call MPI_FINALIZE(ierr)
    call closefiles
    stop
  endif

  call MPI_BCAST( itime_start, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( xl, 1, MPI_REAL, 0, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( yl, 1, MPI_REAL, 0, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( zl, 1, MPI_REAL, 0, MPI_COMM_WORLD,ierr)

  call MPI_BCAST( theta, 1, MPI_REAL, 0, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( eps, 1, MPI_REAL, 0, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( trun, 1, MPI_REAL, 0, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( tlaser, 1, MPI_REAL, 0, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( qe, 1, MPI_REAL, 0, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( qi, 1, MPI_REAL, 0, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( mass_e, 1, MPI_REAL, 0, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( mass_i, 1, MPI_REAL, 0, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( Zion, 1, MPI_REAL, 0, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( a_ii, 1, MPI_REAL, 0, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( Vplas, 1, MPI_REAL, 0, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( Aplas, 1, MPI_REAL, 0, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( Qplas, 1, MPI_REAL, 0, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( plasma_centre, 3, MPI_REAL, 0, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( focus, 3, MPI_REAL, 0, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( ne, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( ni, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( npart_total, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)

  if (ncpu_merge < 0 ) then

     ! n_cpu > # data sets, so carve up restart data

     if (ncpu_merge == -1) then
 ! automatic splitting of single merged data set in dumps/ directory
        ncpu_merge = n_cpu 
        cme = "dumps"
        lendir=5  
 
     else
        ncpu_merge = -ncpu_merge  ! split several sets
     ! Filename (directory) prefix
       me_read = my_rank/ncpu_merge
       cme = "p" &   
       // achar(me_read/1000+48) &
       // achar(mod(me_read/100,10)+48) &
       // achar(mod(me_read/10,10)+48) &
       // achar(mod(me_read,10)+48)  ! Convert 4-digit PE number into character string
       lendir=11
     endif

     if (my_rank==0) write (*,*) 'Splitting sets by 1:',ncpu_merge

     ! get filename suffix from dump counter
     do i=0,4
        cdump(6-i:6-i) =  achar(mod(itime_start/10**i,10) + 48)  
     end do
     cdump(1:1) = achar(itime_start/10**5 + 48)

 !    cfile=cme(1:lendir)//"/parts_info."//cdump(1:6)
     cfile="dumps/info_"//cme//"."//cdump(1:6)

     if (my_rank==0) write (*,'(a,a)') 'Reading info file',cfile
     open (60,file=cfile)    
     read(60,'(2(9x,i8/))')  timestamp,npp_total  ! Find # particles to be read 
     close(60)
     if (ncpu_merge==n_cpu) npp_total = ne+ni  ! Override total # if dataset to be split amoung all PEs
     np_local = npp_total/ncpu_merge
     nrest = mod(npp_total,ncpu_merge)  ! Remainder

     if (mod(my_rank,ncpu_merge)==ncpu_merge-1) then
        nadd = nrest  ! Put remainder on last PE
     else
        nadd = 0
     endif

     if(mod(my_rank,1000)==0) write(*,*) 'PE ',my_rank,': Reading ',np_local+nadd,' particles out of ',npp_total,' from ',cme(1:lendir)
! Record # particles read in openfiles.out
     write(90,*) 'PE ',my_rank,': Reading ',np_local+nadd,' particles out of ',npp_total,' from ',cme(1:lendir)

 !    cfile=cme(1:lendir)//"/parts_dump."//cdump(1:6)
   cfile="dumps/parts_"//cme//"."//cdump(1:6)
    open (60,file=cfile) 

     !  Skip dummy blocks up to previous PEs 
     nslice_e=0
     nslice_i = 0
     nslice = 0
     do j=0,mod(my_rank,ncpu_merge)-1
        if (debug_level>0) write(ipefile,*) 'skip pass ',j
        do i=1,np_local
           read(60,*) x(i),y(i),z(i),ux(i),uy(i),uz(i),q(i),m(i), &
                axdum,aydum,azdum,phidum, idummy,pelabel(i)
           if (beam_config.eq.5 .and. q(i)>0 .and. x(i) < window_min+dt .and. x(i) > window_min) then
  ! create rezoning slice for wakefield mode - first few blocks should be sufficient
              nslice = nslice+1
              xslice(nslice) = x(i)+x_plasma ! Include offset for new slice
              yslice(nslice) = y(i)
              zslice(nslice) = z(i)
           endif
        end do
     end do

     !  Now read in particles to keep
     read(60,*) (x(i),y(i),z(i),ux(i),uy(i),uz(i),q(i),m(i),axdum,aydum,azdum,phidum, &
          idummy,pelabel(i),i=1,np_local+nadd)
     close (60)

!     if (my_rank==n_cpu-1 .and. debug>0)  write(ipefile,'(a,i8/3(f15.5))') 'Slice  ions:', nslice, &
!          (xslice(i),yslice(i),zslice(i),i=1,nslice)
     np_local = np_local+nadd

  else
     ! n_cpu < # data sets, so merge together

     ioffset = 0
     np_local = 0
     if (ncpu_merge>1) write (*,*) 'Merging sets by ',ncpu_merge,' : 1'

     do ipass = 0, ncpu_merge-1
        me_read = my_rank*ncpu_merge + ipass
        ! Filename (directory) prefix
        cme = "p" &   
       // achar(me_read/1000+48) &
       // achar(mod(me_read/100,10)+48) &
       // achar(mod(me_read/10,10)+48) &
       // achar(mod(me_read,10)+48)  ! Convert 4-digit PE number into character string

        write(90,'(a,a,a6,i8)') 'Reading from ',cme,'itime',itime_start
        ! get filename suffix from dump counter
        do i=0,4
           cdump(6-i:6-i) =  achar(mod(itime_start/10**i,10) + 48)  
        end do
        cdump(1:1) = achar(itime_start/10**5 + 48)

!        cfile=cme//"/parts_info."//cdump(1:6)
        cfile="dumps/info_"//cme//"."//cdump(1:6)

        open (60,file=cfile)    
        read(60,'(2(9x,i8/))')  timestamp,npp_partial  ! Find # particles to be read 
        close(60)


!        cfile=cme//"/parts_dump."//cdump(1:6)
       cfile="dumps/parts_"//cme//"."//cdump(1:6)

        open (60,file=cfile) 

        i1 = ioffset + 1
        i2 = ioffset + npp_partial

        !  Initialise particles: read from file
        read(60,*) (x(i),y(i),z(i),ux(i),uy(i),uz(i),q(i),m(i),axdum,aydum,azdum,phidum, &
             idummy,pelabel(i),i=i1,i2)

        close (60)

        ioffset = ioffset + npp_partial
        np_local = np_local + npp_partial
     enddo

  endif
  call MPI_BARRIER( MPI_COMM_WORLD, ierr)   ! Synchronize first
  call MPI_BCAST( nslice, 1, MPI_INTEGER, n_cpu-1, MPI_COMM_WORLD,ierr)

  ! add displacement vector
  !  x(1:np_local) = x(1:np_local) + displace(1)
  !  y(1:np_local) = y(1:np_local) + displace(2)
  !  z(1:np_local) = z(1:np_local) + displace(3)

  pepid(1:np_local) = my_rank                ! processor ID

  Ex(1:np_local) = 0.       ! zero fields until first force comp.
  Ey(1:np_local) = 0.
  Ez(1:np_local) = 0.
  work(1:np_local) = 1.

  ! Rescale velocities if different temperature required 
  if (T_scale /= 1) then
     if (my_rank==0) write(*,*) 'Rescaling temperature by ',T_scale,' to ',Te_keV
     do i=1,np_local
        if (q(i)<0) then
           ux(i) = ux(i)*sqrt(T_scale)
           uy(i) = uy(i)*sqrt(T_scale)
           uz(i) = uz(i)*sqrt(T_scale)
        endif
     end do
  endif

if (my_rank==0) write(*,*) "Finished reading in particle data" 

end subroutine predef_parts


#endif

! ======================
!
!   DUMP_FIELDS
!
!   Write out field data for 
!     postprocessing
!
!
! ======================

subroutine dump_fields(timestamp)


  use module_physvars
  use module_laser

  implicit none 
  include 'mpif.h'

  real, dimension(0:max(nxh+1,ngx+1)) :: phi_pond, ex_pond, ey_pond, ez_pond, azr, epond
  real, dimension(0:ngx+1) :: rhoe_slice, ex_slice, jxe_slice, rhoi_slice
  real, dimension(0:ngx+1,0:ngy+1,0:ngz+1) :: exg, eyg, ezg, jxeg, jyeg, jzeg
  complex :: aph(nxh)
  real :: dx, dz, dy, xd, yd, zd, simtime, epon_x, epon_y, epon_z, phipond
  real :: uxd, tplot, pha, gamma
  real :: Qtot, Qbox, norm, rhonorm, tpon, bx_em, by_em, az_em,ez_em

  character(30) :: cfile

  character(6) :: cdump
  integer, intent(in) :: timestamp
  integer :: i, j, k
  integer :: icall, ierr
  integer :: jfoc, kfoc, ng, nave
  complex :: yi=(0.,1.)

  icall = timestamp/ivis
  simtime = timestamp*dt

  ! Gather locally accumulated averages
  ng = (ngx+2)*(ngy+2)*(ngz+2)                         ! total # gridpoints

  call MPI_ALLREDUCE(rhoe_loc, rhoe, ng, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(rhoi_loc, rhoi, ng, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(jxe_loc, jxeg, ng, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(jye_loc, jyeg, ng, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(jze_loc, jzeg, ng, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(ex_loc, exg, ng, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(ey_loc, eyg, ng, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(ez_loc, ezg, ng, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)

  ! get filename suffix from dump counter
  do i=0,4
     cdump(6-i:6-i) =  achar(mod(timestamp/10**i,10) + 48)  
  end do
  cdump(1:1) = achar(timestamp/10**5 + 48)

  if (my_rank==0) then
     write(*,'(//a/a,f10.1,a1,f10.1,a1,f10.1)') 'Cycle-average field dump:', &
          'writing out densities, fields on grid ',xl,'*',yl,'*',zl
     ! dump electron and ion densities, electron currents, DC fields within xl*yl*zl
     cfile = "fields/"//cdump//".xyz"
     open (62,file=cfile)
     cfile = "fields/laser_"//cdump//".xyz"
     open (63,file=cfile)
     dx = xl/ngx
     dy = yl/ngy
     dz = zl/ngz
     Qbox = 0.
     do k=1,ngz
        do j=1,ngy
           do i=1,ngx
              Qbox = Qbox + rhoe(i,j,k)*dx*dy*dz
              write(62,'(8e13.3)') rhoe(i,j,k)/omega**2,rhoi(i,j,k)/omega**2, &
                   jxeg(i,j,k), jyeg(i,j,k), jzeg(i,j,k), &
                   exg(i,j,k),  eyg(i,j,k), ezg(i,j,k)
!              xd = (i-0.5)*dx - 20.-focus(1) ! position relative to laser focus
              xd = (i-0.5)*dx - focus(1) ! position relative to laser focus
              yd = (j-0.5)*dy - focus(2)
              zd = (k-0.5)*dz - focus(3)

              laser_type: select case(beam_config)

              case(4)
                 call fpond( 1.57/omega, tpulse,sigma,vosc,omega, rho_upper, &
                      xd,yd,zd,epon_x,epon_y,epon_z,phipond)


                 Tpon = min(1.,tlaser/tpulse) * (sin(omega*tlaser))**2
                 !                 mvis(lcount) = epon_x/omega ! Pond field, EM norm

              case(5)  ! propagating fpond
                 call laser_bullet( tpulse,sigma,vosc,xd,yd,zd,epon_x,epon_y,epon_z,phipond)

              case(24) ! Standing wave fpond Ez, By, Az
                 call emobliq(tlaser,tpulse,sigma,vosc,omega,theta_inc, rho_upper,&
                      xd,yd,zd,epon_x,epon_y,epon_z,phipond, ez_em, bx_em, by_em)

              case(6) ! Plane wave
                 call emplane(tlaser,tpulse,sigma,vosc,omega,xd,yd,zd,ez_em,by_em,bx_em,az_em,phipond)

              case default ! Propagating fpond
                 phipond = 0  

              end select laser_type
              write(63,'(4e13.3)') phipond,epon_x,epon_y,epon_z
           end do
        end do
     end do
     Qtot = SUM(rhoe)*dx*dy*dz  ! including ghost cells
!     write(ifile_cpu,'(4(a,f14.5/))') &
!          'Total charge on grid:',Qbox, &
!          '         ghost cells:',Qtot-Qbox, &
!          '                 sum:',Qtot, &
!          'Initial charge Q_s*Ne = rho0*Vplas = ',Vplas*rho0
     write(*,'(4(a,f14.5/))') &
          'Total charge on grid:',Qbox, &
          '         ghost cells:',Qtot-Qbox, &
          '                 sum:',Qtot, &
          'Initial charge Q_s*Ne = rho0*Vplas = ',Vplas*rho0

     close(62)
     close(63)


     ! x-slices along laser axis
     jfoc = focus(2)/dy
     kfoc = focus(3)/dz

     ! density average line-out along laser axis: nave*nave average, converted to n/nc

     rhoe_slice = 0.
     rhoi_slice = 0.
     ex_slice = 0.
     jxe_slice = 0.

     if (ngz<=5) then
        nave=0
     else
        nave = min(3,ngz/2)
     endif

     norm = (2*nave+1)**2
     rhonorm = norm*omega**2

     do k=kfoc-nave,kfoc+nave
        do j=jfoc-nave,jfoc+nave

           rhoe_slice(1:ngx) = rhoe_slice(1:ngx)+rhoe(1:ngx,j,k)/norm  ! density slice along laser axis: 5x5 average
           rhoi_slice(1:ngx) = rhoi_slice(1:ngx)+rhoi(1:ngx,j,k)/norm
           ex_slice(1:ngx) = ex_slice(1:ngx)+exg(1:ngx,j,k)/norm
           jxe_slice(1:ngx) = jxe_slice(1:ngx)+jxeg(1:ngx,j,k)/norm

        end do
     end do

     ! Laser fields on axis
     do i=1,ngx
        xd=i*dx-focus(1)
        !        yd=sigma/2.
        !        zd=sigma/2.
        yd = 0.
        zd = 0.
        call fpond( 1.57/omega,1.0,sigma,vosc,omega,rho_upper,xd,yd,zd,ex_pond(i),ey_pond(i), &
             ez_pond(i), phi_pond(i))
     end do
     ! Renormalise to EM units
     cfile = "fields/xslice."//cdump
     open (62,file=cfile)
     write(62,'(a)') '!   x      rho_e   rho_i  ex, jxe,  phi_p,  ex_p,  ey_p,   ez_p  '
     write(62,'((9(1pe12.4)))') &
          (i*dx,rhoe_slice(i)/omega**2, rhoi_slice(i)/omega**2, ex_slice(i)/omega,&
          jxe_slice(i), phi_pond(i),ex_pond(i)/omega, ey_pond(i)/omega, ez_pond(i)/omega,i=1,ngx)
     close(62)

! Write out time-averaged E-field profiles to file
     write(*,'(a40,a10)') 'Cycle-averaged lineouts at',cdump
     dx = (xgav_end - xgav_start)/ngav ! spacing for time-ave grid
     cfile = "fields/ex_ave."//cdump
     open (60,file=cfile)
     write(60,'(2(a12))') '!   x      ',' ex  '
     write(60,'((2(1pe12.4)))') &
          (i*dx+xgav_start, ex_ave(i), i=0,ngav)
     close(60)
     dy = yl/ngav ! spacing for radial grid
     cfile = "fields/ey_ave."//cdump
     open (60,file=cfile)
     write(60,'(2(a12))') '!   r      ',' er  '
     write(60,'((2(1pe12.4)))') &
          ((i-1)*dy, ey_ave(i,1), i=1,ngav)
     close(60)


! Laser fields on Helmholtz grid
  tplot=tlaser
  pha=omega*tplot
  do i=1,nxh
     aph(i) = Az_helm(i)*cexp(yi*pha)
     Azr(i) = Real(aph(i))
  end do
!  Azr(1:nxh) = Real(Az_helm(1:nxh)*cexp(yi*pha))

  ! pond force - without gamma factor, as in emfield
  do i=1,nxh
     gamma = sqrt(1. + azr(i)**2)

     epond(i) = .5/dxh*azr(i)*( azr(i+1)-azr(i-1) )/gamma
  end do
!     write(*,*) 'On HH grid: '
!     write(*,*) 'start, end, dxh, phase, nxh:',xh_start, xh_end, dxh, pha, nxh
!     write(*,*) 'x, a, a^iph, |a|, real(a), epond'
!     write(*,'((8(f12.4)))') &
!          (i*dxh+xh_start, az_helm(i),aph(i),abs(az_helm(i)),azr(i), epond(i), i=1,nxh)

     do i=0,nxh
        xd=i*dxh+xh_start
        !        yd=sigma/2.
        !        zd=sigma/2.
        yd = 0.
        zd = 0.
        uxd = 0.
        call fpond_helm( tplot, sigma, omega, &
               xd,yd,zd,uxd,Az_helm,nxh,xh_start, xh_end, dxh, focus(1), &
	       ex_pond(i),ey_pond(i),ez_pond(i),phi_pond(i))


     end do
! Write out Helmholtz fields to file
     write(*,'(a40,a10)') 'Helmholtz lineouts at',cdump
     cfile = "fields/helmholtz."//cdump
     open (60,file=cfile)
     write(60,'(2(a12))') '!   x_helm  ',' rho     Azr    Ex    az^2  '
!     dxh = (xh_end-xh_start)/nxh
     write(60,'((5(1pe12.3)))') &
!          (i*dxh+xh_start-xl/2.+x_plasma/2., rho_helm(i), abs(az_helm(i)**2), i=0,nxh)
          (i*dxh+xh_start, rho_helm(i),Azr(i),ex_pond(i),abs(az_helm(i)), i=0,nxh)
     close(60)
  endif

  icall = icall + 1

  ! Rezero local fields
  rhoe_loc = 0.

  ex_ave = 0.
  ex_loc = 0.
  ey_loc = 0.
  ez_loc = 0.
  jxe_loc = 0.
  jye_loc = 0.
  jze_loc = 0.


end subroutine dump_fields


! ======================
!
!   SLICES
!
!   Write out field data for 
!   1D postprocessing

!
!
! ======================

subroutine slices(timestamp)

  use module_physvars
  use module_laser
  implicit none   

  real, dimension(ngx) :: work1, work2
  real, dimension(ngx) :: phi_pond, ex_pond, ey_pond, ez_pond
  real :: dx, dz, dy, xd, yd, zd, simtime
  real :: Qtot, Qbox
  character(30) :: cfile
  character(6) :: cdump
  integer, intent(in) :: timestamp
  integer :: i, j, k
  integer :: icall
  integer :: jfoc, kfoc

  icall = timestamp/ivis
  simtime = timestamp*dt

  ! get filename suffix from dump counter
  do i=0,4
     cdump(6-i:6-i) =  achar(mod(timestamp/10**i,10) + 48)  
  end do
  cdump(1:1) = achar(timestamp/10**5 + 48)

  if (my_rank==0) then
     !  Fields: electron and ion densities within xl*yl*zl
     cfile = "wf_field."//cdump
     open (62,file=cfile)

     dx = xl/ngx
     dy = yl/ngy
     dz = zl/ngz
     Qbox = 0.
     do k=1,ngz
        do j=1,ngy
           do i=1,ngx
              Qbox = Qbox + rhoe(i,j,k)*dx*dy*dz
              write(62,'(3f13.5,2e13.3)') i*dx,j*dy,k*dz,abs(rhoe(i,j,k)),rhoi(i,j,k)
           end do
        end do
     end do
     Qtot = SUM(rhoe)*dx*dy*dz  ! including ghost cells
     if (debug_level>1) write(ipefile,'(4(a,f14.5/))') &
          'Total charge on grid:',Qbox, &
          '         ghost cells:',Qtot-Qbox, &
          '                 sum:',Qtot, &
          'Initial charge Q_s*Ne = rho0*V = ',Vplas*rho0

     close(62)


     cfile = "field_slice."//cdump
     open (62,file=cfile)

     jfoc = focus(2)/dy
     kfoc = focus(3)/dz

     ! density average line-out along laser axis: 5x5 average, converted to n/nc

     work1 = 0.
     work2 = 0.
     do k=kfoc-2,kfoc+2
        do j=jfoc-2,jfoc+2
           work1(1:ngx) = work1(1:ngx)+rhoi(1:ngx,j,k)/25./omega**2  ! density slice along laser axis: 5x5 average
           work2(1:ngx) = work2(1:ngx)+rhoe(1:ngx,j,k)/25./omega**2
        end do
     end do

     do i=1,ngx
        xd=i*dx-focus(1)
        yd=sigma/2.
        zd=sigma/2.
        laser_model: select case(beam_config)

        case(4)  ! standing wave fpond
           call fpond(1.57/omega,1.0,sigma,vosc,omega,rho_upper,xd,yd,zd,ex_pond(i),ey_pond(i), ez_pond(i), phi_pond(i))
        case(5)  ! propagating fpond
           call laser_bullet( tpulse,sigma,vosc, & 
                xd,yd,zd,ex_pond(i),ey_pond(i),ez_pond(i),phi_pond(i))
        case default
           ex_pond(i)=0
           ey_pond(i)=0
           ez_pond(i)=0
           phi_pond(i)=0
        end select laser_model

     end do

     write(62,'(7f13.5)') (i*dx,work1(i),work2(i), &
          phi_pond(i),ex_pond(i), ey_pond(i), ez_pond(i),i=1,ngx)
     close(62)

  endif
  icall = icall + 1

end subroutine slices

end module module_io
