! ==============================================
!
!                PREDEF_PARTS
!
!  Read particle distribution from file
!
! ==============================================

subroutine predef_parts

  use physvars
  use treevars
  implicit none
  include 'mpif.h'

  integer :: i, ipe, idummy, ierr

  character(30) :: cinfile, cdump, cfile
  character(9) :: ct
  character(5) :: cme
  integer :: ner, nir, np_beamr, npartr, iconf, iens, timestamp
  real :: epsr, thetar, xlr, ylr, zlr, boxr
  real :: axdum, aydum, azdum,phidum
  integer :: ioffset, i1, i2, npp_partial, npp_total, ipass, me_read, j, nrest, nadd
  integer :: nslice_e, nslice_i

  if (me == 0) then 


     ! input file for PE me

     cinfile="parts_info.in"
     open(80,file=cinfile)

     ! Root reads info block in run directory to check that run parameters correct

     read(80,'(7(9x,i8/),8(9x,f12.5/))')  &    ! info block - skip variable names
          itime_start, npartr, &
	  ner, nir, np_beamr, iconf, iens, &
	  xlr, ylr, zlr, boxr, &
  	  epsr, thetar, &
          tlaser, trun   
     close(80)
     write(6,'(/a/a/a,i5)') 'RESTART:','Reading run data from info block: parts_info.in','Timestep:',itime_start

     if (me==0) write(*,'(7(a12,i12/),7(a12,f12.5/))')  &    ! info block - skip variable names
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
          'trun: ',trun  

     if (ner /= ne .or. nir /= ni) then
        write(*,*) '*** Particle nos. in input deck do not match those in restart file parts_info.in'
        stop
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

  call MPI_BCAST( xl, 1, MPI_REAL8, 0, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( yl, 1, MPI_REAL8, 0, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( zl, 1, MPI_REAL8, 0, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( boxsize, 1, MPI_REAL8, 0, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( theta, 1, MPI_REAL8, 0, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( eps, 1, MPI_REAL8, 0, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( trun, 1, MPI_REAL8, 0, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( tlaser, 1, MPI_REAL8, 0, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( ne, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( ni, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( npart, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( itime_start, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)

  if (nmerge < 0 ) then

     ! num_pe > # data sets, so carve up restart data

     if (nmerge == -1) then
        nmerge = num_pe  ! automatic splitting of single data set
     else
        nmerge = -nmerge  ! split several sets
     endif

     if (me==0) write (*,*) 'Splitting sets by 1:',nmerge

     ! Filename (directory) prefix
     me_read = me/nmerge
     cme = "pe"// achar(me_read/100+48) // achar(mod(me_read/10,10)+48) // achar(mod(me_read,10)+48)  

     ! get filename suffix from dump counter
     do i=0,4
        cdump(6-i:6-i) =  achar(mod(itime_start/10**i,10) + 48)  
     end do
     cdump(1:1) = achar(itime_start/10**5 + 48)

     cfile=cme//"/parts_info."//cdump(1:6)


     open (60,file=cfile)    
     read(60,'(2(9x,i8/))')  timestamp,npp_total  ! Find # particles to be read 
     close(60)

     npp = npp_total/nmerge
     nrest = mod(npp_total,nmerge)  ! Remainder

     if (mod(me,nmerge)==nmerge-1) then
        nadd = nrest  ! Put remainder on last PE
     else
        nadd = 0
     endif

     write(ipefile,*) 'PE ',me,': Reading ',npp+nadd,' particles out of ',npp_total,' from ',cme
     cfile=cme//"/parts_dump."//cdump(1:6)
     open (60,file=cfile) 



     !  Skip dummy blocks up to previous PEs 
     nslice_e=0
     nslice_i = 0
     nslice = 0
     do j=0,mod(me,nmerge)-1
        write(ipefile,*) 'skip pass ',j
        do i=1,npp
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
          idummy,pelabel(i),i=1,npp+nadd)
     close (60)

!     if (me==num_pe-1)  write(ipefile,'(a,i8/3(f15.5))') 'Slice  ions:', nslice, &
!          (xslice(i),yslice(i),zslice(i),i=1,nslice)
     npp = npp+nadd

  else
     ! num_pe < # data sets, so merge together

     ioffset = 0
     npp = 0
     if (nmerge>1) write (*,*) 'Merging sets by ',nmerge,' : 1'

     do ipass = 0, nmerge-1
        me_read = me*nmerge + ipass
        ! Filename (directory) prefix
        cme = "pe"// achar(me_read/100+48) // achar(mod(me_read/10,10)+48) // achar(mod(me_read,10)+48)  
        write(*,'(a,a)') 'Reading from ',cme
        ! get filename suffix from dump counter
        do i=0,4
           cdump(6-i:6-i) =  achar(mod(itime_start/10**i,10) + 48)  
        end do
        cdump(1:1) = achar(itime_start/10**5 + 48)

        cfile=cme//"/parts_info."//cdump(1:6)


        open (60,file=cfile)    
        read(60,'(2(9x,i8/))')  timestamp,npp_partial  ! Find # particles to be read 
        close(60)


        cfile=cme//"/parts_dump."//cdump(1:6)
        open (60,file=cfile) 

        i1 = ioffset + 1
        i2 = ioffset + npp_partial

        !  Initialise particles: read from file
        read(60,*) (x(i),y(i),z(i),ux(i),uy(i),uz(i),q(i),m(i),axdum,aydum,azdum,phidum, &
             idummy,pelabel(i),i=i1,i2)

        close (60)

        ioffset = ioffset + npp_partial
        npp = npp + npp_partial
     enddo

  endif
  call MPI_BCAST( nslice, 1, MPI_INTEGER, num_pe-1, MPI_COMM_WORLD,ierr)

  ! add displacement vector
  !  x(1:npp) = x(1:npp) + displace(1)
  !  y(1:npp) = y(1:npp) + displace(2)
  !  z(1:npp) = z(1:npp) + displace(3)

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



end subroutine predef_parts



