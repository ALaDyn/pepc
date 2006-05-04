! ==============================================
!
!                PREDEF_PARTS
!
!  Read particle distribution from file
!
! ==============================================

subroutine predef_parts

  use treevars
  use physvars
  use utils

  implicit none
  integer (kind=4) :: i, ipe, idummy

  character(30) :: cinfile, cdump, cfile
  character(9) :: ct
  character(5) :: cme
  integer :: ner, nir, np_beamr, npartr, iconf, iens, timestamp
  real :: epsr, thetar, xlr, ylr, zlr, boxr
  real :: axdum, aydum, azdum,phidum, dum(6)
  integer :: ioffset, i1, i2, npp_partial, npp_total, ipass, me_read, j, nrest, nadd
  integer :: nslice_e, nslice_i

  if (me == 0) then 


     ! input file for PE me

     cinfile="parts_info.in"
     open(80,file=cinfile)

     ! Root reads info block in run directory to check that run parameters correct

     ! Read in info block, skipping variable names

     read(80,*)  itime_start 
     read(80,*)  npartr 
     read(80,*)  ner 
     read(80,*)  nir
     read(80,*)  iens
     read(80,*)  xlr
     read(80,*)  ylr
     read(80,*)  zlr
     read(80,*)  boxsize
     read(80,*)  epsr  
     read(80,*)  thetar  
     write(*,*) "itime ",itime_start,"iens ",iens     
     write(*,*) "npart ",npartr,"  ndust ", ner,"  nstar ",nir
     write(*,'(4(a9,f12.2))') "xlr ",xlr,"  ylr ",ylr,"  zlr ",zlr,"  boxsize ",boxsize
     write(*,'(2(a9,f12.2))') "eps ",epsr,"  theta",thetar
     close(80)



     write(6,'(/a/a/a,i5)') 'RESTART:','Reading run data from info block: parts_info.in','Timestep:',itime_start

     trun = itime_start*dt

     if (ner /= ne .or. nir /= ni) then
        write(*,*) '*** Warning: particle nos. in input deck do not match those in restart file parts_info.in'
        write(*,*) '*** - stopping program -check inputs..'
        stop
     endif

     if (epsr /= eps) then
        write(*,*) '*** Warning: potential cutoff eps changed - check inputs'
     endif

     if (thetar /= theta) then
        write(*,*) '*** Warning: MAC changed - check inputs'
     endif


  endif



  call MPI_BARRIER( MPI_COMM_WORLD, ierr)   ! Synchronize first

  call MPI_BCAST( xl, one, MPI_REAL8, root, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( yl, one, MPI_REAL8, root, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( zl, one, MPI_REAL8, root, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( boxsize, one, MPI_REAL8, root, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( theta, one, MPI_REAL8, root, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( eps, one, MPI_REAL8, root, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( ne, one, MPI_INTEGER, root, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( ni, one, MPI_INTEGER, root, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( npart, one, MPI_INTEGER, root, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( itime_start, one, MPI_INTEGER, root, MPI_COMM_WORLD,ierr)

  if (nmerge < 0 ) then
     ! num_pe > # data sets, so carve up restart data
     if (nmerge == -1) then
        nmerge = num_pe  ! automatic splitting
     else
        nmerge = -nmerge
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
           read(60,*) x(i),y(i),z(i),ux(i),uy(i),uz(i),q(i),axdum, &
                dum(1), dum(2), dum(3), dum(4), dum(5), dum(6)  ! include fdisc dummies (14 reals)
        end do
     end do

     !  Now read in particles to keep
     read(60,*) (x(i),y(i),z(i),ux(i),uy(i),uz(i),q(i),axdum, &
          dum(1),dum(2),dum(3),dum(4),dum(5),dum(6),i=1,npp+nadd)
     close (60)

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
        read(60,*) (x(i),y(i),z(i),ux(i),uy(i),uz(i),q(i),axdum,i=i1,i2)

        close (60)

        ioffset = ioffset + npp_partial
        npp = npp + npp_partial
     enddo

  endif




  pepid(1:npp) = me                ! processor ID
  pelabel(1:npp) = me*npp + (/ (i,i=1,npp) /)     ! Dust labels


  ax(1:npp) = 0.       ! zero accelerations until first force comp.
  ay(1:npp) = 0.
  az(1:npp) = 0.
  work(1:npp) = 1.


!   STARS

  if (me.eq.num_pe-1) then
     npp = npp + ni   ! add stars to last CPU

     ! input file for stars

     cinfile="star_info.in"
     open(81,file=cinfile)

     ! Root reads info block in run directory to check that run parameters correct
     do i=1,ni
        read(81,*) m(ner+i)
     end do
     do i=1,ni
        read(81,*)x(ner+i),y(ner+i),z(ner+i)
     end do
     do i=1,ni
        read(81,*)ux(ner+i),uy(ner+i),uz(ner+i)
     end do
     write(*,'(a14,f12.5,a7,3f12.5)')"star 1 mass: ",m(ner+1),"  pos ",x(ner+1),y(ner+1),z(ner+1)
     write(*,'(a14,f12.5,a7,3f12.5)')"star 2 mass: ",m(ner+2),"  pos ",x(ner+2),y(ner+2),z(ner+2)
     close(81)
     do i=1,ni
        work(ner+i)=0.0
        pelabel(ner+i)=ne+i     
        pepid(ner+i) = me                ! processor ID
        ax(ner+i) = 0.       ! zero accelerations until first force comp.
        ay(ner+i) = 0.
        az(ner+i) = 0.
     end do
  endif

end subroutine predef_parts






