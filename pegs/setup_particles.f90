! ==============================================
!
!                SETUP_PARTICLES
!
!  Read in particle data from configuration file
!  or particle dump
!
! ==============================================

subroutine setup_particles

  use treevars
  use physvars

  implicit none
  integer :: i, ipe, idummy=0, ifile
  real :: t_walk, t_force
  real :: axdum, aydum, azdum,phidum, dum(6)
  integer :: ioffset, i1, i2, npp_partial, npp_total, ipass, me_read, j, nrest, nadd
  integer :: timestamp
  character(30) :: cinfile, cdump, cfile
  character(9) :: ct
  character(11) :: cme



  if (nmerge < 0 ) then
     ! # CPUs is bigger than # data sets, so carve up restart data
     if (nmerge == -1) then
        nmerge = num_pe  ! automatic splitting
     else
        nmerge = -nmerge ! user-defined
     endif

     if (me==0) write (*,*) 'Splitting sets by 1:',nmerge

     ! Filename (directory) prefix
     me_read = me/nmerge
     cme = "data/pe" &   
       // achar(me_read/1000+48) &
       // achar(mod(me_read/100,10)+48) &
       // achar(mod(me_read/10,10)+48) &
       // achar(mod(me_read,10)+48)  ! Convert 4-digit PE number into character string

     ! get filename suffix from dump counter
     do i=0,4
        cdump(6-i:6-i) =  achar(mod(itime_start/10**i,10) + 48)  
     end do
     cdump(1:1) = achar(itime_start/10**5 + 48)

     cfile=cme//"/parts_info."//cdump(1:6)


     open (60,file=cfile)    
     read(60,*)  timestamp
     read(60,*)  npp_total  ! Find # particles to be read from this set
     close(60)

     npp = npp_total/nmerge   ! # particles to be read locally
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

     do j=0,mod(me,nmerge)-1
        write(ipefile,*) 'skip pass ',j
        do i=1,npp
           read(60,*) x(i),y(i),z(i),ux(i),uy(i),uz(i),q(i) & 
           ,axdum, dum(1), dum(2), dum(3), dum(4), dum(5), dum(6) 
 ! include fdisc dummies (14 reals)

        end do
     end do

     !  Now read in particles to keep
     read(60,*) (x(i),y(i),z(i),ux(i),uy(i),uz(i),q(i) & 
         ,axdum, dum(1),dum(2),dum(3),dum(4),dum(5),dum(6) &
          ,i=1,npp+nadd)
     close (60)

     npp = npp+nadd
     nmerge = 1  !  next restart will use same # cpus by default

  else
     ! num_pe < # data sets, so merge together

     ioffset = 0
     npp = 0
     if (nmerge>1) write (*,*) 'Merging sets by ',nmerge,' : 1'

     do ipass = 0, nmerge-1
        me_read = me*nmerge + ipass
        ! Filename (directory) prefix
        cme = "data/pe" &   
       // achar(me_read/1000+48) &
       // achar(mod(me_read/100,10)+48) &
       // achar(mod(me_read/10,10)+48) &
       // achar(mod(me_read,10)+48)  ! Convert 4-digit PE number into character string

        write(*,'(a,a)') 'Reading particles from ',cme
        ! get filename suffix from dump counter
        do i=0,4
           cdump(6-i:6-i) =  achar(mod(itime_start/10**i,10) + 48)  
        end do
        cdump(1:1) = achar(itime_start/10**5 + 48)

        cfile=cme//"/parts_info."//cdump(1:6)


        open (60,file=cfile)    
        read(60,*)  timestamp
        read(60,*)  npp_partial  ! Find # particles to be read 
        close(60)


        cfile=cme//"/parts_dump."//cdump(1:6)
        open (60,file=cfile) 

        i1 = ioffset + 1
        i2 = ioffset + npp_partial

        !  Initialise particles: read from file in PEGS format
        read(60,*) (x(i),y(i),z(i),ux(i),uy(i),uz(i),q(i),m(i), &
             ax(i), ay(i), az(i), &  ! electric field = m.a/q
             pot(i), &  ! potential
             pepid(i), pelabel(i),i=i1,i2)

        close (60)

        ioffset = ioffset + npp_partial
        npp = npp + npp_partial
     enddo
     nmerge = 1  ! next restart will use same # cpus by default

  endif

  !  npp is the local number of dust particles (~ npart/num_pe) 

  !  now initialise/reset remaining particle properties

  pepid(1:npp) = me                ! processor ID
  pelabel(1:npp) = me*npp + (/ (i,i=1,npp) /)     ! Dust labels

  ! reset disc masses according to info file
  q(1:npp)=mdisc/npart
  m(1:npp)=mdisc/npart

  ax(1:npp) = 0.       ! zero accelerations until first force comp.
  ay(1:npp) = 0.
  az(1:npp) = 0.
  work(1:npp) = 1.   ! load-balancing weight

end subroutine setup_particles


