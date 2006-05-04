!  ================================
!
!         SETUP
!
!   $Revision: 1.13 $
!
!     Initialise constants and 
!      simulation variables
!
!  ================================


subroutine setup

  use treevars
  use physvars
  use utils

  implicit none
  include 'mpif.h'

  integer :: ierr,ibig, machinebits, maxleaf, maxtwig,k, npb_pe
  integer :: root=0
  real :: q_factor, x_offset, z_offset, fnn
  character*7 :: ensembles(1:5)=(/ &
       'const U','Te, Ti ','glob Te','loc  Te','Ti only' /)

  integer :: i, ipe, idummy, ifile

  character(30) :: cinfile, cdump, cfile
  character(9) :: ct
  character(5) :: cme
  integer :: timestamp
  real :: axdum, aydum, azdum,phidum, dum(6)
  integer :: ioffset, i1, i2, npp_partial, npp_total, ipass, me_read, j, nrest, nadd



  !  Default input set

  ! switches
  perf_anal = .false.
  tree_debug = .false.
  domain_debug = .false.
  load_balance = .true.
  walk_balance = .true.
  walk_debug = .false.
  coulomb = .true.
  initial_config = 1         ! random sphere

  ! particles
  nep = 0 ! # plasma electrons per PE
  nip = 0
  ne = 100  ! Total # plasma electrons
  ni = 100  ! total # plasma ions
  np_beam = 0  ! Total # beam particles

  xl = 2
  yl = 2
  zl = 2

  ! physics stuff
  force_const = 1./3.
  rho0 = 1.0
  theta = 0.8
  mdisc = 0.1 

  eps = 0.1
  fnn = 5  ! Neighbour search radius multiplier (x a_ii)
  displace(1:3) = (/0.,0.,0./)

  ! control
  nt = 600
  dt = 0.2
  trun = 0.
  vis_on=.true.
  ivis = 1
  ivis_fields = 1
  itime_start = 0
  idens=1
  ngx = 200   ! Grid size for plots
  ngy = 200
  ngz = 10



  ! Get run parameters from 'parts_info.in'

  if (me == 0) then 

     ! input file for PE me

     cinfile="parts_info.in"
     open(80,file=cinfile)

     ! Read in info block, skipping variable names

     read(80,*)  itime_start 
     read(80,*)  npart 
     read(80,*)  ne 
     read(80,*)  ni
     read(80,*)  scheme
     read(80,*)  mdisc
     read(80,*)  xl
     read(80,*)  yl
     read(80,*)  zl
     read(80,*)  boxsize
     read(80,*)  box_centre(1:3)  ! shift vector applied to all particles at start of run
     read(80,*)  eps  
     read(80,*)  theta  
     read(80,*)  nmerge 
     read(80,*)  nt     ! # timesteps
     read(80,*)  dt     ! initial timestep
     read(80,*) ivis    !  particle visualisation frequency 
     read(80,*) ivis_fields ! field vis. freq.
     read(80,*) idump  ! dump frequency
     read(80,*) iprot  !  written o/p freq.


     !     write(*,*) "itime ",itime_start," scheme ",scheme     
     !     write(*,*) "npart ",npart,"  ndust ", ne,"  nstar ",ni
     !     write(*,'(4(a9,f12.2))') "xl ",xl,"  yl ",yl,"  zl ",zl,"  shift ",box_centre
     !    write(*,'(3(a9,f12.2))') "mdisc ",mdisc,"  eps ",eps,"  theta",theta
     close(80)



     write(6,'(/a/a/a,i5)') 'RESTART:','Read run data from info block: parts_info.in','Timestep:',itime_start

     trun = itime_start*dt

  endif



  ! Send global input parameters to all other processes

  call MPI_BARRIER( MPI_COMM_WORLD, ierr)   ! Synchronize first

  call MPI_BCAST( xl, 1, MPI_REAL8, root, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( yl, 1, MPI_REAL8, root, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( zl, 1, MPI_REAL8, root, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( boxsize, 1, MPI_REAL8, root, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( box_centre, 3, MPI_REAL8, root, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( theta, 1, MPI_REAL8, root, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( eps, 1, MPI_REAL8, root, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( mdisc, 1, MPI_REAL8, root, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( dt, 1, MPI_REAL8, root, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( ne, 1, MPI_INTEGER, root, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( ni, 1, MPI_INTEGER, root, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( scheme, 1, MPI_INTEGER, root, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( npart, 1, MPI_INTEGER, root, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( itime_start, 1, MPI_INTEGER, root, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( nt, 1, MPI_INTEGER, root, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( ivis, 1, MPI_INTEGER, root, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( ivis_fields, 1, MPI_INTEGER, root, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( idump, 1, MPI_INTEGER, root, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( iprot, 1, MPI_INTEGER, root, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( nmerge, 1, MPI_INTEGER, root, MPI_COMM_WORLD,ierr)



  !  Array allocation
!  call setup_arrays(npart)
  call pepc_setup(my_rank,n_cpu,npart,theta,debug_tree,np_mult,fetch_mult) 


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
     cme = "pe"// achar(me_read/100+48) // achar(mod(me_read/10,10)+48) // achar(mod(me_read,10)+48)  

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
        cme = "pe"// achar(me_read/100+48) // achar(mod(me_read/10,10)+48) // achar(mod(me_read,10)+48)  
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


  !   STARS - read in using root CPU

  if (me.eq.0) then

     cinfile="star_info."//cdump(1:6)    ! input file for stars
     open(81,file=cinfile)

     do i=1,ni
        read(81,*) m_star(i)
     end do
     do i=1,ni
        read(81,*) x_star(i), y_star(i), z_star(i)
     end do
     do i=1,ni
        read(81,*) ux_star(i), uy_star(i), uz_star(i)
     end do
     write(*,'(a14,f12.5,a7,3f12.5)') &
          "star 1 mass: ",m_star(1),"  pos ",x_star(1),y_star(1),z_star(1)
     write(*,'(a14,f12.5,a7,3f12.5)') &
          "star 2 mass: ",m_star(2),"  pos ",x_star(2),y_star(2),z_star(2)
     close(81)

  endif

  ! Broadcast star data to other CPUs
  call MPI_BCAST( x_star, ni, MPI_REAL8, 0, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( y_star, ni, MPI_REAL8, 0, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( z_star, ni, MPI_REAL8, 0, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( ux_star, ni, MPI_REAL8, 0, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( uy_star, ni, MPI_REAL8, 0, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( uz_star, ni, MPI_REAL8, 0, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( m_star, ni, MPI_REAL8, 0, MPI_COMM_WORLD,ierr)


  !  Derive some constants from input data

!  box_centre =  (/ xl/4., 0., 0. /) ! Centre of box

  qe = mdisc/npart
  qi = 1.0
  mass_e = 1


  if (itime_start.eq.0) then
  !  Translate stars and dust from (0,0,0) to centre of box

  x_star(1:ni) = x_star(1:ni) + box_centre(1)
  y_star(1:ni) = y_star(1:ni) + box_centre(2)
  z_star(1:ni) = z_star(1:ni) + box_centre(3)

  x(1:npp) = x(1:npp) + box_centre(1)
  y(1:npp) = y(1:npp) + box_centre(2)
  z(1:npp) = z(1:npp) + box_centre(3)

     write(ipefile,'(a14,f12.5,a17,3f12.5)') &
          "star 1 mass: ",m_star(1)," shifted pos ",x_star(1),y_star(1),z_star(1)
     write(ipefile,'(a14,f12.5,a17,3f12.5)') &
          "star 2 mass: ",m_star(2)," shifted pos ",x_star(2),y_star(2),z_star(2)
  endif


  ! Write out run info to terminal and .out file

  if (me==0) then
     do ifile = 6,15,9

        write (ifile,*) ' Scheme: ',ensembles(scheme)
        write (ifile,'(a20,1pe12.3)') ' mass of disc particles: ',qe
        write (ifile,'(a20,1pe12.3)') ' Max timestep: ',0.45*sqrt(3.)*eps**2/qe
        write (ifile,'(a20,1pe12.3)') ' Neighbour search radius: ',r_neighbour


        write (ifile,*) ' no.gas/dust: ', ne
        write (ifile,*) ' no.stars: ', ni

        write (ifile,*) ' Particles per PE: ', npp


        write (ifile,'(/a/a)') ' Switches:','--------'
        write (ifile,'(a20,l3)') ' Gravitational forces: ',coulomb

        write (ifile,'(a20,l3)') ' load balance: ',load_balance
        write (ifile,'(a20,l3)') ' walk balance: ',walk_balance
        write (ifile,'(a20,l3)') ' restart: ',restart

        write (ifile,'(a20,l3)') ' domain debug: ',domain_debug
        write (ifile,'(a20,l3)') ' walk debug: ',walk_debug
        write (ifile,'(a20,l3)') ' dump tree: ',dump_tree
        write (ifile,'(a20,l3)') ' performance anal.: ',perf_anal
        write (ifile,'(a20,l3)') ' visit: ',vis_on
        write (ifile,'(a20,l3/)') ' steering: ',steering
        write (ifile,*) 'Max address in #-table: ',2**nbaddr-1
        machinebits = bit_size(1_8)    ! # bits in integer variable (hardware) 
        write (ifile,*) 'Machine bit-size = ',machinebits
        write (ifile,*) 'Max permitted particles / PE: ', nppm
        write (ifile,*) 'Max size of interaction lists: ', nintmax
	write (ifile,*) 'Shortlist length: ',nshortm
        write (ifile,*) 'Memory needed for lists = ',4*nintmax*nshortm*8/2**20,' Mbytes'

        write (ifile,'(/a)') 'Other inputs:'
        !	write(ifile,NML=pepcdata)

     end do
  endif



end subroutine setup

