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

  use physvars
  use utils

  implicit none
  include 'mpif.h'

  integer :: ierr,ibig, machinebits, maxleaf, maxtwig,k, npb_pe
  integer :: root=0
  real :: q_factor, x_offset, z_offset, fnn

  integer :: i, ipe, idummy, ifile

  character(30) :: cinfile, cdump, cfile
  character(9) :: ct
  character(5) :: cme
  integer :: timestamp




  !  Default input set

  ! switches

  debug_tree=2
  initial_config = 1         ! random sphere

  ! particles

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

  if (my_rank == 0) then 

     ! input file for PE me

     cinfile="parts_info.in"
     open(80,file=cinfile)

     ! Read in info block, skipping variable names

     read(80,*)  itime_start 
     read(80,*)  npart_total
     read(80,*)  ndisc 
     read(80,*)  nstar
     read(80,*)  scheme
     read(80,*)  mdisc
     read(80,*)  xl
     read(80,*)  yl
     read(80,*)  zl
     read(80,*)  boxmax
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
  call MPI_BCAST( boxmax, 1, MPI_REAL8, root, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( box_centre, 3, MPI_REAL8, root, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( theta, 1, MPI_REAL8, root, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( eps, 1, MPI_REAL8, root, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( mdisc, 1, MPI_REAL8, root, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( dt, 1, MPI_REAL8, root, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( ndisc, 1, MPI_INTEGER, root, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( nstar, 1, MPI_INTEGER, root, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( scheme, 1, MPI_INTEGER, root, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( npart_total, 1, MPI_INTEGER, root, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( itime_start, 1, MPI_INTEGER, root, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( nt, 1, MPI_INTEGER, root, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( ivis, 1, MPI_INTEGER, root, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( ivis_fields, 1, MPI_INTEGER, root, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( idump, 1, MPI_INTEGER, root, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( iprot, 1, MPI_INTEGER, root, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( nmerge, 1, MPI_INTEGER, root, MPI_COMM_WORLD,ierr)



end subroutine setup

