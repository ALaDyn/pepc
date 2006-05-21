! ==============================================
!
!                SETUP_STARS
!
!  Read in star data from configuration file
!  or particle dump
!
! ==============================================

subroutine setup_stars

  use treevars
  use physvars

  implicit none
  include 'mpif.h'

  integer :: i, ipe, idummy=0, ifile, ierr
  character(30) :: cinfile, cdump, cfile
  character(9) :: ct
  character(5) :: cme


  !   STARS - read in using root CPU

     ! get filename suffix from dump counter
     do i=0,4
        cdump(6-i:6-i) =  achar(mod(itime_start/10**i,10) + 48)  
     end do
     cdump(1:1) = achar(itime_start/10**5 + 48)

  if (my_rank.eq.0) then
     write(*,*) 'Setting up stars ...'
     cinfile="star_info."//cdump(1:6)    ! input file for stars
     open(81,file=cinfile)

     do i=1,nstar
        read(81,*) m_star(i)
     end do
     do i=1,nstar
        read(81,*) x_star(i), y_star(i), z_star(i)
     end do
     do i=1,nstar
        read(81,*) ux_star(i), uy_star(i), uz_star(i)
     end do
     write(*,'(a14,f12.5,a7,3f12.5)') &
          "star 1 mass: ",m_star(1),"  pos ",x_star(1),y_star(1),z_star(1)
     write(*,'(a14,f12.5,a7,3f12.5)') &
          "star 2 mass: ",m_star(2),"  pos ",x_star(2),y_star(2),z_star(2)
     close(81)

  endif

  ! Broadcast star data to other CPUs
  call MPI_BCAST( x_star, nstar, MPI_REAL8, 0, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( y_star, nstar, MPI_REAL8, 0, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( z_star, nstar, MPI_REAL8, 0, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( ux_star, nstar, MPI_REAL8, 0, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( uy_star, nstar, MPI_REAL8, 0, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( uz_star, nstar, MPI_REAL8, 0, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( m_star, nstar, MPI_REAL8, 0, MPI_COMM_WORLD,ierr)


  !  Derive some constants from input data

!  box_centre =  (/ xl/4., 0., 0. /) ! Centre of box

  qe = mdisc/ndisc
  qi = 1.0
  mass_e = 1


  if (itime_start.eq.0) then
  !  Translate stars and dust from (0,0,0) to centre of box

  x_star(1:nstar) = x_star(1:nstar) + box_centre(1)
  y_star(1:nstar) = y_star(1:nstar) + box_centre(2)
  z_star(1:nstar) = z_star(1:nstar) + box_centre(3)

  x(1:npp) = x(1:npp) + box_centre(1)
  y(1:npp) = y(1:npp) + box_centre(2)
  z(1:npp) = z(1:npp) + box_centre(3)

     write(ipefile,'(a14,f12.5,a17,3f12.5)') &
          "star 1 mass: ",m_star(1)," shifted pos ",x_star(1),y_star(1),z_star(1)
     write(ipefile,'(a14,f12.5,a17,3f12.5)') &
          "star 2 mass: ",m_star(2)," shifted pos ",x_star(2),y_star(2),z_star(2)
  endif

  if (my_rank==0) write(*,*) '... done'

end subroutine setup_stars


