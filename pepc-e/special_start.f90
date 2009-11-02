! ==============================================
!
!                SPECIAL_START
!
!  Initialise set of particles with zero velocity
!
! ==============================================


! portable random number generator, see numerical recipes
! check for the random numbers:
! the first numbers should be 0.2853809, 0.2533582 and 0.2533582

subroutine par_rand(res)
  
  use physvars
  implicit none

  integer :: idum, idum2, iy, j, k
  integer :: iv(32)

  real, intent(inout) :: res

  integer :: IM1, IM2, IMM1, IA1, IA2, IQ1, IQ2, IR1, IR2, NTAB, NDIV
  real    :: AM, RNMX

  save

  data idum, idum2 /-1, 123456789/

  IM1 = 2147483563
  IM2 = 2147483399
  AM  = 1.0/IM1
  IMM1 = IM1-1 
  IA1 = 40014
  IA2 = 40692
  IQ1 = 53668
  IQ2 = 52774
  IR1 = 12211
  IR2 = 3791
  NTAB = 32
  NDIV = 1+IMM1/NTAB
  RNMX = 1.0 - 1.2e-7

  if (idum < 0) then
  
     idum = 1
     idum2 = idum

     do j = NTAB+7,0,-1
        k = idum/IQ1
        idum = IA1 * (idum-k*IQ1) - k*IR1
        if (idum < 0 ) idum = idum + IM1
        
        if (j<NTAB) iv(j+1) = idum
        
     end do
     iy = iv(1)
  
  end if

  k = idum/IQ1
  idum = IA1 * (idum-k*IQ1) - k*IR1
  if (idum < 0) idum = idum + IM1

  k = idum2/IQ2
  idum2 = IA2 * (idum2-k*IQ2) - k*IR2
  if (idum2 < 0) idum2 = idum2 + IM2
  
  j = iy/NDIV + 1
  iy = iv(j)-idum2
  iv(j) = idum
  
  if (iy < 1) iy = iy + IMM1
  res = AM*iy
  if (res > RNMX) res = RNMX
  
end subroutine par_rand

subroutine special_start(iconf)

  use physvars
  use utils
  implicit none
  include 'mpif.h'

  integer, intent(in) :: iconf  ! Configuration switch
  integer :: p,i,j,k,l,n1,iseed1,iseed2,iseed3,face_nr,iseed=-17,decomp,max_num,plt, mpi_cnt, ierr
  real*4 :: par_rand_res
  real*8 :: gamma0,yt,zt,xt,qt,mt,c_status, r1
  character(50) :: cinfile, cdump
  character(50) :: dfile
  character(35) :: cme
  integer :: np_local_max

  config: select case(iconf)
  case(1)

     if (my_rank == 0) write(*,*) "Using special start... case 1 (homogeneous distribution)"

     ! get the largest np_local
     call MPI_REDUCE(np_local, np_local_max, 1, MPI_INTEGER, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
     call MPI_BCAST(np_local_max, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

     do mpi_cnt = 0, n_cpu-1
        do p = 1, np_local_max
           
           xt = 0.
           yt = 0.
           zt = 0.
           
           call par_rand(par_rand_res)
           xt = par_rand_res
           call par_rand(par_rand_res)
           yt = par_rand_res
           call par_rand(par_rand_res)
           zt = par_rand_res

           if ( my_rank == mpi_cnt .and. p <= np_local ) then
              
              ux(p) = 0.
              uy(p) = 0.
              uz(p) = 0.
              
              z(p) = zt 
              y(p) = yt
              x(p) = xt
              
           end if
        end do
     end do

  case(2)

     if (my_rank == 0) write(*,*) "Using special start... case 2 (one sphere benchmark)"

     ! get the largest np_local
     call MPI_REDUCE(np_local, np_local_max, 1, MPI_INTEGER, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
     call MPI_BCAST(np_local_max, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

     do mpi_cnt = 0, n_cpu-1
        do p = 1, np_local_max
           
           xt = 1.
           yt = 1.
           zt = 1.
           
           do while ( (xt*xt + yt*yt + zt*zt) > 1)
              call par_rand(par_rand_res)
              xt = -1.0 + 2.*par_rand_res
              call par_rand(par_rand_res)
              yt = -1.0 + 2.*par_rand_res
              call par_rand(par_rand_res)
              zt = -1.0 + 2.*par_rand_res
           end do
           
           if ( my_rank == mpi_cnt .and. p <= np_local ) then

              xt = xt*0.1
              yt = yt*0.1
              zt = zt*0.1
              
              ux(p) = 0.
              uy(p) = 0.
              uz(p) = 0.
              
              z(p) = 0.5 + zt 
              y(p) = 0.5 + yt
              x(p) = 0.5 + xt
              
           end if
        end do
     end do

  case(3)
     if (my_rank == 0) write(*,*) "Using special start... case 3 (two sphere benchmark)"

     ! get the largest np_local
     call MPI_REDUCE(np_local, np_local_max, 1, MPI_INTEGER, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
     call MPI_BCAST(np_local_max, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

     do mpi_cnt = 0, n_cpu-1
        do p = 1, np_local_max
           
           xt = 1.
           yt = 1.
           zt = 1.
           
           do while ( (xt*xt + yt*yt + zt*zt) > 1)
              call par_rand(par_rand_res)
              xt = -1.0 + 2.*par_rand_res
              call par_rand(par_rand_res)
              yt = -1.0 + 2.*par_rand_res
              call par_rand(par_rand_res)
              zt = -1.0 + 2.*par_rand_res
           end do
           
           xt = xt*0.1
           yt = yt*0.1
           zt = zt*0.1
           
           call par_rand(par_rand_res)
           if (par_rand_res .lt. 0.5) then
              xt = xt + 0.2
           else
              xt = xt - 0.2
           end if
           
           if ( my_rank == mpi_cnt .and. p <= np_local ) then
                            
              ux(p) = 0.
              uy(p) = 0.
              uz(p) = 0.
              
              z(p) = 0.5 + zt 
              y(p) = 0.5 + yt
              x(p) = 0.5 + xt
              
           end if
        end do
     end do

  case(4)

     if (my_rank == 0) write(*,*) "Using special start... case 4: Plummer distribution (core cut)"

     ! get the largest np_local
     call MPI_REDUCE(np_local, np_local_max, 1, MPI_INTEGER, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
     call MPI_BCAST(np_local_max, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

     do mpi_cnt = 0, n_cpu-1
        do p = 1, np_local_max
     
           xt = 0.
           yt = 0.
           zt = 0.
           r1 = 0.
           


           do while (.not.(r1 .gt. 0.1 .and. r1 .lt. 3))
              call par_rand(par_rand_res)
              r1 = (par_rand_res**(-0.2D01/0.3D01)-0.1D01)**(-0.5D00)
           end do
           call par_rand(par_rand_res)
           zt = (0.1D01-0.2D01*par_rand_res)*r1
           call par_rand(par_rand_res)
           xt = (r1**0.2D01-zt**0.2D01)**(0.5D00)*cos(0.2D01*pi*par_rand_res)
           yt = (r1**0.2D01-zt**0.2D01)**(0.5D00)*sin(0.2D01*pi*par_rand_res)                     

           if ( my_rank == mpi_cnt .and. p <= np_local) then

              ux(p) = 0.
              uy(p) = 0.
              uz(p) = 0.
              
              z(p) = zt 
              y(p) = yt
              x(p) = xt

           end if
        end do
     end do

  end select config

  q(1:nep)                  = qe        ! plasma electrons
  q(nep + 1:np_local)       = qi        ! plasma ions (need Z* here)
  m(1:nep)                  = mass_e    ! electron mass
  m(nep + 1:np_local)       = mass_i    ! ion mass
  pelabel(1:nep)            = my_rank * nep + (/(i, i = 1, nep)/)      ! Electron labels
  pelabel(nep + 1:np_local) = ne + my_rank * nip + (/(i, i = 1, nip)/) ! Ion labels

112 ex(1:np_local) = 0.
  ey(1:np_local) = 0.
  ez(1:np_local) = 0.
  pot(1:np_local) = 0.

  work(1:np_local) = 1.

end subroutine special_start
