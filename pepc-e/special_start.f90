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
  integer :: p,mpi_cnt, ierr
  real*4 :: par_rand_res,fr,mu
  real*4 :: yt,zt,xt,r1,dx
  integer :: np_local_max
  real*8 :: delta(3)
  integer :: i,j,k,n(3), myidx, globalidx
  real*8 a,b,c,cth,sth,cphi,sphi,s

  integer :: fances(-1:n_cpu-1)

  real*8 zero
  parameter(zero=0.d0)
  real*8 one
  parameter(one=1.d0)
  real*8 mone
  parameter(mone=-one)
  real*8 two
  parameter(two=2.d0)
  real*8 five
  parameter(five=5.d0)
  real*8 nine
  parameter(nine=9.d0)

  ! get the largest np_local
  call MPI_REDUCE(np_local, np_local_max, 1, MPI_INTEGER, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(np_local_max, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  call MPI_SCAN(np_local, fances(my_rank), 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLGATHER(MPI_IN_PLACE, 0, 0, fances(0), 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
  fances(-1) = 0

  config: select case(iconf)
  case(1)

     if (my_rank == 0) write(*,*) "Using special start... case 1 (homogeneous distribution)"

     do mpi_cnt = 0, n_cpu-1
        do p = 1, (fances(mpi_cnt) - fances(mpi_cnt-1))
           
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

     do mpi_cnt = 0, n_cpu-1
        do p = 1, (fances(mpi_cnt) - fances(mpi_cnt-1))
           
           xt = 1.0_8
           yt = 1.0_8
           zt = 1.0_8
           
           do while ( (xt*xt + yt*yt + zt*zt) > 1.0_8)
              call par_rand(par_rand_res)
              xt = -1.0_8 + 2.0_8*par_rand_res
              call par_rand(par_rand_res)
              yt = -1.0_8 + 2.0_8*par_rand_res
              call par_rand(par_rand_res)
              zt = -1.0_8 + 2.0_8*par_rand_res
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

     do mpi_cnt = 0, n_cpu-1
        do p = 1, (fances(mpi_cnt) - fances(mpi_cnt-1))
           
           xt = 1.0_8
           yt = 1.0_8
           zt = 1.0_8
           
           do while ( (xt*xt + yt*yt + zt*zt) > 1.0_8)
              call par_rand(par_rand_res)
              xt = -1.0_8 + 2.0_8*par_rand_res
              call par_rand(par_rand_res)
              yt = -1.0_8 + 2.0_8*par_rand_res
              call par_rand(par_rand_res)
              zt = -1.0_8 + 2.0_8*par_rand_res
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
        do p = 1, fances(mpi_cnt) - fances(mpi_cnt-1)
     
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


  case(5)
       
     if (my_rank == 0) write(*,*) "Using special start... case 5: 2D disc"
     ! get the largest np_local
     call MPI_REDUCE(np_local, np_local_max, 1, MPI_INTEGER, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
     call MPI_BCAST(np_local_max, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

     do mpi_cnt = 0, n_cpu-1
        do p = 1, fances(mpi_cnt) - fances(mpi_cnt-1)
           
           xt = 1.
           yt = 1.
           zt = 0.
           
           do while ( (xt*xt + yt*yt) > 1)
              call par_rand(par_rand_res)
              xt = -1.0 + 2.*par_rand_res
              call par_rand(par_rand_res)
              yt = -1.0 + 2.*par_rand_res
           end do
           
           if ( my_rank == mpi_cnt .and. p <= np_local ) then

              xt = xt*0.1
              yt = yt*0.1
              
              ux(p) = 0.
              uy(p) = 0.
              uz(p) = 0.
              
              z(p) = 0.5
              y(p) = 0.5 + yt
              x(p) = 0.5 + xt
              
           end if
        end do
     end do



  case(6)
     
     if (my_rank == 0) write(*,*) "Using special start... case 6 (2D-homogeneous distribution)"

     ! get the largest np_local
     call MPI_REDUCE(np_local, np_local_max, 1, MPI_INTEGER, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
     call MPI_BCAST(np_local_max, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     
     do mpi_cnt = 0, n_cpu-1
        do p = 1, fances(mpi_cnt) - fances(mpi_cnt-1)
           
           xt = 0.
           yt = 0.
           zt = 0.
           
           call par_rand(par_rand_res)
           xt = par_rand_res
           call par_rand(par_rand_res)
           yt = par_rand_res
           call par_rand(par_rand_res)
           zt = 0
           
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

  case(9)

     if (my_rank == 0) write(*,*) "Using special start... case 9 (ion sphere Coulomb explosion)"

     ! Rescale eps to ensure Coulomb energy correctly computed
     eps = r_sphere/100.

     ! get the largest np_local
     call MPI_REDUCE(np_local, np_local_max, 1, MPI_INTEGER, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
     call MPI_BCAST(np_local_max, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

     do mpi_cnt = 0, n_cpu-1
        do p = 1, fances(mpi_cnt) - fances(mpi_cnt-1)
           
           xt = 1.
           yt = 1.
           zt = 1.
           
           do while ( (xt*xt + yt*yt + zt*zt) > 1)
              call par_rand(par_rand_res)
              xt = -1.+2.*par_rand_res
              call par_rand(par_rand_res)
              yt = -1.+2.*par_rand_res
              call par_rand(par_rand_res)
              zt = -1.+2.*par_rand_res
           end do
           
           if ( my_rank == mpi_cnt .and. p <= np_local ) then

              ux(p) = 0.
              uy(p) = 0.
              uz(p) = 0.
              
              z(p) = zl/2. + r_sphere*xt
              y(p) = yl/2. + r_sphere*yt
              x(p) = xl/2. + r_sphere*zt
              
           end if
        end do
     end do

     q(1:np_local)       = 4*pi/3.*rho0*r_sphere**3/ni       ! ion charge
     m(1:np_local)       = 4*pi/3.*rho0*r_sphere**3/ni   ! ion mass
     ux(1:np_local)      = 0.
     uy(1:np_local)      = 0.
     uz(1:np_local)      = 0.

     pelabel(1:np_local) = my_rank * nip + (/(i, i = 1, nip)/) ! Ion labels

     ex(1:np_local) = 0.
     ey(1:np_local) = 0.
     ez(1:np_local) = 0.
     pot(1:np_local) = 0.

     work(1:np_local) = 1.

     !  on rank 0 redefine the last ngx particles for gridded diagnostics
     dx = xl/ngx
     if (my_rank==0) then
       q(np_local-ngx+1:np_local)=0.
       m(np_local-ngx+1:np_local) = 10000.
       x(np_local-ngx+1:np_local) = xl/2. + (/ (i*dx-dx, i=1,ngx) /)  ! Sample positions
       y(np_local-ngx+1:np_local) = yl/2.
       z(np_local-ngx+1:np_local) = zl/2.
       pelabel(np_local-ngx+1:np_local) = ni + (/(i, i = 1, ngx)/) ! grid particle labels
     endif

     return ! skip rest of routine

  case(10)

     if (my_rank == 0) write(*,*) "Using special start... case 10 (Hollow ion sphere)"

     ! get the largest np_local
     call MPI_REDUCE(np_local, np_local_max, 1, MPI_INTEGER, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
     call MPI_BCAST(np_local_max, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

     a = nine/five*sqrt(dble(npart_total-1)/dble(npart_total))
     j = 0

     if (my_rank==0) npart_total = npart_total-ngx
     do i = 1,npart_total

        if(i.eq.1) then
           cth = mone
           sth = zero
           cphi = one
           sphi = zero
           b = cphi
           c = sphi
        elseif(i.eq.npart_total) then
           cth = one
           sth = zero
           cphi = one
           sphi = zero
        else
           cth = dble(2*i-npart_total-1)/dble(npart_total-1)
           sth = two*(sqrt(dble(i-1)/dble(npart_total-1))*sqrt(dble(npart_total-i)/dble(npart_total-1)))
           s = a*sqrt(dble(npart_total-1)/(dble(i-1)*dble(npart_total-i)))
           cphi = b*cos(s)-c*sin(s)
           sphi = c*cos(s)+b*sin(s)
           b = cphi
           c = sphi
        endif
        if (mod(i+my_rank,n_cpu) == 0) then
           j = j+1
           x(j) = xl/2+r_sphere*sth*cphi
           y(j) = yl/2+r_sphere*sth*sphi
           z(j) = zl/2+r_sphere*cth
        end if
     end do

     if (my_rank==0) npart_total = npart_total+ngx

     q(1:j)       = 4*pi/3.*rho0*r_sphere**3/ni       ! ion charge
     m(1:j)       = 4*pi/3.*rho0*r_sphere**3/ni   ! ion mass

     pelabel(1:j) = my_rank * nip + (/(i, i = 1, nip)/) ! Ion labels

     ex(1:np_local) = 0.
     ey(1:np_local) = 0.
     ez(1:np_local) = 0.
     pot(1:np_local) = 0.

     work(1:np_local) = 1.

!  On task 0 redefine the last ngx particles for gridded diagnostics
     dx = xl/ngx
     if (my_rank==0) then
       q(j+1:np_local) = 0.
       m(j+1:np_local) = 10000.
       x(j+1:np_local) = xl/2. + (/ (i*dx-dx, i=1,np_local-(j+1)) /)  ! Sample positions
       y(j+1:np_local) = yl/2.
       z(j+1:np_local) = zl/2.
       pelabel(j+1:np_local) = ni + (/(i, i = 1,np_local-(j+1))/) ! grid particle labels
     endif

     return ! skip rest of routine

  case(11)
     if (my_rank == 0) write(*,*) "Using special start... case 11 (ion sphere Coulomb explosion with density profile from [Phys. Rev. Letters 91, 143401 (2003)])"

     ! Rescale eps to ensure Coulomb energy correctly computed
     eps = r_sphere/100.

     ! get the largest np_local
     call MPI_REDUCE(np_local, np_local_max, 1, MPI_INTEGER, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
     call MPI_BCAST(np_local_max, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

     do mpi_cnt = 0, n_cpu-1
        do p = 1, fances(mpi_cnt) - fances(mpi_cnt-1)
           
           xt = 1.
           yt = 1.
           zt = 1.
           par_rand_res = 100.
           fr = 0.  
           mu = 1.
           
           do while ( fr < par_rand_res)
              call par_rand(par_rand_res)
              xt = -1.+ 8.*par_rand_res*r_sphere
              call par_rand(par_rand_res)
              yt = -1.+ 8.*par_rand_res*r_sphere
              call par_rand(par_rand_res)
              zt = -1.+ 8.*par_rand_res*r_sphere
              r1 = sqrt(xt*xt + yt*yt + zt*zt) ! normalized radius
              fr = 1/(  (  1 +  r1**(3.*mu) )**(1./mu+1.)  )
              call par_rand(par_rand_res)
           end do
           
           if ( my_rank == mpi_cnt .and. p <= np_local ) then

              ux(p) = 0.
              uy(p) = 0.
              uz(p) = 0.
              
              z(p) = zl/2. + r_sphere*xt
              y(p) = yl/2. + r_sphere*yt
              x(p) = xl/2. + r_sphere*zt
              
           end if
        end do
     end do

     q(1:np_local)       = 4*pi/3.*rho0*r_sphere**3/ni       ! ion charge
     m(1:np_local)       = 4*pi/3.*rho0*r_sphere**3/ni   ! ion mass
     ux(1:np_local)      = 0.
     uy(1:np_local)      = 0.
     uz(1:np_local)      = 0.

     pelabel(1:np_local) = my_rank * nip + (/(i, i = 1, nip)/) ! Ion labels

     ex(1:np_local) = 0.
     ey(1:np_local) = 0.
     ez(1:np_local) = 0.
     pot(1:np_local) = 0.

     work(1:np_local) = 1.

     !  on rank 0 redefine the last ngx particles for gridded diagnostics
     dx = xl/ngx
     if (my_rank==0) then
       q(np_local-ngx+1:np_local)=0.
       m(np_local-ngx+1:np_local) = 10000.
       x(np_local-ngx+1:np_local) = xl/2. + (/ (i*dx-dx, i=1,ngx) /)  ! Sample positions
       y(np_local-ngx+1:np_local) = yl/2.
       z(np_local-ngx+1:np_local) = zl/2. 
       pelabel(np_local-ngx+1:np_local) = ni + (/(i, i = 1, ngx)/) ! grid particle labels
     endif

     return ! skip rest of routine


  end select config

  q(1:nep)                  = qe        ! plasma electrons
  q(nep + 1:np_local)       = qi        ! plasma ions (need Z* here)
  m(1:nep)                  = mass_e    ! electron mass
  m(nep + 1:np_local)       = mass_i    ! ion mass
  pelabel(1:nep)            = fances(my_rank-1) + (/(i, i = 1, nep)/)      ! Electron labels
  pelabel(nep + 1:np_local) = ne + fances(my_rank-1) + nep + (/(i, i = 1, (fances(my_rank) - fances(my_rank-1) - nep))/) ! Ion labels

  ex(1:np_local) = 0.
  ey(1:np_local) = 0.
  ez(1:np_local) = 0.
  pot(1:np_local) = 0.

  work(1:np_local) = 1.

end subroutine special_start
