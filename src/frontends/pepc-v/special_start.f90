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
subroutine par_rand(res, iseed)
  
  use physvars
  implicit none

  integer :: idum, idum2, iy, j, k
  integer :: iv(32)

  real, intent(inout) :: res
  integer, intent(in), optional :: iseed

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
  
     if (present(iseed)) then
       idum = iseed
     else
       idum = 1
     endif

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
     implicit none
     include 'mpif.h'

     integer :: j, k, ind, ind0, i, m
     real*8 :: part_2d, rc, xi1, xi2, xi, part_3d, eta1, eta2, eta, v(3)
     real*8, dimension(3,3) :: D1, D2, D3, D4   !< rotation matrices for ring setup
     real*8, allocatable :: xp(:), yp(:), zp(:), volp(:), wxp(:), wyp(:), wzp(:)  !< helper arrays for ring setups

!     interface
!         subroutine par_rand(res, iseed)
!            real, intent(inout) :: res
!            integer, intent(in), optional :: iseed
!         end subroutine
!     end interface

      ! Set up particle data
     config: select case(ispecial)
     case(1)                               ! Vortex ring setup, side-by-side

     allocate(xp(ns), yp(ns), zp(ns), volp(ns), wxp(ns), wyp(ns), wzp(ns))

     j = 0

     do k = 1, nc
        part_2d = 2*pi/(8*k)
        rc = (1+12*k**2)/(6*k)*rl

        do l = 1,8*k
           j = j+1
           xi1 = part_2d*(l-1)
           xi2 = part_2d*l
           xi = (xi2-xi1)/2+xi1
           xp(j) = rc*cos(xi)
           yp(j) = rc*sin(xi)
           zp(j) = 0
           volp(j) = (2*pi**2*(r_torus+(2*k+1)*rl)*((2*k+1)*rl)**2-2*pi**2*(r_torus+(2*k-1)*rl)*((2*k-1)*rl)**2)/(8*k*Nphi)
           wxp(j) = 0.
           wyp(j) = 0.
           wzp(j) = g*exp(-(rc/rmax)**2)
        end do
     end do

     xp(ns) = 0.
     yp(ns) = 0.
     zp(ns) = 0.
     wxp(ns) = 0.
     wyp(ns) = 0.
     wzp(ns) = g
     volp(ns) = 2*pi**2*(r_torus+rl)*rl**2/Nphi

     j = 0
     ind0 = 0
     ind = 0
     part_3d = 2*pi/Nphi
     do m = 1, Nphi
        eta1 = part_3d*(m-1)
        eta2 = part_3d*m
        eta = (eta2 - eta1)/2+eta1
        v(1) = cos(eta2)
        v(2) = sin(eta2)
        v(3) = 0
        D1 = reshape( (/ -1.0+2*v(1)**2, 2*v(2)*v(1), 0.0D0,  2*v(1)*v(2), -1.0+2*v(2)**2, 0.0D0, 0.0D0, 0.0D0, -1.0D0 /), (/3,3/) )
        v(1) = cos(eta)
        v(2) = sin(eta)
        v(3) = 0
        D2 = reshape( (/ v(1)**2, v(2)*v(1), -v(2),  v(1)*v(2), v(2)**2, v(1), v(2), -v(1), 0.0D0 /), (/3,3/) )
        D3 = reshape( (/ -1.0+2*v(1)**2, 2*v(2)*v(1), 0.0D0,  2*v(1)*v(2), -1.0+2*v(2)**2, 0.0D0, 0.0D0, 0.0D0, -1.0D0 /), (/3,3/) )
        D4 = matmul(D1,D3)
        do i = 1, Ns
           if (m==1) then
              v(1) = xp(i) + (r_torus+rmax)*cos(eta)
              v(2) = yp(i) + (r_torus+rmax)*sin(eta)
              v(3) = zp(i);
              xp(i) = dot_product(v,D2(1:3,1))
              yp(i) = dot_product(v,D2(1:3,2))
              zp(i) = dot_product(v,D2(1:3,3))
              v(1) = wxp(i);
              v(2) = wyp(i);
              v(3) = wzp(i);
              wxp(i) = dot_product(v,D2(1:3,1))
              wyp(i) = dot_product(v,D2(1:3,2))
              wzp(i) = dot_product(v,D2(1:3,3))
           else
              v(1) = xp(i)
              v(2) = yp(i)
              v(3) = zp(i)
              xp(i) = dot_product(v,D4(1:3,1))
              yp(i) = dot_product(v,D4(1:3,2))
              zp(i) = dot_product(v,D4(1:3,3))
              v(1) = wxp(i)
              v(2) = wyp(i)
              v(3) = wzp(i)
              wxp(i) = dot_product(v,D4(1:3,1))
              wyp(i) = dot_product(v,D4(1:3,2))
              wzp(i) = dot_product(v,D4(1:3,3))
           end if
           ind0 = ind0 + 1
           if (mod(ind0-1,n_cpu) == my_rank) then
              ind = ind + 1
              if (ind .gt. nppm-1) then
                 write(*,*) 'something is wrong here: to many particles in init',my_rank,ind,nppm,n
                 call MPI_ABORT(MPI_COMM_WORLD,ierr)
                 stop
              end if
              part(ind)%x(1) = xp(i) - (torus_offset(1)-(rmax-(1+12*nc**2)/(6*nc)*rl))/2.0
              part(ind)%x(2) = yp(i)
              part(ind)%x(3) = zp(i)
              part(ind)%alpha(1) = wxp(i)*volp(i)
              part(ind)%alpha(2) = wyp(i)*volp(i)
              part(ind)%alpha(3) = wzp(i)*volp(i)
              ind = ind + 1
              part(ind)%x(1) = xp(i) + (torus_offset(1)-(rmax-(1+12*nc**2)/(6*nc)*rl))/2.0
              part(ind)%x(2) = yp(i)
              part(ind)%x(3) = zp(i)
              part(ind)%alpha(1) = wxp(i)*volp(i)
              part(ind)%alpha(2) = wyp(i)*volp(i)
              part(ind)%alpha(3) = wzp(i)*volp(i)
           end if
        end do
     end do
     np = ind
     deallocate(xp,yp,zp,volp,wxp,wyp,wzp)

  case(2) ! Vortex ring setup, offset collision

     allocate(xp(ns), yp(ns), zp(ns), volp(ns), wxp(ns), wyp(ns), wzp(ns))

     j = 0

     do k = 1, nc
        part_2d = 2*pi/(8*k)
        rc = (1+12*k**2)/(6*k)*rl

        do l = 1,8*k
           j = j+1
           xi1 = part_2d*(l-1)
           xi2 = part_2d*l
           xi = (xi2-xi1)/2+xi1
           xp(j) = rc*cos(xi)
           yp(j) = rc*sin(xi)
           zp(j) = 0
           volp(j) = (2*pi**2*(r_torus+(2*k+1)*rl)*((2*k+1)*rl)**2-2*pi**2*(r_torus+(2*k-1)*rl)*((2*k-1)*rl)**2)/(8*k*Nphi)
           wxp(j) = 0.
           wyp(j) = 0.
           wzp(j) = g*exp(-(rc/rmax)**2)
        end do
     end do

     xp(ns) = 0.
     yp(ns) = 0.
     zp(ns) = 0.
     wxp(ns) = 0.
     wyp(ns) = 0.
     wzp(ns) = g
     volp(ns) = 2*pi**2*(r_torus+rl)*rl**2/Nphi

     j = 0
     ind0 = 0
     ind = 0
     part_3d = 2*pi/Nphi
     do m = 1, Nphi
        eta1 = part_3d*(m-1)
        eta2 = part_3d*m
        eta = (eta2 - eta1)/2+eta1
        v(1) = cos(eta2)
        v(2) = sin(eta2)
        v(3) = 0
        D1 = reshape( (/ -1.0+2*v(1)**2, 2*v(2)*v(1), 0.0D0,  2*v(1)*v(2), -1.0+2*v(2)**2, 0.0D0, 0.0D0, 0.0D0, -1.0D0 /), (/3,3/) )
        v(1) = cos(eta)
        v(2) = sin(eta)
        v(3) = 0
        D2 = reshape( (/ v(1)**2, v(2)*v(1), -v(2),  v(1)*v(2), v(2)**2, v(1), v(2), -v(1), 0.0D0 /), (/3,3/) )
        D3 = reshape( (/ -1.0+2*v(1)**2, 2*v(2)*v(1), 0.0D0,  2*v(1)*v(2), -1.0+2*v(2)**2, 0.0D0, 0.0D0, 0.0D0, -1.0D0 /), (/3,3/) )
        D4 = matmul(D1,D3)
        do i = 1, Ns
           if (m==1) then
              v(1) = xp(i) + (r_torus+rmax)*cos(eta)
              v(2) = yp(i) + (r_torus+rmax)*sin(eta)
              v(3) = zp(i);
              xp(i) = dot_product(v,D2(1:3,1))
              yp(i) = dot_product(v,D2(1:3,2))
              zp(i) = dot_product(v,D2(1:3,3))
              v(1) = wxp(i);
              v(2) = wyp(i);
              v(3) = wzp(i);
              wxp(i) = dot_product(v,D2(1:3,1))
              wyp(i) = dot_product(v,D2(1:3,2))
              wzp(i) = dot_product(v,D2(1:3,3))
           else
              v(1) = xp(i)
              v(2) = yp(i)
              v(3) = zp(i)
              xp(i) = dot_product(v,D4(1:3,1))
              yp(i) = dot_product(v,D4(1:3,2))
              zp(i) = dot_product(v,D4(1:3,3))
              v(1) = wxp(i)
              v(2) = wyp(i)
              v(3) = wzp(i)
              wxp(i) = dot_product(v,D4(1:3,1))
              wyp(i) = dot_product(v,D4(1:3,2))
              wzp(i) = dot_product(v,D4(1:3,3))
           end if
           ind0 = ind0 + 1
           if (mod(ind0-1,n_cpu) == my_rank) then
              ind = ind + 1
              if (ind .gt. nppm-1) then
                 write(*,*) 'something is wrong here: to many particles in init',my_rank,ind,nppm,n
                 call MPI_ABORT(MPI_COMM_WORLD,ierr)
                 stop
              end if
              part(ind)%x(1) = xp(i) - (torus_offset(1)-(rmax-(1+12*nc**2)/(6*nc)*rl))/2.0
              part(ind)%x(2) = yp(i) - (torus_offset(2)-(rmax-(1+12*nc**2)/(6*nc)*rl))/2.0
              part(ind)%x(3) = zp(i) - (torus_offset(3)-(rmax-(1+12*nc**2)/(6*nc)*rl))/2.0
              part(ind)%alpha(1) = wxp(i)*volp(i)
              part(ind)%alpha(2) = wyp(i)*volp(i)
              part(ind)%alpha(3) = wzp(i)*volp(i)
           end if
        end do
     end do

     j = 0

     do k = 1, nc
        part_2d = 2*pi/(8*k)
        rc = (1+12*k**2)/(6*k)*rl

        do l = 1,8*k
           j = j+1
           xi1 = part_2d*(l-1)
           xi2 = part_2d*l
           xi = (xi2-xi1)/2+xi1
           xp(j) = rc*cos(xi)
           yp(j) = rc*sin(xi)
           zp(j) = 0
           volp(j) = (2*pi**2*(r_torus+(2*k+1)*rl)*((2*k+1)*rl)**2-2*pi**2*(r_torus+(2*k-1)*rl)*((2*k-1)*rl)**2)/(8*k*Nphi)
           wxp(j) = 0.
           wyp(j) = 0.
           wzp(j) = G*exp(-(rc/rmax)**2)
        end do
     end do

     xp(ns) = 0.
     yp(ns) = 0.
     zp(ns) = 0.
     wxp(ns) = 0.
     wyp(ns) = 0.
     wzp(ns) = g
     volp(Ns) = 2*pi**2*(r_torus+rl)*rl**2/Nphi

     j = 0
     ind0 = 0
     part_3d = 2*pi/Nphi
     do m = 1, Nphi
        eta1 = part_3d*(m-1)
        eta2 = part_3d*m
        eta = (eta2 - eta1)/2+eta1
        v(1) = cos(eta2)
        v(2) = sin(eta2)
        v(3) = 0
        D1 = reshape( (/ -1.0+2*v(1)**2, 2*v(2)*v(1), 0.0D0,  2*v(1)*v(2), -1.0+2*v(2)**2, 0.0D0, 0.0D0, 0.0D0, -1.0D0 /), (/3,3/) )
        v(1) = cos(eta)
        v(2) = sin(eta)
        v(3) = 0
        D2 = reshape( (/ v(1)**2, v(2)*v(1), -v(2),  v(1)*v(2), v(2)**2, v(1), v(2), -v(1), 0.0D0 /), (/3,3/) )
        D3 = reshape( (/ -1.0+2*v(1)**2, 2*v(2)*v(1), 0.0D0,  2*v(1)*v(2), -1.0+2*v(2)**2, 0.0D0, 0.0D0, 0.0D0, -1.0D0 /), (/3,3/) )
        D4 = matmul(D1,D3)
        do i = 1, Ns
           if (m==1) then
              v(1) = xp(i) + (r_torus+rmax)*cos(eta)
              v(2) = yp(i) + (r_torus+max)*sin(eta)
              v(3) = zp(i);
              xp(i) = dot_product(v,D2(1:3,1))
              yp(i) = dot_product(v,D2(1:3,2))
              zp(i) = dot_product(v,D2(1:3,3))
              v(1) = wxp(i);
              v(2) = wyp(i);
              v(3) = -wzp(i);
              wxp(i) = dot_product(v,D2(1:3,1))
              wyp(i) = dot_product(v,D2(1:3,2))
              wzp(i) = dot_product(v,D2(1:3,3))
           else
              v(1) = xp(i)
              v(2) = yp(i)
              v(3) = zp(i)
              xp(i) = dot_product(v,D4(1:3,1))
              yp(i) = dot_product(v,D4(1:3,2))
              zp(i) = dot_product(v,D4(1:3,3))
              v(1) = wxp(i)
              v(2) = wyp(i)
              v(3) = -wzp(i)
              wxp(i) = dot_product(v,D4(1:3,1))
              wyp(i) = dot_product(v,D4(1:3,2))
              wzp(i) = dot_product(v,D4(1:3,3))
           end if
           ind0 = ind0 + 1
           if (mod(ind0-1,n_cpu) == my_rank) then
              ind = ind + 1
              if (ind .gt. nppm-1) then
                 write(*,*) 'something is wrong here: to many particles in init',my_rank,ind,nppm,N_all
                 call MPI_ABORT(MPI_COMM_WORLD,ierr)
                 stop
              end if
              part(ind)%x(1) = xp(i) + (torus_offset(1)-(rmax-(1+12*nc**2)/(6*nc)*rl))/2.0
              part(ind)%x(2) = yp(i) + (torus_offset(2)-(rmax-(1+12*nc**2)/(6*nc)*rl))/2.0
              part(ind)%x(3) = zp(i) + (torus_offset(3)-(rmax-(1+12*nc**2)/(6*nc)*rl))/2.0
              part(ind)%alpha(1) = wxp(i)*volp(i)
              part(ind)%alpha(2) = wyp(i)*volp(i)
              part(ind)%alpha(3) = wzp(i)*volp(i)
           end if
        end do
     end do
     np = ind

     deallocate(xp,yp,zp,volp,wxp,wyp,wzp)

     if (dbo_main) write(*,*) 'Using spacings h,sigma,m_h,m_sigma:',h,4*h,m_h,G*m_h**0.95

  end select config

  call MPI_ALLREDUCE(np,n,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)


