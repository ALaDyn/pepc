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
  use module_fmm_framework
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
  use module_fmm_framework
  use module_units
  use module_velocity_setup
  implicit none
  include 'mpif.h'

  integer, intent(in) :: iconf  ! Configuration switch
  integer :: p,mpi_cnt, ierr
  real*4 :: par_rand_res
  real*8 :: yt,zt,xt,r1
  integer :: np_local_max
  real*8 :: delta(3)
  integer :: i,j,k,n(3), myidx, globalidx

  config: select case(iconf)
  case(1)

     if (my_rank == 0) write(*,*) "Using special start... case 1 (homogeneous distribution)"

     ! get the largest np_local
     call MPI_REDUCE(np_local, np_local_max, 1, MPI_INTEGER, MPI_MAX, 0, MPI_COMM_PEPC, ierr)
     call MPI_BCAST(np_local_max, 1, MPI_INTEGER, 0, MPI_COMM_PEPC, ierr)

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
     call MPI_REDUCE(np_local, np_local_max, 1, MPI_INTEGER, MPI_MAX, 0, MPI_COMM_PEPC, ierr)
     call MPI_BCAST(np_local_max, 1, MPI_INTEGER, 0, MPI_COMM_PEPC, ierr)

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
     call MPI_REDUCE(np_local, np_local_max, 1, MPI_INTEGER, MPI_MAX, 0, MPI_COMM_PEPC, ierr)
     call MPI_BCAST(np_local_max, 1, MPI_INTEGER, 0, MPI_COMM_PEPC, ierr)

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
     call MPI_REDUCE(np_local, np_local_max, 1, MPI_INTEGER, MPI_MAX, 0, MPI_COMM_PEPC, ierr)
     call MPI_BCAST(np_local_max, 1, MPI_INTEGER, 0, MPI_COMM_PEPC, ierr)

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

!!$  case(5)
!!$       
!!$     if (my_rank == 0) write(*,*) "Using special start... case 5: Sphere"
!!$
!!$     do i = 1,npart_total
!!$
!!$        if(i.eq.1) then
!!$           cth = mone
!!$           sth = zero
!!$           cphi = one
!!$           sphi = zero
!!$           b = cphi
!!$           c = sphi
!!$        elseif(i.eq.npart_total) then
!!$           cth = one
!!$           sth = zero
!!$           cphi = one
!!$           sphi = zero
!!$        else
!!$           cth = dble(2*i-npart_total-1)/dble(npart_total-1)
!!$           sth = two*(sqrt(dble(i-1)/dble(npart_total-1))*sqrt(dble(npart_total-i)/dble(npart_total-1)))
!!$           s = a*sqrt(dble(npart_total-1)/(dble(i-1)*dble(npart_totaly-i)))
!!$           cphi = b*cos(s)-c*sin(s)
!!$           sphi = c*cos(s)+b*sin(s)
!!$           b = cphi
!!$           c = sphi
!!$        endif
!!$        if (mod(i+my_rank,n_cpu) == 0) then
!!$           j = j+1
!!$           x(i) = sth*cphi 
!!$           y(i) = sth*sphi 
!!$           z(i) = cth
!!$           ux(i) = 0.
!!$           uy(i) = 0.
!!$           uz(i) = 0.
!!$           if (
!!$        end if
!!$     end do

  case(6)
     
     if (my_rank == 0) write(*,*) "Using special start... case 6 (2D-homogeneous distribution)"

     ! get the largest np_local
     call MPI_REDUCE(np_local, np_local_max, 1, MPI_INTEGER, MPI_MAX, 0, MPI_COMM_PEPC, ierr)
     call MPI_BCAST(np_local_max, 1, MPI_INTEGER, 0, MPI_COMM_PEPC, ierr)
     
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

  case(7)

     n(1) = NINT((npart_total/8)**(1./3.))*2
     n(2) = NINT((npart_total/n(1)/4)**(1./2.))*2
     n(3) = npart_total/n(1)/n(2)

     if (my_rank == 0) write(*,*) "Using special start... case 7 (3D Madelung Setup) with", n, "particles per edge"

     delta(1) = sqrt(dot_product(t_lattice_1,t_lattice_1))/real(n(1))
     delta(2) = sqrt(dot_product(t_lattice_2,t_lattice_2))/real(n(2))
     delta(3) = sqrt(dot_product(t_lattice_3,t_lattice_3))/real(n(3))

     myidx     = 0
     globalidx = 0

     do i = 0, n(1)-1
       do j = 0, n(2)-1
         do k = 0, n(3)-1

           globalidx = globalidx + 1

           if ( mod((globalidx-1)/2,n_cpu) == my_rank) then ! distribute pairs of electron and ion, since np_local is constructed a bit weird
             myidx = myidx + 1

             x(myidx)  = (i + 0.5)*delta(1)
             y(myidx)  = (j + 0.5)*delta(2)
             z(myidx)  = (k + 0.5)*delta(3)
             ux(myidx) = 0
             uy(myidx) = 0
             uz(myidx) = 0

             if (mod(i+j+k, 2) == 0) then
               q(myidx) = qe
               m(myidx) = mass_e
             else
               q(myidx) = qi
               m(myidx) = mass_i
             end if

           end if

         end do
       end do
     end do

     pelabel(1:nep)            = my_rank * nep + (/(i, i = 1, nep)/)      ! Electron labels
     pelabel(nep + 1:np_local) = ne + my_rank * nip + (/(i, i = 1, nip)/) ! Ion labels
     ex(1:np_local) = 0.
     ey(1:np_local) = 0.
     ez(1:np_local) = 0.
     pot(1:np_local) = 0.

     work(1:np_local) = 1.

     if (myidx .ne. np_local) write(*,*) "ERROR in special_start(7): PE", my_rank, "set up", myidx, &
           "particles, but np_local=", np_local, "globalidx=", globalidx, "npart_total=",npart_total

     !write(*,*) "Particle positions: "
     !do i=1,np_local
     !  write(*,*) my_rank,i,x(i), y(i), z(i), q(i), pelabel(i)
     !end do

     return ! do not execute the stuff below in this case, since this messes up charges etc.

   case(8)

     if (my_rank == 0) write(*,*) "Using special start... case 8 (homogeneous electron distribution, ionic lattice)"

     n(1) = NINT(ni**(1./3.))
     delta(1) = 1./real(n(1))
     delta(2) = 1./real(n(1))
     delta(3) = 1./real(n(1))

     myidx     = 0
     globalidx = 0

     do i = 0, n(1)-1
       do j = 0, n(1)-1
         do k = 0, n(1)-1

           globalidx = globalidx + 1

           if ( mod((globalidx-1)/2,n_cpu) == my_rank) then ! distribute pairs of electron and ion, since np_local is constructed a bit weird
             myidx = myidx + 1

             x(np_local-myidx+1)  = (i + 0.5)*delta(1)
             y(np_local-myidx+1)  = (j + 0.5)*delta(2)
             z(np_local-myidx+1)  = (k + 0.5)*delta(3)
             ux(np_local-myidx+1) = 0
             uy(np_local-myidx+1) = 0
             uz(np_local-myidx+1) = 0

             q(np_local-myidx+1) = qi
             m(np_local-myidx+1) = mass_i
           end if

         end do
       end do
     end do

     do i = 1,ne
       globalidx = globalidx + 1

       if ( mod((globalidx-1)/2,n_cpu) == my_rank) then ! distribute pairs of electron and ion, since np_local is constructed a bit weird
         myidx = myidx + 1

         call par_rand(par_rand_res)
         x(np_local-myidx+1)  = par_rand_res
         call par_rand(par_rand_res)
         y(np_local-myidx+1)  = par_rand_res
         call par_rand(par_rand_res)
         z(np_local-myidx+1)  = par_rand_res
      
         ux(np_local-myidx+1) = 0
         uy(np_local-myidx+1) = 0
         uz(np_local-myidx+1) = 0

         q(np_local-myidx+1) = qe
         m(np_local-myidx+1) = mass_e
       end if
     end do

     if (myidx .ne. np_local) write(*,*) "ERROR in special_start(8): PE", my_rank, "set up", myidx, &
           "particles, but np_local=", np_local, "globalidx=", globalidx, "npart_total=",npart_total
     if (globalidx .ne. npart_total) write(*,*) "ERROR in special_start(8): PE", my_rank, "set up globalidx=", globalidx, ", but npart_total=",npart_total

  end select config

  q(1:nep)                  = qe        ! plasma electrons
  q(nep + 1:np_local)       = qi        ! plasma ions (need Z* here)
  m(1:nep)                  = mass_e    ! electron mass
  m(nep + 1:np_local)       = mass_i    ! ion mass
  pelabel(1:nep)            = my_rank * nep + (/(i, i = 1, nep)/)      ! Electron labels
  pelabel(nep + 1:np_local) = ne + my_rank * nip + (/(i, i = 1, nip)/) ! Ion labels

  if (my_rank==0) write(*,*) 'Initializing particle velocities to vte =',vte,' vti =',vti
     if (vte > 0) then
        call maxwell3(ux,uy,uz,nppm,1,nep,vte)
     else
        call cold_start(ux,uy,uz,nppm,1,nep)
     endif

     if (vti > 0) then
        call maxwell3(ux,uy,uz,nppm,nep+1,nip,vti)
     else
        call cold_start(ux,uy,uz,nppm,nep+1,nip)
     endif

  ex(1:np_local) = 0.
  ey(1:np_local) = 0.
  ez(1:np_local) = 0.
  pot(1:np_local) = 0.

  work(1:np_local) = 1.


end subroutine special_start
