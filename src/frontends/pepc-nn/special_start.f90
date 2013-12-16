! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2013 Juelich Supercomputing Centre, 
!                         Forschungszentrum Juelich GmbH,
!                         Germany
! 
! PEPC is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! PEPC is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
! 
! You should have received a copy of the GNU Lesser General Public License
! along with PEPC.  If not, see <http://www.gnu.org/licenses/>.
!

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
  use module_mirror_boxes
  implicit none
  include 'mpif.h'

  interface
    subroutine par_rand(res, iseed)
      real, intent(inout) :: res
      integer, intent(in), optional :: iseed
    end subroutine
  end interface

  interface
    function GetSphereCenter(idx)
      integer, intent(in) :: idx
      real*8, dimension(3) :: GetSphereCenter
    end function
  end interface


  integer, intent(in) :: iconf  ! Configuration switch
  integer :: p,mpi_cnt, ierr
  real*4 :: par_rand_res
  real*8 :: yt,zt,xt, r1
  integer :: np_local_max
  real*8 :: delta(3)
  integer :: i,j,k,n(3), myidx, globalidx

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
  call MPI_ALLGATHER(MPI_IN_PLACE, 1, MPI_INTEGER, fances(0), 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
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
              
              particles(p)%x = [xt, yt, zt]
              
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
              
              particles(p)%x = [xt, yt, zt] + 0.5
              
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
                            
              particles(p)%x = [xt, yt, zt] + 0.5
              
           end if
        end do
     end do

  case(342)
     if (my_rank == 0) write(*,*) "Using special start... case 342 (42 spheres benchmark)"

     r_sphere = 0.05

     ! initialize random number generator with some arbitrary seed
     call par_rand(par_rand_res, my_rank + 13)

        do p = 1, (fances(my_rank) - fances(my_rank-1))

           xt = 2.0_8
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

           call par_rand(par_rand_res)

           delta = 1._8*GetSphereCenter(nint(42.*par_rand_res))

           particles(p)%x = [xt, yt, zt] * r_sphere + delta

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

              particles(p)%x = [xt, yt, zt]

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
              
              particles(p)%x = [xt, yt, 0._8] + 0.5
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
              
              particles(p)%x = [xt, yt, zt]
              
           end if
        end do
     end do


  case(7)

    n(1) = NINT((npart_total/8)**(1./3.))*2
    n(2) = NINT((npart_total/n(1)/4)**(1./2.))*2
    n(3) = npart_total/n(1)/n(2)

    if (my_rank == 0) write(*,*) "Using special start... case 7 (3D Madelung Setup) with", n, "particles per edge"

    qe = -1.0
    qi = -qe

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

            particles(myidx)%x = ([i, j, k] + 0.5) * delta

          end if

        end do
      end do
    end do

    particles(1:nep)%label            = my_rank * nep + (/(i, i = 1, nep)/)      ! Electron labels
    particles(nep + 1:np_local)%label = ne + my_rank * nip + (/(i, i = 1, nip)/) ! Ion labels

    if (myidx .ne. np_local) write(*,*) "ERROR in special_start(7): PE", my_rank, "set up", myidx, &
        "particles, but np_local=", np_local, "globalidx=", globalidx, "npart_total=",npart_total

    !write(*,*) "Particle positions: "
    !do i=1,np_local
    !  write(*,*) my_rank,i,x(i), y(i), z(i), q(i), pelabel(i)
    !end do

    return ! do not execute the stuff below in this case, since this messes up charges etc.

  case(8)

     if (my_rank == 0) write(*,*) "Using special start... case 8 (fast homogeneous distribution)"

     ! initialize random number generator with some arbitrary seed
     call par_rand(par_rand_res, my_rank + 13)

     do p = 1, (fances(my_rank) - fances(my_rank-1))
           
        call par_rand(par_rand_res)
        xt = par_rand_res
        call par_rand(par_rand_res)
        yt = par_rand_res
        call par_rand(par_rand_res)
        zt = par_rand_res

        particles(p)%x = [xt, yt, zt]
              
     end do

  end select config

  particles(1:nep)%data%q             = qe        ! plasma electrons
  particles(nep + 1:np_local)%data%q  = qi        ! plasma ions (need Z* here)
  particles(1:nep)%label            = fances(my_rank-1) + (/(i, i = 1, nep)/)      ! Electron labels
  particles(nep + 1:np_local)%label = ne + fances(my_rank-1) + nep + (/(i, i = 1, (fances(my_rank) - fances(my_rank-1) - nep))/) ! Ion labels

end subroutine special_start



function GetSphereCenter(idx)
  implicit none
  integer, intent(in) :: idx
  real*8, dimension(3) :: GetSphereCenter

  INTEGER, DIMENSION(3, 3) :: array = reshape([ 1, 2, 3, 4,&
   5, 6, 7, 8, 9 ], shape(array))


  real*8, dimension(3, 50) :: random_vectors = reshape([  0.7011 ,  0.0942 ,  0.0012, &
    0.6663 ,  0.5985 ,  0.4624,&
    0.5391 ,  0.4709 ,  0.4243,&
    0.6981 ,  0.6959 ,  0.4609,&
    0.6665 ,  0.6999 ,  0.7702,&
    0.1781 ,  0.6385 ,  0.3225,&
    0.1280 ,  0.0336 ,  0.7847,&
    0.9991 ,  0.0688 ,  0.4714,&
    0.1711 ,  0.3196 ,  0.0358,&
    0.0326 ,  0.5309 ,  0.1759,&
    0.5612 ,  0.6544 ,  0.7218,&
    0.8819 ,  0.4076 ,  0.4735,&
    0.6692 ,  0.8200 ,  0.1527,&
    0.1904 ,  0.7184 ,  0.3411,&
    0.3689 ,  0.9686 ,  0.6074,&
    0.4607 ,  0.5313 ,  0.1917,&
    0.9816 ,  0.3251 ,  0.7384,&
    0.1564 ,  0.1056 ,  0.2428,&
    0.8555 ,  0.6110 ,  0.9174,&
    0.6448 ,  0.7788 ,  0.2691,&
    0.3763 ,  0.4235 ,  0.7655,&
    0.1909 ,  0.0908 ,  0.1887,&
    0.4283 ,  0.2665 ,  0.2875,&
    0.4820 ,  0.1537 ,  0.0911,&
    0.1206 ,  0.2810 ,  0.5762,&
    0.5895 ,  0.4401 ,  0.6834,&
    0.2262 ,  0.5271 ,  0.5466,&
    0.3846 ,  0.4574 ,  0.4257,&
    0.5830 ,  0.8754 ,  0.6444,&
    0.2518 ,  0.5181 ,  0.6476,&
    0.2904 ,  0.9436 ,  0.6790,&
    0.6171 ,  0.6377 ,  0.6358,&
    0.2653 ,  0.9577 ,  0.9452,&
    0.8244 ,  0.2407 ,  0.2089,&
    0.9827 ,  0.6761 ,  0.7093,&
    0.7302 ,  0.2891 ,  0.2362,&
    0.3439 ,  0.6718 ,  0.1194,&
    0.5841 ,  0.6951 ,  0.6073,&
    0.1078 ,  0.0680 ,  0.4501,&
    0.9063 ,  0.2548 ,  0.4587,&
    0.8797 ,  0.2240 ,  0.6619,&
    0.8178 ,  0.6678 ,  0.7703,&
    0.2607 ,  0.8444 ,  0.3502,&
    0.5944 ,  0.3445 ,  0.6620,&
    0.0225 ,  0.7805 ,  0.4162,&
    0.4253 ,  0.6753 ,  0.8419,&
    0.3127 ,  0.0067 ,  0.8329,&
    0.1615 ,  0.6022 ,  0.2564,&
    0.1788 ,  0.3868 ,  0.6135,&
    0.4229 ,  0.9160 ,  0.5822], shape(random_vectors))



  GetSphereCenter = 1._8*random_vectors(:, idx + 1)

end function GetSphereCenter

