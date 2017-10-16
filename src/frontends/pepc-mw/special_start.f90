! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2017 Juelich Supercomputing Centre, 
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






subroutine special_start(iconf)

  use physvars
  use module_mirror_boxes
  use module_icosahedron
  use module_diagnostics
  use module_units
  use module_velocity_setup
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
  real*4 :: par_rand_res,mu
  real*8 :: yt,zt,xt, r1,fr,r(3)
  real*4 ::dx
  integer :: currlayer, particletype
  real*8 a,b,c,cth,sth,cphi,sphi,s

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


  config: select case(iconf)



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

     q(1:np_local)       = qi       ! ion charge
     m(1:np_local)       = mass_i   ! ion mass
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

     q(1:j)       = qi       ! ion charge
     m(1:j)       = mass_i   ! ion mass

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

     q(1:np_local)       = qi       ! ion charge
     m(1:np_local)       = mass_i   ! ion mass
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


end subroutine special_start




