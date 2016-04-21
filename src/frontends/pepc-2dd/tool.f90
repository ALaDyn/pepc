! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
!
! Copyright (C) 2002-2014 Juelich Supercomputing Centre,
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

!!!!!!!!!!!!!!!!!!!!
!! module tool
!!!!!!!!!!!!!!!!!!!!



module module_tool

  use module_pepc_kinds
  use module_pepc_types
  use module_interaction_Specific_types
  use module_globals
  use module_shortcut
  implicit none

  contains

  subroutine random(array)
    implicit none
    real*8 :: array(:)
    integer :: i

    do i = 1,size(array)
       array(i) = par_rand()
    end do

  end subroutine random


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>
  !> portable random number generator, see numerical recipes
  !> check for the random numbers:
  !> the first numbers should be 0.2853809, 0.2533582 and 0.0934685
  !> the parameter iseed is optional
  !>
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function par_rand(iseed)
    implicit none
    real :: par_rand
    integer, intent(in), optional :: iseed

    integer, parameter :: IM1  = 2147483563
    integer, parameter :: IM2  = 2147483399
    real,    parameter :: AM   = 1.0/IM1
    integer, parameter :: IMM1 = IM1-1
    integer, parameter :: IA1  = 40014
    integer, parameter :: IA2  = 40692
    integer, parameter :: IQ1  = 53668
    integer, parameter :: IQ2  = 52774
    integer, parameter :: IR1  = 12211
    integer, parameter :: IR2  = 3791
    integer, parameter :: NTAB = 32
    integer, parameter :: NDIV = 1+IMM1/NTAB
    real,    parameter :: eps_ = 1.2e-7 ! epsilon(eps_)
    real,    parameter :: RNMX = 1.0 - eps_

    integer :: j, k
    integer, volatile, save :: idum  = -1
    integer, volatile, save :: idum2 =  123456789
    integer, volatile, save :: iy    =  0
    integer, volatile, save :: iv(NTAB)


    if (idum <=0 .or. present(iseed)) then
       if (present(iseed)) then
          idum = iseed
       else
          if (-idum < 1) then
             idum = 1
          else
             idum = -idum
          endif
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
    par_rand = AM*iy
    if (par_rand > RNMX) par_rand = RNMX

  end function par_rand


  subroutine random_gauss(list)
    implicit none
    real*8, intent(inout) :: list(:)
    real*8  :: v(2), r, p
    integer :: n, i

!    pi = 2.0_8*acos(0.0_8)
    n  = size(list)

    do i=1, n, 2

       call random(v)

       r = sqrt(-2.0_8 * log(v(1)))
       p = 2.0_8*pi*v(2)

       list(i)                = r * sin(p)
       if((i+1)<=n) list(i+1) = r * cos(p)

    end do

  end subroutine


  subroutine random_boltzmann(list,vth,vdrift)
    implicit none
    real(kind_particle), intent(out)   :: list(:)
    real(kind_particle), intent(in)    :: vth(3),vdrift(3)
    real(kind_particle)                :: v(2), r, theta
    integer(kind_particle)             :: i,n

    n = size(list,kind= kind_particle)

    call random(v)
    r         = sqrt(-two*log( one - 0.9999_8*v(1) ) )
    theta     = two*pi*v(2)
    list(1)   = vdrift(1) + vth(1)*r*cos(theta)
    if (n .ge. 2) list(2)   = vdrift(2) + vth(2)*r*sin(theta)
    
    if (n .eq. 3) then 
        call random(v)
        r         = sqrt(-two*log( one - 0.9999_8*v(1) ) )
        theta     = two*pi*v(2)
        list(3)   = vdrift(3) + vth(3)*r*cos(theta)
    endif


  end subroutine



  function cross_product(x,y)
    implicit none
    real(kind_particle)            :: cross_product(3)
    real(kind_particle),intent(in) :: x(3),y(3)


    cross_product(1) = x(2)*y(3) - x(3)*y(2)
    cross_product(2) = x(3)*y(1) - x(1)*y(3)
    cross_product(3) = x(1)*y(2) - x(2)*y(1)


  end function cross_product
  
    function double_cross_product_left(x,y,z)
    implicit none
    real(kind_particle)            :: double_cross_product_left(3)
    real(kind_particle),intent(in) :: x(3),y(3),z(3)


    double_cross_product_left(1:3) = dot_product( x(1:3),z(1:3) )*y(1:3) - dot_product( x(1:3),y(1:3) )*z(1:3)


  end function double_cross_product_left

  

!subroutine change_unit(p,q,flagin,flagout)
!    use module_pepc
!    use module_shortcut,only: pi
!    use module_globals,only: qe,c,me
!    implicit none
!
!    type(t_particle), allocatable, intent(in)    :: p(:)
!    type(t_particle), allocatable, intent(out)   :: q(:)
!    integer(kind_particle)       , intent(in)    :: flagin,flagout
!
!    integer(kind_particle)                       :: ip
!    real(kind_particle)                          :: etilde,phitilde,atilde,btilde,mtilde,qtilde,xtilde,&
!                                                    ttilde,vtilde,wp,ctilde
!
!
!
!    if ( flagin .ne. 0 .or. flagout .ne. 1) then
!        write(*,*) " Wrong choice of flagin in conversion of unit"
!        call exit(1)
!    end if
!    if ( flagout .ne. 0 .or. flagout .ne. 1) then
!        write(*,*) " Wrong choice of flagout in conversion of unit"
!        call exit(1)
!    end if
!
!    select case (flagin)
!          case (0)  !
!
!            etilde          =  one/Volume
!            phitilde        =  one/Volume
!            btilde          =  one/Volume
!            atilde          =  one/Volume
!            jtilde          =  one/Volume
!            rhotilde        =  one/Volume
!
!          case (1)  !
!            etilde          =  one
!            phitilde        =  one
!            btilde          =  one
!            atilde          =  one
!            jtilde          =  one
!            rhotilde        =  one
!
!          case default
!
!
!    end select
!
!end subroutine change_unit

subroutine copy_particle(p,q,np)
    implicit none

    type(t_particle), allocatable, intent(in)    :: p(:)
    type(t_particle), allocatable, intent(out)   :: q(:)
    integer(kind_particle)       , intent(in)    :: np

    integer(kind_particle)                       :: ip,rc

    if (allocated(q)) deallocate(q)
    allocate(q(np),stat=rc )
    do ip = 1,np

        q(ip)%label             = p(ip)%label
        q(ip)%x(1:3)            = p(ip)%x(1:3)
        q(ip)%data%v(1:3)       = p(ip)%data%v(1:3)
        q(ip)%data%q            = p(ip)%data%q
        q(ip)%data%m            = p(ip)%data%m
        q(ip)%results%E(1:3)    = p(ip)%results%E(1:3)
        q(ip)%results%B(1:3)    = p(ip)%results%B(1:3)
        q(ip)%results%A(1:3)    = p(ip)%results%A(1:3)
        q(ip)%results%J(1:3)    = p(ip)%results%J(1:3)
        q(ip)%results%Jirr(1:3) = p(ip)%results%Jirr(1:3)
        q(ip)%results%pot       = p(ip)%results%pot

    enddo

end subroutine copy_particle

function gyrofrequency(p)
    implicit none
    include 'mpif.h'
    type(t_particle), allocatable, intent(in)    :: p(:)
    real(kind_particle)                          :: gyrofrequency,Bnorm,Bnormg
    integer(kind_particle)                       :: ip,rc

    Bnorm  = zero
    Bnormg = zero

    do ip = 1,np
        Bnorm = Bnorm + dot_product(p(ip)%results%B,p(ip)%results%B)
    enddo

    call MPI_ALLREDUCE(Bnorm, Bnormg, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
    Bnormg  = sqrt(Bnormg)/real(tnp,kind=kind_particle)

    gyrofrequency = abs(p(1)%data%q)/p(1)%data%m*Bnormg
    do ip = 1,np
        gyrofrequency = max( gyrofrequency,p(ip)%data%q/p(ip)%data%m*Bnormg )
    enddo

end function gyrofrequency

function inv2x2(p)
    use module_shortcut, only: prec
    implicit none
    real(kind_particle), intent(in)              :: p(1:4)
    real(kind_particle)                          :: inv2x2(1:4),det

    inv2x2    = -p
    inv2x2(1) =  p(4)
    inv2x2(4) =  p(4)
    det      =  ( p(1)*p(4) - p(2)*p(3) )
    if (det .le. prec) write(*,*) " Warning, attempt to invert a singular matrix "
    inv2x2    =  inv2x2/det

end function inv2x2

function inv3x3(p)
    use module_shortcut, only: prec
    implicit none
    real(kind_particle), intent(in)              :: p(1:9)
    real(kind_particle)                          :: inv3x3(1:9),det

    
    det      =  p(1)*( p(5)*p(9) - p(6)*p(8) ) - p(2)*( p(4)*p(9) - p(6)*p(7) ) + p(3)*( p(4)*p(8) - p(5)*p(7) )
    if (det .le. prec) write(*,*) " Warning, attempt to invert a singular matrix "
    
    inv3x3(1)= p(5)*p(9) - p(6)*p(8)
    inv3x3(2)= p(3)*p(8) - p(2)*p(9)
    inv3x3(3)= p(2)*p(6) - p(3)*p(5)
    inv3x3(4)= p(6)*p(7) - p(4)*p(9)
    inv3x3(5)= p(1)*p(9) - p(3)*p(7)
    inv3x3(6)= p(3)*p(4) - p(1)*p(6)
    inv3x3(7)= p(4)*p(8) - p(5)*p(7)
    inv3x3(8)= p(2)*p(7) - p(1)*p(8)
    inv3x3(9)= p(1)*p(3) - p(2)*p(4)
    
    inv3x3   = inv3x3/det

end function inv3x3


function integrate(field)
    use encap
    implicit none
    include 'mpif.h'
    type(field_grid_t), intent(in)               :: field
    real(kind_particle)                          :: integrate_loc,integrate,dx,dy,estimated_loc,estimated
    integer(kind_particle)                       :: ip,jp, nx,ny,np,rc

    nx             = field%n(1)
    ny             = field%n(2)
    dx             = field%dx(1)
    dy             = field%dx(2)

    integrate      = zero
    integrate_loc  = zero
    estimated_loc  = zero
    estimated      = zero

    do jp = 2,ny-1
        do ip = 1,nx-1
            integrate_loc = integrate_loc + quarter*dx*dy*( field%p(ip)%results%J(3) + field%p(ip+1)%results%J(3) + field%p((jp-1)*nx)%results%J(3) + field%p(jp*nx)%results%J(3) )
!            it = it +1
!            write(*,*) it, field%p(ip)%x(1),field%p(ip+1)%x(1),field%p((jp-1)*nx)%x(2),field%p(jp*nx)%x(2)
        enddo
!        integrate_loc = integrate_loc !+ dot_product(p(ip)%results%B,p(ip)%results%B)
    enddo


    do ip = 1,field%nl
        estimated_loc = estimated_loc + field%p(ip)%results%J(3)
    enddo

    call MPI_ALLREDUCE(integrate_loc, integrate, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
    call MPI_ALLREDUCE(estimated_loc, estimated, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)

    if (root)    write(*,*) "current I ===",integrate, estimated/field%nl*(0.2)**2*pi



end function integrate

end module module_tool
