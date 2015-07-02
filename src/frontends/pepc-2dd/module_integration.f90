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
!! Time integration module
!!!!!!!!!!!!!!!!!!!!

module tintegration

  use module_pepc
  use module_pepc_types
  use module_timings
  use module_debug
  use module_pepc_kinds
  use module_interaction_Specific_types
  use module_globals, only: dt,np
  use module_shortcut
  implicit none
  save
  private

    ! shortcut notations
!  real(kind_physics), parameter :: zero             =  0._kind_physics
!  real(kind_physics), parameter :: quarter          =  0.25_kind_physics
!  real(kind_physics), parameter :: half             =  0.5_kind_physics
!  real(kind_physics), parameter :: one              =  1._kind_physics
!  real(kind_physics), parameter :: two              =  2._kind_physics



  public midpoint
  public trapezoidal_rule
  public inv_m
  public update_step
  public leapfrog
  public boris
  public cross_product


  contains

  function cross_product(u,v)
    implicit none
    real(kind_particle), intent(in)     :: u(3),v(3)
    real(kind_particle)                 :: cross_product(3)

    cross_product(1) = u(2)*v(3) - u(3)*v(2)
    cross_product(2) = u(3)*v(1) - u(1)*v(3)
    cross_product(3) = u(1)*v(2) - u(2)*v(1)

  end function cross_product


  subroutine leapfrog(p)
    use module_mirror_boxes
    use module_shortcut
!    use module_globals
    implicit none

    type(t_particle), allocatable, intent(inout) :: p(:)
!    real(kind_particle)          , intent(in)    :: dt
    integer(kind_particle) :: ip
    real*8  :: fact

    fact = half/pi*dt

    do ip=1, np
      p(ip)%data%v = p(ip)%data%v + fact * p(ip)%data%q / p(ip)%data%m * p(ip)%results%e
      p(ip)%x      = p(ip)%x      + dt   * p(ip)%data%v
      if (p(ip)%x(1) .lt. zero) p(ip)%x(1) = four*pi
      if (p(ip)%x(1) .gt. four*pi) p(ip)%x(1) = zero
      p(ip)%data%v(2) = zero
      p(ip)%x(2) = zero
    end do

  end subroutine leapfrog


  subroutine boris(p)
    use module_mirror_boxes
!    use module_globals
    implicit none

    type(t_particle), allocatable, intent(inout) :: p(:)
!    real(kind_particle)          , intent(in)    :: dt
    integer(kind_particle) :: ip
    real(kind_particle)    :: fact,eta,u(3),uAh(3),h(3),s(3)


    do ip=1, np
      eta          = half*p(ip)%data%q*dt/p(ip)%data%m
      u            = p(ip)%data%v + eta*p(ip)%results%E
      h            = eta*p(ip)%results%B
      s            = two*h/(1.0_8 + dot_product(h,h) )
      uAh          = u + cross_product( u,h )
      uAh          = cross_product( uAh , s )

      p(ip)%data%v = u + uAh + eta* p(ip)%results%e
      p(ip)%x      = p(ip)%x      + dt   * p(ip)%data%v
    end do

  end subroutine boris


  subroutine update_step(np,x_tn1)
    implicit none

    integer(kind_particle)       , intent(in)  :: np
    type(t_particle)             , intent(inout)  :: x_tn1(np)

    integer(kind_particle)                     :: ip,jp
    real(kind_particle)                        :: d2,d(3),q,rd2,eps2

    eps2 = 1.0e-4


    do ip = 1,np

        x_tn1(ip)%results%pot = zero
        x_tn1(ip)%results%E   = zero

        do jp = 1,np
            if (jp.ne.ip) then
                q                     = x_tn1(ip)%data%q
                d(1:3)                = x_tn1(ip)%x - x_tn1(jp)%x
                d2                    = dot_product(d,d) + eps2
                rd2                   = one/d2

                x_tn1(ip)%results%pot = x_tn1(ip)%results%pot - q*log(d2)
                x_tn1(ip)%results%E   = x_tn1(ip)%results%E   + two*q*d*rd2
            endif
        enddo


    enddo

  end subroutine update_step

  subroutine midpoint(np,dt,x_tn,x_tn1,res)
    implicit none
!    include 'mpif.h'

    integer(kind_particle)       , intent(in) :: np
!    integer(kind_particle)                    :: n
    real(kind_particle)          , intent(in) :: dt
    type(t_particle)             , intent(in) :: x_tn(np)    ! Actual values of the particles
    real(kind_particle)          , intent(in) :: x_tn1(12*np) ! Next iteration of the particles
    real(kind_particle)          , intent(out):: res(12*np)   ! Residual

    type(t_particle), allocatable             :: r(:)
    integer(kind_particle)                    :: ip,jp,rc,ierr
    real(kind_particle)                       :: v(3),x(3),rot(3),p_tn1(3),p_tn(3),e,m,gamma_tn,gamma_tn1,fact

    !n  = 6*np

    if ( allocated(r) ) deallocate(r)
    allocate( r(np) , stat = rc )

    do ip = 1, np

        jp                = (ip-1)*12
        m                 = x_tn(ip)%data%m
        e                 = x_tn(ip)%data%q

        r(ip)%label       = x_tn(ip)%label
        r(ip)%data%q      = e
        r(ip)%data%m      = m

        x(1:3)            = x_tn1(jp+1:jp+3)
        v(1:3)            = x_tn1(jp+4:jp+6)/m
        gamma_tn          = one!sqrt( one + dot_product(v,v) )
        v(1:3)            = v(1:3)/gamma_tn

        r(ip)%x           = half*( x + x_tn(ip)%x )
        r(ip)%data%v      = half*( v + x_tn(ip)%data%v )

        r(ip)%results%e   = zero
        r(ip)%results%pot = zero
        r(ip)%results%B   = zero
        r(ip)%results%A   = zero
        r(ip)%results%Jirr= zero
        r(ip)%results%J   = zero
        r(ip)%work        = one

    enddo

!    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call pepc_particleresults_clear(r)
    call pepc_grow_tree(r)
    call pepc_traverse_tree(r)
    call pepc_timber_tree()
!    call update_step(np,r)

    fact = half/pi

    do ip = 1,np

        jp                = (ip-1)*12
        m                 = x_tn(ip)%data%m
        e                 = x_tn(ip)%data%q

        gamma_tn          = one!one/sqrt( one - dot_product( x_tn(ip)%data%v, x_tn(ip)%data%v ) )

        p_tn(1:3)         = gamma_tn*m*x_tn(ip)%data%v(1:3)

        p_tn1(1:3)        = x_tn1(jp+4:jp+6)/m!( x_tn1(jp+4:jp+6) - e*x_tn1(jp+7:jp+9) )/m
        gamma_tn1         = one!sqrt( one + dot_product( p_tn1, p_tn1 ) )

        p_tn1(1:3)        = p_tn1(1:3)/gamma_tn1  ! this is the velocity v

        !!! v cross B
        rot(1)            = p_tn1(2)*r(ip)%results%B(3) - p_tn1(3)*r(ip)%results%B(2)
        rot(2)            = p_tn1(3)*r(ip)%results%B(1) - p_tn1(1)*r(ip)%results%B(3)
        rot(3)            = p_tn1(1)*r(ip)%results%B(2) - p_tn1(2)*r(ip)%results%B(1)

        res(jp+4:jp+6)    = x_tn1(jp+4:jp+6)   - p_tn(1:3)          - fact*e*dt*(  r(ip)%results%E(1:3)  )
        res(jp+1:jp+3)    = x_tn1(jp+1:jp+3)   - x_tn(ip)%x(1:3)    - dt/m*( half*( p_tn(1:3) + x_tn1(jp+4:jp+6 ) ) )
        res(jp+7:jp+9)    = x_tn1(jp+7:jp+9)   - r(ip)%results%A
        res(jp+10:jp+12)  = x_tn1(jp+10:jp+12) - r(ip)%results%B

    enddo

    deallocate(r)

  end subroutine midpoint

  subroutine trapezoidal_rule(np,dt,x_tn,x_tn1,res)
    implicit none

    integer(kind_particle)       , intent(in) :: np
!    integer(kind_particle)                    :: n
    real(kind_particle)          , intent(in) :: dt
    type(t_particle)             , intent(in) :: x_tn(np)    ! Actual values of the particles
    real(kind_particle)          , intent(in) :: x_tn1(12*np) ! Next iteration of the particles
    real(kind_particle)          , intent(out):: res(12*np)   ! Residual

    type(t_particle), allocatable             :: r(:)
    integer(kind_particle)                    :: ip,jp,rc
    real(kind_particle)                       :: v(3),x(3),rot(3),p_tn1(3),p_tn(3),e,m,gamma,gamma_tn1

    !n  = 6*np

    if ( allocated(r) ) deallocate(r)
    allocate( r(np) , stat = rc )

    do ip = 1, np

        jp                    = (ip-1)*12

        m                     = x_tn(ip)%data%m
        e                     = x_tn(ip)%data%q

        r(ip)%label           = x_tn(ip)%label
        r(ip)%data%q          = e
        r(ip)%data%m          = m
        r(ip)%x(1:3)          = x_tn1(jp+1:jp+3)
        r(ip)%data%v(1:3)     = x_tn1(jp+4:jp+6)/m


        r(ip)%results%e   = zero
        r(ip)%results%pot = zero
        r(ip)%results%B   = zero
        r(ip)%results%A   = zero
        r(ip)%results%Jirr= zero
        r(ip)%results%J   = zero
        r(ip)%work        = one

    enddo

    call pepc_particleresults_clear(r)
    call pepc_grow_tree(r)
    call pepc_traverse_tree(r)
    call pepc_timber_tree()

    do ip = 1,np

        jp                    = (ip-1)*12

        m                     = x_tn(ip)%data%m
        e                     = x_tn(ip)%data%q
        p_tn(1:3)             = m*x_tn(ip)%data%v(1:3)

        res(jp+4:jp+6)        = x_tn1(jp+4:jp+6) - p_tn(1:3) - e*dt*( half*( x_tn(ip)%results%E(1:3) + r(ip)%results%E(1:3) ) )
        res(jp+1:jp+3)        = x_tn1(jp+1:jp+3) - x_tn(ip)%x(1:3) - dt/m*( half*( p_tn(1:3) + x_tn1(jp+4:jp+6) ) )

    enddo

    deallocate(r)

  end subroutine trapezoidal_rule


  subroutine inv_m(np,dt,x_tn,x_tn1,res)
    implicit none

    integer(kind_particle)       , intent(in) :: np
!    integer(kind_particle)                    :: n
    real(kind_particle)          , intent(in) :: dt
    type(t_particle)             , intent(in) :: x_tn(np)    ! Actual values of the particles
    real(kind_particle)          , intent(in) :: x_tn1(12*np) ! Next iteration of the particles
    real(kind_particle)          , intent(out):: res(12*np)   ! Residual

    type(t_particle), allocatable             :: r(:)
    integer(kind_particle)                    :: ip,jp,rc
    real(kind_particle)                       :: v(3),x(3),rot(3),p_tn1(3),p_tn(3),e,m,gamma,gamma_tn1

    !n  = 6*np

    if ( allocated(r) ) deallocate(r)
    allocate( r(np) , stat = rc )

    do ip = 1, np

        jp             = (ip-1)*12
        res(jp+1)      = exp(-x_tn1(jp+1)) - x_tn1(jp+2)
        res(jp+2)      = (x_tn1(jp+2)) - (x_tn1(jp+1))
        res(jp+3)      = exp(-x_tn1(jp+3)) - (x_tn1(jp+3))

        res(jp+4)      = exp(-x_tn1(jp+4)) - (x_tn1(jp+5))
        res(jp+5)      = (x_tn1(jp+5)) - (x_tn1(jp+4))
        res(jp+6)      = exp(-x_tn1(jp+6)) - (x_tn1(jp+6))

    enddo

  end subroutine inv_m

  end module
