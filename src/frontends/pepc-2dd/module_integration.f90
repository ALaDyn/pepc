! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
!
! Copyright (C) 2002-2016 Juelich Supercomputing Centre,
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

module module_integration

  use module_pepc
  use module_pepc_types
  use module_timings
  use module_debug
  use module_pepc_kinds
  use module_interaction_Specific_types
  use module_globals, only: dt,ixdim,ivdim,vtilde,lorentz_tilde,folder
  use module_shortcut
  implicit none
  save
  private

  public midpoint3D
  public trapezoidal3D
  public midpoint_electrostatic
  public hamiltonian_boris
  public leapfrog
  public midpoint
  public trapezoidal
  public implicit_nkPtr
  public implicit_picardPtr
  public explicit_ptr
  public march
!  public boris


  type pv_particle
    real(kind_physics), dimension(3) :: x, v
    real(kind_physics)               :: gamma
  end type pv_particle

  type pma_particle
    real(kind_physics), dimension(3) :: x, P, A
    real(kind_physics)               :: gamma
  end type pma_particle
  
   abstract interface

      subroutine explicit_ptr(np,dt,p)
      use module_pepc_types
      use module_pepc_kinds
      implicit none
        integer(kind_particle)       , intent(in)    :: np
        real(kind_particle)          , intent(in)    :: dt
        type(t_particle),allocatable , intent(inout) :: p(:)    ! Actual values of the particles

      end subroutine explicit_ptr

      subroutine implicit_nkPtr(np,dt,x_tn,x_tn1,res)
      use module_pepc_types
      use module_pepc_kinds
      implicit none

        integer(kind_particle)       , intent(in)    :: np
        real(kind_particle)          , intent(in)    :: dt
        type(t_particle)             , intent(in)    :: x_tn(np)    ! Actual values of the particles
        real(kind_particle)          , intent(in)    :: x_tn1(9*np) ! Next iteration of the particles
        real(kind_particle)          , intent(out)   :: res(9*np)   ! Residual

      end subroutine implicit_nkPtr

      subroutine implicit_picardPtr(np,dt,p,ierr,itc,tolg)
      use module_pepc_types
      use module_pepc_kinds
      implicit none

        integer(kind_particle)          , intent(in)     :: np
        real(kind_particle)             , intent(in)     :: dt
        type(t_particle), allocatable   , intent(inout)  :: p(:)
        integer(kind_particle)          , intent(out)    :: ierr,itc
        real(kind_particle)             , intent(out)    :: tolg

      end subroutine implicit_picardPtr
   end interface

  contains


  subroutine march(np,dt,p,ischeme,adv)
  use newton_krylov, only: nsolgm
  use module_globals, only: my_rank,root
  implicit none
        integer(kind_particle)       , intent(in)    :: np
        character(255)               , intent(in)    :: ischeme ! numerical scheme
        integer(kind_particle)       , intent(in)    :: adv ! referred to the time integrator: 0 explicit, 1 implicit, newton-kylov, 2 implicit, picard
        real(kind_particle)          , intent(in)    :: dt
        type(t_particle), allocatable, intent(inout) :: p(:)    ! Actual values of the particles

        procedure (implicit_picardPtr), pointer      :: pic_ptr => null ()
        procedure (implicit_nkPtr)    , pointer      :: nk_ptr  => null ()
        procedure (explicit_ptr)      , pointer      :: exp_ptr => null ()
        real(kind_particle)                          :: x_tn1(9*np) ! NK-Iteration
        real(kind_particle)                          :: res(9*np) ,static_gmres(3) ! Residual - NK-Iteration and statistics
        integer(kind_particle)                       :: errmsg,iter,rc=100    ! Picard Iteration
        real(kind_particle)                          :: errnk        ! Global Error - Picard Iteration

        if (adv .eq. 0) then

            select case (ischeme)
!                case ("boris")
!                    if (my_rank.eq. 0) write(*,'(a)')        " ====== Boris Integrator"
!                    exp_ptr => boris
                case ("leapfrog")
                    if (my_rank.eq. 0) write(*,'(a)')        " ====== Leapfrog Integrator"
                    exp_ptr => leapfrog
                case default
                    if (my_rank.eq. 0) write(*,'(a)')        " ====== Leapfrog Integrator"
                    exp_ptr => leapfrog
                end select

            call exp_ptr(np,dt,p)

        elseif (adv.eq.1) then
            select case (ischeme)
                case ("midpoint_picard")
                    if (my_rank.eq. 0) write(*,'(a)')        " ====== Midpoint Picard"
                    pic_ptr => midpoint
                case ("midpoint_picard_electrostatic")
                    if (my_rank.eq. 0) write(*,'(a)')        " ====== Midpoint Picard Electrostatic"
                    pic_ptr => midpoint_electrostatic
                case ("trapezoidal_picard")
                    if (my_rank.eq. 0) write(*,'(a)')        " ====== Trapezoidal_picard"
                    pic_ptr => trapezoidal
                case ("trapezoidal_broyden")
                    if (my_rank.eq. 0) write(*,'(a)')        " ====== Trapezoidal_picard"
                    pic_ptr => trapezoidal_broyden
                case ("hamiltonian_boris")
                    if (my_rank.eq. 0) write(*,'(a)')        " ====== Boris - Hamiltonian Integrator"
                    pic_ptr => hamiltonian_boris 
                case default
                    if (my_rank.eq. 0) write(*,'(a)')        " ====== Boris - Hamiltonian Integrator"
                    pic_ptr => hamiltonian_boris    
                    
                end select
                call pic_ptr(np,dt,p,errmsg,iter,errnk)
                
                if (root) then
                    write(*,'(a,i4,2x,es12.4,i4)') " == Resual : ", iter,errnk,errmsg
                    open(unit=rc,file=trim(folder)//trim("residuo_")//trim(adjustl(ischeme))//".dat",form='formatted',status='unknown',position='append')
                    write(rc,*) iter,errnk,errmsg
                    close (rc )
                endif

        elseif (adv.eq.2) then

            select case (ischeme)
                case ("midpoint3D")
                    nk_ptr => midpoint3D
                case ("trapezoidal3D")
                    nk_ptr => trapezoidal3D
                case default
                    nk_ptr => midpoint3D
                end select
                call nsolgm(np,dt,p,nk_ptr,errmsg,iter,errnk,static_gmres)
                
                if (root) then
                    write(*,'(a,es12.4)') " == Resual : ", iter,errnk,errmsg
                    open(unit=rc,file=trim(folder)//trim("residuo_")//trim(adjustl(ischeme))//".dat",form='formatted',status='unknown',position='append')
                    write(rc,*) errmsg,iter,errnk
                    close (rc )
                endif
                
        elseif (adv.eq.-1) then
            if (my_rank.eq. 0) write(*,'(a)')        " ====== Wrong choice of time integrator"
            call exit(1)

        endif


  end subroutine march



  subroutine leapfrog(np,dt,p)
    use helper         ,only: iperiodic_particles
    use module_shortcut,only: one
    use module_globals ,only: periodicity_particles

    implicit none
    integer(kind_particle)       , intent(in)    :: np
    type(t_particle), allocatable, intent(inout) :: p(:)
    real(kind_particle)          , intent(in)    :: dt

    integer(kind_particle) :: ip
    real(kind_particle)    :: v(3),gamma


    do ip=1, np
      gamma             = one!/sqrt( one - dot_product( p(ip)%data%v/vtilde, p(ip)%data%v/vtilde ) ) 
      v                 = gamma*p(ip)%data%v + dt * p(ip)%data%q / p(ip)%data%m * p(ip)%results%e
      gamma             = sqrt( one + dot_product( v/vtilde, v/vtilde ) )

      p(ip)%data%v      = v/gamma
      p(ip)%x           = p(ip)%x      + dt   * p(ip)%data%v

      p(ip)%x(3)        = zero
      if (periodicity_particles) call iperiodic_particles(ip,np,p)

    end do

  end subroutine leapfrog


!subroutine boris(np,dt,p)
!    use module_tool   , only: cross_product
!    implicit none
!    include 'mpif.h'
!
!    integer(kind_particle)       , intent(in)    :: np
!    real(kind_particle)          , intent(in)    :: dt
!    type(t_particle),allocatable , intent(inout) :: p(:)    ! Actual values of the particles
!
!    integer(kind_particle)                    :: ip,jp,rc,ierr
!    real(kind_particle)                       :: v(3),e,m,gamma,&
!                                                 alpha,Bnorm,B(3),v_(3),vpar(3),vnor(3)
!
!
!    do ip=1, np
!
!      e                 = p(ip)%data%q
!      m                 = p(ip)%data%m
!      B                 = p(ip)%results%B
!      Bnorm             = B(1)**2 + B(2)**2 + B(3)**2
!
!      alpha             = half*e*dt/m
!
!      v_                = p(ip)%data%v + alpha*( p(ip)%results%E )
!      vpar              = alpha*( dot_product(v_,B) )*B
!      vnor              = cross_product(v_,B)
!      v                 = (v_ + alpha*( vpar + vnor ) )/( one + alpha**2*Bnorm )
!
!      p(ip)%data%v      = two*v - p(ip)%data%v
!      p(ip)%x           = p(ip)%x + dt*v
!
!    end do
!
!  end subroutine boris


  subroutine midpoint3D(np,dt,x_tn,x_tn1,res)
    use module_tool     , only: cross_product
    implicit none
    include 'mpif.h'

    integer(kind_particle)       , intent(in) :: np
!    integer(kind_particle)                    :: n
    real(kind_particle)          , intent(in) :: dt
    type(t_particle)             , intent(in) :: x_tn(np)    ! Actual values of the particles
    real(kind_particle)          , intent(in) :: x_tn1(9*np) ! Next iteration of the particles
    real(kind_particle)          , intent(out):: res(9*np)   ! Residual

    type(t_particle), allocatable             :: r(:)
    integer(kind_particle)                    :: ip,jp,rc,ierr
    real(kind_particle)                       :: v(3),x(3),p(3),rot(3),e,m,gamma

    if ( allocated(r) ) deallocate(r)
    allocate( r(size(x_tn, kind=kind_particle)) , stat = rc )
    if(rc.ne.0) write(*,*) " === particle allocation error!"

    do ip = 1, np

        jp                = (ip-1)*9
        m                 = x_tn(ip)%data%m
        e                 = x_tn(ip)%data%q

        r(ip)%label       = x_tn(ip)%label
        r(ip)%data%q      = e
        r(ip)%data%m      = m

        x(1:3)            = x_tn1(jp+1:jp+3)
        v(1:3)            = (x_tn1(jp+4:jp+6) - e/lorentz_tilde*x_tn1(jp+7:jp+9))/m
        gamma             = sqrt( one + dot_product(v/vtilde,v/vtilde) )
        v(1:3)            = v(1:3)/gamma

        r(ip)%x           = half*( x + x_tn(ip)%x )
        r(ip)%data%v      = half*( v + x_tn(ip)%data%v )

        r(ip)%results%e   = zero
        r(ip)%results%pot = zero
        r(ip)%results%B   = zero
        r(ip)%results%A   = zero
        r(ip)%results%dxA = zero
        r(ip)%results%dyA = zero
        r(ip)%results%Jirr= zero
        r(ip)%results%J   = zero
        r(ip)%work        = one

    enddo

    call pepc_particleresults_clear(r)
    call pepc_grow_tree(r)
    call pepc_traverse_tree(r)
    call pepc_restore_particles(r)
    call pepc_timber_tree()

    do ip = 1,np

        jp                = (ip-1)*9
        m                 = x_tn(ip)%data%m
        e                 = x_tn(ip)%data%q

        gamma             = one/sqrt( one - dot_product( x_tn(ip)%data%v/vtilde, x_tn(ip)%data%v/vtilde ) )

        p(1:3)            = gamma*m*x_tn(ip)%data%v(1:3) + e/lorentz_tilde*x_tn(ip)%results%A

        !!! v cross B

        rot               = cross_product( r(ip)%data%v/lorentz_tilde,r(ip)%results%B )

        res(jp+4:jp+6)    = x_tn1(jp+4:jp+6)   - p(1:3)             - e*dt*(  r(ip)%results%E(1:3) + rot  )
        res(jp+1:jp+3)    = x_tn1(jp+1:jp+3)   - x_tn(ip)%x(1:3)    - dt*r(ip)%data%v
        res(jp+7:jp+9)    = x_tn1(jp+7:jp+9)   - r(ip)%results%A

    enddo

    deallocate(r)

  end subroutine midpoint3D

  subroutine trapezoidal3D(np,dt,x_tn,x_tn1,res)
    use module_tool   , only: cross_product
    use module_globals, only: norm_factor
    implicit none
    include 'mpif.h'

    integer(kind_particle)       , intent(in) :: np
!    integer(kind_particle)                    :: n
    real(kind_particle)          , intent(in) :: dt
    type(t_particle)             , intent(in) :: x_tn(np)     ! Actual values of the particles
    real(kind_particle)          , intent(in) :: x_tn1(9*np) ! Next iteration of the particles
    real(kind_particle)          , intent(out):: res(9*np)   ! Residual

    type(t_particle), allocatable             :: r(:)
    integer(kind_particle)                    :: ip,jp,rc,ierr
    real(kind_particle)                       :: v(3),vn(3),x(3),p(3),rot(3),e,m,gamma,B(3),Bn(3)
    !n  = 6*np

    if ( allocated(r) ) deallocate(r)
    allocate( r(size(x_tn, kind=kind_particle)) , stat = rc )
    if(rc.ne.0) write(*,*) " === particle allocation error!"

    do ip = 1, np

        jp                = (ip-1)*9
        m                 = x_tn(ip)%data%m
        e                 = x_tn(ip)%data%q

        r(ip)%label       = x_tn(ip)%label
        r(ip)%data%q      = e
        r(ip)%data%m      = m

        v(1:3)            = (x_tn1(jp+4:jp+6) - e/lorentz_tilde*x_tn1(jp+7:jp+9))/m
        gamma             = one!sqrt( one + dot_product(v,v) )
        v(1:3)            = v(1:3)/gamma

        r(ip)%x           = x_tn1(jp+1:jp+3)
        r(ip)%data%v      = v

        r(ip)%results%e   = zero
        r(ip)%results%pot = zero
        r(ip)%results%B   = zero
        r(ip)%results%A   = zero
        r(ip)%results%dxA = zero
        r(ip)%results%dyA = zero
        r(ip)%results%Jirr= zero
        r(ip)%results%J   = zero
        r(ip)%work        = one

    enddo

    !!! Here we update the fields at x^n+1 and v^n+1

    call pepc_particleresults_clear(r)
    call pepc_grow_tree(r)
    call pepc_traverse_tree(r)
    call pepc_restore_particles(r)
    call pepc_timber_tree()


    do ip = 1,np

        jp                = (ip-1)*9
        m                 = x_tn(ip)%data%m
        e                 = x_tn(ip)%data%q

        x                 = x_tn(ip)%x       ! x^n
        v                 = x_tn(ip)%data%v  ! v^n
        gamma             = one!one/sqrt( 1 - dot_product(v(1:ivdim)/vtilde,v(1:ivdim)/vtilde) )
        p                 = gamma*m*v +  e/lorentz_tilde*x_tn(ip)%results%A ! p^n
        B                 = x_tn(ip)%results%B

        vn                = r(ip)%data%v  ! v^n
        Bn                = r(ip)%results%B
        rot               = cross_product( ( v + vn )/lorentz_tilde , half*( B + Bn ) )



        res(jp+1:jp+3)    = x_tn1(jp+1:jp+3) - x(1:3) - half*dt*( v(1:3) + vn(1:3) )
        res(jp+4:jp+6)    = x_tn1(jp+4:jp+6) - p(1:3) - half*e*dt*(  r(ip)%results%E(1:3) + x_tn(ip)%results%E(1:3) + rot(1:3)  )
        res(jp+7:jp+9)    = x_tn1(jp+7:jp+9) - r(ip)%results%A(1:3)

    enddo

    deallocate(r)

  end subroutine trapezoidal3D

    subroutine midpoint(np,dt,p,ierr,itc,tolg)
    use module_tool   , only: cross_product
    use module_globals, only: root,my_rank,ixdim,ivdim,tnp,lorentz_tilde,vtilde,periodicity_particles
    use helper        , only: iperiodic_particles,inormalize
    use module_pepc

    implicit none
    include 'mpif.h'

    integer(kind_particle)          , intent(in)     :: np
    integer(kind_particle)                           :: n
    real(kind_particle)             , intent(in)     :: dt
    type(t_particle), allocatable   , intent(inout)  :: p(:)
    integer(kind_particle)          , intent(out)    :: ierr,itc
    real(kind_particle)             , intent(out)    :: tolg

    real(kind_particle)                              :: e,m,toll,stop_tol,alpha,rot(3),grad(3),v(3)
    type(t_particle), allocatable                    :: r(:)
    type(pma_particle), allocatable                  :: x(:)
    integer(kind_particle)                           :: rc,ip,maxit

    if (allocated(r)) deallocate(r)
    if (allocated(x)) deallocate(x)
    allocate( r(np), stat=rc )
    allocate( x(np), stat=rc )

    ierr            = 0
    itc             = 0
    maxit           = 100
    tolg            = one
    stop_tol        = prec!tentominusseven!prec


    !!! Initial Guess
    do ip = 1,np

            e            =  p(ip)%data%q
            m            =  p(ip)%data%m
    
            x(ip)%P      =  m*p(ip)%data%v + e/lorentz_tilde*p(ip)%results%A  ! p(ip)%data%v == gamma*v
            x(ip)%x      =  p(ip)%x 
            x(ip)%A      =  p(ip)%results%A

            r(ip)%label       = p(ip)%label
            r(ip)%data%q      = e
            r(ip)%data%m      = m
    enddo
    ! main iteration loop
    !

    do while(tolg .gt. stop_tol .and. itc .lt. maxit)

        itc = itc + 1

            do ip = 1,np

              e                 =  p(ip)%data%q
              m                 =  p(ip)%data%m

              v                 = ( x(ip)%P - e/lorentz_tilde*x(ip)%A )/m
              r(ip)%data%g      = sqrt( one + dot_product(v/vtilde,v/vtilde) ) + p(ip)%data%g
              
              r(ip)%data%v      = (v + p(ip)%data%v)!/( p(ip)%data%g + r(ip)%data%g ) ! p(ip)%data%v == gamma*v

              r(ip)%x           = half*(x(ip)%x + p(ip)%x)
              r(ip)%x(3)        = zero

              r(ip)%label       = p(ip)%label
              r(ip)%data%q      = e
              r(ip)%data%m      = m

              r(ip)%results%e   = zero
              r(ip)%results%pot = zero
              r(ip)%results%B   = zero
              r(ip)%results%A   = zero
              r(ip)%results%dxA = zero
              r(ip)%results%dyA = zero
              r(ip)%results%Jirr= zero
              r(ip)%results%J   = zero
              r(ip)%work        = one
              
              if (periodicity_particles) call iperiodic_particles(ip,np,r)

            enddo


        !!! Update fields

            call pepc_particleresults_clear(r)
            call pepc_grow_tree(r)
            call pepc_traverse_tree(r)
            call pepc_restore_particles(r)
            call pepc_timber_tree()

            toll = zero
            tolg = zero

            do ip = 1,np

              e                 = p(ip)%data%q
              m                 = p(ip)%data%m

              r(ip)%x           = p(ip)%x +      dt*r(ip)%data%v
              r(ip)%x(3)        = zero
              
              call  inormalize(ip,np,r)
              
              rot               = cross_product( r(ip)%data%v, r(ip)%results%B )
              grad(1)           = v(1)*r(ip)%results%dxA(1) + v(2)*r(ip)%results%dyA(1)
              grad(2)           = v(1)*r(ip)%results%dxA(2) + v(2)*r(ip)%results%dyA(2)
              grad(3)           = v(1)*r(ip)%results%dxA(3) + v(2)*r(ip)%results%dyA(3)
              
              rot               = rot/r(ip)%data%g/lorentz_tilde
              grad              = grad/r(ip)%data%g/lorentz_tilde
              
              r(ip)%data%v      = m*p(ip)%data%v + e/lorentz_tilde*p(ip)%results%A + dt*e*( r(ip)%results%E + rot + grad )

              toll              = toll + dot_product( r(ip)%data%v    - x(ip)%P , r(ip)%data%v    - x(ip)%P )         &
                                       + dot_product( r(ip)%x         - x(ip)%x , r(ip)%x         - x(ip)%x )         &
                                       + dot_product( r(ip)%results%A - x(ip)%A , r(ip)%results%A - x(ip)%A )

              if (periodicity_particles) call iperiodic_particles(ip,np,r)                         
              x(ip)%P           =  r(ip)%data%v
              x(ip)%x           =  r(ip)%x
              x(ip)%x(3)        =  zero
              x(ip)%A           =  r(ip)%results%A

              
            enddo

            call MPI_ALLREDUCE(toll, tolg, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
            tolg  = sqrt(tolg)/real(tnp,kind=kind_particle)


    enddo

    if (tolg > stop_tol)  ierr = 1
    

    do ip = 1,np

        e                 = p(ip)%data%q
        m                 = p(ip)%data%m

        p(ip)%data%v      = ( x(ip)%P - e/lorentz_tilde*x(ip)%A )/m
        p(ip)%data%g      = sqrt( one + dot_product(p(ip)%data%v/vtilde,p(ip)%data%v/vtilde) )
        p(ip)%data%v      = p(ip)%data%v
        p(ip)%x           = x(ip)%x + dt*p(ip)%data%v/p(ip)%data%g
        p(ip)%x(3)        = zero
        if (periodicity_particles) call iperiodic_particles(ip,np,p)

    enddo

    deallocate( r,x )


    end subroutine midpoint
    
    
    subroutine midpoint_electrostatic(np,dt,p,ierr,itc,tolg)
    use module_globals, only: root,my_rank,ixdim,ivdim,tnp,periodicity_particles
    use helper        , only: iperiodic_particles,inormalize
    use module_pepc

    implicit none
    include 'mpif.h'

    integer(kind_particle)          , intent(in)     :: np
    integer(kind_particle)                           :: n
    real(kind_particle)             , intent(in)     :: dt
    type(t_particle), allocatable   , intent(inout)  :: p(:)
    integer(kind_particle)          , intent(out)    :: ierr,itc
    real(kind_particle)             , intent(out)    :: tolg

    real(kind_particle)                              :: e,m,toll,stop_tol,v(3)
    type(t_particle), allocatable                    :: r(:)
    type(pma_particle), allocatable                  :: x(:)
    integer(kind_particle)                           :: rc,ip,maxit

    if (allocated(r)) deallocate(r)
    if (allocated(x)) deallocate(x)
    allocate( r(np), stat=rc )
    allocate( x(np), stat=rc )

    ierr            = 0
    itc             = 0
    maxit           = 100
    tolg            = one
    stop_tol        = prec!prec


    !!! Initial Guess
    do ip = 1,np

            e                 =  p(ip)%data%q

            x(ip)%P           =  p(ip)%data%v
            x(ip)%x           =  p(ip)%x 

            r(ip)%label       =  p(ip)%label
            r(ip)%data%q      =  e
            r(ip)%data%m      =  m
            
    enddo
    ! main iteration loop
    !

    do while(tolg .gt. stop_tol .and. itc .lt. maxit)

        itc = itc + 1

            do ip = 1,np

              e                 = p(ip)%data%q
              m                 = p(ip)%data%m 

              v                 = x(ip)%P
              
              r(ip)%data%g      = sqrt(one + dot_product(v/vtilde,v/vtilde) )
              r(ip)%data%v      = v!half*(v       + p(ip)%data%v)
              r(ip)%x           = x(ip)%x!half*(x(ip)%x + p(ip)%x     )
              r(ip)%x(3)        = zero

              r(ip)%label       = p(ip)%label
              r(ip)%data%q      = e
              r(ip)%data%m      = m

              r(ip)%results%E   = zero
              r(ip)%results%pot = zero
              r(ip)%results%B   = zero
              r(ip)%results%A   = zero
              r(ip)%results%dxA = zero
              r(ip)%results%dyA = zero
              r(ip)%results%Jirr= zero
              r(ip)%results%J   = zero
              r(ip)%work        = one
              
              if (periodicity_particles) call iperiodic_particles(ip,np,r)

            enddo


        !!! Update fields

            call pepc_particleresults_clear(r)
            call pepc_grow_tree(r)
            call pepc_traverse_tree(r)
            call pepc_restore_particles(r)
            call pepc_timber_tree()

            toll = zero
            tolg = zero

            do ip = 1,np

              e                 = p(ip)%data%q
              m                 = p(ip)%data%m
              
              call  inormalize(ip,np,r)

              r(ip)%data%v      = p(ip)%data%v  + half*dt*e/m*( r(ip)%results%E +  p(ip)%results%E ) 
              r(ip)%x           = p(ip)%x             + dt*    ( r(ip)%data%v    +  p(ip)%data%v    )/( r(ip)%data%g + p(ip)%data%g  ) 
              r(ip)%x(3)        = zero


              toll              = toll + dot_product(r(ip)%data%v - x(ip)%P, r(ip)%data%v - x(ip)%P)         &
                                       + dot_product(r(ip)%x      - x(ip)%x, r(ip)%x      - x(ip)%x)           
                                       
              if (periodicity_particles) call iperiodic_particles(ip,np,r)
              x(ip)%P           =  r(ip)%data%v
              x(ip)%x           =  r(ip)%x
              x(ip)%x(3)        =  zero

              
            enddo

            call MPI_ALLREDUCE(toll, tolg, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
            tolg  = sqrt(tolg)/real(tnp,kind=kind_particle)

    enddo

    if (tolg > stop_tol)  ierr = 1
    

    do ip = 1,np

        p(ip)%x           = x(ip)%x + dt*( p(ip)%data%v + x(ip)%P )/( p(ip)%data%g + r(ip)%data%g )  
        p(ip)%x(3)        = zero
        p(ip)%data%v      = ( x(ip)%P  )
        p(ip)%data%g      = sqrt( one + dot_product(p(ip)%data%v/vtilde,p(ip)%data%v/vtilde) )
          
        
        
        if (periodicity_particles) call iperiodic_particles(ip,np,p)

    enddo

    deallocate( r,x )


    end subroutine midpoint_electrostatic



    subroutine trapezoidal(np,dt,p,ierr,itc,tolg)
    use module_tool   , only: cross_product
    use module_globals, only: root,my_rank,ixdim,ivdim,tnp,periodicity_particles
    use helper        , only: iperiodic_particles,inormalize
    use module_pepc

    implicit none
    include 'mpif.h'

    integer(kind_particle)          , intent(in)     :: np
    integer(kind_particle)                           :: n
    real(kind_particle)             , intent(in)     :: dt
    type(t_particle), allocatable   , intent(inout)  :: p(:)
    integer(kind_particle)          , intent(out)    :: ierr,itc
    real(kind_particle)             , intent(out)    :: tolg

    real(kind_particle)                              :: e,m,toll,stop_tol,alpha,v(1:3),grad(1:3)!,rot(1:3)
    type(t_particle), allocatable                    :: r(:)
    type(pma_particle), allocatable                  :: x(:)
    integer(kind_particle)                           :: rc,ip,maxit

    if (allocated(r)) deallocate(r)
    if (allocated(x)) deallocate(x)
    allocate( r(np), stat=rc )
    allocate( x(np), stat=rc )

    ierr            = 0
    itc             = 0
    maxit           = 100
    tolg            = one
    stop_tol        = tentominusnine!prec!tentominusseven!prec


    !!! Initial Guess
    do ip = 1,np

            e            =  p(ip)%data%q
            m            =  p(ip)%data%m

            x(ip)%P      =  m*p(ip)%data%v + e/lorentz_tilde*p(ip)%results%A
            x(ip)%x      =  p(ip)%x 
            x(ip)%x(3)   =  zero
            x(ip)%A      =  p(ip)%results%A

            r(ip)%label       = p(ip)%label
            r(ip)%data%q      = e
            r(ip)%data%m      = m
    enddo
    ! main iteration loop
    !

    do while(tolg .gt. stop_tol .and. itc .lt. maxit)

        itc = itc + 1

            do ip = 1,np

              e                 =  p(ip)%data%q
              m                 =  p(ip)%data%m

              r(ip)%data%v      = ( x(ip)%P - e/lorentz_tilde*x(ip)%A )/m
!              if ( r(ip)%label .eq. 1 ) r(ip)%data%v = zero
              r(ip)%data%g      = sqrt( one + dot_product(r(ip)%data%v/vtilde,r(ip)%data%v/vtilde) )
              r(ip)%x           = x(ip)%x
              r(ip)%x(3)        = zero

              r(ip)%label       = p(ip)%label
              r(ip)%results%E   = zero
              r(ip)%results%pot = zero
              r(ip)%results%B   = zero
              r(ip)%results%A   = zero
              r(ip)%results%dxA = zero
              r(ip)%results%dyA = zero
              r(ip)%results%Jirr= zero
              r(ip)%results%J   = zero
              r(ip)%work        = one
              
              if (periodicity_particles) call iperiodic_particles(ip,np,r)

            enddo


        !!! Update fields

            call pepc_particleresults_clear(r)
            call pepc_grow_tree(r)
            call pepc_traverse_tree(r)
            call pepc_restore_particles(r)
            call pepc_timber_tree()

            toll = zero
            tolg = zero

            do ip = 1,np

              e                 = p(ip)%data%q
              m                 = p(ip)%data%m
              
              call  inormalize(ip,np,r)

              v                 = ( r(ip)%data%v + p(ip)%data%v )
              
!              if ( r(ip)%label .eq. 1 ) v = zero
              
!              rot               = cross_product( v , ( r(ip)%results%B + p(ip)%results%B ) )
              
!              grad(1)           = v(1)*r(ip)%results%dxA(1) + v(2)*r(ip)%results%dyA(1)
!              grad(2)           = v(1)*r(ip)%results%dxA(2) + v(2)*r(ip)%results%dyA(2)
!              grad(3)           = v(1)*r(ip)%results%dxA(3) + v(2)*r(ip)%results%dyA(3)
!              
!              grad(1)           = grad(1) + v(1)*p(ip)%results%dxA(1) + v(2)*p(ip)%results%dyA(1)
!              grad(2)           = grad(2) + v(1)*p(ip)%results%dxA(2) + v(2)*p(ip)%results%dyA(2)
!              grad(3)           = grad(3) + v(1)*p(ip)%results%dxA(3) + v(2)*p(ip)%results%dyA(3)
              
              grad              = zero
              
              grad(1)           = v(1)*r(ip)%results%dxA(1) + v(2)*r(ip)%results%dxA(2) + v(3)*r(ip)%results%dxA(3)
              grad(2)           = v(1)*r(ip)%results%dyA(2) + v(2)*r(ip)%results%dyA(2) + v(3)*r(ip)%results%dyA(3)
                            
              grad(1)           = grad(1) + v(1)*p(ip)%results%dxA(1) + v(2)*p(ip)%results%dxA(2) + v(3)*p(ip)%results%dxA(3)
              grad(2)           = grad(2) + v(1)*p(ip)%results%dyA(2) + v(2)*p(ip)%results%dyA(2) + v(3)*p(ip)%results%dyA(3)
              
!              rot               = rot/( r(ip)%data%g + p(ip)%data%g )/lorentz_tilde
              grad              = grad/(r(ip)%data%g + p(ip)%data%g )/lorentz_tilde
              
              
              r(ip)%x           = p(ip)%x + dt*v/( r(ip)%data%g + p(ip)%data%g )
              
              r(ip)%data%v(1:2) = m*p(ip)%data%v(1:2) + e/lorentz_tilde*p(ip)%results%A(1:2)  &
                                  + half*dt*e*( r(ip)%results%E(1:2)  + p(ip)%results%E(1:2)  + grad(1:2) )    ! Canonical Momentum
                                  
              r(ip)%x(3)        = zero                    
              r(ip)%data%v(3)   = m*p(ip)%data%v(3) + e/lorentz_tilde*p(ip)%results%A(3)                    

              toll              = toll + dot_product(r(ip)%data%v(1:2)     - x(ip)%P(1:2) , r(ip)%data%v(1:2)    - x(ip)%P(1:2))   &
                                       + dot_product(r(ip)%x               - x(ip)%x      , r(ip)%x              - x(ip)%x)        &
                                       + dot_product(r(ip)%results%A(1:2)  - x(ip)%A(1:2) , r(ip)%results%A(1:2) - x(ip)%A(1:2))
                                       
              if (periodicity_particles) call iperiodic_particles(ip,np,r)
                                       
              x(ip)%P           =  r(ip)%data%v
              x(ip)%x           =  r(ip)%x
              x(ip)%x(3)        =  zero
              x(ip)%A           =  r(ip)%results%A
              
              

            enddo

            call MPI_ALLREDUCE(toll, tolg, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
            tolg  = sqrt(tolg)/real(tnp,kind=kind_particle)


    enddo

    if (tolg > stop_tol)  ierr = 1
    

    do ip = 1,np

        e                 = p(ip)%data%q
        m                 = p(ip)%data%m

        p(ip)%data%v      = ( x(ip)%P - e/lorentz_tilde*x(ip)%A )/m
!        p(ip)%data%v(3)   =  p(ip)%data%v(3) + e/lorentz_tilde/m*( p(ip)%results%A(3) - x(ip)%A(3) )
!        if ( p(ip)%label .eq. 1 ) p(ip)%data%v = zero
        p(ip)%data%g      = sqrt( one + dot_product( p(ip)%data%v/vtilde ,  p(ip)%data%v/vtilde ) )
        p(ip)%x           = x(ip)%x
        p(ip)%x(3)        = zero
        if (periodicity_particles) call iperiodic_particles(ip,np,p)
        
    enddo

    deallocate( r,x )


    end subroutine trapezoidal

    
    
    subroutine trapezoidal_broyden(np,dt,p,ierr,itc,tolg)
    use module_tool   , only: cross_product
    use module_globals, only: root,my_rank,ixdim,ivdim,tnp,periodicity_particles
    use helper        , only: iperiodic_particles,inormalize
    use module_pepc

    implicit none
    include 'mpif.h'

    integer(kind_particle)          , intent(in)     :: np
    real(kind_particle)             , intent(in)     :: dt
    type(t_particle), allocatable   , intent(inout)  :: p(:)
    integer(kind_particle)          , intent(out)    :: ierr,itc
    real(kind_particle)             , intent(out)    :: tolg

    real(kind_particle)                              :: e,m,toll,stop_tol,alpha,rot(1:3),v(1:3),grad(1:3),z(1:5)
    type(t_particle), allocatable                    :: r(:)
    type(pma_particle), allocatable                  :: x(:)
    integer(kind_particle)                           :: rc,ip,jp,maxit
    real(kind_particle)  , allocatable               :: s(:,:)
 
    ierr            = 0
    itc             = 0
    maxit           = 100
    tolg            = one
    stop_tol        = prec!tentominusseven!prec
    
    if (allocated(r)) deallocate(r)
    if (allocated(s)) deallocate(s)
    if (allocated(x)) deallocate(x)
    allocate( r(np), stat=rc )
    allocate( x(np), stat=rc )
    allocate( s(maxit,5*np), stat=rc )
    
    !!! Initial Guess
    
    
    
    do ip = 1,np
            
        e           =    p(ip)%data%q
        rot         =    cross_product( p(ip)%data%v , p(ip)%results%B )
        grad(1)     =    p(ip)%data%v(1)*p(ip)%results%dxA(1) + p(ip)%data%v(2)*p(ip)%results%dyA(1)
        grad(2)     =    p(ip)%data%v(1)*p(ip)%results%dxA(2) + p(ip)%data%v(2)*p(ip)%results%dyA(2)
        grad(3)     =    p(ip)%data%v(1)*p(ip)%results%dxA(3) + p(ip)%data%v(2)*p(ip)%results%dyA(3)
        
        rot         = rot/lorentz_tilde/p(ip)%data%g
        grad        = rot/lorentz_tilde/p(ip)%data%g
        
        s(1,ip)     =    dt*p(ip)%data%v(1)
        s(1,np+ip)  =    dt*p(ip)%data%v(2)
        s(1,2*np+ip)=    e*dt*(p(ip)%results%E(1) + rot(1) + grad(1) )
        s(1,3*np+ip)=    e*dt*(p(ip)%results%E(2) + rot(2) + grad(2) )   
        s(1,4*np+ip)=    e*dt*(                     rot(3) + grad(3) )   
        
        x(ip)%A(1:3)=    p(ip)%results%A(1:3)
        
        
        r(ip)%x     = p(ip)%x
        r(ip)%x(3)  = zero
        r(ip)%data%v= p(ip)%data%v 
        r(ip)%data%g= p(ip)%data%g
        
            
    enddo
    ! main iteration loop
    !

    do while(tolg .gt. stop_tol .and. itc .lt. maxit)

        itc = itc + 1

            do ip = 1,np

              e                 =  p(ip)%data%q
              m                 =  p(ip)%data%m

              r(ip)%x(1)        = r(ip)%x(1)      + s(itc,ip)
              r(ip)%x(2)        = r(ip)%x(2)      + s(itc,np+ip)
              r(ip)%x(3)        = zero
              r(ip)%data%v(1)   = r(ip)%data%v(1) + ( s(itc,2*np+ip) - e*x(ip)%A(1) )/m
              r(ip)%data%v(2)   = r(ip)%data%v(2) + ( s(itc,3*np+ip) - e*x(ip)%A(2) )/m
              r(ip)%data%v(3)   = r(ip)%data%v(3) + ( s(itc,4*np+ip) - e*x(ip)%A(3) )/m
              
              r(ip)%data%g      = sqrt( one + dot_product(r(ip)%data%v/vtilde,r(ip)%data%v/vtilde) )
              

              r(ip)%label       = p(ip)%label
              r(ip)%results%E   = zero
              r(ip)%results%pot = zero
              r(ip)%results%B   = zero
              r(ip)%results%A   = zero
              r(ip)%results%dxA = zero
              r(ip)%results%dyA = zero
              r(ip)%results%Jirr= zero
              r(ip)%results%J   = zero
              r(ip)%work        = one
              
              if (periodicity_particles) call iperiodic_particles(ip,np,r)

            enddo


        !!! Update fields

            call pepc_particleresults_clear(r)
            call pepc_grow_tree(r)
            call pepc_traverse_tree(r)
            call pepc_restore_particles(r)
            call pepc_timber_tree()

            toll = zero
            tolg = zero
 
            do ip = 1,np

              e                 = p(ip)%data%q
              m                 = p(ip)%data%m
              
              call  inormalize(ip,np,r)

              v                 = ( r(ip)%data%v + p(ip)%data%v )
              
              rot               = cross_product( v , ( r(ip)%results%B + p(ip)%results%B ) )
              
              grad(1)           = v(1)*r(ip)%results%dxA(1) + v(2)*r(ip)%results%dyA(1)
              grad(2)           = v(1)*r(ip)%results%dxA(2) + v(2)*r(ip)%results%dyA(2)
              grad(3)           = v(1)*r(ip)%results%dxA(3) + v(2)*r(ip)%results%dyA(3)
              
              grad(1)           = grad(1) + v(1)*p(ip)%results%dxA(1) + v(2)*p(ip)%results%dyA(1)
              grad(2)           = grad(2) + v(1)*p(ip)%results%dxA(2) + v(2)*p(ip)%results%dyA(2)
              grad(3)           = grad(3) + v(1)*p(ip)%results%dxA(3) + v(2)*p(ip)%results%dyA(3)
              
              rot               = rot/( r(ip)%data%g + p(ip)%data%g )/lorentz_tilde
              grad              = grad/(r(ip)%data%g + p(ip)%data%g )/lorentz_tilde
              
              r(ip)%x           = p(ip)%x + dt*v/( r(ip)%data%g + p(ip)%data%g )
              r(ip)%x(3)        = zero
              r(ip)%data%v      = m*p(ip)%data%v + e/lorentz_tilde*p(ip)%results%A  &
                                  + half*dt*e*( r(ip)%results%E  + p(ip)%results%E + rot + grad)    ! Canonical Momentum
                                  
                                  
              z(1:2)            =  r(ip)%x(1:2)      - x(ip)%x(1:2)
              z(3:5)            =  r(ip)%data%v(1:3) - x(ip)%P(1:3)
              
              do jp = 2,itc
                  
                  alpha = s(jp,ip)*s(jp-1,ip) + s(jp,np+ip)*s(jp-1,np+ip)     + s(jp,2*np+ip)*s(jp-1,2*np+ip) &
                                              + s(jp,3*np+ip)*s(jp-1,3*np+ip) + s(jp,4*np+ip)*s(jp-1,4*np+ip)   
                  alpha = alpha/( s(jp-1,ip)**2 + s(jp-1,np+ip)**2 + s(jp-1,2*np+ip)**2 + s(jp-1,3*np+ip)**2 + s(jp-1,4*np+ip)**2 )                              
                  z     = z + alpha*z
                  
              enddo
              
              alpha             =    s(itc,ip)*z(1) + s(itc,np+ip)*z(2) + s(itc,2*np+ip)*z(3) + s(itc,3*np+ip)*z(4) + s(itc,4*np+ip)*z(5)
              alpha             =    alpha/( s(itc,ip)**2 + s(itc,np+ip)**2 + s(itc,2*np+ip)**2 + s(itc,3*np+ip)**2 + s(itc,4*np+ip)**2 )
              alpha             =    one - alpha
              s(itc+1,ip)       =    z(1)/alpha
              s(itc+1,np+ip)    =    z(2)/alpha
              s(itc+1,2*np+ip)  =    z(3)/alpha
              s(itc+1,3*np+ip)  =    z(4)/alpha 
              s(itc+1,4*np+ip)  =    z(5)/alpha 

              toll              = toll + s(itc+1,ip)**2 + s(itc+1,np+ip)**2 + s(itc+1,2*np+ip)**2 + s(itc+1,3*np+ip)**2 +&
                                         s(itc+1,4*np+ip)**2   
                                         
              if (periodicity_particles) call iperiodic_particles(ip,np,r)
                                       
              x(ip)%x           = r(ip)%x
              x(ip)%x(3)        = zero
              x(ip)%P           = r(ip)%data%v
              x(ip)%A           = r(ip)%results%A              

            enddo

            call MPI_ALLREDUCE(toll, tolg, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
!            tolg  = sqrt(tolg)/real(tnp,kind=kind_particle)


    enddo

    if (tolg > stop_tol)  ierr = 1
    

    do ip = 1,np

        e                 = p(ip)%data%q
        m                 = p(ip)%data%m

        p(ip)%data%v      = ( x(ip)%P - e/lorentz_tilde*x(ip)%A )/m
        p(ip)%data%g      = sqrt( one + dot_product( p(ip)%data%v/vtilde ,  p(ip)%data%v/vtilde ) )
        p(ip)%x           = x(ip)%x
        p(ip)%x(3)        = zero
        if (periodicity_particles) call iperiodic_particles(ip,np,p)
        
    enddo

    deallocate( r,x,s )


    end subroutine trapezoidal_broyden
    
    subroutine res_trapezoidal_broyden(np,dt,p,r,res)
    use module_tool   , only: cross_product
    use module_pepc

    implicit none
    include 'mpif.h'

    integer(kind_particle)          , intent(in)     :: np
    real(kind_particle)             , intent(in)     :: dt
    type(t_particle), allocatable   , intent(in)     :: p(:),r(:)

    real(kind_particle)                              :: e,m,rot(1:3),v(1:3),grad(1:3),Eirr(1:3)
    real(kind_particle), allocatable                 :: res(:,:)
    integer(kind_particle)                           :: ip,rc
    
    if (allocated(res)) deallocate(res)
    allocate( res(5,np), stat=rc )
    
    do ip = 1,np
        
        
        Eirr              = zero
        Eirr(1:2)         = half*( p(ip)%results%E(1:2) + r(ip)%results%E(1:2) )
        
        v                 = p(ip)%data%v + r(ip)%data%v
        rot               = cross_product( v , p(ip)%results%B + r(ip)%results%B )
        
        grad(1)           = v(1)*p(ip)%results%dxA(1) + v(2)*p(ip)%results%dyA(1)
        grad(2)           = v(1)*p(ip)%results%dxA(2) + v(2)*p(ip)%results%dyA(2)
        grad(3)           = v(1)*p(ip)%results%dxA(3) + v(2)*p(ip)%results%dyA(3)
              
        grad(1)           = grad(1) + v(1)*r(ip)%results%dxA(1) + v(2)*r(ip)%results%dyA(1)
        grad(2)           = grad(2) + v(1)*r(ip)%results%dxA(2) + v(2)*r(ip)%results%dyA(2)
        grad(3)           = grad(3) + v(1)*r(ip)%results%dxA(3) + v(2)*r(ip)%results%dyA(3)
        
        rot               = rot/lorentz_tilde/( r(ip)%data%g + p(ip)%data%g )
        grad              = grad/lorentz_tilde/( r(ip)%data%g + p(ip)%data%g )
        
        res(1:2,ip)       = p(ip)%x(1:2)      -  r(ip)%x(1:2) - dt*v(1:2)/( r(ip)%data%g + p(ip)%data%g )
        res(3:5,ip)       = p(ip)%data%v(1:3) -  r(ip)%data%v(1:3) - half*dt*( Eirr + rot + grad  )
                
    enddo
    
    
    end subroutine res_trapezoidal_broyden
    
    subroutine hamiltonian_boris(np,dt,p,ierr,itc,tolg)
    use module_tool   , only: cross_product
    use module_globals, only: root,my_rank,ixdim,ivdim,tnp,periodicity_particles
    use helper        , only: iperiodic_particles,inormalize
    use module_pepc
!!!!  This integrator does not work with relativistic velocities
    implicit none
    include 'mpif.h'

    integer(kind_particle)          , intent(in)     :: np
    integer(kind_particle)                           :: n
    real(kind_particle)             , intent(in)     :: dt
    type(t_particle), allocatable   , intent(inout)  :: p(:)
    integer(kind_particle)          , intent(out)    :: ierr,itc
    real(kind_particle)             , intent(out)    :: tolg

    real(kind_particle)                              :: Pn(1:3),Pold(1:3),P_(1:3),B(1:3),e,m,Bnorm,Ppar(1:3),Pnor(1:3),toll,&
                                                        Ef(1:3),Af(1:3),beta,x_(1:2),stop_tol,alpha,gamma,gamma_old,eta,    &
                                                        gradA(1:3),gradP(1:3),delta
    
    type(t_particle), allocatable                    :: r(:)
    type(pma_particle), allocatable                  :: x(:)
    integer(kind_particle)                           :: rc,ip,maxit

    if (allocated(r)) deallocate(r)
    if (allocated(x)) deallocate(x)
    allocate( r(np), stat=rc )
    allocate( x(np), stat=rc )

    ierr            = 0
    itc             = 0
    maxit           = 50
    tolg            = one
    stop_tol        = prec


    !!! Initial Guess
    do ip = 1,np

            e            =  p(ip)%data%q
            m            =  p(ip)%data%m

            x(ip)%P      =  p(ip)%data%v + half*dt*e/m*p(ip)%results%E !+ dt*e/m*( p(ip)%results%E + p(ip)%results%Et + cross_product(p(ip)%data%v/vtilde,p(ip)%results%B ) )
            x(ip)%x      =  p(ip)%x      
            x(ip)%x(3)   =  zero
            x(ip)%A      =  p(ip)%results%A
            x(ip)%gamma  =  one/sqrt( one - dot_product( x(ip)%P/vtilde, x(ip)%P/vtilde ) ) !one!
                
            r(ip)%label       = p(ip)%label
            r(ip)%data%q      = e
            r(ip)%data%m      = m
    enddo
    ! main iteration loop
    !

    do while(tolg .gt. stop_tol .and. itc .lt. maxit)

        itc = itc + 1

            do ip = 1,np

              r(ip)%data%v      = p(ip)%data%v!half*(x(ip)%v + p(ip)%data%v)!
              r(ip)%x           = half*(x(ip)%x + p(ip)%x)!x(ip)%x
              r(ip)%x(3)        = zero

              r(ip)%results%E   = zero
              r(ip)%results%pot = zero
              r(ip)%results%B   = zero
              r(ip)%results%A   = zero
              r(ip)%results%dxA = zero
              r(ip)%results%dyA = zero
              r(ip)%results%Jirr= zero
              r(ip)%results%J   = zero
              r(ip)%work        = one
              
              if (periodicity_particles) call iperiodic_particles(ip,np,r)

            enddo


        !!! Update fields

            call pepc_particleresults_clear(r)
            call pepc_grow_tree(r)
            call pepc_traverse_tree(r)
            call pepc_restore_particles(r)
            call pepc_timber_tree()

            toll = zero
            tolg = zero

            do ip = 1,np

              e                 = p(ip)%data%q
              m                 = p(ip)%data%m
              call  inormalize(ip,np,r)
              
              B                 = half*( p(ip)%results%B + p(ip)%results%B )!r(ip)%results%B!
              Af(1:3)           = half*( p(ip)%results%A + p(ip)%results%A )
              Ef(1:2)           = half*( p(ip)%results%E(1:2) + p(ip)%results%E(1:2) )
              Ef(3)             = zero
              
              Bnorm             = B(1)**2 + B(2)**2 + B(3)**2

              alpha             = half*e*dt
              beta              = half*e/m/lorentz_tilde/x(ip)%gamma*dt
              eta               = half/m/x(ip)%gamma*(e/lorentz_tilde)**2*dt
              delta             = half/m/x(ip)%gamma*e*dt
              
              gamma_old         = one/sqrt( one - dot_product( p(ip)%data%v/vtilde, p(ip)%data%v/vtilde )  ) !one!
              Pold              = gamma_old*m*p(ip)%data%v + e/lorentz_tilde*p(ip)%results%A
              
              gradA(1)          = p(ip)%results%A(1)*p(ip)%results%dxA(1) + p(ip)%results%A(2)*p(ip)%results%dyA(1) 
              gradA(2)          = p(ip)%results%A(1)*p(ip)%results%dxA(2) + p(ip)%results%A(2)*p(ip)%results%dyA(2)
              gradA(3)          = p(ip)%results%A(1)*p(ip)%results%dxA(3) + p(ip)%results%A(2)*p(ip)%results%dyA(3)
              
              gradP(1)          = Pold(1)*p(ip)%results%dxA(1) + Pold(2)*p(ip)%results%dyA(1) 
              gradP(2)          = Pold(1)*p(ip)%results%dxA(2) + Pold(2)*p(ip)%results%dyA(2)
              gradP(3)          = Pold(1)*p(ip)%results%dxA(3) + Pold(2)*p(ip)%results%dyA(3)
              
              gradP(1)          = gradP(1) + Pold(1)*r(ip)%results%dxA(1) + Pold(2)*r(ip)%results%dyA(1) 
              gradP(2)          = gradP(2) + Pold(1)*r(ip)%results%dxA(2) + Pold(2)*r(ip)%results%dyA(2)
              gradP(3)          = gradP(3) + Pold(1)*r(ip)%results%dxA(3) + Pold(2)*r(ip)%results%dyA(3)
                            
              P_                = Pold + ( alpha*Ef - eta*cross_product(Af,B ) - eta*gradA + half*delta*gradP )
              x_(1:2)           = p(ip)%x(1:2) - beta*Af(1:2) 
              
              Ppar              = delta*( dot_product(P_,B) )*B
              Pnor              = cross_product(P_,B)
              Pn                = (P_ + delta*( Ppar + Pnor ) )/( one + (delta)**2*Bnorm )
!              
!              Pn                = Pold + e/m/lorentz_tilde*dt*( half*gradP - e/lorentz_tilde*gradA  )!two*Pn - Pold
              Af(1:3)           = half*( r(ip)%results%A + p(ip)%results%A )
              
              
              r(ip)%data%v      = ( Pn - e/lorentz_tilde*Af )/m !half*x(ip)%P +  half*( Pn - e/lorentz_tilde*Af )/m
              gamma             = sqrt( one + dot_product( r(ip)%data%v/vtilde, r(ip)%data%v/vtilde ) ) !one!

              r(ip)%data%v      = r(ip)%data%v/gamma
              
              r(ip)%x(1:2)      = two*x_(1:2) - p(ip)%x(1:2) + dt/m/gamma*Pn(1:2)
              r(ip)%x(3)        = zero
              

              toll              = toll + dot_product(r(ip)%data%v-x(ip)%P,r(ip)%data%v-x(ip)%P) 
              toll              = toll + dot_product(r(ip)%x-x(ip)%x,r(ip)%x-x(ip)%x)
              toll              = toll + dot_product(r(ip)%results%A-x(ip)%A,r(ip)%results%A-x(ip)%A)
              
              
              if (periodicity_particles) call iperiodic_particles(ip,np,r)

              x(ip)%P           =  r(ip)%data%v
              x(ip)%x           =  r(ip)%x
              x(ip)%x(3)        =  zero 
              x(ip)%A           =  r(ip)%results%A
              x(ip)%gamma       =  gamma

            enddo

            call MPI_ALLREDUCE(toll, tolg, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
            tolg  = sqrt(tolg)/real(tnp,kind=kind_particle)


    enddo

    if (tolg > stop_tol)  ierr = 1
    

    do ip = 1,np

        p(ip)%data%v      = x(ip)%P
        p(ip)%x           = x(ip)%x
        p(ip)%x(3)        = zero

    enddo

    deallocate( r,x )


    end subroutine hamiltonian_boris



  end module module_integration
