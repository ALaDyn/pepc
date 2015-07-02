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

module newton_krylov

  !use module_pepc
  use module_pepc_types
  !use module_timings
  !use module_debug
  use module_pepc_kinds
  use module_shortcut
  !use module_interaction_Specific_types
  implicit none
  include 'mpif.h'
  save
  private

    ! shortcut notations
!  real(kind_physics), parameter :: zero             =  0._kind_physics
!  real(kind_physics), parameter :: dotone           =  0.1_kind_physics
!  real(kind_physics), parameter :: half             =  0.5_kind_physics
!  real(kind_physics), parameter :: dotnine          =  0.9_kind_physics
!  real(kind_physics), parameter :: one              =  1._kind_physics
!  real(kind_physics), parameter :: tentominusfive   =  1.E-5_kind_physics
!  real(kind_physics), parameter :: tentominusseven  =  1.E-7_kind_physics
!  real(kind_physics), parameter :: tentominuseight  =  1.E-8_kind_physics

  public dirder
  public gmres
  public nsolgm
  public linear


  contains

    subroutine linear(n,dt,x_tn,x_tn1,res)

    implicit none
    integer(kind_particle)         , intent(in)     :: n
    real(kind_particle)            , intent(in)     :: dt
    type(t_particle), allocatable  , intent(in)     :: x_tn(:)
!    real(kind_particle),allocatable, intent(in)     :: x_tn(:)!,x_tn1(:)
!    real(kind_particle),allocatable, intent(out)    :: res(:)

    real(kind_particle),             intent(in)     :: x_tn1(n)
    real(kind_particle),             intent(out)    :: res(n)

    integer(kind_particle)                          :: i,rc

!    if ( allocated(res) )   deallocate(res)
!    allocate(res(n)   , stat=rc )

!    write(*,*) "linear", n, size(x_tn) , size(x_tn1) , size(res)
    do i = 1,n
        res(i) = 2.0*x_tn1(i) + 1.0
    enddo

    end subroutine


    subroutine dirder(np,dt,x_tn,x_tn1,w,funz,f0,z)
    ! Finite difference directional derivative
    ! Approximate f'(x) w
    !
    ! C. T. Kelley, November 25, 1993
    !
    ! This code comes with no guarantee or warranty of any kind.
    !
    ! Translated in Fortran90 by G. Lapenta
    !
    ! function z = dirder(x,w,f,f0)
    !
    ! inputs:
    !           x, w = point and direction
    !           f = function
    !           f0 = f(x), in nonlinear iterations
    !                f(x) has usually been computed
    !                before the call to dirder
    implicit none
    integer(kind_particle)         , intent(in)     :: np
    integer(kind_particle)                          :: n
    real(kind_particle)            , intent(in)     :: dt
    type(t_particle), allocatable  , intent(in)     :: x_tn(:)
!    real(kind_particle),allocatable, intent(in)     :: x_tn(:)!,x_tn1(:),w(:),f0(:)
!    real(kind_particle),allocatable, intent(out)    :: z(:)
    real(kind_particle),             intent(in)     :: x_tn1(12*np),w(12*np),f0(12*np)
    real(kind_particle),             intent(out)    :: z(12*np)
    external                                           funz

!    real(kind_particle),allocatable                 :: f1(:),del(:)
    real(kind_particle)                             :: f1(12*np),del(12*np)
    real(kind_particle)                             :: epsnew,dot
    integer(kind_particle)                          :: rc

!    if ( allocated(z) )     deallocate(z)
!
!    allocate(z(n)   , stat=rc )
!    allocate(f1(n)  , stat=rc )
!    allocate(del(n) , stat=rc )
    ! Hardwired difference increment.
    n      = 12*np
    epsnew = tentominusseven

    ! scale the step
    dot = dot_product(w,w)

    if ( dot .eq. zero )  then
        z(1:n) = zero
        return
    endif

    epsnew = epsnew/dot
    dot = dot_product(x_tn1,x_tn1)
    if ( dot .gt. zero ) epsnew = epsnew*sqrt( dot )
    !
    ! del and f1 could share the same space if storage
    ! is more important than clarity
    !
    del(1:n) = x_tn1(1:n) + epsnew*w(1:n)

    call funz(np,dt,x_tn,del,f1)
    z(1:n) = (f1(1:n) - f0(1:n))/epsnew

!    deallocate(f1)
!    deallocate(del)

    end subroutine dirder


    subroutine gmres(np,dt, res_in, f,x_tn, x_in, errtol,itsub,phi,iter)
    !
    ! nkrylov= dimension of the vector
    ! residu= initial residue
    ! f= external procedure to be computed
    ! x_in= initial iterate
    ! error= error tolerance
    ! itsub= maximum number of iterations
    ! phi=final iterate
    ! iter=number of iterations
    !
    !
      implicit none

      integer(kind_particle)          , intent(in)     :: np,itsub
      integer(kind_particle)                           :: nkrylov
      real(kind_particle)             , intent(in)     :: dt
      type(t_particle), allocatable   , intent(in)     :: x_tn(:)
!      real(kind_particle), allocatable, intent(in)     :: x_tn(:)!,res_in(:),x_in(:)
      real(kind_particle),             intent(in)      :: res_in(12*np),x_in(12*np)
      integer(kind_particle)          , intent(out)    :: iter
!      real(kind_particle), allocatable, intent(out)    :: phi(:)
      real(kind_particle),             intent(out)     :: phi(12*np)
      real(kind_particle)                              :: rnorm,bnorm,errtol,rnorm_loc
!      real(kind_particle), allocatable                 :: residu(:),aq(:),wk1(:)
      real(kind_particle)                              :: residu(12*np),aq(12*np),wk1(12*np)

      real(kind_particle)                              :: q(12*np,itsub),dot,dot_loc,g1,g2,ws
      real(kind_particle)                              :: g(itsub+1,2),s(itsub+1),h(itsub+1,itsub+1)
      integer(kind_particle)                           :: k,ip,ii,jj,iterp1
      integer(kind_particle)                           :: rc = 0
    !
      external f


      nkrylov = 12*np
      bnorm   = sqrt(sum(res_in**2))
      errtol  = errtol*bnorm

      !     ZERO WORKING ARRAYS

      aq(1:nkrylov)           = zero
      q(1:nkrylov,1:itsub)    = zero
      phi(1:nkrylov)          = zero
      s(1:itsub+1)            = zero
      g(1:itsub+1,1:2)        = zero
      h(1:itsub+1,1:itsub+1)  = zero


      !     CALCULATE THE INITIAL RESIDUAL ERROR
      call dirder(np,dt,x_tn,x_in, aq, f, res_in,residu)
      residu(1:nkrylov) = residu(1:nkrylov) - res_in(1:nkrylov)

      dot_loc = dot_product(residu,residu) ! changed
      call MPI_ALLREDUCE(dot_loc           , dot       , 1, MPI_KIND_PARTICLE, MPI_SUM, MPI_COMM_WORLD, rc) ! changed

      dot    = sqrt(dot)  ! changed

      if( dot<errtol )  return

      s(1)           = dot
      q(1:nkrylov,1) = residu(1:nkrylov)/dot

    !     ****************************************************************

    !     begin gmres

    !     ****************************************************************


      iter  = 0
      rnorm = dot

!      do while ( ( iter .le. itsub ).and.( rnorm .gt. errtol ) )
      do iter = 1,itsub

!         iter = iter + 1
         iterp1         = iter + 1

         !     normalize direction
         aq(1:nkrylov) = q(1:nkrylov,iter)
         !     compute A times preconditioned q
         call dirder(np,dt,x_tn,x_in, aq, f, res_in,wk1)
         aq(1:nkrylov) = wk1(1:nkrylov)
         !    orthogonalize:

         do  k = 1,iter
            dot = dot_product( aq(1:nkrylov) , q(1:nkrylov,k) )
            h(k,iter) = dot
            aq(1:nkrylov) = aq(1:nkrylov) - dot*q(1:nkrylov,k)
         enddo

         dot_loc        = dot_product(aq,aq) ! changed
         call MPI_ALLREDUCE(dot_loc           , dot       , 1, MPI_KIND_PARTICLE, MPI_SUM, MPI_COMM_WORLD, rc) ! changed
         dot            = sqrt(dot)
         h(iterp1,iter) = dot

         !     apply previous Givens rotations to h (update QR factorization)

         do  k = 1,iter-1
            ws          = g(k,1)*h(k,iter) - g(k,2)*h(k+1,iter)
            h(k+1,iter) = g(k,2)*h(k,iter) + g(k,1)*h(k+1,iter)
            h(k,iter)   = ws
         enddo

         !     compute next Givens rotation
         g1         =  h(iter,iter)
         g2         =  h(iterp1,iter)
         ws         =  sqrt(g1*g1+g2*g2)
         g1         =  g1/ws
         g2         = -g2/ws
         g(iter,1)  = g1
         g(iter,2)  = g2

         !     apply g to h
         h(iter,iter)   = g1*h(iter,iter) - g2*h(iterp1,iter)
         h(iterp1,iter) = zero

         !     apply g to s
         ws             = g1*s(iter) - g2*s(iterp1)
         s(iterp1)      = g2*s(iter) + g1*s(iterp1)
         s(iter)        = ws

         !     |s(iter+1)| is the norm of the current residual
         !     check for convergence
         rnorm_loc      = abs(s(iterp1))                                                                           ! changed
         call MPI_ALLREDUCE(rnorm_loc           , rnorm       , 1, MPI_KIND_PARTICLE, MPI_SUM, MPI_COMM_WORLD, rc) ! changed
         if ( (rnorm .le.errtol ).or.(iter .eq. itsub) ) go to 2

         !     normalize next q

         q(1:nkrylov,iterp1) = aq(1:nkrylov)/dot

      enddo

2     continue
      !    update the solution
      !    solve h*y=sbar  (sbar contains elements 1,..,iter of s)
      !    store y in s
!      q(1:nkrylov,iterp1) = q(1:nkrylov,iterp1)*dot
      s(iter) = s(iter)/h(iter,iter)

      do ii=iter-1,1,-1
         ws = zero
         do jj=ii+1,iter
            ws = ws + h(ii,jj)*s(jj)
         enddo
         s(ii) = ( s(ii) - ws )/h(ii,ii)
      enddo

      !    compute new phi

      wk1 = zero

      do  ii=1,iter
         do  ip=1,nkrylov
            wk1(ip) = wk1(ip) + s(ii)*q(ip,ii)
         enddo
      enddo

      phi(1:nkrylov)  =   phi(1:nkrylov) + wk1(1:nkrylov)

    end subroutine gmres


    subroutine nsolgm(np,dt,x_tn,f,sol,ierr,itc,fnrm,static_gmres)

    ! Newton-GMRES locally convergent solver for f(x) = 0
    !
    ! Uses Eisenstat-Walker forcing term
    !
    ! C. T. Kelley, July 1, 1994
    !
    ! Translated into Fortran 90, G. Lapenta
    !
    ! This code comes with no guarantee or warranty of any kind.
    !
    ! function [sol, it_hist, ierr] = nsolgm(x,f,tol,parms)
    !
    ! inputs:
    !        initial iterate = x
    !        function = f
    !        tol = [atol, rtol] relative/absolute
    !            error tolerances for the nonlinear iteration
    !        parms = [maxit, maxitl, etamax]
    !            maxit = maxmium number of nonlinear iterations
    !                default = 40
    !            maxitl = maximum number of inner iterations
    !                default = 40
    !            |etamax| = Maximum error tolerance for residual in inner
    !                iteration. The inner iteration terminates
    !                when the relative linear residual is
    !                smaller than eta*| F(x_c) |. eta is determined
    !                by the modified Eisenstat-Walker formula if etamax > 0.
    !                If etamax < 0, then eta = |etamax| for the entire
    !                iteration.
    !                default: etamax=.9
    !
    ! output:
    !        sol = solution
    !        it_hist(maxit,3) = scaled l2 norms of nonlinear residuals
    !            for the iteration, number function evaluations,
    !            and number of steplength reductions
    !        ierr = 0 upon successful termination
    !        ierr = 1 if either after maxit iterations
    !             the termination criterion is not satsified
    !             or the ratio of successive nonlinear residuals
    !             exceeds 1. In this latter case, the iteration
    !             is terminted.
    !
    ! internal parameters:
    !       debug = turns on/off iteration statistics display as
    !               the iteration progresses
    !
    ! Requires fdgmres.m and givapp.m
    !
!    use module_global:
    implicit none
!    include 'mpif.h'

    integer(kind_particle)          , intent(in)     :: np
    integer(kind_particle)                           :: n
    real(kind_particle)             , intent(in)     :: dt
    type(t_particle), allocatable   , intent(in)     :: x_tn(:)
!    real(kind_particle), allocatable, intent(in)     :: x_tn(:)
    integer(kind_particle)          , intent(out)    :: ierr,itc
    real(kind_particle)             , intent(out)    :: fnrm
    real(kind_particle)             , intent(out)    :: sol(12*np)
!    real(kind_particle), allocatable, intent(out)    :: sol(:)

!    real(kind_particle), allocatable                 :: step(:),f0(:)
    real(kind_particle)                              :: step(12*np),f0(12*np)
    integer(kind_particle)                           :: maxit,lmaxit,debug,                 &
                                                        inner_it_count,gmkmax
    integer(kind_particle)                           :: rc,ip,jp
    real(kind_particle)                              :: atol,rtol,etamax,gamma,gmerrtol,        &
                                                        fnrmo,rat,etaold,etanew,stop_tol,       &
                                                        lf,static_gmres(3),tmp_min,tmp_max,     &
                                                        tmp_mean,fnrm_loc
!    character(255)                                   :: str_proc
    external f


    !
    ! set the debug parameter, 1 turns display on, otherwise off
    !
    n     = 12*np
    debug = 0
    !
    ! initialize it_hist, ierr, and set the iteration parameters
    !
    gamma   = dotnine
    atol    = tentominuseight
    rtol    = tentominuseight
    maxit   = 100
    lmaxit  = 100
    etamax  = tentominuseight

    ierr    = 0

    gmerrtol= abs(etamax)
    gmkmax  = lmaxit
    itc     = 0
    !
    ! evaluate f at the initial iterate
    ! compute the stop tolerance


    sol(1:n) = zero
!    do ip = 1,np
!        jp = (ip-1)*12
!        sol(jp+1:jp+3)      = x_tn(ip)%x
!        lf                  = one!sqrt( 1 - dot_product(x_tn(ip)%data%v,x_tn(ip)%data%v) )
!        sol(jp+4:jp+6)      = x_tn(ip)%data%m*x_tn(ip)%data%v/lf
!        sol(jp+7:jp+9)      = x_tn(ip)%results%A
!        sol(jp+10:jp+12)    = x_tn(ip)%results%B
!    enddo

    call f(np,dt,x_tn,sol,f0)

    fnrm_loc        =  sum(f0**2)/dble(n)                                                                       ! changed
    call MPI_ALLREDUCE(fnrm_loc           , fnrm       , 1, MPI_KIND_PARTICLE, MPI_SUM, MPI_COMM_WORLD, rc)     ! changed

    fnrm            = sqrt(fnrm)
    fnrmo           = one
    stop_tol        = atol + rtol*fnrm
    static_gmres(1) = lmaxit
    static_gmres(2) = zero
    static_gmres(3) = zero
    !
    ! main iteration loop
    !

    do while(fnrm .gt. stop_tol .and. itc .lt. maxit)
    !
    ! keep track of the ratio (rat = fnrm/frnmo)
    ! of successive residual norms and
    ! the iteration counter (itc)
    !
        rat     = fnrm/fnrmo
        fnrmo   = fnrm
        itc     = itc + 1


    ! compute the step using a GMRES routine especially designed
    ! for this purpose
    !
    !	call gmres(n, f0, f, x, &
    !     gmerrtol,gmkmax,step,inner_it_count)

    ! to use  Kelley's gmres uncomment
    !	call fdgmres(n,f0,f,x,gmerrtol,gmkmax,step,error,inner_it_count)
    ! to use my usual gmres uncomment
        call gmres(np,dt,f0,f,x_tn,sol,gmerrtol,gmkmax,step,inner_it_count)
!        write( str_proc , '(i10)' ) my_rank
!        open (unit=my_rank,file=trim("nk/nk_")//trim("_")//trim(adjustl(str_proc))//".dat",action="write",status="replace")
!        write (rc,*) itc,inner_it_count

        static_gmres(1) =  min( static_gmres(1) , real( inner_it_count , kind = kind_particle ) )
        static_gmres(2) =  static_gmres(2) + real( inner_it_count , kind = kind_particle )
        static_gmres(3) =  max( static_gmres(3) , real( inner_it_count , kind = kind_particle ) )


        sol(1:n)        =    sol(1:n) + step(1:n)

        call f(np,dt,x_tn,sol,f0)

        fnrm_loc = sum(f0**2) / dble(n)                                                                             ! changed
        call MPI_ALLREDUCE(fnrm_loc           , fnrm       , 1, MPI_KIND_PARTICLE, MPI_SUM, MPI_COMM_WORLD, rc)     ! changed

        fnrm = sqrt(fnrm)                                                                                           ! changed

        rat  =fnrm/fnrmo

    !   adjust eta
    !
        if ( etamax > zero ) then

            etaold = gmerrtol
            etanew = gamma*rat**2

            if ( gamma*etaold**2 > dotone)  etanew = max( etanew,gamma*etaold**2 )
            gmerrtol = min(etanew,etamax)
            gmerrtol = max(gmerrtol,half*stop_tol/fnrm)

        endif
        write(*,*) itc, fnrm , stop_tol

    enddo

    ! on failure, set the error flag
    !
    if (fnrm > stop_tol)  ierr = 1
    static_gmres(2) = static_gmres(2)/itc


    end subroutine nsolgm


end module
