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
   use module_globals, only: vtilde
   !use module_interaction_Specific_types
   use mpi
   implicit none
   save
   private

   public dirder
   public gmres
   public nsolgm

contains

   subroutine dirder(np, dt, x_tn, x_tn1, w, funz, f0, z)
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
      integer(kind_particle), intent(in)        :: np
      real(kind_particle), intent(in)           :: dt
      type(t_particle), allocatable, intent(in) :: x_tn(:)
!    real(kind_particle),allocatable, intent(in)     :: x_tn(:)!,x_tn1(:),w(:),f0(:)
!    real(kind_particle),allocatable, intent(out)    :: z(:)
      real(kind_particle), intent(in)           :: x_tn1(9 * np), w(9 * np), f0(9 * np)
      real(kind_particle), intent(out)          :: z(9 * np)
      external funz

      real(kind_particle)                       :: del(9 * np), epsnew, dot
      integer(kind_particle)                    :: n

      n = 9 * np
      epsnew = tentominusseven

      ! scale the step
      dot = dot_product(w, w)
      z(1:n) = zero

      if (dot .eq. zero) return

      epsnew = epsnew / dot
      dot = dot_product(x_tn1, x_tn1)
      if (dot .gt. zero) epsnew = epsnew * sqrt(dot)
      !
      ! del and f1 could share the same space if storage
      ! is more important than clarity
      !
      del(1:n) = x_tn1(1:n) + epsnew * w(1:n)

      call funz(np, dt, x_tn, del, z)
      z(1:n) = (z(1:n) - f0(1:n)) / epsnew

   end subroutine dirder

   subroutine gmres(np, dt, res_in, f, x_tn, x_in, errtol, itsub, phi, iter)
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
      use module_globals, only: my_rank
      implicit none

      integer(kind_particle), intent(in)        :: np, itsub
      integer(kind_particle)                    :: nkrylov
      real(kind_particle), intent(in)           :: dt
      type(t_particle), allocatable, intent(in) :: x_tn(:)
!      real(kind_particle), allocatable, intent(in)     :: x_tn(:)!,res_in(:),x_in(:)
      real(kind_particle), intent(in)           :: res_in(9 * np), x_in(9 * np)
      integer(kind_particle), intent(out)       :: iter
!      real(kind_particle), allocatable, intent(out)    :: phi(:)
      real(kind_particle), intent(out)          :: phi(9 * np)
      real(kind_particle)                       :: rnorm, bnorm, errtol, rnorm_loc
!      real(kind_particle), allocatable                 :: residu(:),aq(:),wk1(:)
      real(kind_particle)                       :: residu(9 * np), aq(9 * np), wk1(9 * np)

      real(kind_particle)                       :: q(9 * np, itsub), dot, dot_loc, g1, g2, ws
      real(kind_particle)                       :: g(itsub + 1, 2), s(itsub + 1), h(itsub + 1, itsub + 1)
      integer(kind_particle)                    :: k, ip, ii, jj, iterp1, ierr_MPI_B
      integer(kind_particle)                    :: rc = 0
      !
      external f

      nkrylov = 9 * np
      bnorm = sqrt(sum(res_in**2))
      errtol = errtol * bnorm

      !     ZERO WORKING ARRAYS

      aq(1:nkrylov) = zero
      q(1:nkrylov, 1:itsub) = zero
      phi(1:nkrylov) = zero
      s(1:itsub + 1) = zero
      g(1:itsub + 1, 1:2) = zero
      h(1:itsub + 1, 1:itsub + 1) = zero

      !     CALCULATE THE INITIAL RESIDUAL ERROR
      call dirder(np, dt, x_tn, x_in, aq, f, res_in, residu)
      residu(1:nkrylov) = residu(1:nkrylov) - res_in(1:nkrylov)

      dot_loc = dot_product(residu, residu) ! changed
      dot = zero
      call MPI_ALLREDUCE(dot_loc, dot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc) ! changed

      dot = sqrt(dot)  ! changed

      if (dot .lt. errtol) return

      s(1) = dot
      q(1:nkrylov, 1) = residu(1:nkrylov) / dot

      !     ****************************************************************

      !     begin gmres

      !     ****************************************************************

      iter = 0
      rnorm = dot

!      do while ( ( iter .le. itsub ).and.( rnorm .gt. errtol ) )
      do iter = 1, itsub

!         iter = iter + 1
         iterp1 = iter + 1

         !     normalize direction
         aq(1:nkrylov) = q(1:nkrylov, iter)
         !     compute A times preconditioned q
         call dirder(np, dt, x_tn, x_in, aq, f, res_in, wk1)
         aq(1:nkrylov) = wk1(1:nkrylov)
         !    orthogonalize:

         do k = 1, iter
            dot = dot_product(aq(1:nkrylov), q(1:nkrylov, k))
            h(k, iter) = dot
            aq(1:nkrylov) = aq(1:nkrylov) - dot * q(1:nkrylov, k)
         end do

!         dot        = dot_product(aq,aq)

         dot_loc = dot_product(aq, aq) ! changed
         call MPI_ALLREDUCE(dot_loc, dot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc) ! changed
         dot = sqrt(dot)
         h(iterp1, iter) = dot

         !     apply previous Givens rotations to h (update QR factorization)

         do k = 1, iter - 1
            ws = g(k, 1) * h(k, iter) - g(k, 2) * h(k + 1, iter)
            h(k + 1, iter) = g(k, 2) * h(k, iter) + g(k, 1) * h(k + 1, iter)
            h(k, iter) = ws
         end do

         !     compute next Givens rotation
         g1 = h(iter, iter)
         g2 = h(iterp1, iter)
         ws = sqrt(g1 * g1 + g2 * g2)
         g1 = g1 / ws
         g2 = -g2 / ws
         g(iter, 1) = g1
         g(iter, 2) = g2

         !     apply g to h
         h(iter, iter) = g1 * h(iter, iter) - g2 * h(iterp1, iter)
         h(iterp1, iter) = zero

         !     apply g to s
         ws = g1 * s(iter) - g2 * s(iterp1)
         s(iterp1) = g2 * s(iter) + g1 * s(iterp1)
         s(iter) = ws

         !     |s(iter+1)| is the norm of the current residual
         !     check for convergence
!         rnorm      = abs(s(iterp1))
         rnorm_loc = abs(s(iterp1))                                                                           ! changed
         call MPI_ALLREDUCE(rnorm_loc, rnorm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc) ! changed
         if ((rnorm .le. errtol) .or. (iter .eq. itsub)) go to 2

         !     normalize next q

         q(1:nkrylov, iterp1) = aq(1:nkrylov) / dot

      end do

2     continue
      !    update the solution
      !    solve h*y=sbar  (sbar contains elements 1,..,iter of s)
      !    store y in s
!      q(1:nkrylov,iterp1) = q(1:nkrylov,iterp1)*dot
      s(iter) = s(iter) / h(iter, iter)

      do ii = iter - 1, 1, -1
         ws = zero
         do jj = ii + 1, iter
            ws = ws + h(ii, jj) * s(jj)
         end do
         s(ii) = (s(ii) - ws) / h(ii, ii)
      end do

      !    compute new phi

      wk1 = zero

      do ii = 1, iter
         do ip = 1, nkrylov
            wk1(ip) = wk1(ip) + s(ii) * q(ip, ii)
         end do
      end do

      phi(1:nkrylov) = phi(1:nkrylov) + wk1(1:nkrylov)

   end subroutine gmres

   subroutine nsolgm(np, dt, x_tn, f, ierr, itc, fnrm, static_gmres)

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
      use module_tool, only: cross_product
      use module_globals, only: root, my_rank, ixdim, ivdim
      use mpi
      implicit none

      integer(kind_particle), intent(in)           :: np
      integer(kind_particle)                       :: n
      real(kind_particle), intent(in)              :: dt
      type(t_particle), allocatable, intent(inout) :: x_tn(:)
!    real(kind_particle), allocatable, intent(in)     :: x_tn(:)
      integer(kind_particle), intent(out)          :: ierr, itc
      real(kind_particle), intent(out)             :: fnrm
!    type(t_particle), allocatable   , intent(out)    :: particle(:)
      real(kind_particle)                          :: sol(9 * np)
!    real(kind_particle), allocatable, intent(out)    :: sol(:)

!    real(kind_particle), allocatable                 :: step(:),f0(:)
      real(kind_particle)                          :: step(9 * np), f0(9 * np)
      integer(kind_particle)                       :: maxit, lmaxit, debug, &
                                                      inner_it_count, gmkmax
      integer(kind_particle)                       :: rc, ip, jp, ierr_MPI_B
      real(kind_particle)                          :: atol, rtol, etamax, gamma, gmerrtol, &
                                                      fnrmo, rat, etaold, etanew, stop_tol, &
                                                      lf, static_gmres(3), tmp_min, tmp_max, &
                                                      tmp_mean, fnrm_loc, m, P(3), rot(3), e, v(3), &
                                                      B(3), A(3), Eirr(3), x(3), &
                                                      eta, u(3), uAh(3), h(3), s(3)
!    character(255)                                   :: str_proc
      external f

      !
      ! set the debug parameter, 1 turns display on, otherwise off
      !
      n = 9 * np
      debug = 0
      !
      ! initialize it_hist, ierr, and set the iteration parameters
      !
      gamma = dotnine
      atol = tentominuseight!tentominuseight
      rtol = tentominuseight!tentominuseight
      maxit = 100
      lmaxit = 100
      etamax = tentominuseight!tentominuseight

      ierr = 0

      gmerrtol = abs(etamax)
      gmkmax = lmaxit
      itc = 0
      !
      ! evaluate f at the initial iterate
      ! compute the stop tolerance

!    sol(1:n) = zero
      do ip = 1, np
         jp = (ip - 1) * 9

         rot = zero
         x = zero
         v = zero
         A = zero
         Eirr = zero
         B = zero

         x(1:ixdim) = x_tn(ip)%x(1:ixdim)
         Eirr(1:ixdim) = x_tn(ip)%results%E(1:ixdim)
         v(1:ivdim) = x_tn(ip)%data%v(1:ivdim)
         e = x_tn(ip)%data%q
         m = x_tn(ip)%data%m
         A = dotone * x_tn(ip)%results%A
         B = x_tn(ip)%results%B

         lf = one!one/sqrt( one - dot_product( v(1:ivdim) , v(1:ivdim) ) )
         rot = cross_product(v / vtilde, B)

         if (ivdim .eq. 1) then
            rot = zero
            A = zero
         end if

         eta = half * e * dt / m
         u = v + dotone * eta * x_tn(ip)%results%E
         h = dotone * eta * B
         s = two * h / (1.0_8 + dot_product(h, h))
         uAh = u + cross_product(u, h)
         uAh = cross_product(uAh, s)

         v = u + uAh + eta * x_tn(ip)%results%e
         x = x_tn(ip)%x + dotone * dt * x_tn(ip)%data%v

         sol(jp + 1:jp + 3) = x !+ dotone*dt*v
         sol(jp + 4:jp + 6) = m * v * lf !+ A !+ dotone*e*( Eirr +  rot )
         sol(jp + 7:jp + 9) = zero!(one+dotone)*A

      end do

      call f(np, dt, x_tn, sol, f0)

      fnrm_loc = sum(f0**2) / dble(n)                                                                       ! changed
      fnrm = zero
      call MPI_ALLREDUCE(fnrm_loc, fnrm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)     ! changed

      fnrm = sqrt(fnrm)
      fnrmo = one
      stop_tol = atol + rtol * fnrm
      static_gmres(1) = lmaxit
      static_gmres(2) = zero
      static_gmres(3) = zero
      !
      ! main iteration loop
      !

      do while (fnrm .gt. stop_tol .and. itc .lt. maxit)
         !
         ! keep track of the ratio (rat = fnrm/frnmo)
         ! of successive residual norms and
         ! the iteration counter (itc)
         !
!        write(*,*) "Inside Newton", fnrm, itc, my_rank
         rat = fnrm / fnrmo
         fnrmo = fnrm
         itc = itc + 1

         ! compute the step using a GMRES routine especially designed
         ! for this purpose
         !
         !        call gmres(n, f0, f, x, &
         !     gmerrtol,gmkmax,step,inner_it_count)

         ! to use  Kelley's gmres uncomment
         !        call fdgmres(n,f0,f,x,gmerrtol,gmkmax,step,error,inner_it_count)
         ! to use my usual gmres uncomment
!        call MPI_Barrier(MPI_COMM_WORLD, ierr_MPI_B)
         call gmres(np, dt, f0, f, x_tn, sol, gmerrtol, gmkmax, step, inner_it_count)
!        write( str_proc , '(i10)' ) my_rank
!        open (unit=my_rank,file=trim("nk/nk_")//trim("_")//trim(adjustl(str_proc))//".dat",action="write",status="replace")
!        write (rc,*) itc,inner_it_count

         static_gmres(1) = min(static_gmres(1), real(inner_it_count, kind=kind_particle))
         static_gmres(2) = static_gmres(2) + real(inner_it_count, kind=kind_particle)
         static_gmres(3) = max(static_gmres(3), real(inner_it_count, kind=kind_particle))

         sol(1:n) = sol(1:n) + step(1:n)

         call f(np, dt, x_tn, sol, f0)

!        fnrm = sum(f0**2) / dble(n)
         fnrm_loc = sum(f0**2) / dble(n)                                                                             ! changed
         call MPI_ALLREDUCE(fnrm_loc, fnrm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)     ! changed

         fnrm = sqrt(fnrm)                                                                                           ! changed

         rat = fnrm / fnrmo

         !   adjust eta
         !
         if (etamax .gt. zero) then

            etaold = gmerrtol
            etanew = gamma * rat**2

            if (gamma * etaold**2 .gt. dotone) etanew = max(etanew, gamma * etaold**2)
            gmerrtol = min(etanew, etamax)
            gmerrtol = max(gmerrtol, half * stop_tol / fnrm)

         end if

      end do

      ! on failure, set the error flag
      !
      if (fnrm .gt. stop_tol) ierr = 1
      static_gmres(2) = static_gmres(2) / itc

      if (root) write (*, *) itc, fnrm, stop_tol, ierr

      do ip = 1, np

         jp = (ip - 1) * 9

         x = zero
         v = zero
         Eirr = zero
         A = zero

         x(1:ixdim) = sol(jp + 1:jp + ixdim)
         P(1:ivdim) = sol(jp + 4:jp + 3 + ivdim)
         A(1:ivdim) = sol(jp + 7:jp + 6 + ivdim)

         e = x_tn(ip)%data%q
         m = x_tn(ip)%data%m

         if (ivdim .eq. 1) A = zero

!        P                          = ( P - e/vtilde*A )/m
         P = (P) / m

         lf = one!sqrt( one + dot_product( P(1:ivdim)/vtilde,P(1:ivdim)/vtilde ) )
         x_tn(ip)%data%v = P / lf
         x_tn(ip)%x = x !+ dt*( x_tn(ip)%data%v )

      end do

   end subroutine nsolgm

end module
