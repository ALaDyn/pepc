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
   !use module_interaction_Specific_types
   implicit none
   save
   private

   ! shortcut notations
   real(kind_physics), parameter :: zero = 0._kind_physics
   real(kind_physics), parameter :: dotone = 0.1_kind_physics
   real(kind_physics), parameter :: half = 0.5_kind_physics
   real(kind_physics), parameter :: dotnine = 0.9_kind_physics
   real(kind_physics), parameter :: one = 1._kind_physics
   real(kind_physics), parameter :: tentominusthree = 1.E-3_kind_physics
   real(kind_physics), parameter :: tentominusfive = 1.E-5_kind_physics
   real(kind_physics), parameter :: tentominusseven = 1.E-7_kind_physics
   real(kind_physics), parameter :: tentominuseight = 1.E-8_kind_physics

   public dirder
   public gmres
   public nsolgm
   public givapp
   public linear

contains

   subroutine linear(n, dt, x_tn, x_tn1, res)

      implicit none
      integer(kind_particle), intent(in)        :: n
      real(kind_particle), intent(in)           :: dt
      type(t_particle), allocatable, intent(in) :: x_tn(:)
!    real(kind_particle),allocatable, intent(in)     :: x_tn(:)!,x_tn1(:)
!    real(kind_particle),allocatable, intent(out)    :: res(:)

      real(kind_particle), intent(in)           :: x_tn1(n)
      real(kind_particle), intent(out)          :: res(n)

      integer(kind_particle)                    :: i, rc

!    if ( allocated(res) )   deallocate(res)
!    allocate(res(n)   , stat=rc )

!    write(*,*) "linear", n, size(x_tn) , size(x_tn1) , size(res)
      do i = 1, n
         res(i) = 2.0 * x_tn1(i) + 1.0
      end do

   end subroutine

   subroutine givapp(dx, dy, cs, sn)
      implicit none
      real(kind_particle), intent(inout)     :: dx, dy
      real(kind_particle), intent(in)        :: cs, sn
      real(kind_particle)                    :: tmp

      tmp = cs * dx + sn * dy
      dy = -sn * dx + cs * dy
      dx = tmp

   end subroutine

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
      integer(kind_particle)                    :: n
      real(kind_particle), intent(in)           :: dt
      type(t_particle), allocatable, intent(in) :: x_tn(:)
!    real(kind_particle),allocatable, intent(in)     :: x_tn(:)!,x_tn1(:),w(:),f0(:)
!    real(kind_particle),allocatable, intent(out)    :: z(:)
      real(kind_particle), intent(in)           :: x_tn1(6 * np), w(6 * np), f0(6 * np)
      real(kind_particle), intent(out)          :: z(6 * np)
      external funz

!    real(kind_particle),allocatable                 :: f1(:),del(:)
      real(kind_particle)                       :: f1(6 * np), del(6 * np)
      real(kind_particle)                       :: epsnew, normx, normw, xs
      integer(kind_particle)                    :: rc

!    if ( allocated(z) )     deallocate(z)
!
!    allocate(z(n)   , stat=rc )
!    allocate(f1(n)  , stat=rc )
!    allocate(del(n) , stat=rc )
      ! Hardwired difference increment.
      n = 6 * np
      epsnew = tentominusseven

      ! scale the step
      normw = sqrt(dot_product(w, w))

      if (normw .eq. zero) then
         z(1:n) = zero
         return
      end if

      normx = sqrt(dot_product(x_tn1, x_tn1))
      xs = dot_product(w, x_tn1) / normw

      if (xs .ne. zero) epsnew = epsnew * max(abs(xs), one) * abs(xs) / xs
      epsnew = epsnew / normw
      if (normx .gt. zero) epsnew = epsnew * normx
      !
      ! del and f1 could share the same space if storage
      ! is more important than clarity
      !
      del(1:n) = x_tn1(1:n) + epsnew * w(1:n)

      call funz(np, dt, x_tn, del, f1)
      z(1:n) = (f1(1:n) - f0(1:n)) / epsnew

   end subroutine dirder

   subroutine gmres(np, dt, f0, f, x_tn, x_init, errtol, maxiter, phi, tot_iter)
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

      integer(kind_particle), intent(in)        :: np, maxiter
      integer(kind_particle)                    :: nkrylov
      real(kind_particle), intent(in)           :: dt
      type(t_particle), allocatable, intent(in) :: x_tn(:)
!      real(kind_particle), allocatable, intent(in)     :: x_tn(:)!,res_in(:),x_in(:)
      real(kind_particle), intent(in)           :: f0(6 * np), x_init(6 * np)
      integer(kind_particle), intent(out)       :: tot_iter
!      real(kind_particle), allocatable, intent(out)    :: phi(:)
      real(kind_particle), intent(out)          :: phi(6 * np)
      real(kind_particle)                       :: normb, errtol, hr, nu
!      real(kind_particle), allocatable                 :: residu(:),aq(:),wk1(:)
      real(kind_particle)                       :: r(6 * np), vv(6 * np), w(6 * np)

      real(kind_particle)                       :: vM(6 * np, maxiter + 1), rho
      real(kind_particle)                       :: s(maxiter + 1), h(maxiter + 1, maxiter), &
                                                   cn(maxiter + 1), sn(maxiter + 1), y(maxiter + 1)
      integer(kind_particle)                    :: ii, jj, iterp1, k
!      integer(kind_particle)                           :: rc = 0
      !
      external f

      !    Initialization

      nkrylov = 6 * np
      r(1:nkrylov) = -f0(1:nkrylov)
      vv(1:nkrylov) = zero
      w(1:nkrylov) = zero
      vM(1:nkrylov, 1:maxiter + 1) = zero
      s(1:maxiter + 1) = zero
      h(maxiter + 1, 1:maxiter) = zero
      cn(1:maxiter + 1) = zero
      sn(1:maxiter + 1) = zero
      y(1:maxiter + 1) = zero

      rho = sqrt(dot_product(f0, f0))
      normb = rho
      s(1) = rho
      if (normb .eq. zero) normb = one
      errtol = errtol * normb
      k = 1

      if (rho .lt. errtol) return
      vM(1:nkrylov, 1) = r(1:nkrylov) / rho

      do while ((rho .gt. errtol) .and. (k .lt. maxiter))

!!!        CALL DIRECTIONAL DERIVATIVE function
         vv(1:nkrylov) = vM(1:nkrylov, k)
         call dirder(np, dt, x_tn, x_init, vv, f, f0, w)
         vM(1:nkrylov, k + 1) = w(1:nkrylov)

         do jj = 1, k + 1

            h(jj, k) = zero

            do ii = 1, nkrylov
               h(jj, k) = h(jj, k) + vM(ii, jj) * vM(ii, k + 1)
            end do

            vM(1:nkrylov, k + 1) = vM(1:nkrylov, k + 1) - h(jj, k) * vM(1:nkrylov, jj)

         end do

         h(k + 1, k) = sqrt(dot_product(vM(1:nkrylov, k + 1), vM(1:nkrylov, k + 1)))

!!!      reorthogonalize??

         if (tentominusthree * h(k + 1, k) .le. tentominuseight) then

            w(1:nkrylov) = vM(1:nkrylov, k + 1)

            do jj = 1, k + 1
               hr = zero
               do ii = 1, nkrylov
                  hr = hr + vM(ii, jj) * vM(ii, k + 1)
               end do

               h(jj, k) = h(jj, k) + hr
               w(1:nkrylov) = w(1:nkrylov) - hr * vM(1:nkrylov, jj)

            end do

            vM(1:nkrylov, k + 1) = w(1:nkrylov)
            hr = sqrt(dot_product(w, w))
            h(k + 1, k) = hr
         end if

!!!!     WATCH OUT FOR HAPPY BREAKDOWN

         if (h(k + 1, k) .ne. zero) vM(1:nkrylov, k + 1) = vM(1:nkrylov, k + 1) / h(k + 1, k)

         if (k .gt. 1) then

            do jj = 1, k
               call givapp(h(jj + 1, k), h(jj, k), cn(jj), sn(jj))
            end do
            y(1:k + 1) = h(1:k + 1, k)

         end if

         nu = sqrt(h(k, k)**2 + h(k + 1, k)**2)

         if (nu .ne. zero) then

            cn(k) = h(k, k) / nu
            sn(k) = -h(k + 1, k) / nu
            h(k, k) = cn(k) * h(k, k) - sn(k) * h(k + 1, k)
            h(k + 1, k) = zero

            call givapp(s(k + 1), s(k), cn(k), sn(k))

         end if

!!!     UPDATE RESIDUAL NORM

         rho = abs(s(k))
         k = k + 1

      end do
      write (*, *) h
      k = k - 1
      y(k) = s(k) / h(k, k)

      do ii = k - 1, 1, -1
         hr = zero
         do jj = ii + 1, k + 1
            hr = hr + h(ii, jj) * y(jj)
         end do
         y(ii) = (s(ii) - hr) / h(ii, ii)
      end do

      !    compute new phi

      phi(1:nkrylov) = zero

      do ii = 1, k
         do jj = 1, nkrylov
            phi(jj) = phi(jj) + y(ii) * vM(jj, ii)
         end do
      end do

      tot_iter = k

   end subroutine gmres

   subroutine nsolgm(np, dt, x_tn, f, sol, ierr, itc, fnrm, static_gmres)

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
      implicit none

      integer(kind_particle), intent(in)        :: np
      integer(kind_particle)                    :: n
      real(kind_particle), intent(in)           :: dt
      type(t_particle), allocatable, intent(in) :: x_tn(:)
!    real(kind_particle), allocatable, intent(in)     :: x_tn(:)
      integer(kind_particle), intent(out)       :: ierr, itc
      real(kind_particle), intent(out)          :: fnrm
      real(kind_particle), intent(out)          :: sol(6 * np)
!    real(kind_particle), allocatable, intent(out)    :: sol(:)

!    real(kind_particle), allocatable                 :: step(:),f0(:)
      real(kind_particle)                       :: step(6 * np), f0(6 * np)
      integer(kind_particle)                    :: maxit, lmaxit, &
                                                   inner_it_count, gmkmax
      integer(kind_particle)                    :: rc, ip, jp
      real(kind_particle)                       :: atol, rtol, etamax, gamma, gmerrtol, &
                                                   fnrmo, rat, etaold, etanew, stop_tol, &
                                                   lf, static_gmres(3), tmp_min, tmp_max, tmp_mean

      external f

      !
      ! set the debug parameter, 1 turns display on, otherwise off
      !
      n = 6 * np
      !
      ! initialize it_hist, ierr, and set the iteration parameters
      !
      gamma = dotnine
      atol = tentominusfive
      rtol = tentominusfive
      maxit = 100
      lmaxit = 100
      etamax = tentominusfive

      ierr = 0

      gmerrtol = abs(etamax)
      gmkmax = lmaxit
      itc = 0
      !
      ! evaluate f at the initial iterate
      ! compute the stop tolerance

!    sol(1:n) = zero
      do ip = 1, np
         jp = (ip - 1) * 6
         sol(jp + 1:jp + 3) = x_tn(ip)%x
         lf = sqrt(1 - dot_product(x_tn(ip)%data%v, x_tn(ip)%data%v))
         sol(jp + 4:jp + 6) = x_tn(ip)%data%m * x_tn(ip)%data%v / lf
      end do

      call f(np, dt, x_tn, sol, f0)

      fnrm = sqrt(sum(f0**2) / dble(n))

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

         call gmres(np, dt, f0, f, x_tn, sol, gmerrtol, gmkmax, step, inner_it_count)
         write (*, *) " FATTO FATTO FATTO FATTO FATTO"

         static_gmres(1) = min(static_gmres(1), real(inner_it_count, kind=kind_particle))
         static_gmres(2) = static_gmres(2) + real(inner_it_count, kind=kind_particle)
         static_gmres(3) = max(static_gmres(3), real(inner_it_count, kind=kind_particle))
         sol(1:n) = sol(1:n) + step(1:n)

         call f(np, dt, x_tn, sol, f0)

         fnrm = sqrt(sum(f0**2) / dble(n))
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

   end subroutine nsolgm

end module
