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
   use module_globals, only: vtilde, lorentz_tilde, folder
   use module_shortcut
   implicit none
   save
   private

   !  public iNKmidpoint_2D3V
   !  public NKmidpoint_2D3V
   !  public iNKtrapezoidal_2D3V
   !  public NKtrapezoidal_2D3V
   public trapezoidal_electrostatic
   !  public boris
   public leapfrog
   !  public midpoint
   public trapezoidal
   public implicit_nkPtr
   public implicit_picardPtr
   public explicit_ptr
   public march
   !  public inverse
   !  public boris

   type pv_particle
      real(kind_physics), dimension(1:3) :: x, v
      real(kind_physics)                 :: gamma
   end type pv_particle

   type pma_particle
      real(kind_physics), dimension(1:3) :: x, P, A, dxA, dyA
      real(kind_physics)                 :: gamma
   end type pma_particle

   type broyden_particle
      real(kind_physics), dimension(1:3) :: x, P, A, dxA, dyA
      real(kind_physics), dimension(1:5) :: res, res_old
      real(kind_physics)                 :: gamma, norm_res, error
   end type broyden_particle

   abstract interface

      subroutine explicit_ptr(np, dt, p)
         use module_pepc_types
         use module_pepc_kinds
         implicit none
         integer(kind_particle), intent(in) :: np
         real(kind_particle), intent(in)    :: dt
         type(t_particle), allocatable, intent(inout) :: p(:)    ! Actual values of the particles

      end subroutine explicit_ptr

      subroutine implicit_nkPtr(n, dt, p, r, s, res)
         use module_pepc_types
         use module_pepc_kinds
         implicit none

         integer(kind_particle), intent(in) :: n
         real(kind_particle), intent(in)    :: dt
         type(t_particle), intent(in)       :: p, r    ! Actual values of the particles
         real(kind_particle), intent(in)    :: s(n)
         real(kind_particle), intent(out)   :: res(n)       ! Residual

      end subroutine implicit_nkPtr

      subroutine implicit_nkPtr_update(np, n, dt, p, s, res)
         use module_pepc_types
         use module_pepc_kinds
         implicit none

         integer(kind_particle), intent(in) :: n, np
         real(kind_particle), intent(in)    :: dt
         type(t_particle), intent(in)       :: p(1:np)
         real(kind_particle), intent(in)    :: s(n, np)
         real(kind_particle), intent(out)   :: res(n, np)

      end subroutine implicit_nkPtr_update

      subroutine implicit_nkPtr_scheme(n, p, step, r)
         use module_pepc_types
         use module_pepc_kinds
         implicit none

         integer(kind_particle), intent(in) :: n
         type(t_particle), intent(in)       :: p
         real(kind_particle), intent(in)    :: step(1:n)
         type(t_particle), intent(out)      :: r    ! Actual values of the particles

      end subroutine implicit_nkPtr_scheme

      subroutine implicit_picardPtr(np, dt, p, itc, tolg)
         use module_pepc_types
         use module_pepc_kinds
         implicit none

         integer(kind_particle), intent(in)           :: np
         real(kind_particle), intent(in)              :: dt
         type(t_particle), allocatable, intent(inout) :: p(:)
         integer(kind_particle), intent(out)          :: itc
         real(kind_particle), intent(out)             :: tolg

      end subroutine implicit_picardPtr
   end interface

contains

   subroutine march(np, dt, p, ischeme, adv)
      !  use newton_krylov, only: nsolgm
      use module_globals, only: root
      implicit none
      integer(kind_particle), intent(in)           :: np
      character(255), intent(in)                   :: ischeme ! numerical scheme
      integer(kind_particle), intent(in)           :: adv ! referred to the time integrator: 0 explicit, 1 implicit, newton-kylov, 2 implicit, picard
      real(kind_particle), intent(in)              :: dt
      type(t_particle), allocatable, intent(inout) :: p(:)    ! Actual values of the particles

      procedure(implicit_picardPtr), pointer       :: pic_ptr => null()
      procedure(implicit_nkPtr), pointer           :: nk_ptr => null()
      procedure(implicit_nkPtr_update), pointer    :: nk_ptr_update => null()
      procedure(implicit_nkPtr_scheme), pointer    :: nk_ptr_scheme => null()
      procedure(explicit_ptr), pointer             :: exp_ptr => null()
      integer(kind_particle)                       :: errmsg_nk, iter                          ! Picard Iteration
      integer                                      :: rc
      real(kind_particle)                          :: errnk, static_gmres!,errmsg(1:3,1:5)        ! Global Error - Picard Iteration

      if (adv .eq. 0) then

         select case (ischeme)
         case ("leapfrog")
            if (root) write (*, '(a)') " ====== Leapfrog Integrator"
            exp_ptr => leapfrog
         case ("explicit")
            if (root) write (*, '(a)') " ====== Explicit Integrator"
            exp_ptr => explicit
         case ("euler_method")
            if (root) write (*, '(a)') " ====== Euler Method"
            exp_ptr => euler_method
         case ("euler_method3d")
            if (root) write (*, '(a)') " ====== Euler Method 3D"
            exp_ptr => euler_method3d
         case ("explicit_verlet")
            if (root) write (*, '(a)') " ====== Explicit Verlet"
            exp_ptr => explicit_verlet

         case default
            if (root) write (*, '(a)') " ====== Leapfrog Integrator"
            exp_ptr => leapfrog
         end select

         call exp_ptr(np, dt, p)

      elseif (adv .eq. 1) then
         select case (ischeme)
            !                case ("midpoint_picard")
            !                    if (root) write(*,'(a)')        " ====== Midpoint Picard"
            !                    pic_ptr => midpoint
         case ("trapezoidal_picard_electrostatic")
            if (root) write (*, '(a)') " ====== Trapezoidal Picard Electrostatic"
            pic_ptr => trapezoidal_electrostatic
         case ("trapezoidal_picard")
            if (root) write (*, '(a)') " ====== Trapezoidal_picard"
            pic_ptr => trapezoidal
            !                case ("trapezoidal_broyden2D3V")
            !                    if (root) write(*,'(a)')        " ====== Trapezoidal_Broyden"
            !                    pic_ptr => broyden_trapezoidal2D3V
            !                case ("mfnk_trapezoidal2D3V")
            !                    if (root) write(*,'(a)')        " ====== Modified Newton-Gmres"
            !                    pic_ptr => mfnk_trapezoidal2D3V
            !
            !                case ("boris")
            !                    if (root) write(*,'(a)')        " ====== Boris  Integrator"
            !                    pic_ptr => boris
            !                case ("hamiltonian_boris")
            !                    if (root) write(*,'(a)')        " ====== Boris - Hamiltonian - Integrator"
            !                    pic_ptr => hamiltonian_boris
            !                case default
            !                    if (root) write(*,'(a)')        " ====== Boris - Hamiltonian Integrator"
            !                    pic_ptr => hamiltonian_boris

         end select
         call pic_ptr(np, dt, p, iter, errnk)

         if (root) then
            !                    write(*,'(a,i4,2x,es12.4,es12.4,es12.4,es12.4,es12.4,es12.4)') " == Resual : ", iter,errnk,errmsg
            write (*, '(a,i4,2x,es12.4)') " == Resual : ", iter, errnk
            open (newunit=rc, file=trim(folder)//trim("residuo_")//trim(adjustl(ischeme))//".dat", form='formatted', status='unknown', position='append')
            write (rc, *) iter, errnk!(:,1),errmsg(:,2),errmsg(:,3),errmsg(:,4),errmsg(:,5)
            close (rc)
         end if

      elseif (adv .eq. 2) then

         if (root) write (*, '(a)') " ====== NK Solver is deactivated"
         call exit(1)
         !            select case (ischeme)
         !                case ("midpoint2VD")
         !                    if (root) write(*,'(a)')        " ====== Midpoint_Rule with NK solver"
        !!                    nk_ptr        => iNKmidpoint_2D3V
        !!                    nk_ptr_update =>  NKmidpoint_2D3V
        !!                    nk_ptr_scheme =>  NKmidpoint
         !                case ("trapezoidal2V3D")
         !                    if (root) write(*,'(a)')        " ====== Trapezoidal_Rule with NK solver"
         !                    nk_ptr_update =>  NKtrapezoidal_2D3V
         !                case default
        !!                    nk_ptr        => iNKmidpoint_2D3V
        !!                    nk_ptr_update =>  NKmidpoint_2D3V
        !!                    nk_ptr_scheme =>  NKmidpoint
         !                end select
         !
         !                call nsolgm(np,dt,p,nk_ptr_update,errmsg_nk,iter,errnk,static_gmres)
         !
         !                if (root) then
         !                    write(*,'(a,i4,2x,es12.4)') " == Resual : ", iter,errnk
         !                    open(unit=rc,file=trim(folder)//trim("residuo_")//trim(adjustl(ischeme))//".dat",form='formatted',status='unknown',position='append')
         !                    write(rc,*) errmsg_nk,iter,errnk
         !                    close (rc )
         !                endif

      elseif (adv .eq. -1) then
         if (root) write (*, '(a)') " ====== Wrong choice of time integrator"
         call exit(1)

      end if

   end subroutine march

   subroutine leapfrog(np, dt, p)
      use module_tool, only: cross_product
      use helper, only: iperiodic_particles
      use module_shortcut, only: one
      use module_globals, only: periodicity_particles, B0

      implicit none

      !------------------------------------------------------------------------------------------------------------------------------
      !
      !    This subroutine is meant to compute the advance in time with the well known leap-frog integration.
      !    Notice this integrator is very efficient for electrostatic problems, but it fails for Electro-Darwin formulation.
      !
      !    x^(t+1)      = x^(t) +   dt*v^(t+1/2)
      !    v^(t+1)      = v^(t) + q*dt*E^(t)
      !
      !    np represents the number of particles of the current process
      !
      !    dt is the time step
      !
      !    p is the array of particles at the old time step, "t".
      !
      !------------------------------------------------------------------------------------------------------------------------------

      integer(kind_particle), intent(in)    :: np
      type(t_particle), allocatable, intent(inout) :: p(:)
      real(kind_particle), intent(in)    :: dt

      integer(kind_particle) :: ip
      real(kind_particle)    :: vxb(1:3), v(1:3)

      do ip = 1, np

         v = p(ip)%data%v / p(ip)%data%g
         vxb = cross_product(v, B0)

         p(ip)%data%v = p(ip)%data%v + dt * p(ip)%data%q / p(ip)%data%m * (p(ip)%results%e + vxb)
         !      p(ip)%data%v      = p(ip)%data%v + dt * p(ip)%data%q / p(ip)%data%m * ( vxb )
         p(ip)%data%g = sqrt(one + sum((p(ip)%data%v / vtilde)**2))

         v = p(ip)%data%v / p(ip)%data%g

         p(ip)%x = p(ip)%x + dt * v

         p(ip)%x(3) = zero
         if (periodicity_particles) call iperiodic_particles(p(ip))

      end do

   end subroutine leapfrog

   subroutine explicit(np, dt, p)
      use helper, only: iperiodic_particles
      use module_shortcut, only: one
      use module_globals, only: periodicity_particles, pold

      implicit none

      !------------------------------------------------------------------------------------------------------------------------------
      !
      !    This subroutine is meant to compute the advance in time with the well known leap-frog integration.
      !    Notice this integrator is very efficient for electrostatic problems, but it fails for Electro-Darwin formulation.
      !
      !    x^(t+1)      = x^(t) +   dt*v^(t+1/2)
      !    v^(t+1)      = v^(t) + q*dt*E^(t)
      !
      !    np represents the number of particles of the current process
      !
      !    dt is the time step
      !
      !    p is the array of particles at the old time step, "t".
      !
      !------------------------------------------------------------------------------------------------------------------------------

      integer(kind_particle), intent(in)    :: np
      type(t_particle), allocatable, intent(inout) :: p(:)
      real(kind_particle), intent(in)    :: dt

      integer(kind_particle) :: ip
      real(kind_particle)    :: grad(1:3), v(1:3), m, e

      do ip = 1, np

         v = p(ip)%data%v
         grad = zero
         grad(1) = v(1) * p(ip)%results%dxA(1) + v(2) * p(ip)%results%dxA(2) + v(3) * p(ip)%results%dxA(3)
         grad(2) = v(1) * p(ip)%results%dyA(2) + v(2) * p(ip)%results%dyA(2) + v(3) * p(ip)%results%dyA(3)
         grad = grad / lorentz_tilde

         m = p(ip)%data%m
         e = p(ip)%data%q

         p(ip)%data%v = m * p(ip)%data%v + e / lorentz_tilde * pold(ip)%results%A + dt * e * (p(ip)%results%e + grad)

         p(ip)%data%v = (m * p(ip)%data%v - e / lorentz_tilde * p(ip)%results%A) / m

         p(ip)%data%g = sqrt(one + sum((p(ip)%data%v / vtilde)**2))

         p(ip)%x = p(ip)%x + dt * p(ip)%data%v / p(ip)%data%g
         p(ip)%x(3) = zero
         if (periodicity_particles) call iperiodic_particles(p(ip))

      end do

   end subroutine explicit

   subroutine trapezoidal_electrostatic(np, dt, p, itc, toll)
      use module_tool, only: cross_product
      use module_globals, only: periodicity_particles, ischeme, root, my_rank, folder, step, pold, newmark_x, newmark_v, newmark_Es, &
                                newmark_Ei, newmark_B, newmark_g, dA_1, dA__1, dA_0, poldold, B0
      use helper, only: iperiodic_particles, inormalize, write_particles_ascii, write_particles_vtk
      use module_pepc
      use mpi
      implicit none

      !------------------------------------------------------------------------------------------------------------------------------
      !
      !    Ref.: Leimkuhler B., Reich S. -Simulating Hamiltonian dynamics
      !
      !    This subroutine is meant to push in time the particles. We implement here  the time integration scheme (Trapezoidal rule):
      !
      !    x^(t+1)      = x^(t) +   dt*v^(t+1/2)
      !    P^(t+1)      = P^(t) + q*dt*( E^(t+1/2) + grad( A^(t+1/2) . v^(t+1/2) ) )
      !
      !    Since the Trapezoidal rule is an implicit scheme, the nonlinearity is resolved with a Picard iteration.
      !    Picard method is very simple to implement, compared to Newton-Krylov, but it is not guaranteed the method converge toward
      !    the solution. In fact, it is, required the function  (the map which advance in time) is a contraction in the phase space.
      !
      !    np represents the number of particles in the current process
      !
      !    dt is the time step
      !
      !    p is the array of particles at the old time step, "t".
      !
      !
      !------------------------------------------------------------------------------------------------------------------------------
      integer(kind_particle), intent(in)     :: np
      real(kind_particle), intent(in)     :: dt
      type(t_particle), allocatable, intent(inout)  :: p(:)
      integer(kind_particle), intent(out)    :: itc
      real(kind_particle), intent(out)    :: toll

      real(kind_particle)                              :: e, m, stop_tol, v(1:3), grad(1:3), vnablaA(1:3), Eirr(1:3), &
                                                          coef1, coef2, vxb(1:3)

      !    real(kind_particle)                              :: err_A(1:3),err_dxA(1:3),err_dyA(1:3),err_x(1:3)
      real(kind_particle)                              :: err_x, err_v, rat, toll_old
      type(t_particle), allocatable                    :: r(:)!,error(:)
      type(pma_particle), allocatable                  :: s(:)
      integer(kind_particle)                           :: ip, j, maxit
      integer                                          :: rc, filehandle
      character(100)                                   :: filename
      integer                                          :: iteration

      if (allocated(r)) deallocate (r)
      if (allocated(s)) deallocate (s)
      !    if (allocated(error)) deallocate(error)
      allocate (r(np), stat=rc)
      allocate (s(np), stat=rc)
      !    allocate( error(np), stat=rc )

      itc = 0
      maxit = 10
      toll = one
      stop_tol = tentominusnine!prec!tentominusseven!prec

      !    alpha           = theta_Ar      ! theta parameter - A        - in the rhs: grad(A*v) 1/2 works ok. if alpha = 0 => explicit
      !    beta            = theta_Al      ! theta parameter - A        - in the canonical momentum. For large A beta > 0 is recommended.
      !    eta             = theta_vr  ! theta parameter - velocity - in the rhs: grad(A*v) 1/2 works ok. if eta = 0 => explicit
      !    theta           = theta_vl  ! theta parameter - velocity - in the canonical momentum.

      !------------------------------------------------------------------------------------------------------------------------------
      !        Here, we define the initial guess, all iterative methods need an initial guess, supposed to be close enough to the solution.
      !        We assume the guess solution is the configuration at old time t.
      !------------------------------------------------------------------------------------------------------------------------------

      do ip = 1, np

         e = p(ip)%data%q
         m = p(ip)%data%m

         s(ip)%x = p(ip)%x
         s(ip)%x(3) = zero
         s(ip)%P = p(ip)%data%v

         r(ip)%label = p(ip)%label
         r(ip)%data%q = e
         r(ip)%data%m = m
      end do

      !------------------------------------------------------------------------------------------------------------------------------
      !
      !                            main iteration loop
      !
      !------------------------------------------------------------------------------------------------------------------------------

      do while (toll .gt. stop_tol .and. itc .lt. maxit)

         itc = itc + 1

         do ip = 1, np
            !------------------------------------------------------------------------------------------------------------------------------
            !        At each iteration we update the canonical momentums, according to the Trapezoidal rule.
            !        This method assumes the fields are computed at the future values:
            !
            !                x^(t+1)
            !                P^(t+1) = m*gamma*v^(t+1) + q/c*A^(t+1)
            !
            !        We slightly modify the definition of velocity, we replace the natural definition of   v^(t+1/2) = ( v^(t) + v^(t+1) )/2
            !        with  v^(t+1/2) = ( v^(t) + v^(t+1) )/( gamma^(t) + gamma^(t+1) ).
            !
            !        It turns out that for non relativistic velocities the two definitions are equivalent.
            !
            !        So, the fields are evaluated as average values, for instance E^(t+1/2) = ( E( x^(t+1) ) + E( x^(t) ) )/2, etc.
            !        Notice:  s(ip)%P is the canonical momentum at time t+1, but regarding the previous iteration of the Picard method.
            !        In this piece of code s(ip) contains the information of the particle at time time t+1.
            !        We are redefining the particle position r(ip) as explained before.
            !        Be aware the fields are calculated at the phase space point x, gamma*v
            !------------------------------------------------------------------------------------------------------------------------------

            e = p(ip)%data%q
            m = p(ip)%data%m

            r(ip)%data%v = newmark_v * s(ip)%P + (one - newmark_v) * p(ip)%data%v

            r(ip)%data%g = sqrt(one + sum((r(ip)%data%v / vtilde)**2))
            r(ip)%x = newmark_x * s(ip)%x + (one - newmark_x) * p(ip)%x
            r(ip)%x(3) = zero

            r(ip)%label = p(ip)%label
            r(ip)%results%E = zero
            r(ip)%results%pot = zero
            r(ip)%results%B = zero
            r(ip)%results%A = zero
            r(ip)%results%dxA = zero
            r(ip)%results%dyA = zero
            r(ip)%results%Jirr = zero
            r(ip)%results%J = zero
            r(ip)%work = one

            if (periodicity_particles) call iperiodic_particles(r(ip))

         end do

         !------------------------------------------------------------------------------------------------------------------------------
         !            The fields are, as written above, evaluated at the future predicted phase space point, t+1.
         !------------------------------------------------------------------------------------------------------------------------------

         call pepc_particleresults_clear(r)
         call pepc_grow_tree(r)
         call pepc_traverse_tree(r)
         call pepc_restore_particles(r)
         call pepc_timber_tree()

         !            do ip = 1,np
         !                call  inormalize(r(ip))
         !            enddo
         !
         !            iteration = step*10+itc
         !            call write_particles_ascii(iteration, r)

         toll_old = toll
         toll = zero
         err_x = zero
         err_v = zero

         do ip = 1, np

            !------------------------------------------------------------------------------------------------------------------------------
            !            Once the fields are updated, we push the particle in time according to the Trapezoidal rule:
            !
            !              x^(t+1) = x^(t) +   dt*v^(t+1/2)
            !              P^(t+1) = m*gamma*v^(t+1) + q/c*A^(t+1)
            !              P^(t+1) = P^(t) + q*dt*(E^(t+1/2) + q/c*grad( A^(t+1/2).v^(t+1/2) )
            !------------------------------------------------------------------------------------------------------------------------------

            e = p(ip)%data%q
            m = p(ip)%data%m

            call inormalize(r(ip))

            v = (newmark_g * r(ip)%data%v + (one - newmark_g) * p(ip)%data%v) / (newmark_g * r(ip)%data%g + (one - newmark_g) * p(ip)%data%g)

            grad = zero
            !              vnablaA           = zero
            !              gradE             = zero
            !              nabla2A           = zero
            Eirr = zero

            !              gradE(1)          = v(1)*r(ip)%results%dxE(1) + v(2)*r(ip)%results%dxE(2)
            !              gradE(2)          = v(1)*r(ip)%results%dyE(1) + v(2)*r(ip)%results%dyE(2)

            !              gradE             = gradE/lorentz_tilde

            Eirr = newmark_Ei * r(ip)%results%E + (one - newmark_Ei) * p(ip)%results%E

            !              coef1             = -    dot_product( Eirr , v )/p(ip)%data%g**2/lorentz_tilde**2
            !              coef2             = -two*dot_product( grad , v )/p(ip)%data%g**2/lorentz_tilde**2

            vxb = cross_product(v, B0)
            !

            r(ip)%x = p(ip)%x + dt * v            !                                                       &
            !                                + half*dt**2*e*( m*p(ip)%data%v + e/lorentz_tilde*p(ip)%results%A )*coef1                 &
            !                                + half*dt**2*e*( m*p(ip)%data%v + e/lorentz_tilde*p(ip)%results%A )*coef2                 &
            !                                + half*dt**2*e*( grad + vnablaA + Eirr  )/m/p(ip)%data%g

            r(ip)%data%v = p(ip)%data%v + dt * e / m * (Eirr + vxb) !                     + & ! Canonical Momentum
            !                                       - half*dt**2*e/m*( gradE + nabla2A     )

            r(ip)%x(3) = zero

            if (periodicity_particles) call iperiodic_particles(r(ip))

            err_x = sum((s(ip)%x(1:2) - r(ip)%x(1:2))**2)
            err_v = sum((s(ip)%P(1:3) - r(ip)%data%v(1:3))**2)

            toll = max(toll, sqrt(err_x + err_v))

         end do

         call MPI_ALLREDUCE(MPI_IN_PLACE, toll, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, rc)
         rat = toll / toll_old

         if (.false. .and. (rat .gt. one) .and. (itc .eq. maxit)) then !( itc .eq. maxit ) then !

            if (root) then
               write (*, '(a,es12.4,es12.4)') " == Picard iteration is not converging : ", toll, toll_old
            end if

            write (filename, '(a,"errors_",i6.6,".dat")') trim(folder), my_rank
            open (newunit=filehandle, file=trim(filename), STATUS='REPLACE')

            do ip = 1, np

               err_x = sum((s(ip)%x(1:2) - r(ip)%x(1:2))**2)
               err_v = sum((s(ip)%P(1:3) - r(ip)%data%v(1:3))**2)

               write (filehandle, *) ip, err_x, err_v

            end do

            close (filehandle)
         end if

         !            endif

         do ip = 1, np

            s(ip)%P = r(ip)%data%v
            s(ip)%x = r(ip)%x
            s(ip)%x(3) = zero

         end do

      end do

      do ip = 1, np

         e = p(ip)%data%q
         m = p(ip)%data%m

         v = p(ip)%data%v
         p(ip)%data%v = s(ip)%P

         !        if ( p(ip)%label .eq. 1 ) p(ip)%data%v = zero
         p(ip)%data%g = sqrt(one + sum((p(ip)%data%v / vtilde)**2))
         p(ip)%x = s(ip)%x
         p(ip)%x(3) = zero
         if (periodicity_particles) call iperiodic_particles(p(ip))

      end do

      deallocate (r, s)!,error )

   end subroutine trapezoidal_electrostatic

   subroutine trapezoidal(np, dt, p, itc, toll)
      use module_tool, only: cross_product
      use module_globals, only: periodicity_particles, ischeme, root, my_rank, folder, step, pold, newmark_x, newmark_v, newmark_Es, &
                                newmark_Ei, newmark_B, newmark_g, dA_1, dA__1, dA_0, poldold, B0
      use helper, only: iperiodic_particles, inormalize, write_particles_ascii, write_particles_vtk
      use module_pepc
      use mpi
      implicit none

      !------------------------------------------------------------------------------------------------------------------------------
      !
      !    Ref.: Leimkuhler B., Reich S. -Simulating Hamiltonian dynamics
      !
      !    This subroutine is meant to push in time the particles. We implement here  the time integration scheme (Trapezoidal rule):
      !
      !    x^(t+1)      = x^(t) +   dt*v^(t+1/2)
      !    P^(t+1)      = P^(t) + q*dt*( E^(t+1/2) + grad( A^(t+1/2) . v^(t+1/2) ) )
      !
      !    Since the Trapezoidal rule is an implicit scheme, the nonlinearity is resolved with a Picard iteration.
      !    Picard method is very simple to implement, compared to Newton-Krylov, but it is not guaranteed the method converge toward
      !    the solution. In fact, it is, required the function  (the map which advance in time) is a contraction in the phase space.
      !
      !    np represents the number of particles in the current process
      !
      !    dt is the time step
      !
      !    p is the array of particles at the old time step, "t".
      !
      !
      !------------------------------------------------------------------------------------------------------------------------------
      integer(kind_particle), intent(in)           :: np
      real(kind_particle), intent(in)              :: dt
      type(t_particle), allocatable, intent(inout) :: p(:)
      integer(kind_particle), intent(out)          :: itc
      real(kind_particle), intent(out)             :: toll

      real(kind_particle)                          :: e, m, stop_tol, v(1:3), grad(1:3), vnablaA(1:3), Aold(1:3), Anew(1:3), Eirr(1:3), &
                                                      Esol(1:3), gradE(1:3), nabla2A(1:3), vxb(1:3), coef1, coef2

      !    real(kind_particle)                          :: err_A(1:3),err_dxA(1:3),err_dyA(1:3),err_x(1:3)
      real(kind_particle)                          :: err_A, err_dxA, err_dyA, err_x, err_v, rat, toll_old
      type(t_particle), allocatable                :: r(:)!,error(:)
      type(pma_particle), allocatable              :: s(:)
      integer(kind_particle)                       :: ip, j, maxit
      integer                                      :: rc, filehandle
      character(100)                               :: filename
      integer                                      :: iteration

      if (allocated(r)) deallocate (r)
      if (allocated(s)) deallocate (s)
      !    if (allocated(error)) deallocate(error)
      allocate (r(np), stat=rc)
      allocate (s(np), stat=rc)
      !    allocate( error(np), stat=rc )

      itc = 0
      maxit = 10
      toll = one
      stop_tol = tentominusnine!prec!tentominusseven!prec

      !    alpha           = theta_Ar      ! theta parameter - A        - in the rhs: grad(A*v) 1/2 works ok. if alpha = 0 => explicit
      !    beta            = theta_Al      ! theta parameter - A        - in the canonical momentum. For large A beta > 0 is recommended.
      !    eta             = theta_vr  ! theta parameter - velocity - in the rhs: grad(A*v) 1/2 works ok. if eta = 0 => explicit
      !    theta           = theta_vl  ! theta parameter - velocity - in the canonical momentum.

      !------------------------------------------------------------------------------------------------------------------------------
      !        Here, we define the initial guess, all iterative methods need an initial guess, supposed to be close enough to the solution.
      !        We assume the guess solution is the configuration at old time t.
      !------------------------------------------------------------------------------------------------------------------------------

      do ip = 1, np

         e = p(ip)%data%q
         m = p(ip)%data%m

         s(ip)%x = p(ip)%x
         !            s(ip)%x(3)   =  zero
         s(ip)%A = p(ip)%results%A
         s(ip)%P = p(ip)%data%v
         s(ip)%dxA = p(ip)%results%dxA
         s(ip)%dyA = p(ip)%results%dyA

         r(ip)%label = p(ip)%label
         r(ip)%data%q = e
         r(ip)%data%m = m
      end do

      !------------------------------------------------------------------------------------------------------------------------------
      !
      !                            main iteration loop
      !
      !------------------------------------------------------------------------------------------------------------------------------

      do while (toll .gt. stop_tol .and. itc .lt. maxit)

         itc = itc + 1

         do ip = 1, np
            !------------------------------------------------------------------------------------------------------------------------------
            !        At each iteration we update the canonical momentums, according to the Trapezoidal rule.
            !        This method assumes the fields are computed at the future values:
            !
            !                x^(t+1)
            !                P^(t+1) = m*gamma*v^(t+1) + q/c*A^(t+1)
            !
            !        We slightly modify the definition of velocity, we replace the natural definition of   v^(t+1/2) = ( v^(t) + v^(t+1) )/2
            !        with  v^(t+1/2) = ( v^(t) + v^(t+1) )/( gamma^(t) + gamma^(t+1) ).
            !
            !        It turns out that for non relativistic velocities the two definitions are equivalent.
            !
            !        So, the fields are evaluated as average values, for instance E^(t+1/2) = ( E( x^(t+1) ) + E( x^(t) ) )/2, etc.
            !        Notice:  s(ip)%P is the canonical momentum at time t+1, but regarding the previous iteration of the Picard method.
            !        In this piece of code s(ip) contains the information of the particle at time time t+1.
            !        We are redefining the particle position r(ip) as explained before.
            !        Be aware the fields are calculated at the phase space point x, gamma*v
            !------------------------------------------------------------------------------------------------------------------------------

            e = p(ip)%data%q
            m = p(ip)%data%m

            r(ip)%data%v = newmark_v * s(ip)%P + (one - newmark_v) * p(ip)%data%v

            r(ip)%data%g = sqrt(one + sum((r(ip)%data%v / vtilde)**2))
            r(ip)%x = newmark_x * s(ip)%x + (one - newmark_x) * p(ip)%x
            !              r(ip)%x(3)        = zero

            r(ip)%label = p(ip)%label
            r(ip)%results%E = zero
            r(ip)%results%pot = zero
            r(ip)%results%B = zero
            r(ip)%results%A = zero
            r(ip)%results%dxA = zero
            r(ip)%results%dyA = zero
            r(ip)%results%Jirr = zero
            r(ip)%results%J = zero
            r(ip)%work = one

            if (periodicity_particles) call iperiodic_particles(r(ip))

         end do

         !------------------------------------------------------------------------------------------------------------------------------
         !            The fields are, as written above, evaluated at the future predicted phase space point, t+1.
         !------------------------------------------------------------------------------------------------------------------------------

         call pepc_particleresults_clear(r)
         call pepc_grow_tree(r)
         call pepc_traverse_tree(r)
         call pepc_restore_particles(r)
         call pepc_timber_tree()

         !            do ip = 1,np
         !                call  inormalize(r(ip))
         !            enddo
         !
         !            iteration = step*10+itc
         !            call write_particles_ascii(iteration, r)

         toll_old = toll
         toll = zero
         err_x = zero
         err_v = zero
         err_A = zero
         err_dxA = zero
         err_dyA = zero

         do ip = 1, np

            !------------------------------------------------------------------------------------------------------------------------------
            !            Once the fields are updated, we push the particle in time according to the Trapezoidal rule:
            !
            !              x^(t+1) = x^(t) +   dt*v^(t+1/2)
            !              P^(t+1) = m*gamma*v^(t+1) + q/c*A^(t+1)
            !              P^(t+1) = P^(t) + q*dt*(E^(t+1/2) + q/c*grad( A^(t+1/2).v^(t+1/2) )
            !------------------------------------------------------------------------------------------------------------------------------

            e = p(ip)%data%q
            m = p(ip)%data%m

            call inormalize(r(ip))

            v = (newmark_g * r(ip)%data%v + (one - newmark_g) * p(ip)%data%v) / (newmark_g * r(ip)%data%g + (one - newmark_g) * p(ip)%data%g)

            grad = zero
            !              vnablaA           = zero
            !              gradE             = zero
            !              nabla2A           = zero
            Eirr = zero
            Esol = zero

            grad(1) = v(1) * r(ip)%results%dxA(1) + v(2) * r(ip)%results%dxA(2) + v(3) * r(ip)%results%dxA(3)
            grad(2) = v(1) * r(ip)%results%dyA(1) + v(2) * r(ip)%results%dyA(2) + v(3) * r(ip)%results%dyA(3)
            grad(2) = v(1) * r(ip)%results%dzA(1) + v(2) * r(ip)%results%dzA(2) + v(3) * r(ip)%results%dzA(3)

            grad = newmark_B * grad
            grad(1) = grad(1) + (one - newmark_B) * (v(1) * p(ip)%results%dxA(1) + v(2) * p(ip)%results%dxA(2) + v(3) * p(ip)%results%dxA(3))
            grad(2) = grad(2) + (one - newmark_B) * (v(1) * p(ip)%results%dyA(1) + v(2) * p(ip)%results%dyA(2) + v(3) * p(ip)%results%dyA(3))
            grad(3) = grad(3) + (one - newmark_B) * (v(1) * p(ip)%results%dzA(1) + v(2) * p(ip)%results%dzA(2) + v(3) * p(ip)%results%dzA(3))

            !              vnablaA(1)        = v(1)*r(ip)%results%dxA(1) + v(2)*r(ip)%results%dyA(1)
            !              vnablaA(2)        = v(1)*r(ip)%results%dxA(2) + v(2)*r(ip)%results%dyA(2)
            !
            !              nabla2A(1)        = ( v(1)*r(ip)%results%dxxA(1) + v(2)*r(ip)%results%dxxA(2) + v(3)*r(ip)%results%dxxA(3) )*v(1) +&
            !                                  ( v(1)*r(ip)%results%dxyA(1) + v(2)*r(ip)%results%dxyA(2) + v(3)*r(ip)%results%dxyA(3) )*v(2)
            !
            !              nabla2A(2)        = ( v(1)*r(ip)%results%dxyA(1) + v(2)*r(ip)%results%dxyA(2) + v(3)*r(ip)%results%dxyA(3) )*v(1) +&
            !                                  ( v(1)*r(ip)%results%dyyA(1) + v(2)*r(ip)%results%dyyA(2) + v(3)*r(ip)%results%dyyA(3) )*v(2)
            !
            !              gradE(1)          = v(1)*r(ip)%results%dxE(1) + v(2)*r(ip)%results%dxE(2)
            !              gradE(2)          = v(1)*r(ip)%results%dyE(1) + v(2)*r(ip)%results%dyE(2)

            grad = grad / lorentz_tilde
            !              vnablaA           = vnablaA/lorentz_tilde
            !              gradE             = gradE/lorentz_tilde
            !              nabla2A           = nabla2A/lorentz_tilde
            !
            Eirr = newmark_Ei * r(ip)%results%E + (one - newmark_Ei) * p(ip)%results%E

            !              coef1             = -    dot_product( Eirr , v )/p(ip)%data%g**2/lorentz_tilde**2
            !              coef2             = -two*dot_product( grad , v )/p(ip)%data%g**2/lorentz_tilde**2

            !
            Anew = (dA_1 * r(ip)%results%A + dA_0 * p(ip)%results%A + dA__1 * pold(ip)%results%A) / lorentz_tilde / dt
            Aold = (dA_1 * p(ip)%results%A + dA_0 * pold(ip)%results%A + dA__1 * poldold(ip)%results%A) / lorentz_tilde / dt
            Esol = newmark_Es * Anew + (one - newmark_Es) * Aold

            vxb = cross_product(v, B0)

            r(ip)%x = p(ip)%x + dt * v            !                                                       &
            !                                + half*dt**2*e*( m*p(ip)%data%v + e/lorentz_tilde*p(ip)%results%A )*coef1                 &
            !                                + half*dt**2*e*( m*p(ip)%data%v + e/lorentz_tilde*p(ip)%results%A )*coef2                 &
            !                                + half*dt**2*e*( grad + vnablaA + Eirr  )/m/p(ip)%data%g

            r(ip)%data%v = p(ip)%data%v + dt * e / m * (Eirr - Esol + grad + vxb)!                     + & ! Canonical Momentum
            !                                       - half*dt**2*e/m*( gradE + nabla2A     )

            !              r(ip)%x(3)        = zero

            if (periodicity_particles) call iperiodic_particles(r(ip))

            err_x = sum((s(ip)%x(1:2) - r(ip)%x(1:2))**2)
            err_v = sum((s(ip)%P(1:3) - r(ip)%data%v(1:3))**2)

            toll = max(toll, sqrt(err_x + err_v))

         end do

         call MPI_ALLREDUCE(MPI_IN_PLACE, toll, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, rc)
         rat = toll / toll_old

         if (.false. .and. (rat .gt. one) .and. (itc .eq. maxit)) then !( itc .eq. maxit ) then !

            if (root) then
               write (*, '(a,es12.4,es12.4)') " == Picard iteration is not converging : ", toll, toll_old
            end if

            write (filename, '(a,"errors_",i6.6,".dat")') trim(folder), my_rank
            open (newunit=filehandle, file=trim(filename), STATUS='REPLACE')

            do ip = 1, np

               err_x = sum((s(ip)%x(1:2) - r(ip)%x(1:2))**2)
               err_v = sum((s(ip)%P(1:3) - r(ip)%data%v(1:3))**2)
               err_A = sum((s(ip)%A(1:3) - r(ip)%results%A(1:3))**2)
               err_dxA = sum((s(ip)%dxA(1:3) - r(ip)%results%dxA(1:3))**2)
               err_dyA = sum((s(ip)%dyA(1:3) - r(ip)%results%dyA(1:3))**2)

               write (filehandle, *) ip, err_x, err_v, err_A, err_dxA, err_dyA

            end do

            close (filehandle)
         end if

         !            endif

         do ip = 1, np

            s(ip)%P = r(ip)%data%v
            s(ip)%x = r(ip)%x
            !              s(ip)%x(3)        =  zero
            s(ip)%A = r(ip)%results%A

            s(ip)%dxA(1:3) = r(ip)%results%dxA(1:3)
            s(ip)%dyA(1:3) = r(ip)%results%dyA(1:3)

         end do

      end do

      do ip = 1, np

         e = p(ip)%data%q
         m = p(ip)%data%m

         v = p(ip)%data%v
         p(ip)%data%v = s(ip)%P

         !        if ( p(ip)%label .eq. 1 ) p(ip)%data%v = zero
         p(ip)%data%g = sqrt(one + sum((p(ip)%data%v / vtilde)**2))
         p(ip)%x = s(ip)%x
         !        p(ip)%x(3)        = zero
         if (periodicity_particles) call iperiodic_particles(p(ip))

      end do

      deallocate (r, s)!,error )

   end subroutine trapezoidal

   subroutine euler_method(np, dt, p)
      use module_tool, only: cross_product
      use module_globals, only: periodicity_particles, ischeme, root, my_rank, folder, step, B0
      use helper, only: iperiodic_particles, inormalize, write_particles_ascii, write_particles_vtk
      use module_pepc
      use mpi
      implicit none

      !------------------------------------------------------------------------------------------------------------------------------
      !
      !    Ref.: Leimkuhler B., Reich S. -Simulating Hamiltonian dynamics
      !
      !    This subroutine is meant to push in time the particles. We implement here  the time integration scheme (Trapezoidal rule):
      !
      !    x^(t+1)      = x^(t) +   dt*v^(t+1/2)
      !    P^(t+1)      = P^(t) + q*dt*( E^(t+1/2) + grad( A^(t+1/2) . v^(t+1/2) ) )
      !
      !    Since the Trapezoidal rule is an implicit scheme, the nonlinearity is resolved with a Picard iteration.
      !    Picard method is very simple to implement, compared to Newton-Krylov, but it is not guaranteed the method converge toward
      !    the solution. In fact, it is, required the function  (the map which advance in time) is a contraction in the phase space.
      !
      !    np represents the number of particles in the current process
      !
      !    dt is the time step
      !
      !    p is the array of particles at the old time step, "t".
      !
      !
      !------------------------------------------------------------------------------------------------------------------------------
      integer(kind_particle), intent(in)           :: np
      real(kind_particle), intent(in)              :: dt
      type(t_particle), allocatable, intent(inout) :: p(:)
      !    integer(kind_particle)          , intent(out)    :: itc
      !    real(kind_particle)             , intent(out)    :: toll

      real(kind_particle)                          :: e, m, v(1:3), grad(1:3), vxb(1:3)

      type(t_particle), allocatable                :: r(:)!,error(:)
      integer(kind_particle)                       :: ip
      integer                                      :: rc

      if (allocated(r)) deallocate (r)
      allocate (r(np), stat=rc)

      !------------------------------------------------------------------------------------------------------------------------------
      !        Here, we define the initial guess, all iterative methods need an initial guess, supposed to be close enough to the solution.
      !        We assume the guess solution is the configuration at old time t.
      !------------------------------------------------------------------------------------------------------------------------------

      do ip = 1, np
         !------------------------------------------------------------------------------------------------------------------------------
         !        At each iteration we update the canonical momentums, according to the Trapezoidal rule.
         !        This method assumes the fields are computed at the future values:
         !
         !                x^(t+1)
         !                P^(t+1) = m*gamma*v^(t+1) + q/c*A^(t+1)
         !
         !        We slightly modify the definition of velocity, we replace the natural definition of   v^(t+1/2) = ( v^(t) + v^(t+1) )/2
         !        with  v^(t+1/2) = ( v^(t) + v^(t+1) )/( gamma^(t) + gamma^(t+1) ).
         !
         !        It turns out that for non relativistic velocities the two definitions are equivalent.
         !
         !        So, the fields are evaluated as average values, for instance E^(t+1/2) = ( E( x^(t+1) ) + E( x^(t) ) )/2, etc.
         !        Notice:  s(ip)%P is the canonical momentum at time t+1, but regarding the previous iteration of the Picard method.
         !        In this piece of code s(ip) contains the information of the particle at time time t+1.
         !        We are redefining the particle position r(ip) as explained before.
         !        Be aware the fields are calculated at the phase space point x, gamma*v
         !------------------------------------------------------------------------------------------------------------------------------

         r(ip)%data%q = p(ip)%data%q
         r(ip)%data%m = p(ip)%data%m

         r(ip)%data%v = p(ip)%data%v
         r(ip)%data%g = sqrt(one + sum((r(ip)%data%v / vtilde)**2))
         r(ip)%x = p(ip)%x + dt * p(ip)%data%v / p(ip)%data%g
         r(ip)%x(3) = zero

         r(ip)%label = p(ip)%label
         r(ip)%results%E = zero
         r(ip)%results%pot = zero
         r(ip)%results%B = zero
         r(ip)%results%A = zero
         r(ip)%results%dxA = zero
         r(ip)%results%dyA = zero
         r(ip)%results%Jirr = zero
         r(ip)%results%J = zero
         r(ip)%work = one

         if (periodicity_particles) call iperiodic_particles(r(ip))

      end do

      !------------------------------------------------------------------------------------------------------------------------------
      !            The fields are, as written above, evaluated at the future predicted phase space point, t+1.
      !------------------------------------------------------------------------------------------------------------------------------

      call pepc_particleresults_clear(r)
      call pepc_grow_tree(r)
      call pepc_traverse_tree(r)
      call pepc_restore_particles(r)
      call pepc_timber_tree()

      do ip = 1, np

         !------------------------------------------------------------------------------------------------------------------------------
         !            Once the fields are updated, we push the particle in time according to the Trapezoidal rule:
         !
         !              x^(t+1) = x^(t) +   dt*v^(t)
         !              P^(t+1) = m*gamma*v^(t+1) + q/c*A^(t+1)
         !              P^(t+1) = P^(t) + q*dt*(E^(t+1) + q/c*grad( A^(t+1).v^(t) )
         !------------------------------------------------------------------------------------------------------------------------------

         e = p(ip)%data%q
         m = p(ip)%data%m

         call inormalize(r(ip))

         v = p(ip)%data%v / p(ip)%data%g

         grad = zero
         grad(1) = v(1) * r(ip)%results%dxA(1) + v(2) * r(ip)%results%dxA(2) + v(3) * r(ip)%results%dxA(3)
         grad(2) = v(1) * r(ip)%results%dyA(1) + v(2) * r(ip)%results%dyA(2) + v(3) * r(ip)%results%dyA(3)
         grad = grad / lorentz_tilde

         vxb = cross_product(v, B0)

         !              p(ip)%data%v      = p(ip)%data%v   + dt*e/m*( r(ip)%results%E )
         !              p(ip)%data%v      = p(ip)%data%v   + dt/m*e*( vxb )
         p(ip)%data%v = m * p(ip)%data%v + e / lorentz_tilde * p(ip)%results%A + dt * e * (r(ip)%results%E + grad + vxb)
         p(ip)%data%v = (p(ip)%data%v - e / lorentz_tilde * r(ip)%results%A) / m
         p(ip)%data%g = sqrt(one + sum((p(ip)%data%v / vtilde)**2))

         p(ip)%x = r(ip)%x
         !              p(ip)%x           = p(ip)%x        + half*dt*( p(ip)%data%v/p(ip)%data%g + v )
         p(ip)%x(3) = zero

         if (periodicity_particles) call iperiodic_particles(p(ip))

      end do

      deallocate (r)

   end subroutine euler_method

   subroutine euler_method3d(np, dt, p)
      use module_tool, only: cross_product
      use module_globals, only: periodicity_particles, B0
      use helper, only: iperiodic_particles, inormalize, write_particles_ascii, write_particles_vtk
      use module_pepc
      use mpi
      implicit none

      !------------------------------------------------------------------------------------------------------------------------------
      !
      !    Ref.: Leimkuhler B., Reich S. -Simulating Hamiltonian dynamics
      !
      !    This subroutine is meant to push in time the particles. We implement here  the time integration scheme (Trapezoidal rule):
      !
      !    x^(t+1)      = x^(t) +   dt*v^(t+1/2)
      !    P^(t+1)      = P^(t) + q*dt*( E^(t+1/2) + grad( A^(t+1/2) . v^(t+1/2) ) )
      !
      !    Since the Trapezoidal rule is an implicit scheme, the nonlinearity is resolved with a Picard iteration.
      !    Picard method is very simple to implement, compared to Newton-Krylov, but it is not guaranteed the method converge toward
      !    the solution. In fact, it is, required the function  (the map which advance in time) is a contraction in the phase space.
      !
      !    np represents the number of particles in the current process
      !
      !    dt is the time step
      !
      !    p is the array of particles at the old time step, "t".
      !
      !
      !------------------------------------------------------------------------------------------------------------------------------
      integer(kind_particle), intent(in)           :: np
      real(kind_particle), intent(in)              :: dt
      type(t_particle), allocatable, intent(inout) :: p(:)
      !    integer(kind_particle)          , intent(out)    :: itc
      !    real(kind_particle)             , intent(out)    :: toll

      real(kind_particle)                          :: e, m, v(1:3), grad(1:3), vxb(1:3)

      type(t_particle), allocatable                :: r(:)!,error(:)
      integer(kind_particle)                       :: ip
      integer                                      :: rc

      if (allocated(r)) deallocate (r)
      allocate (r(np), stat=rc)

      !------------------------------------------------------------------------------------------------------------------------------
      !        Here, we define the initial guess, all iterative methods need an initial guess, supposed to be close enough to the solution.
      !        We assume the guess solution is the configuration at old time t.
      !------------------------------------------------------------------------------------------------------------------------------

      do ip = 1, np
         !------------------------------------------------------------------------------------------------------------------------------
         !        At each iteration we update the canonical momentums, according to the Trapezoidal rule.
         !        This method assumes the fields are computed at the future values:
         !
         !                x^(t+1)
         !                P^(t+1) = m*gamma*v^(t+1) + q/c*A^(t+1)
         !
         !        We slightly modify the definition of velocity, we replace the natural definition of   v^(t+1/2) = ( v^(t) + v^(t+1) )/2
         !        with  v^(t+1/2) = ( v^(t) + v^(t+1) )/( gamma^(t) + gamma^(t+1) ).
         !
         !        It turns out that for non relativistic velocities the two definitions are equivalent.
         !
         !        So, the fields are evaluated as average values, for instance E^(t+1/2) = ( E( x^(t+1) ) + E( x^(t) ) )/2, etc.
         !        Notice:  s(ip)%P is the canonical momentum at time t+1, but regarding the previous iteration of the Picard method.
         !        In this piece of code s(ip) contains the information of the particle at time time t+1.
         !        We are redefining the particle position r(ip) as explained before.
         !        Be aware the fields are calculated at the phase space point x, gamma*v
         !------------------------------------------------------------------------------------------------------------------------------

         r(ip)%data%q = p(ip)%data%q
         r(ip)%data%m = p(ip)%data%m

         r(ip)%data%v = p(ip)%data%v
         r(ip)%data%g = sqrt(one + sum((r(ip)%data%v / vtilde)**2))
         r(ip)%x = p(ip)%x + dt * p(ip)%data%v / p(ip)%data%g

         r(ip)%label = p(ip)%label
         r(ip)%results%E = zero
         r(ip)%results%pot = zero
         r(ip)%results%B = zero
         r(ip)%results%A = zero
         r(ip)%results%dxA = zero
         r(ip)%results%dyA = zero
         r(ip)%results%dzA = zero
         r(ip)%results%Jirr = zero
         r(ip)%results%J = zero
         r(ip)%work = one

         if (periodicity_particles) call iperiodic_particles(r(ip))

      end do

      !------------------------------------------------------------------------------------------------------------------------------
      !            The fields are, as written above, evaluated at the future predicted phase space point, t+1.
      !------------------------------------------------------------------------------------------------------------------------------

      call pepc_particleresults_clear(r)
      call pepc_grow_tree(r)
      call pepc_traverse_tree(r)
      call pepc_restore_particles(r)
      call pepc_timber_tree()

      do ip = 1, np

         !------------------------------------------------------------------------------------------------------------------------------
         !            Once the fields are updated, we push the particle in time according to the Trapezoidal rule:
         !
         !              x^(t+1) = x^(t) +   dt*v^(t)
         !              P^(t+1) = m*gamma*v^(t+1) + q/c*A^(t+1)
         !              P^(t+1) = P^(t) + q*dt*(E^(t+1) + q/c*grad( A^(t+1).v^(t) )
         !------------------------------------------------------------------------------------------------------------------------------

         e = p(ip)%data%q
         m = p(ip)%data%m

         call inormalize(r(ip))

         v = p(ip)%data%v / p(ip)%data%g

         grad = zero
         grad(1) = v(1) * r(ip)%results%dxA(1) + v(2) * r(ip)%results%dxA(2) + v(3) * r(ip)%results%dxA(3)
         grad(2) = v(1) * r(ip)%results%dyA(1) + v(2) * r(ip)%results%dyA(2) + v(3) * r(ip)%results%dyA(3)
         grad(3) = v(1) * r(ip)%results%dzA(1) + v(2) * r(ip)%results%dzA(2) + v(3) * r(ip)%results%dzA(3)
         grad = grad / lorentz_tilde

         vxb = cross_product(v, B0)

         p(ip)%x = p(ip)%x + dt * v
         p(ip)%data%v = m * p(ip)%data%v + e / lorentz_tilde * p(ip)%results%A + dt * e * (r(ip)%results%E + grad + vxb)
         p(ip)%data%v = (p(ip)%data%v - e / lorentz_tilde * r(ip)%results%A) / m
         p(ip)%data%g = sqrt(one + sum((p(ip)%data%v / vtilde)**2))

         if (periodicity_particles) call iperiodic_particles(p(ip))

      end do

      deallocate (r)

   end subroutine euler_method3d

   subroutine explicit_verlet(np, dt, p)
      use module_tool, only: cross_product
      use module_globals, only: periodicity_particles, ischeme, root, my_rank, folder, step, B0
      use helper, only: iperiodic_particles, inormalize, write_particles_ascii, write_particles_vtk
      use module_pepc
      use mpi
      implicit none

      !------------------------------------------------------------------------------------------------------------------------------
      !
      !    Ref.: Leimkuhler B., Reich S. -Simulating Hamiltonian dynamics
      !
      !    This subroutine is meant to push in time the particles. We implement here  the time integration scheme (Trapezoidal rule):
      !
      !    x^(t+1)      = x^(t) +   dt*v^(t+1/2)
      !    P^(t+1)      = P^(t) + q*dt*( E^(t+1/2) + grad( A^(t+1/2) . v^(t+1/2) ) )
      !
      !    Since the Trapezoidal rule is an implicit scheme, the nonlinearity is resolved with a Picard iteration.
      !    Picard method is very simple to implement, compared to Newton-Krylov, but it is not guaranteed the method converge toward
      !    the solution. In fact, it is, required the function  (the map which advance in time) is a contraction in the phase space.
      !
      !    np represents the number of particles in the current process
      !
      !    dt is the time step
      !
      !    p is the array of particles at the old time step, "t".
      !
      !
      !------------------------------------------------------------------------------------------------------------------------------
      integer(kind_particle), intent(in)           :: np
      real(kind_particle), intent(in)              :: dt
      type(t_particle), allocatable, intent(inout) :: p(:)
      !    integer(kind_particle)          , intent(out)    :: itc
      !    real(kind_particle)             , intent(out)    :: toll

      real(kind_particle)                          :: e, m, v(1:3), grad(1:3), vxb(1:3)

      type(t_particle), allocatable                :: r(:), s(:)!,error(:)
      integer(kind_particle)                       :: ip
      integer                                      :: rc

      if (allocated(r)) deallocate (r)
      allocate (r(np), stat=rc)
      if (allocated(s)) deallocate (s)
      allocate (s(np), stat=rc)

      !------------------------------------------------------------------------------------------------------------------------------
      !        Here, we define the initial guess, all iterative methods need an initial guess, supposed to be close enough to the solution.
      !        We assume the guess solution is the configuration at old time t.
      !------------------------------------------------------------------------------------------------------------------------------

      do ip = 1, np
         !------------------------------------------------------------------------------------------------------------------------------
         !        At each iteration we update the canonical momentums, according to the Trapezoidal rule.
         !        This method assumes the fields are computed at the future values:
         !
         !                x^(t+1)
         !                P^(t+1) = m*gamma*v^(t+1) + q/c*A^(t+1)
         !
         !        We slightly modify the definition of velocity, we replace the natural definition of   v^(t+1/2) = ( v^(t) + v^(t+1) )/2
         !        with  v^(t+1/2) = ( v^(t) + v^(t+1) )/( gamma^(t) + gamma^(t+1) ).
         !
         !        It turns out that for non relativistic velocities the two definitions are equivalent.
         !
         !        So, the fields are evaluated as average values, for instance E^(t+1/2) = ( E( x^(t+1) ) + E( x^(t) ) )/2, etc.
         !        Notice:  s(ip)%P is the canonical momentum at time t+1, but regarding the previous iteration of the Picard method.
         !        In this piece of code s(ip) contains the information of the particle at time time t+1.
         !        We are redefining the particle position r(ip) as explained before.
         !        Be aware the fields are calculated at the phase space point x, gamma*v
         !------------------------------------------------------------------------------------------------------------------------------

         r(ip)%data%q = p(ip)%data%q
         r(ip)%data%m = p(ip)%data%m

         r(ip)%data%v = p(ip)%data%v
         r(ip)%data%g = p(ip)%data%g
         r(ip)%x = p(ip)%x + half * dt * p(ip)%data%v / p(ip)%data%g
         r(ip)%x(3) = zero

         r(ip)%label = p(ip)%label
         r(ip)%results%E = zero
         r(ip)%results%pot = zero
         r(ip)%results%B = zero
         r(ip)%results%A = zero
         r(ip)%results%dxA = zero
         r(ip)%results%dyA = zero
         r(ip)%results%Jirr = zero
         r(ip)%results%J = zero
         r(ip)%work = one

         if (periodicity_particles) call iperiodic_particles(r(ip))

      end do

      !------------------------------------------------------------------------------------------------------------------------------
      !            The fields are, as written above, evaluated at the future predicted phase space point, t+1.
      !------------------------------------------------------------------------------------------------------------------------------

      call pepc_particleresults_clear(r)
      call pepc_grow_tree(r)
      call pepc_traverse_tree(r)
      call pepc_restore_particles(r)
      call pepc_timber_tree()

      do ip = 1, np

         !------------------------------------------------------------------------------------------------------------------------------
         !            Once the fields are updated, we push the particle in time according to the Trapezoidal rule:
         !
         !              x^(t+1) = x^(t) +   dt*v^(t)
         !              P^(t+1) = m*gamma*v^(t+1) + q/c*A^(t+1)
         !              P^(t+1) = P^(t) + q*dt*(E^(t+1) + q/c*grad( A^(t+1).v^(t) )
         !------------------------------------------------------------------------------------------------------------------------------

         call inormalize(r(ip))

         e = p(ip)%data%q
         m = p(ip)%data%m

         v = p(ip)%data%v / p(ip)%data%g

         grad = zero
         grad(1) = v(1) * r(ip)%results%dxA(1) + v(2) * r(ip)%results%dxA(2) + v(3) * r(ip)%results%dxA(3)
         grad(2) = v(1) * r(ip)%results%dyA(1) + v(2) * r(ip)%results%dyA(2) + v(3) * r(ip)%results%dyA(3)
         grad = grad / lorentz_tilde

         vxb = cross_product(v, B0)

         p(ip)%data%v = m * p(ip)%data%v + e / lorentz_tilde * p(ip)%results%A + half * dt * e * (r(ip)%results%E + grad + vxb)
         p(ip)%data%v = (p(ip)%data%v - e / lorentz_tilde * r(ip)%results%A) / m
         p(ip)%data%v = two * p(ip)%data%v - p(ip)%data%g * v
         p(ip)%data%g = sqrt(one + sum((p(ip)%data%v / vtilde)**2))

         p(ip)%x = r(ip)%x + half * dt * p(ip)%data%v / p(ip)%data%g
         p(ip)%x(3) = zero

         !              r(ip)%x           = r(ip)%x + half*dt*s(ip)%data%v/s(ip)%data%g
         !              r(ip)%x(3)        = zero

         !
         !              s(ip)%data%v      = m*p(ip)%data%v + e/lorentz_tilde *p(ip)%results%A + half*dt*e*( r(ip)%results%E + grad + vxb )
         !              s(ip)%data%v      = ( s(ip)%data%v - e/lorentz_tilde *r(ip)%results%A )/m
         !              s(ip)%data%g      = sqrt( one + sum( ( s(ip)%data%v/vtilde )**2 ) )
         !
         !
         !              s(ip)%x           = r(ip)%x + half*dt*s(ip)%data%v/s(ip)%data%g
         !              s(ip)%x(3)        = zero
         !
         !              s(ip)%label       = p(ip)%label
         !              s(ip)%results%E   = zero
         !              s(ip)%results%pot = zero
         !              s(ip)%results%B   = zero
         !              s(ip)%results%A   = zero
         !              s(ip)%results%dxA = zero
         !              s(ip)%results%dyA = zero
         !              s(ip)%results%Jirr= zero
         !              s(ip)%results%J   = zero
         !              s(ip)%work        = one

      end do

      !    call pepc_particleresults_clear(s)
      !    call pepc_grow_tree(s)
      !    call pepc_traverse_tree(s)
      !    call pepc_restore_particles(s)
      !    call pepc_timber_tree()
      !
      !
      !
      !    do ip = 1,np
      !
      !              call  inormalize(s(ip))
      !
      !              e                 = p(ip)%data%q
      !              m                 = p(ip)%data%m
      !
      !              v                 = s(ip)%data%v/s(ip)%data%g
      !
      !              grad              = zero
      !              grad(1)           = v(1)*s(ip)%results%dxA(1) + v(2)*s(ip)%results%dxA(2) + v(3)*s(ip)%results%dxA(3)
      !              grad(2)           = v(1)*s(ip)%results%dyA(1) + v(2)*s(ip)%results%dyA(2) + v(3)*s(ip)%results%dyA(3)
      !              grad              = grad/lorentz_tilde
      !
      !              vxb               = cross_product(v,B0)
      !
      !
      !              p(ip)%data%v      = m*r(ip)%data%v + e/lorentz_tilde *r(ip)%results%A + dt*e*( r(ip)%results%E + grad + vxb )
      !              p(ip)%data%v      = ( p(ip)%data%v - e/lorentz_tilde *s(ip)%results%A )/m
      !              p(ip)%data%g      = sqrt( one + sum( ( p(ip)%data%v/vtilde )**2 ) )
      !
      !
      !              p(ip)%data%v      = s(ip)%data%v
      !              p(ip)%data%g      = s(ip)%data%g
      !
      !
      !              p(ip)%x           = r(ip)%x + half*dt*r(ip)%data%v/r(ip)%data%g
      !              p(ip)%x(3)        = zero
      !
      !
      !
      !    enddo

      deallocate (r, s)

   end subroutine explicit_verlet

end module module_integration
