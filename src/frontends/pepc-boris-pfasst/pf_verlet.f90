module pf_mod_verlet
  use pf_mod_dtype
  use pfm_encap
  implicit none
  integer, parameter, private :: npieces = 1

  interface
     subroutine pf_calc_Efield_p(x, t, level, ctx, E)
       import c_ptr, c_int, pfdp
       type(c_ptr),    intent(in), value :: x, E, ctx
       real(pfdp),     intent(in)        :: t
       integer(c_int), intent(in)        :: level
     end subroutine pf_calc_Efield_p
  end interface

  interface
     subroutine pf_build_rhs_p(E, Q, level, ctx, rhs)
       import c_ptr, c_int, pfdp
       type(c_ptr),    intent(in), value :: E, Q, rhs, ctx
       integer(c_int), intent(in)        :: level
     end subroutine pf_build_rhs_p
  end interface

  interface
     subroutine pf_impl_solver_p(v, level, ctx, v_old, E_old, E_new, SDCint, dt)
       import c_ptr, c_int, pfdp
       type(c_ptr),    intent(in), value :: v, v_old, E_old, E_new, SDCint, ctx
       integer(c_int), intent(in)        :: level
       real(pfdp),     intent(in)        :: dt
     end subroutine pf_impl_solver_p
  end interface

  type :: pf_verlet_t
     procedure(pf_calc_Efield_p), pointer, nopass :: calc_Efield
     procedure(pf_build_rhs_p), pointer, nopass   :: build_rhs
     procedure(pf_impl_solver_p), pointer, nopass :: impl_solver
     real(pfdp), allocatable :: smat(:,:,:)
  end type pf_verlet_t

contains

  !
  ! Perform one SDC sweep on level F and set qend appropriately
  !
  subroutine verlet_sweep(pf, F, t0, dt)
    use pf_mod_timer

    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp),        intent(in)    :: dt, t0
    type(pf_level_t),  intent(inout) :: F

    integer     :: m, n
    real(pfdp)  :: t,dtsq
    real(pfdp)  :: dtsdc(1:F%nnodes-1)
    type(c_ptr) :: fq_int, rhs

    type(pf_verlet_t), pointer :: verlet

    call c_f_pointer(F%sweeper%sweeperctx, verlet)
    call start_timer(pf, TLEVEL+F%level-1)

    !
    ! compute integrals and add fas correction
    !

    call F%encap%create(fq_int, F%level, SDC_KIND_INTEGRAL, F%nvars, F%shape, F%levelctx, F%encap%encapctx)
    call F%encap%setval(fq_int, 0.0_pfdp, SETVAL_XV)

    ! setup temporary variable to contain full right-hand side (based on E-field in F)
    call F%encap%create(rhs, F%level, SDC_KIND_FEVAL, F%nvars, F%shape, F%levelctx, F%encap%encapctx)

    dtsq = dt*dt
    do m = 1, F%nnodes-1
       call F%encap%setval(F%S(m), 0.0_pfdp, SETVAL_XV)
       do n = 1, F%nnodes
          call verlet%build_rhs(F%F(n,1), F%Q(n), F%level, F%levelctx, rhs)
          call F%encap%axpy(F%S(m), dt  *F%smat(m,n,1), rhs, AXPY_aVpV)
          call F%encap%axpy(F%S(m), dtsq*F%smat(m,n,2), rhs, AXPY_aXpX)
       end do
       if (associated(F%tau)) then
          call F%encap%axpy(F%S(m), 1.0_pfdp, F%tau(m), AXPY_aVpV_aXpX)
          ! tau is 0-to-node in our case, we need node-to-node instead
          if (m > 1) then
            call F%encap%axpy(F%S(m), -1.0_pfdp, F%tau(m-1), AXPY_aVpV_aXpX)
          endif
       end if
    end do

    !
    ! do time-stepping
    !

    call F%encap%unpack(F%Q(1), F%q0)

    call verlet%calc_Efield(F%Q(1),t0,F%level, F%levelctx, F%F(1,1))

    t = t0
    dtsdc = dt * (F%nodes(2:F%nnodes) - F%nodes(1:F%nnodes-1))
    do m = 1, F%nnodes-1
       t = t + dtsdc(m)
       call F%encap%setval(fq_int, 0.0_pfdp, SETVAL_XV)

       !  Lower triangular verlet to old new
       do n = 1, m
          call verlet%build_rhs(F%F(n,1), F%Q(n), F%level, F%levelctx, rhs)
          call F%encap%axpy(fq_int, dtsq*F%smat(m,n,3), rhs, AXPY_aXpX)
       end do

       !  Update position term (trapezoid rule)
       call F%encap%copy(F%Q(m+1), fq_int, COPY_X)
       call F%encap%axpy(F%Q(m+1), dtsdc(m), F%Q(1), AXPY_aVpX) !  Add the dt*v_0 term
       call F%encap%axpy(F%Q(m+1), 1.0_pfdp, F%S(m), AXPY_aXpX)  !  Add integration term for p
       call F%encap%axpy(F%Q(m+1), 1.0_pfdp, F%Q(m), AXPY_aXpX)  !  Start m+1 with value from m

       call verlet%calc_Efield(F%Q(m+1), t, F%level, F%levelctx, F%F(m+1,1))

       ! Call implicit solver for updated velocity term
       call verlet%impl_solver(F%Q(m+1), &                ! Output: updated velocity at m+1
                               F%level, &                 ! current level
                               F%levelctx, &              ! level context
                               F%Q(m), &                  ! old velocity at previous node m
                               F%F(m+1,1), &              ! current E-field at m+1, using already updated positions
                               F%F(m,1), &                ! current E-field at previous node m
                               F%S(m),dtsdc(m))           ! mod. right-hand side (i.e. integral for q) + \Delta\tau

    end do

    call F%encap%copy(F%qend, F%Q(F%nnodes), COPY_XV)
    call F%encap%destroy(fq_int)
    call F%encap%destroy(rhs)

    call end_timer(pf, TLEVEL+F%level-1)
  end subroutine verlet_sweep


  !
  ! Evaluate function values
  !
  subroutine verlet_evaluate(F, t, m)
    real(pfdp),       intent(in)    :: t
    integer,          intent(in)    :: m
    type(pf_level_t), intent(inout) :: F
    type(pf_verlet_t), pointer :: verlet
    call c_f_pointer(F%sweeper%sweeperctx, verlet)
    ! Note: this is the E-field, not the full RHS and this is ok (will need this function for spread, interpolate and restrict)
    call verlet%calc_Efield(F%Q(m), t, F%level, F%levelctx, F%F(m,1))
  end subroutine verlet_evaluate


  !
  ! Initialize integration matrices
  !
  subroutine verlet_initialize(F)
    type(pf_level_t), intent(inout) :: F

    real(pfdp) :: dsdc(F%nnodes-1)
    real(pfdp) :: FE(F%nnodes,F%nnodes)
    real(pfdp) :: Trap(F%nnodes,F%nnodes)
    real(pfdp) :: Abar(F%nnodes,F%nnodes)
    real(pfdp) :: Abartil(F%nnodes,F%nnodes)
    real(pfdp) :: Q(F%nnodes,F%nnodes)
    real(pfdp) :: D(F%nnodes,F%nnodes)
    real(pfdp) :: Dtrap(F%nnodes,F%nnodes)

    real(pfdp) :: qtmp(F%nnodes,F%nnodes)

    integer :: m,i,j

    F%smat = 0.0_pfdp

    dsdc = F%nodes(2:F%nnodes) - F%nodes(1:F%nnodes-1)

    !  Build Q from qmat
    Q=0.0_pfdp

    !  Build Abar from Q
    do i = 2,F%nnodes
       Q(i,:) = Q(i-1,:) + F%s0mat(i-1,:)
    end do

    !  Build Abartil
    do i = 1,F%nnodes
       do j = 1,F%nnodes
          Abar(i,j) =  Q(i,j)*(F%nodes(i)-F%nodes(j))
       end do
    end do

    !  Set up Forward Euler and Trapezoid matrices
    FE=0.0_pfdp
    Trap=0.0_pfdp
    do i = 2,F%nnodes
       do j = 1,i-1
          FE(i,j) =  dsdc(j)
       end do
    end do
    do i = 2,F%nnodes
       do j = 1,i-1
          Trap(i,j) = Trap(i,j)+ 0.5d0*dsdc(j)
          Trap(i,j+1) = Trap(i,j+1)+ 0.5d0*dsdc(j)
       end do
    end do

    Abartil = matmul(FE,Trap) + 0.5d0*FE*FE

    !  Make differences
    D = Abar-Abartil
    Dtrap = Q-Trap

    !  Put the differnce increment in Smat
    do m = 1, F%nnodes-1
       F%smat(m,:,1) = (Dtrap(m+1,:)-Dtrap(m,:))
       F%smat(m,:,2) = (D(m+1,:)-D(m,:))
       F%smat(m,:,3) = Abartil(m+1,:)-Abartil(m,:)
       F%smat(m,:,4) = Abar(m+1,:)-Abar(m,:)
    end do

    ! populate qqmat (i.e. Q*Q without first line)
    qtmp(1 ,:) = 0.
    qtmp(2:,:) = F%qmat
    qtmp(:, :) = matmul(qtmp, qtmp)
    F%qqmat(:, :) = qtmp(2:,:)
  end subroutine verlet_initialize


  !
  ! Integrate (0-to-node)
  !
  subroutine verlet_integrate(F, qSDC, fSDC, dt, fintSDC)
    type(pf_level_t), intent(in) :: F
    type(c_ptr),      intent(in) :: qSDC(:)    !< unknowns at SDC nodes
    type(c_ptr),      intent(in) :: fSDC(:, :) !< function values at SDC nodes
    real(pfdp),       intent(in) :: dt
    type(c_ptr),      intent(inout) :: fintSDC(:) !< '0-to-node' integral

    real(pfdp) :: dtsq
    integer :: n, m
    type(c_ptr) :: rhs
    type(pf_verlet_t), pointer :: verlet

    call c_f_pointer(F%sweeper%sweeperctx, verlet)

    ! setup temporary variable to contain full right-hand side (based on E-field in F)
    call F%encap%create(rhs, F%level, SDC_KIND_FEVAL, F%nvars, F%shape, F%levelctx, F%encap%encapctx)

    dtsq = dt*dt

    do n = 1, F%nnodes-1
       call F%encap%setval(fintSDC(n), 0.0_pfdp, SETVAL_XV)

       do m = 1, F%nnodes
          call verlet%build_rhs(fSDC(m,1), qSDC(m), F%level, F%levelctx, rhs)
          ! NOW: rhs%x(1:dim, i) = f ; rhs%v(1:dim, i) = f
#define INTEGRATE_MM
#ifndef INTEGRATE_MM
          ! Qpic := dt(QQ \otimes Ix) + (Q \otimes Iv)
          ! I    := dt*[(Qpic \otimes Id)fSDC](qSDC)
          call F%encap%axpy(fintSDC(n), dtsq*F%qqmat(n,m), rhs, AXPY_aXpX)
          call F%encap%axpy(fintSDC(n), dt  *F%qmat( n,m), rhs, AXPY_aVpV)
#else
          call F%encap%axpy(fintSDC(n), dt *F%s0mat(n,m),   rhs, AXPY_aVpV)
          call F%encap%axpy(fintSDC(n), dtsq*F%smat(n,m,4), rhs, AXPY_aXpX)
#endif
       end do
    end do

    call F%encap%destroy(rhs)

#ifdef INTEGRATE_MM
    ! build 0-to-node from node-to-node integral
    do n = 2, F%nnodes-1
      call F%encap%axpy(fintSDC(n), 1.0_pfdp, fintSDC(n-1), AXPY_aVpV_aXpX)
    end do
#endif
  end subroutine verlet_integrate


  !
  ! Create Verlet sweeper
  !
  subroutine pf_verlet_create(sweeper, calc_Efield, build_rhs, impl_solver)
    type(pf_sweeper_t), intent(inout) :: sweeper
    procedure(pf_calc_Efield_p)       :: calc_Efield
    procedure(pf_build_rhs_p)         :: build_rhs
    procedure(pf_impl_solver_p)       :: impl_solver
    type(pf_verlet_t), pointer        :: verlet
    allocate(verlet)
    verlet%calc_Efield => calc_Efield
    verlet%build_rhs   => build_rhs
    verlet%impl_solver => impl_solver
    sweeper%npieces = 1
    sweeper%sweep       => verlet_sweep
    sweeper%evaluate    => verlet_evaluate
    sweeper%initialize  => verlet_initialize
    sweeper%integrate   => verlet_integrate
    sweeper%destroy     => null() ! we care for destroying the sweeper by ourselves
    sweeper%sweeperctx  =  c_loc(verlet)
  end subroutine pf_verlet_create

  subroutine pf_verlet_destroy(sweeper)
    type(pf_sweeper_t), intent(inout) :: sweeper
    type(pf_verlet_t), pointer :: verlet
    call c_f_pointer(sweeper%sweeperctx, verlet)
    deallocate(verlet)
  end subroutine pf_verlet_destroy

end module pf_mod_verlet

