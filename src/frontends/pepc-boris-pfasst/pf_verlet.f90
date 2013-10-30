module pf_mod_verlet
  use pf_mod_dtype
  implicit none
  integer, parameter, private :: npieces = 1

  interface
     subroutine pf_calc_Efield_p(x, t, level, ctx, a)
       import c_ptr, c_int, pfdp
       type(c_ptr),    intent(in), value :: x, a, ctx
       real(pfdp),     intent(in)        :: t
       integer(c_int), intent(in)        :: level
     end subroutine pf_calc_Efield_p
  end interface

  type :: pf_verlet_t
     procedure(pf_calc_Efield_p), pointer, nopass :: calc_Efield
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
    real(pfdp)  :: t,dtmhalf,dtsq
    real(pfdp)  :: dtsdc(1:F%nnodes-1)
    type(c_ptr) :: fq_int

    type(pf_verlet_t), pointer :: verlet

    call c_f_pointer(F%sweeper%sweeperctx, verlet)
    call start_timer(pf, TLEVEL+F%level-1)

    !
    ! compute integrals and add fas correction
    !

    call F%encap%create(fq_int, F%level, SDC_KIND_INTEGRAL, F%nvars, F%shape, F%levelctx, F%encap%encapctx)
    call F%encap%setval(fq_int, 0.0_pfdp)

    ! FIXME: need to substitute F (which is now the Efield) with the full rhs
    dtsq = dt*dt
    do m = 1, F%nnodes-1
       call F%encap%setval(F%S(m), 0.0_pfdp)
       do n = 1, F%nnodes
          call F%encap%axpy(F%S(m), dt*F%smat(m,n,1), F%F(n,1),1)
          call F%encap%axpy(F%S(m), dtsq*F%smat(m,n,2), F%F(n,1),2)
       end do
       if (associated(F%tau)) then
          call F%encap%axpy(F%S(m), 1.0_pfdp, F%tau(m))
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
       dtmhalf = 0.5d0*dtsdc(m)

       call F%encap%setval(fq_int, 0.0_pfdp)

       !  Lower triangular verlet to old new
       do n = 1, m
          call F%encap%axpy(fq_int, dtsq*F%smat(m,n,3), F%F(n,1), 2)
       end do

       !  Update position term (trapezoid rule)
       call F%encap%copy(F%Q(m+1), fq_int,2)
       call F%encap%axpy(F%Q(m+1), dtsdc(m), F%Q(1), 12) !  Add the dt*v_0 term
       call F%encap%axpy(F%Q(m+1), 1.0_pfdp, F%S(m), 2)  !  Add integration term for p
       call F%encap%axpy(F%Q(m+1), 1.0_pfdp, F%Q(m), 2)  !  Start m+1 with value from m

       call verlet%calc_Efield(F%Q(m+1), t, F%level, F%levelctx, F%F(m+1,1))

       ! Call boris solver for updated velocity term
       call F%encap%boris(F%Q(m+1),                 ! Output: updated velocity at m+1 
                          F%Q(m),                   ! old velocity at previous node m
                          F%F(m+1,1),               ! current E-field at m+1, using already updated positions
                          F%F(m,1),                 ! current E-field at previous node m
                          F%S(m),dtsdc(m))          ! mod. right-hand side (i.e. integral for q) + \Delta\tau

    end do

    call F%encap%copy(F%qend, F%Q(F%nnodes))
    call F%encap%destroy(fq_int)

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
    call verlet%calc_Efield(F%Q(m), t, F%level, F%levelctx, F%F(m,1))
    ! Note: this is the E-field, not the full RHS and this is ok (will need this function for spread, interpolate and restrict)
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

    integer :: m,i,j

    allocate(F%smat(F%nnodes-1,F%nnodes,4))

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
  end subroutine verlet_initialize


  !
  ! Integrate (node-to-node)
  !
  subroutine verlet_integrate(F, qSDC, fSDC, dt, fintSDC)
    type(pf_level_t), intent(in) :: F
    type(c_ptr),      intent(in) :: qSDC(:), fSDC(:, :)
    real(pfdp),       intent(in) :: dt
    type(c_ptr),      intent(inout) :: fintSDC(:)

    real(pfdp) :: dtsdc(1:F%nnodes-1)
    integer :: n, m
    ! FIXME: need to use the full RHS here, not the E-field
    dtsdc = dt * (F%nodes(2:F%nnodes) - F%nodes(1:F%nnodes-1))
    do n = 1, F%nnodes-1
       call F%encap%setval(fintSDC(n), 0.0_pfdp)
       do m = 1, F%nnodes
          call F%encap%axpy(fintSDC(n), dt*F%s0mat(n,m), fSDC(m,1),1)
          call F%encap%axpy(fintSDC(n), dt*dt*F%smat(n,m,4), fSDC(m,1),2)
       end do
    end do
  end subroutine verlet_integrate


  !
  ! Create Verlet sweeper
  !
  subroutine pf_verlet_create(sweeper, calc_Efield)
    type(pf_sweeper_t), intent(inout) :: sweeper
    procedure(pf_calc_Efield_p)      :: calc_Efield
    type(pf_verlet_t), pointer :: verlet
    allocate(verlet)
    verlet%calc_Efield => calc_Efield
    sweeper%npieces = 1
    sweeper%sweep       => verlet_sweep
    sweeper%evaluate    => verlet_evaluate
    sweeper%initialize  => verlet_initialize
    sweeper%integrate   => verlet_integrate
    sweeper%sweeperctx  =  c_loc(verlet)
  end subroutine pf_verlet_create

  subroutine pf_verlet_destroy(sweeper)
    type(pf_sweeper_t), intent(inout) :: sweeper
    type(pf_verlet_t), pointer :: verlet
    call c_f_pointer(sweeper%sweeperctx, verlet)
    deallocate(verlet)
  end subroutine pf_verlet_destroy

end module pf_mod_verlet

