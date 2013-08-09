!
! This is just a placeholder
!

module pf_mod_verlet
  use pf_mod_dtype
  implicit none
  integer, parameter, private :: npieces = 2
  
  interface
     subroutine pf_acceleration_p(x, t, level, ctx, a)
       import c_ptr, c_int, pfdp
       type(c_ptr),    intent(in), value :: x, a, ctx
       real(pfdp),     intent(in)        :: t
       integer(c_int), intent(in)        :: level
     end subroutine pf_acceleration_p
  end interface
  
  type :: pf_verlet_t
     procedure(pf_acceleration_p), pointer, nopass :: acceleration
  end type pf_verlet_t

contains

  ! Perform on SDC sweep on level F and set qend appropriately.
  subroutine verlet_sweep(pf, F, t0, dt)
    use pf_mod_timer

    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp),        intent(in)    :: dt, t0
    type(pf_level_t),  intent(inout) :: F

    integer    :: m, n
    real(pfdp) :: t
    real(pfdp) :: dtsdc(1:F%nnodes-1)
    type(c_ptr) :: rhs

    type(pf_verlet_t), pointer :: verlet
    call c_f_pointer(F%sweeper%ctx, verlet)

    call start_timer(pf, TLEVEL+F%level-1)

    ! FIXME: fill with content here

    call end_timer(pf, TLEVEL+F%level-1)
  end subroutine verlet_sweep

  ! Evaluate function values
  subroutine verlet_evaluate(F, t, m)
    real(pfdp),       intent(in)    :: t
    integer,          intent(in)    :: m
    type(pf_level_t), intent(inout) :: F

    ! FIXME: fill with content here

  end subroutine verlet_evaluate

  ! Initialize smats
  subroutine verlet_initialize(F)
    use pf_mod_dtype
    type(pf_level_t), intent(inout) :: F
    real(pfdp) :: dsdc(F%nnodes-1)

    integer :: m

    allocate(F%smat(F%nnodes-1,F%nnodes,npieces))

    ! FIXME: fill with content here

  end subroutine verlet_initialize

  ! Compute SDC integral
  subroutine verlet_integrate(F, qSDC, fSDC, dt, fintSDC)
    type(pf_level_t), intent(in)    :: F
    type(c_ptr),      intent(in)    :: qSDC(:), fSDC(:, :)
    real(pfdp),       intent(in)    :: dt
    type(c_ptr),      intent(inout) :: fintSDC(:)
    
    ! FIXME: fill with content here

  end subroutine verlet_integrate

  ! Create/destroy VERLET sweeper
  subroutine pf_verlet_create(sweeper, acceleration)
    type(pf_sweeper_t), intent(inout) :: sweeper
    procedure(pf_acceleration_p) :: acceleration

    type(pf_verlet_t), pointer :: verlet

    allocate(verlet)
    verlet%acceleration => acceleration

    sweeper%npieces    =  npieces
    sweeper%sweep      => verlet_sweep
    sweeper%evaluate   => verlet_evaluate
    sweeper%initialize => verlet_initialize
    sweeper%integrate  => verlet_integrate

    sweeper%ctx = c_loc(verlet)
  end subroutine pf_verlet_create

  subroutine pf_verlet_destroy(sweeper)
    type(pf_sweeper_t), intent(inout) :: sweeper

    type(pf_verlet_t), pointer :: verlet
    call c_f_pointer(sweeper%ctx, verlet)

    deallocate(verlet)
  end subroutine pf_verlet_destroy

end module pf_mod_verlet

