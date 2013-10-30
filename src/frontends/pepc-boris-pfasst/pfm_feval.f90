module pfm_feval
    use pfm_encap
    use pf_mod_mpi
    use module_debug
    implicit none
    private

    public feval_init
    public feval_finalize
    public calc_Efield
    public eval_force
    public build_rhs
    public impl_solver

  contains

    !> Initialize feval, i.e. transfer initial u to PFASST data object
    subroutine feval_init(y0, yend, nlevels, levelctx, encapctx, particles) ! FIXME: shouldnt we call this function once per level?
      use iso_c_binding
      use module_pepc_types, only: t_particle
      implicit none

      type(app_data_t), pointer, intent(inout) :: y0, yend
      integer, intent(in) :: nlevels
      type(c_ptr), intent(inout) :: levelctx, encapctx
      type(t_particle), intent(in), pointer :: particles(:)

      type(c_ptr) :: y0_C, yend_c

      type(level_params_t), pointer :: params

      call pepc_status('|------> feval_init()')

      call c_f_pointer(levelctx, params)

      call encap_create(  y0_c, nlevels, -1, (2*params%dim)*params%nparts, [-1], levelctx, encapctx) ! dim*(coordinates and momenta) per particle
      call encap_create(yend_c, nlevels, -1, (2*params%dim)*params%nparts, [-1], levelctx, encapctx) ! dim*(coordinates and momenta) per particle

      call c_f_pointer(  y0_c, y0)
      call c_f_pointer(yend_c, yend)

      call particles_to_encap(y0_c, particles)

    end subroutine feval_init


    !> End feval
    subroutine feval_finalize(y0, yend)
      implicit none

      type(app_data_t), pointer, intent(inout) :: y0, yend

      call pepc_status('|------> feval_finalize()')

      call encap_destroy(c_loc(y0))
      call encap_destroy(c_loc(yend))
    end subroutine feval_finalize


    !> Use pepc to compute accelerations
    !> note here that we're assuming that acceleration
    !> is only a function of x, and that \ddot{x} = a.
    subroutine calc_Efield(xptr, t, level, ctx, Eptr)
      use module_pepc
      use iso_c_binding, only : c_ptr, c_int, c_f_pointer
      use module_debug, only : dbg, DBG_STATS, pepc_status
      implicit none
      type(c_ptr),    intent(in), value :: xptr, Eptr, ctx
      real(pfdp),     intent(in)        :: t
      integer(c_int), intent(in)        :: level

      type(t_particle), allocatable :: particles(:)
      type(app_data_t), pointer :: E,x
      type(level_params_t), pointer :: levelctx
      integer(kind_particle) :: i
      integer, save :: step =0

      call pepc_status('|------> calc_Efield()')
      call ptr_print('x', xptr)
      call ptr_print('E', Eptr)
      call ptr_print('ctx', ctx)

      call c_f_pointer(Eptr, E)
      call c_f_pointer(xptr, x)
      call c_f_pointer(ctx, levelctx)

      step = step + 1

      DEBUG_ASSERT(E%params%nparts==x%params%nparts)

      ! prepare and run PEPC
      allocate(particles(E%params%nparts))
      call encap_to_particles(particles, xptr, ctx)
      call eval_force(particles, levelctx, step, levelctx%comm, clearresults=.true.)

      ! just copy results to E variable (yes, it's weird that we have x- and v-components.. don't ask, won't tell)
      do i=1,size(particles, kind=kind(i))
        E%x(1:E%params%dim, i) = particles(i)%results%e(1:E%params%dim)
        E%x(E%params%dim+1:,i) = 0
        E%v(:,:) = E%x(:,:)
      end do

      deallocate(particles)
    end subroutine calc_Efield


    !> invoke pepc or direct sum
    subroutine eval_force(particles, level_params, step, comm, clearresults)
      use module_debug
      use module_pepc
      use module_interaction_specific, only : intspec_theta2 => theta2
      implicit none
      type(t_particle), allocatable, target, intent(inout) :: particles(:)
      type(level_params_t), intent(in) :: level_params
      integer, intent(in) :: step
      integer(kind_default), intent(in) :: comm
      logical, intent(in) :: clearresults

      if (clearresults) then
          call pepc_particleresults_clear(particles)
      endif

      if (level_params%directforce) then
        call compute_force_direct(particles, comm)
      else
        intspec_theta2 = level_params%theta * level_params%theta
        call pepc_grow_tree(particles)
        call pepc_traverse_tree(particles)
        call pepc_restore_particles(particles)
        if (dbg(DBG_STATS)) call pepc_statistics(step)
        call pepc_timber_tree()
      endif
    end subroutine

    !> compile full right-hand side using the E-field and particle data
    subroutine build_rhs(Eptr, level, ctx, rhsptr)
      use module_pepc
      use iso_c_binding, only : c_ptr, c_int, c_f_pointer
      use module_debug, only : dbg, DBG_STATS, pepc_status
      implicit none
      type(c_ptr),    intent(in), value :: Eptr, rhsptr, ctx
      integer(c_int), intent(in)        :: level
      ! FIXME
    end subroutine

    !> solve for the updated velocity, e.g. using the Boris solver
    subroutine impl_solver(v, level, ctx, v_old, E_old, E_new, SDCint, dt)
      use module_pepc
      use iso_c_binding, only : c_ptr, c_int, c_f_pointer
      use module_debug, only : dbg, DBG_STATS, pepc_status
      implicit none
      type(c_ptr),    intent(in), value :: v, v_old, E_old, E_new, SDCint, ctx
      integer(c_int), intent(in)        :: level
      real(pfdp),     intent(in)        :: dt
      ! FIXME
    end subroutine


    subroutine compute_force_direct(particles, comm)
      use module_pepc_types
      use module_timings
      use module_directsum
      use module_debug, only : pepc_status
      use module_interaction_specific_types, only: t_particle_results
      implicit none
      type(t_particle), intent(inout) :: particles(:) !< input particle data, initializes %x, %data appropriately (and optionally set %label) before calling this function
      integer(kind_default), intent(in) :: comm

      integer(kind_particle) :: i
      type(t_particle_results), allocatable :: directresults(:)

      call pepc_status('DIRECTSUM')

      call timer_start(t_all)
      call directforce(particles, [(i,i=1,size(particles,kind=kind_particle))], size(particles,kind=kind_particle), directresults, comm)
      particles(:)%results = directresults(:)
      deallocate(directresults)
      call timer_stop(t_all)
    end subroutine

end module pfm_feval
