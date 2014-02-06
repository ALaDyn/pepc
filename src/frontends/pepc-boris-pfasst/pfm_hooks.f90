module pfm_hooks
  implicit none
  save
  private

  public dump_particles_hook

contains

  !> particle output
  subroutine dump_particles_hook(pf, level, state, ctx)
    use pf_mod_dtype
    use pf_mod_hooks
    use pfm_encap
    use iso_c_binding
    use pfm_feval
    use module_pepc_types
    use module_debug
    use pepcboris_helper
    use pepcboris_diagnostics
    use module_vtk, only: VTK_STEP_NORMAL
    implicit none
    type(pf_pfasst_t), intent(inout) :: pf
    type(pf_level_t),  intent(inout) :: level
    type(pf_state_t),  intent(in)    :: state
    type(c_ptr),       intent(in)    :: ctx

    type(t_particle), allocatable, target :: particles(:)
    type(level_params_t), pointer :: levelctx
    real(pfdp) ::  t
    integer :: step

    logical, save :: did_prestep = .false.

    type(c_ptr) :: encaptmp

    call pepc_status('------------- dump_particles_hook')

    t = state%t0+state%dt

    call c_f_pointer(ctx, levelctx)

    allocate(particles(levelctx%nparts))

    select case (state%hook)
      case (PF_POST_STEP)
        call encap_to_particles(particles, level%qend, ctx)
        t = state%t0 + state%dt ! yes, this is OK, no multiplication with step as t0 is automatically updated during each step
        step = state%step + 1
      case (PF_PRE_STEP)
        if (did_prestep) return

        did_prestep=.true.
        call encap_create(encaptmp, level%level, -1, level%nvars, level%shape, ctx, level%encap%encapctx)
        call encap_unpack(encaptmp, level%q0)
        call encap_to_particles(particles, encaptmp, ctx)
        call encap_destroy(encaptmp)
        t = state%t0 ! yes, this is OK, t0 is automatically updated during each step
        step = state%step
      case default
        DEBUG_ERROR(*,'wrong hook')
    end select

    if (levelctx%root_stdio) then
      if (state%step == 0) write(*,*)
      write(*,'(a1, a,"| step: ",i0,"/",i0," t=", es10.3, " iter: ",i10, " residual: ", g15.6)',advance='no') &
        char(13), hook_names(state%hook), step, state%nsteps, t, state%iter, level%residual
    endif

    if (levelctx%root_file) then
      ! do diagnostics etc here
      ! FIXME: this output has to synchronized somehow between the different time ranks
      call dump_particles(VTK_STEP_NORMAL, step, state%dt, particles, levelctx%comm, do_average=.false.)
      call dump_energy(t, particles, levelctx, levelctx%comm, do_average=.false.)
      call dump_iterations(step, state%dt, state%hook, state%iter, level%residual)
   endif

  end subroutine
end module pfm_hooks
