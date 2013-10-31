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
    implicit none
    type(pf_pfasst_t), intent(inout) :: pf
    type(pf_level_t),  intent(inout) :: level
    type(pf_state_t),  intent(in)    :: state
    type(c_ptr),       intent(in)    :: ctx

    type(t_particle), allocatable, target :: particles(:)
    type(level_params_t), pointer :: levelctx
    real(pfdp) ::  t
    integer(kind_particle) :: p

    call pepc_status('------------- dump_particles_hook')

    t = state%t0+state%dt

    call c_f_pointer(ctx, levelctx)

    allocate(particles(levelctx%nparts))

    select case (state%hook)
      case (PF_POST_STEP)
        call encap_to_particles(particles, level%qend, ctx)
        t = state%t0 + state%dt ! yes, this is OK, no multiplication with step as t0 is automatically updated during each step
      case default
        DEBUG_ERROR(*,'wrong hook')
    end select

    if (levelctx%root) then
      write(*,'(a,"| step: ",i5," t=", es10.3, " iter: ",i3, " === dumping particles if requested...")') hook_names(state%hook), state%step+1, t,state%iter
    endif

    ! do diagnostics etc here
    do p=1,size(particles,kind=kind(p))
      write(48,*) t, p, particles(p)%x, particles(p)%data%v
    end do
  end subroutine
end module pfm_hooks
