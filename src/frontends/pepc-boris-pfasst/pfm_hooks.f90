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
    implicit none
    type(pf_pfasst_t), intent(inout) :: pf
    type(pf_level_t),  intent(inout) :: level
    type(pf_state_t),  intent(in)    :: state
    type(c_ptr),       intent(in)    :: ctx

    type(t_particle), allocatable, target :: particles(:)
    type(level_params_t), pointer :: levelctx
    real(pfdp) ::  t

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
      case (PF_PRE_STEP)
        if (did_prestep) return

        did_prestep=.true.
        call encap_create(encaptmp, level%level, -1, level%nvars, level%shape, ctx, level%encap%encapctx)
        call encap_unpack(encaptmp, level%q0)
        call encap_to_particles(particles, encaptmp, ctx)
        call encap_destroy(encaptmp)
        t = state%t0 ! yes, this is OK, t0 is automatically updated during each step
      case default
        DEBUG_ERROR(*,'wrong hook')
    end select

    if (levelctx%root) then
      if (state%step == 0) write(*,*)
      write(*,'(a1, a,"| step: ",i5," t=", es10.3, " iter: ",i3, " === dumping particles if requested...")',advance='no') &
        char(13), hook_names(state%hook), state%step+1, t, state%iter
   endif

    ! do diagnostics etc here
    call dump_particles(t, particles, pepcboris_nml%workingmode + IFILE_SUMMAND, do_average=.false.)
    call dump_energy(t, particles, pepcboris_nml%workingmode + IFILE_SUMMAND_ENERGY, &
      levelctx, pepcboris_nml, levelctx%comm, do_average=.false.)

  end subroutine
end module pfm_hooks