module pfm_hooks
  implicit none
  save
  private

  public dump_particles_hook
  public constrain_particles_hook
  public dump_iterstate_hook

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
    use pepc_helper
    use dump_helper
    use module_vtk, only: VTK_STEP_NORMAL
    implicit none
    type(pf_pfasst_t), intent(inout) :: pf
    type(pf_level_t),  intent(inout) :: level
    type(pf_state_t),  intent(in)    :: state
    type(c_ptr),       intent(in)    :: ctx

    type(t_particle), allocatable, target :: particles(:)
    type(level_params_t), pointer :: levelctx
     integer :: step

    logical, save :: did_prestep = .false.

    type(c_ptr) :: encaptmp

    call pepc_status('------------- dump_particles_hook')

    call c_f_pointer(level%levelctx, levelctx)

    allocate(particles(levelctx%nparts))

    select case (state%hook)
      case (PF_POST_STEP)
        call encap_to_particles(particles, level%qend, level%levelctx)
        step = state%step + 1
      case (PF_PRE_STEP)
        if (did_prestep) return

        did_prestep=.true.
        call encap_create(encaptmp, level%level, -1, level%nvars, level%shape, level%levelctx, level%encap%encapctx)
        call encap_unpack(encaptmp, level%q0)
        call encap_to_particles(particles, encaptmp, level%levelctx)
        call encap_destroy(encaptmp)
        step = state%step
      case default
        DEBUG_ERROR(*,'wrong hook')
    end select

    ! do diagnostics etc here
    call perform_all_dumps(step, levelctx%pepc_pars, levelctx%physics_pars, levelctx%time_pars, levelctx%field_grid, particles)

    if (levelctx%pepc_pars%pepc_comm%root_file) then
      call dump_iterations(step, state%dt, state%hook, state%iter, level%residual, levelctx%pepc_pars)
   endif

   deallocate(particles)

  end subroutine



  !> output step, resiude, etc.
  subroutine dump_iterstate_hook(pf, level, state, ctx)
    use pf_mod_dtype
    use pf_mod_hooks
    use pfm_encap
    use iso_c_binding
    use pfm_feval
    use module_pepc_types
    use module_debug
    use pepc_helper
    use dump_helper
    implicit none
    type(pf_pfasst_t), intent(inout) :: pf
    type(pf_level_t),  intent(inout) :: level
    type(pf_state_t),  intent(in)    :: state
    type(c_ptr),       intent(in)    :: ctx

    type(level_params_t), pointer :: levelctx
    real(pfdp) ::  t
    integer :: step

    call pepc_status('------------- dump_iterstate_hook')

    call c_f_pointer(level%levelctx, levelctx)

    select case (state%hook)
      case (PF_POST_STEP, PF_POST_ITERATION)
        t = state%t0 + state%dt ! yes, this is OK, no multiplication with step as t0 is automatically updated during each step
        step = state%step + 1
      case (PF_PRE_STEP)
        t = state%t0+ state%dt ! yes, this is OK, t0 is automatically updated during each step
        step = state%step + 1
        ! FIXME: the residual still contains the value from the prev iteration
        level%residual = -1.0
      case default
        DEBUG_ERROR(*,'wrong hook')
    end select

    if (levelctx%pepc_pars%pepc_comm%root_stdio) then
      call debug_print_timestamp(6)
      write(*,'(x, a,"| step: ",i0,"/",i0," t=", es10.3, " iter: ",i10, " residual: ", g15.6)',advance='yes') &
        hook_names(state%hook), step, state%nsteps, t, state%iter, level%residual
    endif

  end subroutine


  !> folding particles back into the periodic simulation box
  subroutine constrain_particles_hook(pf, level, state, ctx)
    use pf_mod_dtype
    use pf_mod_hooks
    use pfm_encap
    use iso_c_binding
    use module_pepc_types
    use module_debug
    use pepc_helper
    use time_helper
    implicit none
    type(pf_pfasst_t), intent(inout) :: pf
    type(pf_level_t),  intent(inout) :: level
    type(pf_state_t),  intent(in)    :: state
    type(c_ptr),       intent(in)    :: ctx

    type(t_particle), allocatable, target :: particles(:)
    type(level_params_t), pointer :: levelctx

    call pepc_status('------------- constrain_particles_hook')

    select case (state%hook)
      case (PF_POST_STEP)
        call c_f_pointer(level%levelctx, levelctx)
        allocate(particles(levelctx%nparts))
        call encap_to_particles(particles, level%qend, level%levelctx)
        call constrain_particles(levelctx%physics_pars, particles)
        call particles_to_encap(level%qend, particles)
        deallocate(particles)
      case default
        DEBUG_ERROR(*,'wrong hook')
    end select
  end subroutine


  subroutine dump_iterations(step, dt, hook, niter, residual, pepc_pars)
    use encap
    use pepcboris_paralleldump
    implicit none
    integer, intent(in) :: step, niter, hook
    real*8, intent(in) :: residual
    real*8, intent(in) :: dt
    type(pepc_pars_t), intent(in) :: pepc_pars

    character(len=PARALLELDUMP_MAXLEN) :: line

    write(line,*) step*dt, step, hook, niter, residual
    call paralleldump_dump(pepc_pars%workingmode + IFILE_SUMMAND_NITER, line)
  end subroutine

end module pfm_hooks