module pfm_hooks
  implicit none
  save
  private
  
  public track_energy_hook
  public compare_checkpoint_hook
  public dump_particles_hook
  
contains

  !> track kinetic, potential, total energy
  subroutine track_energy_hook(pf, level, state, ctx)
    use pf_mod_dtype
    use pfm_encap
    use iso_c_binding
    use pfm_feval
    use module_pepc_types
    use pepca_diagnostics
    use module_debug
    implicit none
    type(pf_pfasst_t), intent(inout) :: pf
    type(pf_level_t),  intent(inout) :: level
    type(pf_state_t),  intent(in)    :: state
    type(c_ptr),       intent(in)    :: ctx

    type(t_particle), allocatable, target :: particles(:)
    type(app_params_t), pointer :: levelctx
    real*8 :: energies(E_MAXIDX), delen
    real(pfdp) ::  t

    call pepc_status('------------- track_energy_hook')
    
    t = state%t0+state%dt ! yes, this is OK, no multiplication with step as t0 is automatically updated during each step
    
    call c_f_pointer(ctx, levelctx)
    
    allocate(particles(levelctx%nparts))
    call encap_to_particles(particles, level%qend, ctx)
    
    ! for potential energy we will need the potential, i.e. have to call pepc again...
    ! FIXME: this is not very fortunate, should join this with acceleration function in pfm_feval
    call eval_force(particles, levelctx%directforce, state%step, levelctx%comm, clearresults=.true.)
    
    call diagnose_energy(particles, energies, state%step+1, t, levelctx%comm, levelctx%root)
    
    delen = abs(energies(E_TOT)-levelctx%initial_energies(E_TOT))/abs(levelctx%initial_energies(E_TOT))
    
    if (levelctx%root) then
      write(*, '(" hook: ",i3," step: ",i5, " t=", es10.3," iter: ",i3," dH: ",es14.7)') state%hook, state%step+1, t,state%iter, delen
    endif

  end subroutine

    
  !> compare particle set to some checkpoint produced with the classical verlet integrator
  subroutine compare_checkpoint_hook(pf, level, state, ctx)
    use pf_mod_dtype
    use pf_mod_hooks
    use pfm_encap
    use iso_c_binding
    use pfm_feval
    use module_pepc_types
    use pepca_diagnostics
    use module_debug
    use pepca_units, only: unit_time_as_per_simunit
    implicit none
    type(pf_pfasst_t), intent(inout) :: pf
    type(pf_level_t),  intent(inout) :: level
    type(pf_state_t),  intent(in)    :: state
    type(c_ptr),       intent(in)    :: ctx

    type(t_particle), allocatable, target :: particles(:)
    type(app_params_t), pointer :: levelctx
    real*8 :: xerr, verr
    real(pfdp) ::  t
    integer(kind_default) :: itime_in

    call pepc_status('------------- track_energy_hook')
    
    call c_f_pointer(ctx, levelctx)
    
    allocate(particles(levelctx%nparts))
    
    ! compare particle data to checkpoint of some previous run - we use the actual physical simulation time to identify the appropriate checkpoint
    select case (state%hook)
      case (PF_PRE_ITERATION)
        call encap_to_particles(particles, level%qend, ctx)
        t = state%t0 ! yes, this is OK, no multiplication with step as t0 is automatically updated during each step
      case (PF_POST_ITERATION)
        call encap_to_particles(particles, level%qend, ctx)
        t = state%t0+state%dt ! yes, this is OK, no multiplication with step as t0 is automatically updated during each step
      case default
        DEBUG_ERROR(*,'wrong hook')
    end select
    
    itime_in = aint(t*unit_time_as_per_simunit) ! compare pepc.f90, line 144: dumpstep = ...
    call compare_particles_to_checkpoint(particles, itime_in, levelctx%comm, xerr, verr)

    if (levelctx%root) then
      write(*,'(" hook: ",i3," step: ",i5," t=", es10.3, " iter: ",i3," Ex: ",es14.7," Ev: ",es14.7)') &
               state%hook, state%step+1, t,state%iter, xerr, verr
    endif
    
  end subroutine

    
  !> particle output
  subroutine dump_particles_hook(pf, level, state, ctx)
    use pf_mod_dtype
    use pf_mod_hooks
    use pfm_encap
    use iso_c_binding
    use pfm_feval
    use module_pepc_types
    use pepca_diagnostics
    use module_debug
    use pepca_units
    implicit none
    type(pf_pfasst_t), intent(inout) :: pf
    type(pf_level_t),  intent(inout) :: level
    type(pf_state_t),  intent(in)    :: state
    type(c_ptr),       intent(in)    :: ctx

    type(t_particle), allocatable, target :: particles(:)
    type(app_params_t), pointer :: levelctx
    real*8 :: xerr, verr
    real(pfdp) ::  t
    integer(kind_default) :: itime

    call pepc_status('------------- track_energy_hook')
    
    t = state%t0+state%dt
    
    call c_f_pointer(ctx, levelctx)
    
    allocate(particles(levelctx%nparts))

    select case (state%hook)
      case (PF_PRE_ITERATION)
        call encap_to_particles(particles, level%qend, ctx)
        itime = (level%nnodes-1)*(state%step  )
      case (PF_POST_ITERATION)
        call encap_to_particles(particles, level%qend, ctx)
        itime = (level%nnodes-1)*(state%step+1)
      case default
        DEBUG_ERROR(*,'wrong hook')
    end select
    
    ! do diagnostics etc here
    call write_particles_vtk(particles, itime, (level%nnodes-1)*state%nsteps, t*unit_time_fs_per_simunit, levelctx%comm)
    !if (dumpnow(pepca_nml%output_interval(OI_PARTICLES_ASC), step, nt)) call write_particles_ascii(pepca_nml%rank, step, particles, filename)
    !if (dumpnow(pepca_nml%output_interval(OI_PARTICLES_MPI), step, nt)) call write_particles_mpiio(MPI_COMM_SPACE, step, 2*pepca_nml%numparts_total, particles, filename)
    
  end subroutine  
end module pfm_hooks