module pfm_hooks
  implicit none
  save
  private
  
  public track_energy_hook
  
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
    integer(kind_default) :: itime_in

    call pepc_status('------------- track_energy_hook')
    
    t = state%t0+state%dt
    
    call c_f_pointer(ctx, levelctx)
    
    allocate(particles(levelctx%nparts))
    call encap_to_particles(particles, level%qend, ctx)
    
    ! for potential energy we will need the potential, i.e. have to call pepc again...
    ! FIXME: this is not very fortunate, should join this with acceleration function in pfm_feval
    call eval_force(particles, levelctx%directforce, state%step, levelctx%comm, clearresults=.true.)
    
    call diagnose_energy(particles, energies, state%step+1, t, levelctx%comm, levelctx%root)
    
    delen = abs(energies(E_TOT)-levelctx%initial_energies(E_TOT))/abs(levelctx%initial_energies(E_TOT))
    
    if (levelctx%root) then
      write(*, '(" step: ",i5," t=", es10.3," iter: ",i3," dH: ",es14.7)') state%step+1, t,state%iter, delen
    endif
    
    ! compare particle data to checkpoint of some previous run
    ! currently we are more or less aligning leap-frog steps to sdc nodes (at least wrt their number)
    itime_in = (level%nnodes-1)*state%step
    call compare_particles_to_checkpoint(particles, itime_in, levelctx%comm)

    
  end subroutine 
 
end module pfm_hooks
