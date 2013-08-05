module pfm_feval
    use pfm_encap
    use pf_mod_mpi
    use module_debug
    implicit none
    private
    
    public feval_init
    public feval_finalize
    public eval_acceleration

contains

    !> Initialize feval, i.e. transfer initial u to PFASST data object
    subroutine feval_init(y0, yend, nlevels, levelctx, encapctx)
      use iso_c_binding
      use pepca_helper, only: read_particles
      implicit none

      type(app_data_t), pointer, intent(inout) :: y0, yend
      integer, intent(in) :: nlevels
      type(c_ptr), intent(in) :: levelctx, encapctx

      type(c_ptr) :: y0_C, yend_c

      type(app_params_t), pointer :: params

      call pepc_status('|------> feval_init()')

      call c_f_pointer(levelctx, params)

      call encap_create(  y0_c, nlevels, -1, 3*(params%n_el+params%n_ion), [-1], levelctx, encapctx) ! 3 coordinates per particle
      call encap_create(yend_c, nlevels, -1, 3*(params%n_el+params%n_ion), [-1], levelctx, encapctx) ! 3 coordinates per particle

      call c_f_pointer(  y0_c, y0)
      call c_f_pointer(yend_c, yend)

      ! set initial values for particle positions and velocities in y0
       call read_particles(y0%particles, 'E_phase_space.dat', params%n_el, 'I_phase_space.dat', params%n_ion)

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
    subroutine eval_acceleration(xptr, t, level, ctx, aptr)
      use module_pepc
      use iso_c_binding, only : c_ptr, c_int, c_f_pointer
      use module_debug, only : dbg, DBG_STATS, pepc_status
      use pepca_units
      implicit none
      type(c_ptr),    intent(in), value :: xptr, aptr, ctx
      real(pfdp),     intent(in)        :: t
      integer(c_int), intent(in)        :: level
      
      type(app_data_t), pointer :: x, a
      integer(kind_particle) :: i
      
      call pepc_status('|------> eval_acceleration()')
      
      call c_f_pointer(xptr, x)
      call c_f_pointer(aptr, a)

      ! FIXME: do we really store the accelerations in an app_data_t? - this appears to be strange
      call pepc_particleresults_clear(x%particles)
      call pepc_grow_tree(x%particles)
      call pepc_traverse_tree(x%particles)
      
      call pepc_restore_particles(x%particles) ! FIXME: can we avoid this?

      !FIXME if (dbg(DBG_STATS)) call pepc_statistics(step)

      call pepc_timber_tree()
      
      ! compute accelerations from fields
      DEBUG_ASSERT(size(x%particles) == size(a%particles))

      do i=1,size(x%particles, kind=kind(i)) ! FIXME: tis is far too ugly
        a%particles(i)      = x%particles(i)
        ! force from field
        a%particles(i)%results%e = a%particles(i)%data%q/unit_4piepsilon0 * a%particles(i)%results%e
        ! acceleration from force
        a%particles(i)%results%e = a%particles(i)%results%e / ( a%particles(i)%data%m * unit_c )
      end do

    end subroutine eval_acceleration
       

end module pfm_feval
