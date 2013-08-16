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
    subroutine feval_init(y0, yend, nlevels, levelctx, encapctx) ! TODO: shouldnt we call this function once per level?
      use iso_c_binding
      use module_pepc_types, only: t_particle
      implicit none

      type(app_data_t), pointer, intent(inout) :: y0, yend
      integer, intent(in) :: nlevels
      type(c_ptr), intent(inout) :: levelctx, encapctx
 
      type(c_ptr) :: y0_C, yend_c

      type(app_params_t), pointer :: params

      call pepc_status('|------> feval_init()')

      call c_f_pointer(levelctx, params)

      call encap_create(  y0_c, nlevels, -1, (2*params%dim)*(params%n_el+params%n_ion), [-1], levelctx, encapctx) ! dim*(coordinates and momenta) per particle
      call encap_create(yend_c, nlevels, -1, (2*params%dim)*(params%n_el+params%n_ion), [-1], levelctx, encapctx) ! dim*(coordinates and momenta) per particle

      call c_f_pointer(  y0_c, y0)
      call c_f_pointer(yend_c, yend)
     
      call particles_to_encap(y0_c, params%particles)

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
      
      type(t_particle), allocatable :: particles(:)
      type(app_data_t), pointer :: a
      integer(kind_particle) :: i
      
      call pepc_status('|------> eval_acceleration()')
      call ptr_print('x', xptr)
      call ptr_print('a', aptr)
     
      call c_f_pointer(aptr, a)
      
      ! prepare and run PEPC
      allocate(particles(a%params%n_el+a%params%n_ion))
      call encap_to_particles(particles, xptr, ctx)   ! TODO: ctx might be a level context instead of an encap context...
      call pepc_particleresults_clear(particles)
      call pepc_grow_tree(particles)
      call pepc_traverse_tree(particles)
      call pepc_restore_particles(particles)
      !FIXME if (dbg(DBG_STATS)) call pepc_statistics(step)
      call pepc_timber_tree()
      
      ! compute accelerations from fields
      do i=1,size(particles, kind=kind(i))
        ! acceleration from force from field
        a%x(1:a%params%dim, i) = particles(i)%data%q/unit_4piepsilon0 * particles(i)%results%e(1:a%params%dim) / ( particles(i)%data%m * unit_c )
        a%x(a%params%dim+1:,i) = 0
        a%v(:,:) = a%x(:,:)
      end do
      
      deallocate(particles)

    end subroutine eval_acceleration
       

end module pfm_feval
