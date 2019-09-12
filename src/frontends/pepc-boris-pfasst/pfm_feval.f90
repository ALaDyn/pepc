module pfm_feval
    use module_pepc_kinds
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
    public dump_nfeval

    integer, parameter, public :: FEVAL_MODE_NO_EXTERNAL_FIELD   = 0 * 10
    integer, parameter, public :: FEVAL_MODE_FULL_EXTERNAL_FIELD = 1 * 10

    integer, parameter, public :: FEVAL_MODE_NO_INTERNAL_FIELD   = 0 * 1
    integer, parameter, public :: FEVAL_MODE_FULL_INTERNAL_FIELD = 1 * 1

    integer :: num_feval_internal = 0
    integer :: num_feval_external = 0

  contains

    !> dump number of rhs evaluations
    subroutine dump_nfeval(rank, comm)
      use module_pepc_types
      use pepcboris_helper, only : IFILE_SUMMAND_NFEVAL, pepcboris_nml
      use mpi
      implicit none
      integer(kind_pe), intent(in) :: rank
      integer(kind_default), intent(in) :: comm
      integer :: istream
      integer(kind_default) :: mpi_err
      integer :: num_feval_internal_global, num_feval_external_global

      call MPI_REDUCE(num_feval_internal, num_feval_internal_global, 1, MPI_INTEGER, MPI_SUM, 0, comm, mpi_err)
      call MPI_REDUCE(num_feval_external, num_feval_external_global, 1, MPI_INTEGER, MPI_SUM, 0, comm, mpi_err)

      istream = pepcboris_nml%workingmode + IFILE_SUMMAND_NFEVAL

      if (rank == 0) then
        write(istream,'("# overall total number of rhs evaluations (internal and external forces)",/,i0,x,i0)') num_feval_internal, num_feval_external
        write(*,'(/,"overall total num_feval_internal / external: ",i0," / ",i0)') num_feval_internal, num_feval_external
      endif
    end subroutine


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
      use pepcboris_helper
      implicit none
      type(c_ptr),    intent(in), value :: xptr, Eptr, ctx
      real(pfdp),     intent(in)        :: t
      integer(c_int), intent(in)        :: level

      type(t_particle), allocatable :: particles(:)
      type(app_data_t), pointer :: E,x
      type(level_params_t), pointer :: levelctx
      integer(kind_particle) :: i

      call pepc_status('|------> calc_Efield()')
      call c_f_pointer(Eptr, E)
      call c_f_pointer(xptr, x)
      call c_f_pointer(ctx, levelctx)

      DEBUG_ASSERT(E%params%nparts==x%params%nparts)

      ! prepare and run PEPC
      allocate(particles(E%params%nparts))
      call encap_to_particles(particles, xptr, ctx)
                                           ! FIXME: using pepcboris_nml as global variable is not really nice
      call eval_force(particles, levelctx, pepcboris_nml, levelctx%comm, clearresults=.true.)

      ! just copy results to E variable (yes, it's weird that we have x- and v-components.. don't ask, won't tell)
      do i=1,size(particles, kind=kind(i))
        E%x(1:E%params%dim, i) = particles(i)%results%e(1:E%params%dim)
        E%x(E%params%dim+1:,i) = 0
        E%v(:,i) = E%x(:,i)
      end do

      deallocate(particles)
    end subroutine calc_Efield


    !> invoke pepc or direct sum
    subroutine eval_force(particles, level_params, nml, comm, clearresults)
      use module_debug
      use module_pepc
      use module_interaction_specific, only : intspec_theta2 => theta2
      use pepcboris_helper
      implicit none
      type(t_particle), allocatable, target, intent(inout) :: particles(:)
      type(level_params_t), intent(in) :: level_params
      type(pepcboris_nml_t), intent(in) :: nml
      integer(kind_default), intent(in) :: comm
      logical, intent(in) :: clearresults

      integer :: feval_mode_internal, feval_mode_external

      feval_mode_internal = modulo(level_params%feval_mode, 10)
      feval_mode_external = int(level_params%feval_mode/10)*10

      if (clearresults) then
          call pepc_particleresults_clear(particles)
      endif

      select case (feval_mode_internal)
        case (FEVAL_MODE_FULL_INTERNAL_FIELD)
          num_feval_internal = num_feval_internal + 1

          if (level_params%directforce) then
            call compute_force_direct(particles, comm)
          else
            intspec_theta2 = level_params%theta * level_params%theta
            call pepc_grow_tree(particles)
            call pepc_traverse_tree(particles)
            call pepc_restore_particles(particles)
            call pepc_timber_tree()
          endif
        case (FEVAL_MODE_NO_INTERNAL_FIELD)
          ! simply do nothing
        case default
          DEBUG_ERROR(*,'feval_mode_internal invalid: ', feval_mode_internal)
      end select

      select case (feval_mode_external)
        case (FEVAL_MODE_FULL_EXTERNAL_FIELD)
          num_feval_external = num_feval_external + 1

          call apply_external_field(nml, particles)
        case (FEVAL_MODE_NO_INTERNAL_FIELD)
          ! simply do nothing
        case default
          DEBUG_ERROR(*,'feval_mode_external invalid: ', feval_mode_external)
      end select

    end subroutine

    !> compile full right-hand side using the E-field and particle data
    subroutine build_rhs(Eptr, Qptr, level, ctx, rhsptr)
      use module_pepc
      use module_pepc_kinds
      use iso_c_binding, only : c_ptr, c_int, c_f_pointer
      use module_debug, only : dbg, DBG_STATS, pepc_status
      use pepcboris_helper
      implicit none
      type(c_ptr),    intent(in), value :: Eptr, rhsptr, ctx, Qptr
      integer(c_int), intent(in)        :: level
      real(kind_physics) :: B0(3)

      type(t_particle), allocatable :: particles(:)
      type(app_data_t), pointer :: E,rhs,Q
      type(level_params_t), pointer :: levelctx
      integer(kind_particle) :: i

      call pepc_status('|------> build_rhs()')
      call c_f_pointer(Eptr, E)
      call c_f_pointer(Qptr, Q)
      call c_f_pointer(rhsptr, rhs)
      call c_f_pointer(ctx, levelctx)

      DEBUG_ASSERT(E%params%nparts==rhs%params%nparts)
      DEBUG_ASSERT(E%params%dim==rhs%params%dim)

      allocate(particles(E%params%nparts))
      call encap_to_particles(particles, Eptr, ctx)

      B0 = get_magnetic_field()

      do i=1,size(particles, kind=kind(i))
        rhs%x(1:rhs%params%dim, i) = particles(i)%data%q*cross_prod_plus(Q%v(1:Q%params%dim,i),B0,E%x(1:E%params%dim,i)) / particles(i)%data%m
        rhs%x(rhs%params%dim+1:,i) = 0
        rhs%v(:,i) = rhs%x(:,i)
      end do

      deallocate(particles)

    end subroutine

    !> solve for the updated velocity, e.g. using the Boris solver
    subroutine impl_solver(vptr, level, ctx, v_oldptr, E_oldptr, E_newptr, SDCintptr, dt)
      use module_pepc
      use module_pepc_kinds
      use iso_c_binding, only : c_ptr, c_int, c_f_pointer
      use module_debug, only : dbg, DBG_STATS, pepc_status
      use pepcboris_helper
      implicit none
      type(c_ptr),    intent(in), value :: vptr, v_oldptr, E_oldptr, E_newptr, SDCintptr, ctx
      integer(c_int), intent(in)        :: level
      real(pfdp),     intent(in)        :: dt

      real(kind_physics)               :: beta, gam
      real(kind_physics), dimension(3) :: B0, uminus, uprime, uplus, t, s, Emean

      type(t_particle), allocatable :: particles(:)
      type(app_data_t), pointer :: v, v_old, E_new, E_old, SDCint
      type(level_params_t), pointer :: levelctx
      integer(kind_particle) :: i

      call pepc_status('|------> impl_solver()')
      call c_f_pointer(vptr, v)
      call c_f_pointer(v_oldptr, v_old)
      call c_f_pointer(E_newptr, E_new)
      call c_f_pointer(E_oldptr, E_old)
      call c_f_pointer(SDCintptr, SDCint)
      call c_f_pointer(ctx, levelctx)

      B0 = get_magnetic_field()

      allocate(particles(E_new%params%nparts))
      call encap_to_particles(particles, E_newptr, ctx)

      do i = 1,size(particles, kind=kind(i))
         ! charge/mass*time-constant
        beta   = particles(i)%data%q / (2._kind_physics * particles(i)%data%m) * dt
        Emean(:)  = (E_old%v(:,i) + E_new%v(:,i)) / 2._kind_physics
        ! first half step with electric field
        uminus(:) = v_old%v(:,i) + beta * Emean(:) + SDCint%v(:,i) / 2._kind_physics
        ! gamma factor
        !gam    = sqrt( 1._kind_physics + ( dot_product(uminus, uminus) ) / unit_c2 )
        gam    = 1._kind_physics
        ! rotation with magnetic field
        t      = beta/gam * B0
        uprime = cross_prod_plus(uminus, t, uminus)
        s      = 2._kind_physics * t / (1._kind_physics + dot_product(t, t))
        uplus  = cross_prod_plus(uprime, s, uminus)
        ! second half step with electric field
        v%v(:,i) = uplus(:) + beta * Emean(:) + SDCint%v(:,i) / 2._kind_physics
      end do

      deallocate(particles)

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
      allocate(directresults(size(particles,kind=kind_particle)))
      call directforce(particles, [(i,i=1,size(particles,kind=kind_particle))], size(particles,kind=kind_particle), directresults, comm)
      particles(:)%results = directresults(:)
      deallocate(directresults)
      call timer_stop(t_all)
    end subroutine

    subroutine apply_external_field(nml, particles)
      use module_pepc_kinds
      use module_pepc_types
      use pepcboris_helper
      implicit none
      type(t_particle), target, intent(inout) :: particles(:)
      type(pepcboris_nml_t), intent(in) :: nml
      integer(kind_particle) :: p
      real(kind_physics), parameter :: w(3) = [1._kind_physics, 1._kind_physics, -2._kind_physics]
      real(kind_physics) :: w2

      w2 = nml%setup_params(PARAMS_OMEGAE)**2

      do p=1,size(particles,kind=kind(p))
        associate(pt => particles(p))
          pt%results%e   = pt%results%e   + w2 * w*pt%x
          pt%results%pot = pt%results%pot - w2*pt%data%m/pt%data%q/2._kind_physics * dot_product(w,pt%x*pt%x)
        end associate
      end do
    end subroutine

end module pfm_feval
