module pfm_encap
  use iso_c_binding
  use module_pepc_types
  use module_debug
  use pf_mod_dtype
  implicit none

  !> data type for level-dependent application parameters
  !> if any parameters are added here, they also have to be included in #pfm_setup_solver_level_params
  type :: level_params_t
    integer :: level
    integer :: nparts
    integer :: dim
    real*8 :: theta
    logical :: directforce
    integer :: feval_mode
    integer(kind_default) :: comm
  end type

  !> Data encapsulation: data and parameters which will be filled in encap_create using ctx
  type :: app_data_t
    real*8, dimension(:,:), allocatable :: x
    real*8, dimension(:,:), allocatable :: v
    type(level_params_t) :: params
    type(t_particle), pointer :: particles(:) !< this is a pointer to blueprint particles: x and v are stored in app_data_t structure, all other properties are stores here, see encap_to_particles() for details
  end type app_data_t

  integer, parameter :: AXPY_aVpV_aXpX = 0
  integer, parameter :: AXPY_aVpV      = 1
  integer, parameter :: AXPY_aXpX      = 2
  integer, parameter :: AXPY_aVpX      = 12

  integer, parameter :: COPY_XV = 0
  integer, parameter :: COPY_V  = 1
  integer, parameter :: COPY_X  = 2

  integer, parameter :: SETVAL_XV = 0
  integer, parameter :: SETVAL_V  = 1
  integer, parameter :: SETVAL_X  = 2

contains

  subroutine particles_to_encap(ptrenc, p)
    use module_pepc_types, only: t_particle
    implicit none
    type(c_ptr), intent(in) :: ptrenc
    type(t_particle), intent(in) :: p(:)

    integer(kind_particle) :: i
    type(app_data_t), pointer :: enc

    call pepc_status('|----> particles_to_encap()')

    call c_f_pointer(ptrenc, enc)

    DEBUG_ASSERT(enc%params%nparts==size(p))

    do i=1,size(p)
      enc%x(1:enc%params%dim, i) = p(i)%x(     1:enc%params%dim)
      enc%v(1:enc%params%dim, i) = p(i)%data%v(1:enc%params%dim)
    end do
  end subroutine

  subroutine encap_to_particles(p, ptrenc, levelctx)
    use module_pepc_types, only: t_particle
    implicit none
    type(t_particle), intent(out) :: p(:)
    type(c_ptr), intent(in) :: ptrenc
    type(c_ptr), intent(in) :: levelctx

    integer(kind_particle) :: i
    type(app_data_t), pointer :: enc
    type(level_params_t), pointer :: levelparams

    call pepc_status('|----> encap_to_particles()')

    call c_f_pointer(ptrenc, enc)
    call c_f_pointer(levelctx, levelparams)

    DEBUG_ASSERT(enc%params%nparts==size(p))
    DEBUG_ASSERT(enc%params%nparts==size(enc%particles))

    do i=1,size(p)
      p(i)                          = enc%particles(i) ! FIXME: here we set coordinates and velocities but overwrite them again in the next lines
      p(i)%x(     1:enc%params%dim) = enc%x(1:enc%params%dim, i)
      p(i)%x(enc%params%dim+1:)     = 0.
      p(i)%data%v(1:enc%params%dim) = enc%v(1:enc%params%dim, i)
      p(i)%x(enc%params%dim+1:)     = 0.
    end do
  end subroutine

  !> internal helper for getting a pointer value (i.e. the address of the pointee)
  function ptr_val(ptr)
    use iso_c_binding, only : c_ptr
    implicit none
    type(c_ptr), intent(in) :: ptr
    integer :: ptr_val

    ptr_val = transfer (ptr, ptr_val)
  end function


  !> Fill pf_encap_t with pointers to encapsulation functions
  subroutine pfm_encap_init(encap, particles)
    use module_pepc_types, only: t_particle
    implicit none
    type(pf_encap_t), intent(out) :: encap
    type(t_particle), target, intent(inout) :: particles(:)

    encap%encapctx  = c_loc(particles(1)) !< this is a pointer to blueprint particles: x and v are stored in app_data_t structure, all other properties are stores here, see encap_to_particles() for details
    encap%create    => encap_create
    encap%destroy   => encap_destroy
    encap%setval    => encap_setval
    encap%copy      => encap_copy
    encap%pack      => encap_pack
    encap%unpack    => encap_unpack
    encap%axpy      => encap_axpy
    encap%norm      => encap_norm
  end subroutine pfm_encap_init


  !> Allocate/create solution (spatial data set) for the given level.
  !> This is called for each SDC node.
  subroutine encap_create(sol, level, kind, nvars, shape, levelctx, encapctx)
    implicit none
    type(c_ptr),       intent(inout)     :: sol
    integer,           intent(in)        :: level, nvars, shape(:)
    integer,           intent(in)        :: kind
    type(c_ptr),       intent(in), value :: levelctx, encapctx

    type(app_data_t), pointer :: q
    type(level_params_t), pointer :: p

    !call pepc_status('|----> encap_create()')
    call c_f_pointer(levelctx, p)

    DEBUG_ASSERT(nvars==(2*p%dim)*(p%nparts)) ! dim*(coordinates and momenta) per particle

    allocate(q)
    q%params = p
    call c_f_pointer(encapctx, q%particles, shape=[q%params%nparts])

    allocate(q%x(q%params%dim,q%params%nparts))
    allocate(q%v(q%params%dim,q%params%nparts))

    sol = c_loc(q)

  end subroutine encap_create


  !> Deallocate/destroy solution.
  subroutine encap_destroy(ptr)
    implicit none
    type(c_ptr), intent(in), value :: ptr

    type(app_data_t), pointer :: q

    !call pepc_status('|----> encap_destroy()')
    call c_f_pointer(ptr, q)

    deallocate(q%x)
    deallocate(q%v)
    deallocate(q)

  end subroutine encap_destroy


  !> Set solution value.
  subroutine encap_setval(ptr, val, flags)
    implicit none
    type(c_ptr), intent(in), value    :: ptr
    real(pfdp),  intent(in)           :: val
    integer,     intent(in), optional :: flags

    type(app_data_t), pointer :: q
     integer :: which

    !call pepc_status('|----> encap_setval()')
    call c_f_pointer(ptr, q)

    which = 0
    if (present(flags)) which = flags

    select case (which)
    case (SETVAL_XV)
      ! set value for x and v
      q%x = val
      q%v = val
    case (SETVAL_V)
      ! set value for v
      q%v = val
    case (SETVAL_X)
      ! set value for x
      q%x = val
    case default
       DEBUG_ERROR(*, 'Invalid flags')
    end select

  end subroutine encap_setval


  !> Copy solution value on one level.
  subroutine encap_copy(dstptr, srcptr, flags)
    implicit none
    type(c_ptr), intent(in), value    :: dstptr, srcptr
    integer,     intent(in), optional :: flags

    type(app_data_t), pointer :: dst, src
    integer :: which

    !call pepc_status('|----> encap_copy()')
    call c_f_pointer(dstptr,dst)
    call c_f_pointer(srcptr,src)

    DEBUG_ASSERT(src%params%nparts==dst%params%nparts)

    which = COPY_XV
    if (present(flags)) which = flags

    select case (which)
    case (COPY_XV)
      ! copy value for x and v
      dst%x(:,:) = src%x
      dst%v(:,:) = src%v
    case (COPY_V)
      ! copy value for v
      dst%v(:,:) = src%v
    case (COPY_X)
      ! copy value for x
      dst%x(:,:) = src%x
    case default
       DEBUG_ERROR(*, 'Invalid flags')
    end select

  end subroutine encap_copy


  !> Pack solution q into a flat array.
  subroutine encap_pack(z, ptr)
    implicit none
    type(c_ptr), intent(in), value  :: ptr
    real(pfdp),  intent(out)        :: z(:)  !FIXME: an alternative with a frontend-defined MPI type would be very nice - this will avoid packing completely :-)

    type(app_data_t), pointer :: q
    integer(kind_particle) :: i, j

    call pepc_status('|----> encap_pack()')
    call c_f_pointer(ptr, q)

    j = 1
    do i=1,q%params%nparts
      z(j:j+q%params%dim-1) = q%x(1:q%params%dim, i)
      j = j+q%params%dim
      z(j:j+q%params%dim-1) = q%v(1:q%params%dim, i)
      j = j+q%params%dim
    end do

    DEBUG_ASSERT(j==size(z)+1)

  end subroutine encap_pack


  !> Unpack solution from a flat array.
  subroutine encap_unpack(ptr, z)
    implicit none
    type(c_ptr), intent(in), value :: ptr
    real(pfdp),  intent(in)        :: z(:)  !FIXME: an alternative with a frontend-defined MPI type would be very nice - this will avoid packing completely :-)

    type(app_data_t), pointer :: q
    integer(kind_particle) :: i, j

    call pepc_status('|----> encap_unpack()')
    call c_f_pointer(ptr, q)

    j = 1
    do i=1,q%params%nparts
      q%x(1:q%params%dim, i) = z(j:j+q%params%dim-1)
      j = j+q%params%dim
      q%v(1:q%params%dim, i) = z(j:j+q%params%dim-1)
      j = j+q%params%dim
    end do

    DEBUG_ASSERT(j==size(z)+1)

  end subroutine encap_unpack


  !> Compute y = a x + y where a is a scalar and x and y are solutions on the same level.
  subroutine encap_axpy(yptr, a, xptr, flags)
    implicit none
    type(c_ptr), intent(in), value    :: xptr, yptr
    real(pfdp),  intent(in)           :: a
    integer,     intent(in), optional :: flags

    type(app_data_t), pointer :: x, y
    integer :: which

    call pepc_status('|----> encap_axpy()')
    call c_f_pointer(xptr, x)
    call c_f_pointer(yptr, y)

    DEBUG_ASSERT(x%params%nparts==y%params%nparts)

    which = AXPY_aVpV_aXpX
    if (present(flags)) which = flags
    select case (which)
    case (AXPY_aVpV_aXpX)
        y%x(:,:) = a * x%x + y%x
        y%v(:,:) = a * x%v + y%v
    case (AXPY_aVpV)
        y%v(:,:) = a * x%v + y%v
    case (AXPY_aXpX)
        y%x(:,:) = a * x%x + y%x
    case (AXPY_aVpX)
        y%x(:,:) = a * x%v + y%x
    case default
       DEBUG_ERROR(*, 'Invalid flags')
    end select

  end subroutine encap_axpy


  !> Compute norm of solution
  function encap_norm(ptr) result (norm)
    use pf_mod_mpi
    implicit none
    type(c_ptr), intent(in), value :: ptr
    real(pfdp) :: norm, norm_loc

    type(app_data_t), pointer :: q
    integer(kind_default) :: ierr

    call pepc_status('|----> encap_norm()')
    call c_f_pointer(ptr, q)

    norm_loc = max(maxval(abs(q%x)), maxval(abs(q%v)))

    call MPI_ALLREDUCE( norm_loc, norm, 1, MPI_DOUBLE_PRECISION, MPI_MAX, q%params%comm, ierr )

  end function encap_norm

end module pfm_encap
