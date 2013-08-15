module pfm_encap
  use iso_c_binding
  use module_pepc_types
  use module_debug
  use pf_mod_dtype
  implicit none

  !> data type for level-dependent application parameters
  !> if any parameters are added here, they also have to be included in #pfm_setup_solver_level_params
  type :: app_params_t
    integer :: n_el, n_ion
    integer :: dim
    real*8 :: theta
    integer(kind_default) :: comm
  end type

  !> Data encapsulation: data and parameters which will be filled in encap_create using ctx
  type :: app_data_t
    type(t_particle), dimension(:), allocatable :: particles
    type(app_params_t) :: params
  end type app_data_t

contains

  !> internal helper for getting a pointer value (i.e. the address of the pointee)
  function ptr_val(ptr)
    use iso_c_binding, only : c_ptr
    implicit none
    type(c_ptr), intent(in) :: ptr
    integer :: ptr_val
    
    ptr_val = transfer (ptr, ptr_val)
  end function
  
  !> internal helper for printing a pointer value (i.e. the address of the pointee)
  subroutine ptr_print(name, ptr)
    use iso_c_binding, only : c_ptr
    implicit none
    character(*), intent(in) :: name
    type(c_ptr), intent(in) :: ptr
    
    if (dbg(DBG_STATUS)) write(*,'("|------------------------------- &", a,"=0x",Z8.8)') trim(name), ptr_val(ptr)
  end subroutine


  !> Fill pf_encap_t with pointers to encapsulation functions
  subroutine pfm_encap_init(encap)
    type(pf_encap_t), intent(out) :: encap

    encap%create  => encap_create
    encap%destroy => encap_destroy
    encap%setval  => encap_setval
    encap%copy    => encap_copy
    encap%pack    => encap_pack
    encap%unpack  => encap_unpack
    encap%axpy    => encap_axpy
    encap%norm    => encap_norm
  end subroutine pfm_encap_init


  !> Allocate/create solution (spatial data set) for the given level.
  !> This is called for each SDC node.
  subroutine encap_create(sol, level, kind, nvars, shape, levelctx, encapctx)
    type(c_ptr),       intent(inout)     :: sol
    integer,           intent(in)        :: level, nvars, shape(:)
    integer,           intent(in)        :: kind
    type(c_ptr),       intent(in), value :: levelctx, encapctx

    type(app_data_t), pointer :: q
    type(app_params_t), pointer :: p
    
    call pepc_status('|----> encap_create()')

    call c_f_pointer(levelctx, p)

    DEBUG_ASSERT(nvars==2*p%dim*(p%n_el+p%n_ion)) ! dim coordinates and momenta per particle

    allocate(q)
    q%params = p
    
    allocate(q%particles(q%params%n_el+q%params%n_ion))
    sol = c_loc(q)
    
    call ptr_print('sol', sol)

  end subroutine encap_create


  !> Deallocate/destroy solution.
  subroutine encap_destroy(ptr)
    type(c_ptr), intent(in), value :: ptr

    type(app_data_t), pointer :: q

    call pepc_status('|----> encap_destroy()')
    call ptr_print('ptr', ptr)
    call c_f_pointer(ptr, q)

    deallocate(q%particles)
    deallocate(q)

  end subroutine encap_destroy


  !> Set solution value.
  subroutine encap_setval(ptr, val, flags)
    type(c_ptr), intent(in), value    :: ptr
    real(pfdp),  intent(in)           :: val
    integer,     intent(in), optional :: flags

    type(app_data_t), pointer :: q
    integer(kind_particle) :: i
    integer :: which

    call pepc_status('|----> encap_setval()')
    call ptr_print('ptr', ptr)
    call c_f_pointer(ptr, q)

    which = 0
    if (present(flags)) which = flags

    select case (which)
    case (0)
      ! set value for x and v
      do i=1,q%params%n_el+q%params%n_ion
        q%particles(i)%x(:)      = val
        q%particles(i)%data%v(:) = val
      end do
    case (1)
      ! set value for v
      do i=1,q%params%n_el+q%params%n_ion
        q%particles(i)%data%v(:) = val
      end do
    case (2)
      ! set value for x
      do i=1,q%params%n_el+q%params%n_ion
        q%particles(i)%x(:)      = val
      end do
    case default
       DEBUG_ERROR(*, 'Invalid flags')
    end select

  end subroutine encap_setval


  !> Copy solution value on one level.
  subroutine encap_copy(dstptr, srcptr, flags)
    type(c_ptr), intent(in), value    :: dstptr, srcptr
    integer,     intent(in), optional :: flags

    type(app_data_t), pointer :: dst, src

    call pepc_status('|----> encap_copy()')
    call ptr_print('src', dstptr)
    call ptr_print('dst', srcptr)
    call c_f_pointer(dstptr,dst)
    call c_f_pointer(srcptr,src)

    dst = src

  end subroutine encap_copy


  !> Pack solution q into a flat array.
  subroutine encap_pack(z, ptr)
    type(c_ptr), intent(in), value  :: ptr
    real(pfdp),  intent(out)        :: z(:)  !FIXME: an alternative with a frontend-defined MPI type would be very nice - this will avoid packing completely :-)

    type(app_data_t), pointer :: q
    integer(kind_particle) :: i, j

    call pepc_status('|----> encap_pack()')
    call ptr_print('ptr', ptr)
    call c_f_pointer(ptr, q)

    j = 1
    do i=1,q%params%n_el+q%params%n_ion
      z(j:j+q%params%dim-1) = q%particles(i)%x(1:q%params%dim)
      j = j+q%params%dim
      z(j:j+q%params%dim-1) = q%particles(i)%data%v(1:q%params%dim)
      j = j+q%params%dim
    end do

    DEBUG_ASSERT(j==size(z)+1)

  end subroutine encap_pack


  !> Unpack solution from a flat array.
  subroutine encap_unpack(ptr, z)
    type(c_ptr), intent(in), value :: ptr
    real(pfdp),  intent(in)        :: z(:)  !FIXME: an alternative with a frontend-defined MPI type would be very nice - this will avoid packing completely :-)

    type(app_data_t), pointer :: q
    integer(kind_particle) :: i, j

    call pepc_status('|----> encap_unpack()')
    call ptr_print('ptr', ptr)
    call c_f_pointer(ptr, q)

    j = 1
    do i=1,q%params%n_el+q%params%n_ion
      q%particles(i)%x(1:q%params%dim) = z(j:j+q%params%dim-1)
      j = j+q%params%dim
      q%particles(i)%data%v(1:q%params%dim) = z(j:j+q%params%dim-1)
      j = j+q%params%dim
    end do
    
    DEBUG_ASSERT(j==size(z)+1)

  end subroutine encap_unpack


  !> Compute y = a x + y where a is a scalar and x and y are solutions on the same level.
  subroutine encap_axpy(yptr, a, xptr, flags)
    type(c_ptr), intent(in), value    :: xptr, yptr
    real(pfdp),  intent(in)           :: a
    integer,     intent(in), optional :: flags

    type(app_data_t), pointer :: x, y
    integer(kind_particle) :: i
    integer :: which

    call pepc_status('|----> encap_axpy()')
    call ptr_print('xptr', xptr)
    call ptr_print('yptr', yptr)
    call c_f_pointer(xptr, x)
    call c_f_pointer(yptr, y)
    
    DEBUG_ASSERT(x%params%n_el ==y%params%n_el)
    DEBUG_ASSERT(x%params%n_ion==y%params%n_ion)

    which = 0
    if (present(flags)) which = flags
    select case (which)
    case (0)
      do i=1,x%params%n_el+x%params%n_ion
        y%particles(i)%x(:)      = a * x%particles(i)%x(:)      + y%particles(i)%x(:)
        y%particles(i)%data%v(:) = a * x%particles(i)%data%v(:) + y%particles(i)%data%v(:)
      end do
    case (1)
      do i=1,x%params%n_el+x%params%n_ion
        y%particles(i)%data%v(:) = a * x%particles(i)%data%v(:) + y%particles(i)%data%v(:)
      end do
    case (2)
      do i=1,x%params%n_el+x%params%n_ion
        y%particles(i)%x(:)      = a * x%particles(i)%x(:)      + y%particles(i)%x(:)
      end do
    case (12)
      do i=1,x%params%n_el+x%params%n_ion
        y%particles(i)%x(:)      = a * x%particles(i)%data%v(:) + y%particles(i)%x(:)
      end do
    case default
       DEBUG_ERROR(*, 'Invalid flags')
    end select
    
  end subroutine encap_axpy


  !> Compute norm of solution
  function encap_norm(ptr) result (norm)
    use pf_mod_mpi
    type(c_ptr), intent(in), value :: ptr
    real(pfdp) :: norm, norm_loc

    type(app_data_t), pointer :: q
    integer(kind_particle) :: i
    integer(kind_Default) :: ierr

    call pepc_status('|----> encap_norm()')
    call ptr_print('ptr', ptr)
    call c_f_pointer(ptr, q)

    ! TODO
    norm_loc = 0.
    
    do i=1,q%params%n_el+q%params%n_ion
      norm_loc = maxval([norm_loc, maxval(q%particles(i)%x(:)), maxval(q%particles(i)%data%v(:))])
    end do

    call MPI_ALLREDUCE( norm_loc, norm, 1, MPI_DOUBLE_PRECISION, MPI_MAX, q%params%comm, ierr )

  end function encap_norm

end module pfm_encap
