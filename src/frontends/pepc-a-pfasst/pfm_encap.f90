module pfm_encap
  use iso_c_binding
  use pfasst
  use module_pepc_types
  use module_debug
  implicit none

  !> data type for level-dependent application parameters
  type :: app_params_t
   integer :: npart
   real*8 :: theta
  end type

  ! Data encapsulation: 3d array for pmg + solver parameters which will be filled in encap_create using ctx
  type :: app_data_t
   ! TODO: fill with content here (particles and parameters)
   type(t_particle), dimension(:), pointer :: particles
   type(app_params_t) :: params
  end type app_data_t

contains

  !> Fill pf_encap_t with pointers to encapsulation functions
  subroutine pfm_encap_create(encap)
    type(pf_encap_t), intent(out) :: encap

    encap%create  => encap_create
    encap%destroy => encap_destroy
    encap%setval  => encap_setval
    encap%copy    => encap_copy
    encap%pack    => encap_pack
    encap%unpack  => encap_unpack
    encap%axpy    => encap_axpy
    encap%norm    => encap_norm
  end subroutine pfm_encap_create


  !> Allocate/create solution (spatial data set) for the given level.
  !> This is called for each SDC node.
  subroutine encap_create(sol, level, kind, nvars, shape, levelctx, encapctx) ! FIXME: the type of encapctx could be pf_encap_t, shouldn't it?
    type(c_ptr),       intent(inout)     :: sol
    integer,           intent(in)        :: level, nvars, shape(:)
    integer,           intent(in)        :: kind
    type(c_ptr),       intent(in), value :: levelctx, encapctx

    type(app_data_t), pointer :: q
    type(app_params_t), pointer :: p
    
    call pepc_status('|----> encap_create()')

    call c_f_pointer(levelctx, p)

    DEBUG_ASSERT(nvars==2*3*p%npart)

    allocate(q)
    q%params = p
    
    allocate(q%particles(q%params%npart))
    sol = c_loc(q)

  end subroutine encap_create


  !> Deallocate/destroy solution.
  subroutine encap_destroy(ptr)
    type(c_ptr), intent(in), value :: ptr

    type(app_data_t), pointer :: q

    call pepc_status('|----> encap_destroy()')
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

    call pepc_status('|----> encap_setval()')
    call c_f_pointer(ptr,q)

    do i=1,q%params%npart
      q%particles(i)%x(:)      = val
      q%particles(i)%data%v(:) = val
    end do

  end subroutine encap_setval


  !> Copy solution value on one level.
  subroutine encap_copy(dstptr, srcptr, flags)
    type(c_ptr), intent(in), value    :: dstptr, srcptr
    integer,     intent(in), optional :: flags

    type(app_data_t), pointer :: dst, src

    call pepc_status('|----> encap_copy()')
    call c_f_pointer(dstptr,dst) ! FIXME: do we actually need this? direct c_ptr copy? - No: shared pointers would be deallocated twice then, This would only work for moving data
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
    call c_f_pointer(ptr, q)

    j = 1
    do i=1,q%params%npart
      z(j:j+2) = q%particles(i)%x(1:3)
      j = j+3
      z(j:j+2) = q%particles(i)%data%v(1:3)
      j = j+3
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
    call c_f_pointer(ptr, q)

    j = 1
    do i=1,q%params%npart
      q%particles(i)%x(1:3)      = z(j:j+2)
      j = j+3
      q%particles(i)%data%v(1:3) = z(j:j+2)
      j = j+3
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

    call pepc_status('|----> encap_axpy()')
    call c_f_pointer(xptr, x)
    call c_f_pointer(yptr, y)
    
    DEBUG_ASSERT(x%params%npart==y%params%npart)

    ! FIXME: how does this look like for velocities and positions separately?
    do i=1,x%params%npart
      y%particles(i)%x(:)      = a * x%particles(i)%x(:)      + y%particles(i)%x(:)
      y%particles(i)%data%v(:) = a * x%particles(i)%data%v(:) + y%particles(i)%data%v(:)
    end do

  end subroutine encap_axpy


  !> Compute norm of solution
  function encap_norm(ptr) result (norm)
    use pf_mod_mpi
    type(c_ptr), intent(in), value :: ptr
    real(pfdp) :: norm

    type(app_data_t), pointer :: q

    call pepc_status('|----> encap_norm()')
    call c_f_pointer(ptr, q)

    ! TODO
    !norm_loc = maxval(abs(q%array(F%m_start:F%m_end,F%n_start:F%n_end,F%o_start:F%o_end)))
    !call MPI_ALLREDUCE( norm_loc, norm, 1, MPI_DOUBLE_PRECISION, MPI_MAX, F%pmg_comm%mpi_comm, ierr )
    norm = 1.

  end function encap_norm

end module pfm_encap
