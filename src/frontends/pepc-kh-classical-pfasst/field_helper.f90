module field_helper
   implicit none

  character(*), parameter, private :: field_dir = "fields"

  type field_grid_nml_t
    integer, dimension(2) :: n
    real(kind = 8), dimension(2) :: offset, extent
  end type field_grid_nml_t

contains

  subroutine setup_field_grid(field_grid, pepc_comm)
    use pepc_helper, only: para_file_available, para_file_name
    use module_utils, only: create_directory
    use encap
    use module_debug
    implicit none

    type(field_grid_t), intent(out) :: field_grid
    type(pepc_comm_t), intent(in) :: pepc_comm

    type(field_grid_nml_t) :: field_grid_nml
    integer(kind_default) :: mpi_size, mpi_rank
    integer(kind_particle) :: ipl, ipg, ix, iy

    call pepc_status("SETUP FIELD GRID")

    call create_directory(field_dir)

    call read_in_field_grid_params(field_grid_nml, para_file_available, para_file_name)

    mpi_rank = int(pepc_comm%rank_space,  kind = kind(mpi_rank))
    mpi_size = int(pepc_comm%nrank_space, kind = kind(mpi_size))

    field_grid%n = field_grid_nml%n
    field_grid%ntot = product(field_grid%n)
    field_grid%nl = field_grid%ntot / mpi_size
    if (mpi_rank < mod(field_grid%ntot, int(mpi_size, kind=kind_particle))) &
      field_grid%nl = field_grid%nl + 1
    field_grid%offset = field_grid_nml%offset
    field_grid%extent = field_grid_nml%extent
    field_grid%dx = field_grid%extent / field_grid%n

    if (mpi_rank == 0) then
      print *, "== [setup_field_grid]"
      print *, "   arranging ", field_grid%ntot, " field probes on a grid."
      print *, "   nx x ny = ", field_grid%n(1), " x ", field_grid%n(2)
      print *, "   x0 x y0 = ", field_grid%offset(1), " x ", field_grid%offset(2)
      print *, "   Lx x Ly = ", field_grid%extent(1), " x ", field_grid%extent(2)
      print *, "   dx x dy = ", field_grid%dx(1), " x ", field_grid%dx(2)
      print *, ""
    end if

    allocate(field_grid%p(field_grid%nl), &
      field_grid%ne(field_grid%n(1), field_grid%n(2)), &
      field_grid%ni(field_grid%n(1), field_grid%n(2)), &
      field_grid%ne_from_left(field_grid%n(1), field_grid%n(2)), &
      field_grid%ni_from_left(field_grid%n(1), field_grid%n(2)), &
      field_grid%vex(field_grid%n(1), field_grid%n(2)), &
      field_grid%vey(field_grid%n(1), field_grid%n(2)), &
      field_grid%vix(field_grid%n(1), field_grid%n(2)), &
      field_grid%viy(field_grid%n(1), field_grid%n(2)))

    do ipl = 1, field_grid%nl
      ipg = ipl + min(int(mpi_rank, kind=kind_particle), mod(field_grid%ntot, int(mpi_size, kind=kind_particle))) * &
        (field_grid%ntot / mpi_size + 1) + &
        max(0_kind_particle, mpi_rank - mod(field_grid%ntot, int(mpi_size, kind=kind_particle))) * &
        (field_grid%ntot / mpi_size)

      ix = mod(ipg - 1, field_grid%n(1)) + 1
      iy = (ipg - 1) / field_grid%n(1) + 1

      field_grid%p(ipl)%x(1) = (ix - 0.5D0) * field_grid%dx(1) + field_grid%offset(1)
      field_grid%p(ipl)%x(2) = (iy - 0.5D0) * field_grid%dx(2) + field_grid%offset(2)
      field_grid%p(ipl)%x(3) = 0.0D0

      field_grid%p(ipl)%label = ipg
    end do

  end subroutine setup_field_grid


   subroutine read_in_field_grid_params(field_grid_namelist, file_available, file_name)
      use mpi
      implicit none

      type(field_grid_nml_t), intent(out) :: field_grid_namelist
      logical, intent(in) :: file_available
      character(len = 255), intent(in) :: file_name

      integer, dimension(2) :: n = [ 0, 0 ]
      real(kind=8), dimension(2) :: offset = [ 0.0D0, 0.0D0 ]
      real(kind=8), dimension(2) :: extent = [ 0.0D0, 0.0D0 ]

      namelist /field_grid_nml/ n, offset, extent

      integer, parameter :: para_file_id = 10

      if (file_available) then
         open(para_file_id,file=trim(file_name),action='read')
         rewind(para_file_id)
         read(para_file_id, NML=field_grid_nml)
         close(para_file_id)

      end if

      field_grid_namelist%n = n
      field_grid_namelist%offset = offset
      field_grid_namelist%extent = extent

   end subroutine read_in_field_grid_params


  subroutine write_field_grid_params(field_grid, file_name)
    use encap
    implicit none

    type(field_grid_t), intent(in) :: field_grid
    character(len = 255), intent(in) :: file_name

    integer, parameter :: para_file_id = 10

    integer(kind_particle), dimension(2) :: n
    real(kind=8), dimension(2) :: offset
    real(kind=8), dimension(2) :: extent

    namelist /field_grid_nml/ n, offset, extent

    n = field_grid%n
    offset = field_grid%offset
    extent = field_grid%extent

    open(para_file_id, file = trim(file_name), status = 'old', position = &
      'append', action = 'write')
    write(para_file_id, nml = field_grid_nml)
    close(para_file_id)

  end subroutine write_field_grid_params


  subroutine compute_field(pepc_pars, field_grid, p)
    use mpi
    use module_pepc
    use module_tree, only : tree_allocated
    use encap
    use physics_helper
    implicit none

    type(pepc_pars_t), intent(in) :: pepc_pars
    type(field_grid_t), intent(inout) :: field_grid
    type(t_particle), dimension(:), allocatable, intent(inout) :: p

    integer(kind_particle) :: ipl
    integer(kind_default) :: mpi_err
    integer, dimension(2) :: ic
    real(kind=8) :: da, rda

    call pepc_particleresults_clear(field_grid%p)

    if (.not. tree_allocated(global_tree)) then
      call pepc_grow_and_traverse_for_others(p, field_grid%p, 0)
    else
      call pepc_traverse_tree(field_grid%p)
    endif

    field_grid%p(:)%results%e(1) = field_grid%p(:)%results%e(1) * force_const
    field_grid%p(:)%results%e(2) = field_grid%p(:)%results%e(2) * force_const
    field_grid%p(:)%results%pot  = field_grid%p(:)%results%pot  * force_const

    field_grid%ne = 0.0D0
    field_grid%ni = 0.0D0
    field_grid%vex = 0.0D0
    field_grid%vey = 0.0D0
    field_grid%vix = 0.0D0
    field_grid%viy = 0.0D0
    da = product(field_grid%dx)
    rda = 1.0D0 / da

    ! make a spatial histogram of local particle numbers and velocities per species
    do ipl = 1, size(p)
      ! do not count particles outside the grid
      if (any(p(ipl)%x(1:2) < field_grid%offset) .or. &
        any((p(ipl)%x(1:2) - field_grid%offset) >= field_grid%extent)) cycle

      ic = ceiling((p(ipl)%x(1:2) - field_grid%offset) / field_grid%dx)

      if (p(ipl)%data%q < 0.0D0) then ! electron
        field_grid%ne( ic(1), ic(2)) = field_grid%ne( ic(1), ic(2)) + 1
        field_grid%vex(ic(1), ic(2)) = field_grid%vex(ic(1), ic(2)) +  p(ipl)%data%v(1)
        field_grid%vey(ic(1), ic(2)) = field_grid%vey(ic(1), ic(2)) +  p(ipl)%data%v(2)
        if (p(ipl)%label == LABEL_ELECTRON_LEFT) then
          field_grid%ne_from_left( ic(1), ic(2)) = field_grid%ne_from_left( ic(1), ic(2)) + 1
        endif
      else
        field_grid%ni( ic(1), ic(2)) = field_grid%ni( ic(1), ic(2)) + 1
        field_grid%vix(ic(1), ic(2)) = field_grid%vix(ic(1), ic(2)) +  p(ipl)%data%v(1)
        field_grid%viy(ic(1), ic(2)) = field_grid%viy(ic(1), ic(2)) +  p(ipl)%data%v(2)
        if (p(ipl)%label == LABEL_ION_LEFT) then
          field_grid%ni_from_left( ic(1), ic(2)) = field_grid%ni_from_left( ic(1), ic(2)) + 1
        endif
      end if
    end do

    ! accumulate the histograms
    call mpi_allreduce(MPI_IN_PLACE, field_grid%ne, int(field_grid%ntot, kind = kind_default), MPI_REAL8, &
      MPI_SUM, pepc_pars%pepc_comm%comm_space, mpi_err)
    call mpi_allreduce(MPI_IN_PLACE, field_grid%ni, int(field_grid%ntot, kind = kind_default), MPI_REAL8, &
      MPI_SUM, pepc_pars%pepc_comm%comm_space, mpi_err)
    call mpi_allreduce(MPI_IN_PLACE, field_grid%ne_from_left, int(field_grid%ntot, kind = kind_default), MPI_REAL8, &
      MPI_SUM, pepc_pars%pepc_comm%comm_space, mpi_err)
    call mpi_allreduce(MPI_IN_PLACE, field_grid%ni_from_left, int(field_grid%ntot, kind = kind_default), MPI_REAL8, &
      MPI_SUM, pepc_pars%pepc_comm%comm_space, mpi_err)

    call mpi_allreduce(MPI_IN_PLACE, field_grid%vex, int(field_grid%ntot, kind = kind_default), MPI_REAL8, &
      MPI_SUM, pepc_pars%pepc_comm%comm_space, mpi_err)
    call mpi_allreduce(MPI_IN_PLACE, field_grid%vey, int(field_grid%ntot, kind = kind_default), MPI_REAL8, &
      MPI_SUM, pepc_pars%pepc_comm%comm_space, mpi_err)
    call mpi_allreduce(MPI_IN_PLACE, field_grid%vix, int(field_grid%ntot, kind = kind_default), MPI_REAL8, &
      MPI_SUM, pepc_pars%pepc_comm%comm_space, mpi_err)
    call mpi_allreduce(MPI_IN_PLACE, field_grid%viy, int(field_grid%ntot, kind = kind_default), MPI_REAL8, &
      MPI_SUM, pepc_pars%pepc_comm%comm_space, mpi_err)

    ! normalize to fluid quantities
    ! total velocity -> mean velocity
    where (field_grid%ne > 0.0D0)
      field_grid%vex = field_grid%vex / field_grid%ne
      field_grid%vey = field_grid%vey / field_grid%ne
    elsewhere
      field_grid%vex = 0.0D0
      field_grid%vey = 0.0D0
    end where

    where (field_grid%ni > 0.0D0)
      field_grid%vix = field_grid%vix / field_grid%ni
      field_grid%viy = field_grid%viy / field_grid%ni
    elsewhere
      field_grid%vix = 0.0D0
      field_grid%viy = 0.0D0
    end where
    ! particle number -> particle density
    field_grid%ne(:,:) = field_grid%ne * rda
    field_grid%ni(:,:) = field_grid%ni * rda
    field_grid%ne_from_left(:,:) = field_grid%ne_from_left * rda
    field_grid%ni_from_left(:,:) = field_grid%ni_from_left * rda

  end subroutine compute_field


  subroutine write_field_on_grid(pepc_comm, time_pars, step, physics_pars, field_grid)
    use mpi
    use encap
    implicit none

    type(pepc_comm_t), intent(in) :: pepc_comm
    type(time_pars_t), intent(in) :: time_pars
    integer, intent(in) :: step
    type(physics_pars_t), intent(in) :: physics_pars
    type(field_grid_t), intent(in) :: field_grid

    integer(kind_particle) :: my_count, fflatshape(1)
    integer(kind = MPI_OFFSET_KIND) :: my_offset
    real(kind = 8), dimension(:), allocatable :: fflat

    allocate(fflat(field_grid%ntot))
    fflatshape(1) = field_grid%ntot

    my_count = field_grid%nl
    my_offset = min(int(pepc_comm%rank_space, kind=kind_particle), &
      mod(field_grid%ntot, int(pepc_comm%nrank_space, kind=kind_particle))) &
      * (field_grid%ntot / pepc_comm%nrank_space + 1) + &
      max(0_kind_particle, pepc_comm%rank_space - mod(field_grid%ntot, int(pepc_comm%nrank_space, kind=kind_particle))) * &
      (field_grid%ntot / pepc_comm%nrank_space)

    call write_quantity_on_grid("potential", field_grid%p(:)%results%pot)
    call write_quantity_on_grid("ex", field_grid%p(:)%results%e(1))
    call write_quantity_on_grid("ey", field_grid%p(:)%results%e(2))
    fflat(:) = reshape(field_grid%ne, fflatshape)
    call write_quantity_on_grid("ne", fflat(my_offset + 1:))
    fflat(:) = reshape(field_grid%ni, fflatshape)
    call write_quantity_on_grid("ni", fflat(my_offset + 1:))
    fflat(:) = reshape(field_grid%ne_from_left, fflatshape)
    call write_quantity_on_grid("nefromleft", fflat(my_offset + 1:))
    fflat(:) = reshape(field_grid%ni_from_left, fflatshape)
    call write_quantity_on_grid("nifromleft", fflat(my_offset + 1:))
    fflat(:) = reshape(field_grid%vex, fflatshape)
    call write_quantity_on_grid("vex", fflat(my_offset + 1:))
    fflat(:) = reshape(field_grid%vey, fflatshape)
    call write_quantity_on_grid("vey", fflat(my_offset + 1:))
    fflat(:) = reshape(field_grid%vix, fflatshape)
    call write_quantity_on_grid("vix", fflat(my_offset + 1:))
    fflat(:) = reshape(field_grid%viy, fflatshape)
    call write_quantity_on_grid("viy", fflat(my_offset + 1:))

    deallocate(fflat)

    contains

    subroutine write_quantity_on_grid(yname, y)
      use mpi
      use module_debug
      implicit none

      character(*), intent(in) :: yname
      real(kind = 8), dimension(:), intent(in) :: y

      integer(kind = MPI_OFFSET_KIND), parameter :: field_header_size = 512
      integer(kind = 4), parameter :: field_header_magic = 1

      character(1024) :: filename
      integer :: fh, mpi_err, nx_, ny_
      integer(kind = MPI_OFFSET_KIND) :: mpi_disp
      integer, dimension(MPI_STATUS_SIZE) :: mpi_stat

      write(filename, '(a,"/",a,"_",i6.6,".bin")') field_dir, yname, step

      call mpi_file_open(pepc_comm%comm_space, filename, &
        ior(MPI_MODE_RDWR, MPI_MODE_CREATE), MPI_INFO_NULL, fh, mpi_err)

      if (mpi_err .ne. MPI_SUCCESS) then
        DEBUG_ERROR(*, 'write_field_on_grid(): I/O error for ', filename)
      end if

      nx_ = int(field_grid%n(1))
      ny_ = int(field_grid%n(2))

      call mpi_file_set_view(fh, 0_MPI_OFFSET_KIND, MPI_BYTE, MPI_BYTE, &
        'native', MPI_INFO_NULL, mpi_err)

      ! write file header on rank 0
      ! int32   magic
      ! int32   nx
      ! int32   ny
      ! float64 offset_x
      ! float64 offset_y
      ! float64 extent_x
      ! float64 extent_y
      ! int32   step
      ! float64 t
      ! float64 B0
      ! float64 vte
      ! float64 vti
      ! float64 qe
      ! float64 qi
      ! float64 me
      ! float64 mi
      ! float64 min
      ! float64 max
      if (pepc_comm%rank_space == 0) then
        call mpi_file_write(fh, field_header_magic, 1, MPI_INTEGER4, mpi_stat, mpi_err)
        call mpi_file_write(fh, nx_, 1, MPI_INTEGER4, mpi_stat, mpi_err)
        call mpi_file_write(fh, ny_, 1, MPI_INTEGER4, mpi_stat, mpi_err)
        call mpi_file_write(fh, field_grid%offset(1), 1, MPI_REAL8, mpi_stat, mpi_err)
        call mpi_file_write(fh, field_grid%offset(2), 1, MPI_REAL8, mpi_stat, mpi_err)
        call mpi_file_write(fh, field_grid%extent(1), 1, MPI_REAL8, mpi_stat, mpi_err)
        call mpi_file_write(fh, field_grid%extent(2), 1, MPI_REAL8, mpi_stat, mpi_err)
        call mpi_file_write(fh, step, 1, MPI_INTEGER4, mpi_stat, mpi_err)
        call mpi_file_write(fh, step * time_pars%dt, 1, MPI_REAL8, mpi_stat, mpi_err)
        call mpi_file_write(fh, physics_pars%B0, 1, MPI_REAL8, mpi_stat, mpi_err)
        call mpi_file_write(fh, physics_pars%vte, 1, MPI_REAL8, mpi_stat, mpi_err)
        call mpi_file_write(fh, physics_pars%vti, 1, MPI_REAL8, mpi_stat, mpi_err)
        call mpi_file_write(fh, physics_pars%qe, 1, MPI_REAL8, mpi_stat, mpi_err)
        call mpi_file_write(fh, physics_pars%qi, 1, MPI_REAL8, mpi_stat, mpi_err)
        call mpi_file_write(fh, physics_pars%me, 1, MPI_REAL8, mpi_stat, mpi_err)
        call mpi_file_write(fh, physics_pars%mi, 1, MPI_REAL8, mpi_stat, mpi_err)
        call mpi_file_write(fh, minval(y), 1, MPI_REAL8, mpi_stat, mpi_err)
        call mpi_file_write(fh, maxval(y), 1, MPI_REAL8, mpi_stat, mpi_err)

        call mpi_file_get_position(fh, mpi_disp, mpi_err)
        if (mpi_disp > field_header_size) then
          DEBUG_ERROR(*, "write_field_on_grid(): field_header_size is too small: ", field_header_size, " < ", mpi_disp)
        end if
      end if

      call mpi_file_set_view(fh, field_header_size, MPI_REAL8, MPI_REAL8, &
        'native', MPI_INFO_NULL, mpi_err)

      call mpi_file_write_at_all(fh, my_offset, y(:), int(field_grid%nl, kind = kind_default), MPI_REAL8, mpi_stat, mpi_err)

      call mpi_file_sync(fh, mpi_err)
      call mpi_file_close(fh, mpi_err)

    end subroutine

  end subroutine

end module field_helper
