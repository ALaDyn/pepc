module pepc_helper

   implicit none

   type pepc_nml_t
      integer :: np = 0
      integer :: pdump = 0
      integer :: fdump = 0
      integer :: cdump = 0
   end type pepc_nml_t

   logical :: para_file_available
   character(len = 255) :: para_file_name

contains

  !> computes a x b
  pure function cross_prod(a, b)
    implicit none

    real*8, dimension(3), intent(in) :: a, b
    real*8, dimension(3) :: cross_prod

    cross_prod(1) = a(2) * b(3) - a(3) * b(2)
    cross_prod(2) = a(3) * b(1) - a(1) * b(3)
    cross_prod(3) = a(1) * b(2) - b(2) * a(1)
  end function cross_prod

  !> computes a x b + c
  pure function cross_prod_plus(a, b, c)
    implicit none

    real*8, dimension(3), intent(in) :: a, b, c
    real*8, dimension(3) :: cross_prod_plus

    cross_prod_plus(1) = a(2) * b(3) - a(3) * b(2) + c(1)
    cross_prod_plus(2) = a(3) * b(1) - a(1) * b(3) + c(2)
    cross_prod_plus(3) = a(1) * b(2) - b(2) * a(1) + c(3)
  end function cross_prod_plus

  subroutine get_mpi_rank(comm, rank, size)
    use module_pepc_types
    implicit none
    include 'mpif.h'
    integer(kind_default), intent(in) :: comm
    integer(kind_pe), intent(out) :: rank, size
    integer(kind_default) :: mpi_err

    call MPI_COMM_RANK( comm, rank, mpi_err )
    call MPI_COMM_SIZE( comm, size, mpi_err )

  end subroutine

   subroutine pepc_setup(pepc_pars)
      use constants
      use encap
      use module_pepc
      use module_pepc_types, only: t_particle, kind_dim, kind_particle
      implicit none

      type(pepc_pars_t), intent(out) :: pepc_pars

      type(pepc_nml_t) :: pepc_nml

      call pepc_read_parameters_from_first_argument(para_file_available, para_file_name)
      call read_in_params(pepc_nml, para_file_available, para_file_name)

      ! Pass MPI stuff to parameters
      pepc_pars%np = pepc_nml%np

      pepc_pars%pdump = pepc_nml%pdump
      pepc_pars%fdump = pepc_nml%fdump
      pepc_pars%cdump = pepc_nml%cdump

   end subroutine pepc_setup


   subroutine read_in_params(pepc_namelist, file_available, file_name)
      use encap
      use module_mirror_boxes, only: mirror_box_layers
      use module_fmm_periodicity, only: do_extrinsic_correction
      implicit none

      type(pepc_nml_t), intent(out) :: pepc_namelist
      logical, intent(in) :: file_available
      character(len = 255), intent(in) :: file_name

      integer, parameter :: para_file_id = 10

      ! variables for pepc namelist
      integer :: np = 0
      integer :: pdump = 0
      integer :: fdump = 0
      integer :: cdump = 0

      namelist /pepc_nml/ np, pdump, fdump, cdump, mirror_box_layers, do_extrinsic_correction

      if (file_available) then
        open(para_file_id,file=trim(file_name),action='read')
        rewind(para_file_id)
        read(para_file_id, NML=pepc_nml)
        close(para_file_id)
      end if

      pepc_namelist%np = np
      pepc_namelist%pdump = pdump
      pepc_namelist%fdump = fdump
      pepc_namelist%cdump = cdump

   end subroutine read_in_params


   subroutine write_params(pepc_pars, file_name)
      use encap
      use module_mirror_boxes, only: mirror_box_layers
      use module_fmm_periodicity, only: do_extrinsic_correction
      implicit none

      type(pepc_pars_t), intent(in) :: pepc_pars
      character(len = 255), intent(in) :: file_name

      integer, parameter :: para_file_id = 10

      ! variables for pepc namelist
      integer(kind_particle) :: np = 0
      integer :: pdump = 0
      integer :: fdump = 0
      integer :: cdump = 0

      namelist /pepc_nml/ np, pdump, fdump, cdump, mirror_box_layers, do_extrinsic_correction

      np = pepc_pars%np
      pdump = pepc_pars%pdump
      fdump = pepc_pars%fdump
      cdump = pepc_pars%cdump

      open(para_file_id, file = trim(file_name), status = 'old', position = &
        'append', action = 'write')
      write(para_file_id, NML=pepc_nml)
      close(para_file_id)

   end subroutine write_params


  function get_time()
    use mpi
    implicit none

    real(kind = 8) :: get_time

    get_time = MPI_WTIME()
  end function


  ! FIXME: this routine will not in time-parallel mode
  subroutine write_particles(pepc_pars, time_pars, step, p)
    use module_pepc_types
    use module_vtk
    use module_vtk_helpers
    use encap
    implicit none

    type(pepc_pars_t), intent(in) :: pepc_pars
    type(time_pars_t), intent(in) :: time_pars
    type(t_particle), allocatable, intent(in) :: p(:)
    integer, intent(in) :: step

    integer :: vtk_step
    real*8 :: time
    real*8 :: ta, tb

    ta = get_time()

    time = time_pars%dt * step

    if (step .eq. 0) then
      vtk_step = VTK_STEP_FIRST
    else if (step .eq. time_pars%nsteps) then
      vtk_step = VTK_STEP_LAST
    else
      vtk_step = VTK_STEP_NORMAL
    endif

    call vtk_write_particles("particles", pepc_pars%pepc_comm%rank_space, step, time, vtk_step, p)

    tb = get_time()

    if(pepc_pars%pepc_comm%rank_space == 0) write(*,'(a,es12.4)') " == [write particles] time in vtk output [s]      : ", tb - ta

  end subroutine write_particles


  subroutine write_domain(time_pars, step, p)
    use module_pepc_types
    use module_vtk
    use module_vtk_helpers
    use module_pepc, only : global_tree

    use encap
    implicit none

    type(time_pars_t), intent(in) :: time_pars
    integer, intent(in) :: step
    type(t_particle), allocatable, intent(in) :: p(:)

    integer :: vtk_step

    ! output of tree diagnostics
    if (step .eq. 0) then
      vtk_step = VTK_STEP_FIRST
    else if (step .eq. time_pars%nsteps) then
      vtk_step = VTK_STEP_LAST
    else
      vtk_step = VTK_STEP_NORMAL
    endif

    call vtk_write_branches(step,  time_pars%dt * step, vtk_step, global_tree)
    call vtk_write_leaves(step, time_pars%dt * step, vtk_step, global_tree)
    call vtk_write_spacecurve(step, time_pars%dt * step, vtk_step, p)

  end subroutine write_domain

end module pepc_helper
