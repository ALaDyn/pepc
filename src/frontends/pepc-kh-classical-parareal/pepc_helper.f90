module pepc_helper
  use module_pepc_kinds
  implicit none

contains

  subroutine pepc_setup(file_name, pepc_pars)
    use constants
    use encap
    use module_pepc
    use module_pepc_types, only: t_particle
    implicit none

    character(*), intent(in) :: file_name
    type(pepc_pars_t), intent(out) :: pepc_pars

    call pepc_initialize(FRONTEND_NAME, pepc_pars%pepc_comm%mpi_rank, &
      pepc_pars%pepc_comm%mpi_size, .true., comm = pepc_pars%pepc_comm%mpi_comm)

    call pepc_read_parameters_from_file_name(file_name)
    call read_in_params(file_name, pepc_pars)
  end subroutine pepc_setup


  subroutine read_in_params(file_name, pepc_pars)
    use encap
    use module_mirror_boxes, only: mirror_box_layers
    use module_fmm_periodicity, only: do_extrinsic_correction
    implicit none

    character(*), intent(in) :: file_name
    type(pepc_pars_t), intent(out) :: pepc_pars

    integer, parameter :: param_file_id = 10

    ! variables for pepc namelist
    integer :: np = 0
    integer :: pdump = 0
    integer :: fdump = 0
    integer :: cdump = 0

    namelist /pepc_nml/ np, pdump, fdump, cdump, mirror_box_layers, do_extrinsic_correction

    open(param_file_id,file=trim(file_name),action='read')
    rewind(param_file_id)
    read(param_file_id, NML=pepc_nml)
    close(param_file_id)

    pepc_pars%np = np
    pepc_pars%pdump = pdump
    pepc_pars%fdump = fdump
    pepc_pars%cdump = cdump
  end subroutine read_in_params


  subroutine write_params(file_name, pepc_pars)
    use encap
    use module_mirror_boxes, only: mirror_box_layers
    use module_fmm_periodicity, only: do_extrinsic_correction
    implicit none

    character(*), intent(in) :: file_name
    type(pepc_pars_t), intent(in) :: pepc_pars

    integer, parameter :: param_file_id = 10

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

    open(param_file_id, file = trim(file_name), status = 'old', position = &
      'append', action = 'write')
    write(param_file_id, NML=pepc_nml)
    close(param_file_id)
  end subroutine write_params


  function get_time()
    use mpi
    implicit none

    real(kind = 8) :: get_time

    get_time = MPI_WTIME()
  end function


  pure function vtk_step_of_step(step, time_pars) result(vtk_step)
    use module_vtk
    use encap
    implicit none

    integer, intent(in) :: step
    type(time_pars_t), intent(in) :: time_pars

    integer :: vtk_step

    if (step .eq. 0) then
      vtk_step = VTK_STEP_FIRST
    else if (step .eq. time_pars%nsteps) then
      vtk_step = VTK_STEP_LAST
    else
      vtk_step = VTK_STEP_NORMAL
    endif
  end function


  subroutine write_particles(pepc_pars, time_pars, step, p)
    use module_pepc_types
    use module_vtk
    use module_vtk_helpers
    use encap
    use time_helper
    implicit none

    type(pepc_pars_t), intent(in) :: pepc_pars
    type(time_pars_t), intent(in) :: time_pars
    type(t_particle), allocatable, intent(in) :: p(:)
    integer, intent(in) :: step

    integer :: vtk_step
    real*8 :: time
    real*8 :: ta, tb

    ta = get_time()

    time = time_of_step(step, time_pars)
    vtk_step = vtk_step_of_step(step, time_pars)

    call vtk_write_particles("particles", pepc_pars%pepc_comm%mpi_comm, step, time, vtk_step, p)

    tb = get_time()

    if(pepc_pars%pepc_comm%mpi_rank == 0) write(*,'(a,es12.4)') " == [write particles] time in vtk output [s]      : ", tb - ta
  end subroutine write_particles


  subroutine write_domain(time_pars, step, p)
    use module_pepc_types
    use module_vtk
    use module_vtk_helpers
    use module_pepc, only : global_tree

    use encap
    use time_helper
    implicit none

    type(time_pars_t), intent(in) :: time_pars
    integer, intent(in) :: step
    type(t_particle), allocatable, intent(in) :: p(:)

    integer :: vtk_step
    real*8 :: time

    time = time_of_step(step, time_pars)
    vtk_step = vtk_step_of_step(step, time_pars)

    call vtk_write_branches(step,  time, vtk_step, global_tree)
    call vtk_write_leaves(step, time, vtk_step, global_tree)
    call vtk_write_spacecurve(step, time, vtk_step, p)
  end subroutine write_domain
end module pepc_helper
