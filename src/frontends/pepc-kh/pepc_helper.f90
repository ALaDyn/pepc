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


   subroutine pepc_setup(p, pepc_pars)
      use encap
      use module_pepc
      use module_pepc_types, only: t_particle, kind_dim, kind_particle
      implicit none

      type(t_particle), dimension(:), allocatable, intent(out) :: p
      type(pepc_pars_t), intent(out) :: pepc_pars

      type(pepc_nml_t) :: pepc_nml

      call pepc_initialize("pepc-kh", pepc_pars%pepc_comm%mpi_rank, &
        pepc_pars%pepc_comm%mpi_size, .true., comm = pepc_pars%pepc_comm%mpi_comm)

      call pepc_read_parameters_from_first_argument(para_file_available, para_file_name)
      call read_in_params(pepc_nml, para_file_available, para_file_name)

      ! Pass MPI stuff to parameters
      pepc_pars%np = pepc_nml%np
      pepc_pars%npp = pepc_pars%np / pepc_pars%pepc_comm%mpi_size
      if (pepc_pars%pepc_comm%mpi_rank < mod(pepc_pars%np, int(pepc_pars%pepc_comm%mpi_size, kind=kind_particle))) then
        pepc_pars%npp = pepc_pars%npp + 1
      end if
      pepc_pars%pdump = pepc_nml%pdump
      pepc_pars%fdump = pepc_nml%fdump
      pepc_pars%cdump = pepc_nml%cdump

      allocate(p(1:pepc_pars%npp))

   end subroutine pepc_setup


   subroutine read_in_params(pepc_namelist, file_available, file_name)
      use encap
      use mpi
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
      use module_mirror_boxes, only: t_lattice_1, t_lattice_2, t_lattice_3, &
        periodicity, mirror_box_layers
      implicit none

      type(pepc_pars_t), intent(in) :: pepc_pars
      character(len = 255), intent(in) :: file_name

      integer, parameter :: para_file_id = 10

      ! variables for pepc namelist
      integer(kind_particle) :: np = 0
      integer :: pdump = 0
      integer :: fdump = 0
      integer :: cdump = 0

      namelist /pepc_nml/ np, pdump, fdump, cdump,&
        t_lattice_1, t_lattice_2, t_lattice_3, periodicity, mirror_box_layers

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


  subroutine write_particles(pepc_pars, time_pars, step, p)
    use module_pepc_types
    use module_vtk
    use encap
    implicit none
    
    type(pepc_pars_t), intent(in) :: pepc_pars
    type(time_pars_t), intent(in) :: time_pars
    type(t_particle), allocatable, intent(in) :: p(:)
    integer, intent(in) :: step

    integer(kind_particle) :: i
    type(vtkfile_unstructured_grid) :: vtk
    integer :: vtk_step
    real*8 :: time
    real*8 :: ta, tb
    real*8, allocatable :: dummy_ez(:)
    
    ta = get_time()
    
    allocate(dummy_ez(pepc_pars%npp))
    dummy_ez = 0

    time = time_pars%dt * step

    if (step .eq. 0) then
      vtk_step = VTK_STEP_FIRST
    else if (step .eq. time_pars%nsteps) then
      vtk_step = VTK_STEP_LAST
    else
      vtk_step = VTK_STEP_NORMAL
    endif

    call vtk%create_parallel("particles", step, pepc_pars%pepc_comm%mpi_rank, &
      pepc_pars%pepc_comm%mpi_size, time, vtk_step)
    call vtk%write_headers(pepc_pars%npp, 0_kind_particle)
    call vtk%startpoints()
    call vtk%write_data_array("xyz", p(:)%x(1), p(:)%x(2), p(:)%x(3))
    call vtk%finishpoints()
    call vtk%startpointdata()
    call vtk%write_data_array("velocity", p(:)%data%v(1), p(:)%data%v(2), p(:)%data%v(3))
    call vtk%write_data_array("el_field", p(:)%results%e(1), &
                              p(:)%results%e(2), dummy_ez)
    call vtk%write_data_array("el_pot", p(:)%results%pot)
    call vtk%write_data_array("charge", p(:)%data%q)
    call vtk%write_data_array("mass", p(:)%data%m)
    call vtk%write_data_array("work", p(:)%work)
    call vtk%write_data_array("pelabel", p(:)%label)
    call vtk%write_data_array("local index", [(i,i=1,pepc_pars%npp)])
    call vtk%write_data_array("processor", p(:)%pid)
    call vtk%finishpointdata()
    call vtk%dont_write_cells()
    call vtk%write_final()
    call vtk%close()

    deallocate(dummy_ez)

    tb = get_time()

    if(pepc_pars%pepc_comm%mpi_rank == 0) write(*,'(a,es12.4)') " == [write particles] time in vtk output [s]      : ", tb - ta

  end subroutine write_particles


  subroutine write_domain(time_pars, step, p)
    use module_pepc_types
    use module_vtk
    use module_treediags

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

    call write_branches_to_vtk(step,  time_pars%dt * step, vtk_step)
    call write_leaves_to_vtk(step, time_pars%dt * step, vtk_step)
    call write_spacecurve_to_vtk(step, time_pars%dt * step, vtk_step, p)
    
  end subroutine write_domain

end module pepc_helper
