module pepc_helper

   implicit none

   type pepc_nml_t
      integer :: np = 0
      integer :: pdump = 0
      integer :: fdump = 0
      integer :: cdump = 0
      !real(kind=8) :: theta = 0.3D0
   end type pepc_nml_t

   logical :: para_file_available
   character(len = 255) :: para_file_name

contains

   subroutine init_pepc(pepc_comm, pepc_nml, MPI_COMM_SPACE)
      use mpi
      use encap
      use module_pepc
      implicit none

      type(pepc_comm_t), intent(out) :: pepc_comm
      type(pepc_nml_t), intent(out) :: pepc_nml
      integer, intent(in) :: MPI_COMM_SPACE

      !logical, dimension(1:3) :: mpi_periods
      integer :: mpi_err

      pepc_comm%mpi_comm = MPI_COMM_SPACE

      call MPI_COMM_RANK(pepc_comm%mpi_comm, pepc_comm%mpi_rank, mpi_err)
      call MPI_COMM_SIZE(pepc_comm%mpi_comm, pepc_comm%mpi_size, mpi_err)

      call pepc_initialize("pepc-kh" ,pepc_comm%mpi_rank, pepc_comm%mpi_size, .false., 0, pepc_comm%mpi_comm, idim = 2)
      call pepc_read_parameters_from_first_argument(para_file_available, para_file_name)

      call read_in_params(pepc_nml, para_file_available, para_file_name)

   end subroutine init_pepc


   subroutine pepc_setup(p, pepc_pars, pepc_comm, pepc_nml)
      use encap
      use module_pepc, only: pepc_prepare
      use module_pepc_types, only: t_particle
      use module_interaction_specific, only: particleresults_clear
      implicit none

      type(t_particle), dimension(:), allocatable, intent(out) :: p
      type(pepc_pars_t), intent(out) :: pepc_pars
      type(pepc_comm_t), intent(in)  :: pepc_comm
      type(pepc_nml_t),  intent(in)  :: pepc_nml

      ! Pass MPI stuff to parameters
      pepc_pars%pepc_comm = pepc_comm
      pepc_pars%np = pepc_nml%np
      pepc_pars%npp = pepc_pars%np / pepc_comm%mpi_size
      if (pepc_pars%pepc_comm%mpi_rank < mod(pepc_pars%np, pepc_pars%pepc_comm%mpi_size)) then
        pepc_pars%npp = pepc_pars%npp + 1
      end if
      !pepc_pars%theta = pepc_nml%theta
      pepc_pars%pdump = pepc_nml%pdump
      pepc_pars%fdump = pepc_nml%fdump
      pepc_pars%cdump = pepc_nml%cdump

      allocate(p(1:pepc_pars%npp))

      call pepc_prepare(2)

   end subroutine pepc_setup


   subroutine read_in_params(pepc_namelist, file_available, file_name)
      use encap
      use mpi
      use module_mirror_boxes, only: t_lattice_1, t_lattice_2, t_lattice_3, &
        periodicity, mirror_box_layers
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
      !real(kind=8) :: theta = 0.3D0

      namelist /pepc_nml/ np, pdump, fdump, cdump, &
        t_lattice_1, t_lattice_2, t_lattice_3, periodicity, mirror_box_layers ! ,theta

      if (file_available) then
        open(para_file_id,file=trim(file_name),action='read')
        rewind(para_file_id)
        read(para_file_id, NML=pepc_nml)
        close(para_file_id)
      end if

      pepc_namelist%np = np
      !pepc_namelist%theta = theta
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
      integer :: np = 0
      integer :: pdump = 0
      integer :: fdump = 0
      integer :: cdump = 0
      !real(kind=8) :: theta = 0.3D0

      namelist /pepc_nml/ np, pdump, fdump, cdump,&
        t_lattice_1, t_lattice_2, t_lattice_3, periodicity, mirror_box_layers ! ,theta

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

    integer :: i
    type(vtkfile_unstructured_grid) :: vtk
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

    call vtk%create_parallel("particles", step, pepc_pars%pepc_comm%mpi_rank, &
      pepc_pars%pepc_comm%mpi_size, time, vtk_step)
    call vtk%write_headers(pepc_pars%npp, 0)
    call vtk%startpoints()
    call vtk%write_data_array("xyz", pepc_pars%npp, p(:)%x(1), p(:)%x(2), p(:)%x(3))
    call vtk%finishpoints()
    call vtk%startpointdata()
    call vtk%write_data_array("velocity", pepc_pars%npp, p(:)%data%v(1), p(:)%data%v(2), p(:)%data%v(3))
    call vtk%write_data_array("el_field", pepc_pars%npp, p(:)%results%e(1), \
                              p(:)%results%e(2), p(:)%results%e(3))
    call vtk%write_data_array("el_pot", pepc_pars%npp, p(:)%results%pot)
    call vtk%write_data_array("charge", pepc_pars%npp, p(:)%data%q)
    call vtk%write_data_array("mass", pepc_pars%npp, p(:)%data%m)
    call vtk%write_data_array("work", pepc_pars%npp, p(:)%work)
    call vtk%write_data_array("pelabel", pepc_pars%npp, p(:)%label)
    call vtk%write_data_array("local index", pepc_pars%npp, [(i,i=1,pepc_pars%npp)])
    call vtk%write_data_array("processor", pepc_pars%npp, p(:)%pid)
    call vtk%finishpointdata()
    call vtk%dont_write_cells()
    call vtk%write_final()
    call vtk%close()

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
    call write_spacecurve_to_vtk(step, time_pars%dt * step, vtk_step, p)
    
  end subroutine write_domain

end module pepc_helper
