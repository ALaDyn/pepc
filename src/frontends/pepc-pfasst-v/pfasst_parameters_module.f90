module pfasst_parameters_module

  implicit none

  integer, save :: &
       NvarF = -1, &           ! number of fine variables in system
       NvarG = -1, &           ! number of coarse variables in system
       Nodes = 1, &            ! type of nodes: 1 GL, 2: CC
       NnodesF = -1, &         ! number of fine time nodes
       NnodesG = -1, &         ! number of coarse time nodes
       Niter = 7               ! number of pfasst iterations
  integer :: num_space_instances

  namelist /pfasst/ NvarF, NvarG
  namelist /pfasst/ Nodes, NnodesF, NnodesG
  namelist /pfasst/ Niter, num_space_instances

  ! flags etc
  integer, save :: &
       parallel = 1, &             ! 1: pfasst run
                                   ! 0: serial run, stop iterating according to serial_residual_tol
                                   ! -1: serial run, fixed number of iterations
                                   ! -2: serial run, using Runge-Kutta-3 instead of SDC
       echo_timings      = 0, &    ! echo timings if non-zero
       echo_errors       = 0, &    ! echo errors if non-zero
       num_coarse_sweeps = 3, &    ! number of coarse SDC sweeps to perform
       num_fine_sweeps   = 1, &    ! number of fine SDC sweeps to perform
       use_coarse_nodes  = 0, &    ! use "coarse" nodes fine level?
       use_fas = 1                 ! apply fas corrections?

  namelist /pfasst/ parallel, echo_timings, echo_errors
  namelist /pfasst/ num_coarse_sweeps, num_fine_sweeps
  namelist /pfasst/ use_fas, use_coarse_nodes

  real(kind=8), save :: &
       serial_residual_tol = 1e-10, &  ! for serial runs: skip to next
                                       ! time step when the residual is
                                       ! less than the tolerance
       pfasst_residual_tol = 1e-10, &  ! for pfasst runs: send same
                                       ! yend forward if residual is
                                       ! less than this tolerance
       pfasst_u0_tol = 1e-10           ! for pfasst runs: don't sweep
                                       ! if delta u0 is less than this
                                       ! tolerance

  namelist /pfasst/ serial_residual_tol, pfasst_residual_tol, pfasst_u0_tol

  ! input/output
  integer, save :: dump_interval

  character(len=128) :: output = ''  ! output file name
  character(len=128) :: initial = '' ! XXX: initial condition file name

  namelist /pfasst/ initial, dump_interval

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_pfasst_parameters(rank, comm)

    use module_pepc, only : pepc_get_para_file
    implicit none

    integer, intent(in) :: rank, comm

    character(len=255) :: filename
    logical :: read_param_file

    call pepc_get_para_file(read_param_file, filename, rank, comm)

    open(7, FILE=filename)
    read(UNIT=7, NML=pfasst)
    close(7)

  end subroutine init_pfasst_parameters

end module pfasst_parameters_module
