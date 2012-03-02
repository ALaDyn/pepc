module pfasst_module

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

  integer, save :: MPI_COMM_pfasst    !  MPI communicator
  integer, save :: n_cpu_pfasst       !  Total number of processors
  integer, save :: my_rank_pfasst     !  Number of current processor

  ! Important: If you change the number of timers, also adapt the 'ntimings'
  ! parameter in the dumper module. Yes, I know it is ugly...
  integer, parameter :: &
       TTOTAL       = 1, &
       TPREDICTOR   = 2, &
       TITERATION   = 3, &
       TFINE        = 4, &
       TCOARSE      = 5, &
       TINTERPOLATE = 6, &
       TRESTRICT    = 7, &
       TRECEIVE     = 8, &
       TSEND        = 9, &
       TIO          = 10,&
       MAXTIMERS    = 11


  character(len=14), parameter :: timer_names(10) = (/ &
       'total      ', &
       'predictor  ', &
       'iteration  ', &
       'fine       ', &
       'coarse     ', &
       'interpolate', &
       'restrict   ', &
       'receive    ', &
       'send       ', &
       'io         '/)

  real(kind=8), save :: timers(MAXTIMERS)   = 0.0
  real(kind=8), save :: runtimes(MAXTIMERS) = 0.0

  integer, parameter :: pepc_to_pfasst_attributes = 6
  real(kind=8), allocatable :: y0(:)

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


  subroutine init_pfasst_comm(rank, nprocs, comm)
    implicit none
    include 'mpif.h'

    integer, intent(in) :: rank, nprocs, comm

    MPI_COMM_PFASST = comm
    my_rank_pfasst = rank
    n_cpu_pfasst = nprocs

  end subroutine init_pfasst_comm


  subroutine init_pfasst(np)
    implicit none

    integer, intent(in) :: np

    allocate(y0(pepc_to_pfasst_attributes*np))


  end subroutine init_pfasst


  subroutine finish_pfasst()
    implicit none

    deallocate(y0)

  end subroutine finish_pfasst


  subroutine receive(y, nvar, dest, tag)
    implicit none
    include 'mpif.h'

    real(kind=8), intent(out) :: y(nvar)
    integer,      intent(in)  :: nvar, dest, tag
    integer :: ierror, recvstat(MPI_STATUS_SIZE)

    if (my_rank_pfasst > 0) &
         call mpi_recv(y, nvar, MPI_REAL8, dest, tag, MPI_COMM_PFASST, recvstat, ierror)

  end subroutine receive


  subroutine send(y, nvar, dest, tag)
    implicit none
    include 'mpif.h'

    real(kind=8), intent(in) :: y(nvar)
    integer,      intent(in) :: nvar, dest, tag
    integer :: ierror, recvstat(MPI_STATUS_SIZE)

    if (my_rank_pfasst < n_cpu_pfasst-1) &
         call mpi_send(y, nvar, MPI_REAL8, dest, tag, MPI_COMM_PFASST, recvstat, ierror)

  end subroutine send

  subroutine broadcast(y, nvar, root)
    implicit none
    include 'mpif.h'

    real(kind=8), intent(in) :: y(nvar)
    integer,      intent(in) :: nvar, root
    integer :: ierror

    call mpi_bcast(y, nvar, MPI_REAL8, root, MPI_COMM_PFASST, ierror)

  end subroutine broadcast

  subroutine barrier()
    implicit none
    include 'mpif.h'

    integer :: ierror

    call mpi_barrier(MPI_COMM_PFASST, ierror)

  end subroutine barrier

  subroutine start_timer(timer)
    implicit none
    include 'mpif.h'

    integer, intent(in) :: timer

    timers(timer) = mpi_wtime()

  end subroutine start_timer

  subroutine end_timer(timer, iteration_in, echo_timings)
    implicit none
    include 'mpif.h'

    integer, intent(in) :: timer, echo_timings
    integer, intent(in), optional :: iteration_in
    integer :: iteration
    real(kind=8) :: t

    iteration = -1
    if (present(iteration_in)) iteration = iteration_in

    t = mpi_wtime()
    runtimes(timer) = runtimes(timer) + t - timers(timer)

    if (echo_timings > 0) then
       print '("timer:",a16," my_rank_pfasst: ",i3," iter: ",i3," time: ",es14.4,es14.4)', &
            timer_names(timer), my_rank_pfasst, iteration, t-timers(timer), t-timers(TTOTAL)
    end if

  end subroutine end_timer

end module pfasst_module
