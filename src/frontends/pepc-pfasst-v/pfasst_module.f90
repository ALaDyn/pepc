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

  real(8), allocatable, save :: qnodesF(:)
  real(8), allocatable, save :: qnodesG(:)
  real(8), allocatable, save :: SmatF(:,:)       ! F Integration table
  real(8), allocatable, save :: SmatG(:,:)       ! G Integration table
  real(8), allocatable, save :: StilLF(:,:)      ! F Approx Integration table
  real(8), allocatable, save :: StilLG(:,:)      ! G Approx Integration table
  real(8), allocatable, save :: StilnoLF(:,:)    ! F Approx Integration table
  real(8), allocatable, save :: StilnoLG(:,:)    ! G Approx Integration table
  real(8), allocatable, save :: dtfac_sdcF(:)
  real(8), allocatable, save :: dtfac_sdcG(:)
  real(8), allocatable, save :: InterpMat(:,:)    ! Interpolates G to F

  integer, parameter :: pepc_to_pfasst_attributes = 6
  real(kind=8), allocatable :: y0(:)

  ! solutions (flattened)
  real(kind=8), save, allocatable, dimension(:) :: &
       y0F, &                   !  Fine initial condition
       y0G, &                   !  Coarse initial condition
       y0newF, &                !  Parareal initial condition
       y0newG, &                !  Parareal initial condition
       yendG, &                 !  Coarse value to be sent forward
       yendF                    !  Fine value to be sent forward

  ! sdc, indexed by: (i, sdc node)
  real(kind=8), save, allocatable, dimension(:,:) :: &
       ySDC_F, &                !  All y on fine nodes
       ySDC_G, &                !  All y on coarse nodes
       tau_G                    !  Correction

  ! imex sdc, indexed by (i, sdc node, implicit/explicit)
  real(kind=8), save, allocatable, dimension(:,:,:) :: &
       fSDC_F, &                !  All function values (fine)
       fSDC_G                   !  All function values (coarse)

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

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_pfasst_comm(rank, nprocs, comm)
    implicit none
    include 'mpif.h'

    integer, intent(in) :: rank, nprocs, comm

    MPI_COMM_PFASST = comm
    my_rank_pfasst = rank
    n_cpu_pfasst = nprocs

  end subroutine init_pfasst_comm

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_pfasst(np)
    implicit none

    integer, intent(in) :: np
    integer :: nvars

    nvars = pepc_to_pfasst_attributes*np

    allocate(y0(nvars))

    NvarF = nvars
    NvarG = nvars
    NnodesF = 5
    NnodesG = 3

    call init_quadrature()

    if ((n_cpu_pfasst == 1) .and. (parallel == 1)) then
       parallel = 0
    end if

    !  Allocate
    allocate(y0F(NvarF))
    allocate(y0newF(NvarF))
    allocate(yendF(NvarF))

    allocate(ySDC_F(NvarF,NnodesF))
    allocate(fSDC_F(NvarF,NnodesF,2))

    if (parallel > 0) then
       allocate(y0G(NvarG))
       allocate(y0newG(NvarG))
       allocate(yendG(NvarG))

       allocate(ySDC_G(NvarG,NnodesG))
       allocate(fSDC_G(NvarG,NnodesG,2))
       allocate(tau_G(NvarG,NnodesG-1))
    end if

  end subroutine init_pfasst

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine finish_pfasst()
    implicit none

    deallocate(y0F, y0newF, yendF, ySDC_F, fSDC_F)

    if (parallel > 0) then
        deallocate( y0G, y0newG, yendG, ySDC_G, fSDC_G, tau_G)
    end if

    call quadrature_close()

    deallocate(y0)

  end subroutine finish_pfasst

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine receive(y, nvar, dest, tag)
    implicit none
    include 'mpif.h'

    real(kind=8), intent(out) :: y(nvar)
    integer,      intent(in)  :: nvar, dest, tag
    integer :: ierror, recvstat(MPI_STATUS_SIZE)

    if (my_rank_pfasst > 0) &
         call mpi_recv(y, nvar, MPI_REAL8, dest, tag, MPI_COMM_PFASST, recvstat, ierror)

  end subroutine receive

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine send(y, nvar, dest, tag)
    implicit none
    include 'mpif.h'

    real(kind=8), intent(in) :: y(nvar)
    integer,      intent(in) :: nvar, dest, tag
    integer :: ierror, recvstat(MPI_STATUS_SIZE)

    if (my_rank_pfasst < n_cpu_pfasst-1) &
         call mpi_send(y, nvar, MPI_REAL8, dest, tag, MPI_COMM_PFASST, recvstat, ierror)

  end subroutine send

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine broadcast(y, nvar, root)
    implicit none
    include 'mpif.h'

    real(kind=8), intent(in) :: y(nvar)
    integer,      intent(in) :: nvar, root
    integer :: ierror

    call mpi_bcast(y, nvar, MPI_REAL8, root, MPI_COMM_PFASST, ierror)

  end subroutine broadcast

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine barrier()
    implicit none
    include 'mpif.h'

    integer :: ierror

    call mpi_barrier(MPI_COMM_PFASST, ierror)

  end subroutine barrier

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine start_timer(timer)
    implicit none
    include 'mpif.h'

    integer, intent(in) :: timer

    timers(timer) = mpi_wtime()

  end subroutine start_timer

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_quadrature()

    use quadrature_smat

    implicit none

    integer :: k, kk, j, m, Ndiff
    real(8) :: den, num, trat

    allocate(qnodesF(1:NnodesF))
    allocate(qnodesG(1:NnodesG))

    allocate(SmatF(1:NnodesF-1,1:NnodesF))
    allocate(SmatG(1:NnodesG-1,1:NnodesG))

    allocate(StilLF(1:NnodesF-1,1:NnodesF))
    allocate(StilLG(1:NnodesG-1,1:NnodesG))
    allocate(StilnoLF(1:NnodesF-1,1:NnodesF))
    allocate(StilnoLG(1:NnodesG-1,1:NnodesG))

    allocate(dtfac_sdcF(1:NnodesF-1))
    allocate(dtfac_sdcG(1:NnodesG-1))

    Ndiff = NnodesF-NnodesG
    if (Ndiff > 0) then
       allocate(InterpMat(1:Ndiff,1:NnodesG))
    else
       allocate(InterpMat(1,1:NnodesG))
    end if

    !  Nodes and integration matrices
    select case (Nodes)
    case (1)                    ! GL

       if (use_coarse_nodes == 1) then
          call sdcquadGL(NnodesF, 2, SmatF, qnodesF)
       else
          call sdcquadGL(NnodesF, 1, SmatF, qnodesF)
          call sdcquadGL(NnodesF, 2, SmatG, qnodesG)
       end if

    case (2)                    ! CC

       if (use_coarse_nodes == 1) then
          call sdcquadCC(NnodesF, 2, SmatF, qnodesF)
       else
          call sdcquadCC(NnodesF, 1, SmatF, qnodesF)
          call sdcquadCC(NnodesF, 2, SmatG, qnodesG)
       end if

    case (3)

       if (use_coarse_nodes == 1) then
          call sdcquadGR(NnodesF, 2, SmatF, qnodesF)
       else
          call sdcquadGR(NnodesF, 1, SmatF, qnodesF)
          call sdcquadGR(NnodesF, 2, SmatG, qnodesG)
       end if

    case default
       print *, 'Bad case in quadrature_init, Nodes=', Nodes

    end select

    !  Make the Stil matrices
    dtfac_sdcF = qnodesF(2:NnodesF)-qnodesF(1:NnodesF-1)
    dtfac_sdcG = qnodesG(2:NnodesG)-qnodesG(1:NnodesG-1)

    StilLG = 0.0d0
    StilnoLG = 0.0d0
    do m = 1,NnodesG-1
       StilLG(m,m)=dtfac_sdcG(m)
       StilnoLG(m,m+1)=dtfac_sdcG(m)
    end do

    StilLF = 0.0d0
    StilnoLF = 0.0d0
    do m = 1,NnodesF-1
       StilLF(m,m)=dtfac_sdcF(m)
       StilnoLF(m,m+1)=dtfac_sdcF(m)
    end do

    ! Make an interpolation matrix
    ! XXX: there are more loops than necessary here...
    InterpMat=0.0d0
    trat = (NnodesF-1)/(NnodesG-1)
    if (trat == 2) then
       do j = 1,NnodesG
          den = 1.0d0
          do k = 1,NnodesG
             if (j /= k) then
                den = den*(qnodesG(j)-qnodesG(k))
             end if
             do  m = 1,Ndiff
                num = 1.0d0
                do kk = 1,NnodesG
                   if (j /= kk) then
                      num = num*(qnodesF(2*m)-qnodesG(kk))
                   end if
                end do
                InterpMat(m,j)=num/den
             end do
          end do
       end do
    end if

  end subroutine init_quadrature

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine quadrature_close()

    deallocate(dtfac_sdcF)
    deallocate(dtfac_sdcG)
    deallocate(SmatF)
    deallocate(SmatG)
    deallocate(StilLF)
    deallocate(StilnoLF)
    deallocate(StilLG)
    deallocate(StilnoLG)
    deallocate(qnodesF)
    deallocate(qnodesG)
    deallocate(InterpMat)

  end subroutine quadrature_close

end module pfasst_module
