module pfm_helper
  use module_debug
  implicit none
  private
  save

    public pf_nml_t
    public pfm_init_pfasst
    public pfm_fill_pfasst_object
    public pfm_setup_solver_level_params
    public pfm_finalize_solver_level_params

    integer, parameter :: max_nlevels = 20

    !> pfasst parameter collection
    type pf_nml_t
        integer :: niter                = 1
        integer :: nlevels              = 1!max_nlevels
        integer :: num_space_instances  = 1
        logical :: color_space_div      = .true.
        logical :: color_time_div       = .true.
        logical :: echo_errors          = .true.
        logical :: echo_timings         = .false.
        real(kind=8) :: tend            = 1.
        real(kind=8) :: res_tol         = 0d0
        integer :: nsteps               = 1
        integer, dimension(max_nlevels) :: nsweeps = 1
        integer, dimension(max_nlevels) :: nnodes  = 3
        real*8, dimension(max_nlevels)  :: theta = 0.3
        logical, dimension(max_nlevels) :: directforce = .true. ! PEPC does not support single-particle-invocation
        integer, dimension(max_nlevels) :: feval_mode = 11
    end type pf_nml_t

  contains


    !> Read params from file, init MPI communication, split into TIME and SPACE communicators
    subroutine pfm_init_pfasst(pf_nml, MPI_COMM_SPACE, MPI_COMM_TIME)
        use module_pepc_types
        use pf_mod_mpi
        use module_debug
        implicit none

        ! "Multiple threads may call MPI, with no restrictions." - MPI-2.2, p. 385
        integer(kind_default), parameter :: MPI_THREAD_LEVEL = MPI_THREAD_MULTIPLE

        type(pf_nml_t), intent(inout) :: pf_nml
        integer(kind_default), intent(out) :: MPI_COMM_SPACE, MPI_COMM_TIME
        integer(kind_pe) :: rank_world, nrank_world

        integer(kind_default) :: provided
        integer :: mpi_err, color
        integer(kind_pe) :: mpi_size_space

        call pepc_status('|--> init_pfasst()')

        ! Global MPI initialization
        call MPI_INIT_THREAD(MPI_THREAD_LEVEL, provided, mpi_err)
        call MPI_COMM_SIZE( MPI_COMM_WORLD, nrank_world, mpi_err )
        call MPI_COMM_RANK( MPI_COMM_WORLD, rank_world, mpi_err )

        if ((provided < MPI_THREAD_LEVEL) .and. (rank_world == 0)) then
          !inform the user about possible issues concerning MPI thread safety
          write(*,'("Call to MPI_INIT_THREAD failed. Requested/provided level of multithreading:", I2, "/" ,I2)') &
                         MPI_THREAD_LEVEL, provided
          write(*,'(a/)') 'Initializing with provided level of multithreading. This can lead to incorrect results or crashes.'
        end if

        call read_in_pf_params(pf_nml, rank_world, MPI_COMM_WORLD)

        if (mod(nrank_world,pf_nml%num_space_instances) .ne. 0) then
          DEBUG_ERROR(*, 'Well, this is not going to work, num_space_instances must be a factor of nrank_world:', nrank_world, pf_nml%num_space_instances)
        end if

        ! Generate spatial MPI communicator, depending on num_space_instances
        if (pf_nml%color_space_div) then
          color = rank_world/(nrank_world/pf_nml%num_space_instances)
        else
          color = mod(rank_world,pf_nml%num_space_instances)
        endif

        call MPI_COMM_SPLIT(MPI_COMM_WORLD, color, rank_world, MPI_COMM_SPACE, mpi_err)
        call MPI_COMM_SIZE(MPI_COMM_SPACE, mpi_size_space, mpi_err)

        ! Generate temporal MPI communicator
        if (pf_nml%color_space_div) then
          color = rank_world/pf_nml%num_space_instances
        else
          color = mod(rank_world, mpi_size_space)
        endif

        call MPI_COMM_SPLIT(MPI_COMM_WORLD, color, rank_world, MPI_COMM_TIME, mpi_err)

        if (rank_world == 0) write(*,*) 'All right, I can use ', nrank_world,&
                                        ' processors in total. These will be split up into ',pf_nml%num_space_instances,&
                                        ' instances with ',mpi_size_space,' processors each.'
    end subroutine pfm_init_pfasst


    !> Sets up parameters on each level and defines initial RHS (as "old" u value)
    subroutine pfm_setup_solver_level_params(particles, level_params, pf_nml, dim)
        use pf_mod_mpi
        use pfm_encap
        use pepcboris_helper, only: pepcboris_nml
        implicit none

        type(pf_nml_t), intent(in) :: pf_nml
        type(level_params_t), pointer, intent(inout) :: level_params(:)
        integer(kind_dim), intent(in) :: dim
        type(t_particle), allocatable, target, intent(in) :: particles(:)

        integer :: i

        allocate(level_params(pf_nml%nlevels))

        do i = 1, pf_nml%nlevels
          associate (lp=>level_params(i))
            ! add any parameters from level_params_t here
            lp%level  = i
            lp%nparts = size(particles)
            lp%theta  = pf_nml%theta(i)
            lp%directforce = pf_nml%directforce(i)
            lp%feval_mode  = pf_nml%feval_mode(i)
            lp%dim         = dim
            lp%comm        = pepcboris_nml%comm_space
          end associate
        end do
    end subroutine pfm_setup_solver_level_params


    !> Finalizes solver parameters data structure
    subroutine pfm_finalize_solver_level_params(level_params, nlevels)
        use pfm_encap
        implicit none

        type(level_params_t), intent(inout), pointer :: level_params(:)
        integer, intent(in) :: nlevels

        integer :: i

        ! not really anything to do here for now
        do i = 1, nlevels
          associate (lp=>level_params(i))
            ! add any parameters from level_params_t here
            lp%nparts = -1
            lp%theta  = -1
            lp%directforce = .false.
            lp%feval_mode  = -1
            lp%dim    = -1
            lp%comm   = -1
          end associate
        end do

        deallocate(level_params)
    end subroutine pfm_finalize_solver_level_params


    !> Fill PFASST object pf (generated earlier), mostly with read-in or derived parameters
    subroutine pfm_fill_pfasst_object(pf, encap, sweeper, pf_nml, level_params)
        use pf_mod_dtype, only: pf_pfasst_t, pf_sweeper_t, PF_WINDOW_BLOCK, SDC_GAUSS_LOBATTO
        use pfm_encap, only : pf_encap_t, level_params_t
        use pfm_transfer, only : interpolate, restrict
        use pepcboris_helper, only: pepcboris_nml
        use iso_c_binding !, only : c_loc, c_null_ptr ! had to remove this as compile fix for Michael
        implicit none

        type(pf_pfasst_t), intent(inout) :: pf
        type(pf_nml_t), intent(inout) :: pf_nml
        type(pf_encap_t), target, intent(in) :: encap
        type(level_params_t), pointer, intent(in) :: level_params(:)
        type(pf_sweeper_t), target, intent(in) :: sweeper

        integer :: i

        call pepc_status('|--> fill_pfasst_object()')

        do i = 1, pf%nlevels

            pf%levels(i)%levelctx    = c_loc(level_params(i))

            pf%levels(i)%nnodes      = pf_nml%nnodes(i)
            pf%levels(i)%nsweeps     = pf_nml%nsweeps(i)
            pf%levels(i)%nvars       = (2*level_params(i)%dim)*level_params(i)%nparts ! dim*(coordinates and momenta) per particle

            pf%levels(i)%interpolate => interpolate
            pf%levels(i)%restrict    => restrict
            pf%levels(i)%encap       => encap
            pf%levels(i)%sweeper     => sweeper

            pf%levels(i)%Finterp     = .true.

        end do

        pf%niters       = pf_nml%niter
        pf%qtype        = SDC_GAUSS_LOBATTO
        pf%echo_timings = pf_nml%echo_timings .and. pepcboris_nml%rank_space == 0
        pf%window       = PF_WINDOW_BLOCK
        pf%rel_res_tol  = pf_nml%res_tol
        pf%abs_res_tol  = pf_nml%res_tol
    end subroutine pfm_fill_pfasst_object


    !> Simple read-in routine, using the first command line argument
    subroutine read_in_pf_params(pf_namelist, rank, comm)
        use pf_mod_mpi
        implicit none

        type(pf_nml_t), intent(out) :: pf_namelist
        integer, intent(in) :: rank, comm

        integer :: niter
        integer :: nlevels
        integer :: num_space_instances
        logical :: echo_errors
        logical :: echo_timings
        logical :: color_space_div
        logical :: color_time_div
        real(kind=8) :: tend
        real(kind=8) :: res_tol
        integer :: nsteps
        integer, dimension(max_nlevels) :: nsweeps
        integer, dimension(max_nlevels) :: nnodes
        real*8, dimension(max_nlevels) :: theta
        logical, dimension(max_nlevels) :: directforce
        integer, dimension(max_nlevels) :: feval_mode

        namelist /pfasst/ niter, num_space_instances, nlevels, nnodes, echo_timings, echo_errors, tend, nsteps, nsweeps, res_tol, color_space_div, color_time_div, theta, directforce, fevaL_mode

        logical :: available
        character(len=255) :: file_name
        integer :: ierr
        integer, parameter :: para_file_id = 10

        niter               = pf_namelist%niter
        nlevels             = pf_namelist%nlevels
        num_space_instances = pf_namelist%num_space_instances
        echo_errors         = pf_namelist%echo_errors
        echo_timings        = pf_namelist%echo_timings
        color_space_div     = pf_namelist%color_space_div
        color_time_div      = pf_namelist%color_time_div
        tend                = pf_namelist%tend
        res_tol             = pf_namelist%res_tol
        nsteps              = pf_namelist%nsteps
        ! storage order in params is reverse compared to in the code
        nsweeps(1:max_nlevels)     = pf_namelist%nsweeps(max_nlevels:1:-1)
        nnodes(1:max_nlevels)      = pf_namelist%nnodes(max_nlevels:1:-1)
        theta(1:max_nlevels)       = pf_namelist%theta(max_nlevels:1:-1)
        directforce(1:max_nlevels) = pf_namelist%directforce(max_nlevels:1:-1)
        fevaL_mode(1:max_nlevels)  = pf_namelist%feval_mode(max_nlevels:1:-1)

        call pepc_status('|--> read_in_pf_params()')

        ! rank 0 reads in first command line argument
        available = .false.
        if (rank .eq. 0) then
            if( COMMAND_ARGUMENT_COUNT() .ne. 0 ) then
                call GET_COMMAND_ARGUMENT(1, file_name)
                available = .true.
                if(rank .eq. 0) write(*,*) 'found parameter file: ', file_name
            end if
        end if

        ! broadcast file name, read actual inputs from namelist file
        call MPI_BCAST( available, 1, MPI_LOGICAL, 0, comm, ierr )
        if (available) then

            call MPI_BCAST( file_name, 255, MPI_CHARACTER, 0, comm, ierr )
            open(para_file_id,file=trim(file_name),action='read')
            rewind(para_file_id)
            read(para_file_id, NML=pfasst)
            close(para_file_id)

            pf_namelist%nlevels = nlevels
            if (nlevels .gt. max_nlevels) then
                write(*,*) 'Error, nlevels is too large, must be lower than ', max_nlevels
                call MPI_ABORT(MPI_COMM_WORLD, 0, ierr)
            end if

        end if

        pf_namelist%niter                  = niter
        pf_namelist%nlevels                = nlevels
        pf_namelist%num_space_instances    = num_space_instances
        pf_namelist%color_space_div        = color_space_div
        pf_namelist%color_time_div         = color_time_div
        pf_namelist%echo_errors            = echo_errors
        pf_namelist%echo_timings           = echo_timings
        pf_namelist%tend                   = tend
        pf_namelist%res_tol                = res_tol
        pf_namelist%nsteps                 = nsteps
        ! storage order in params is reverse compared to in the code
        pf_namelist%nsweeps(1:nlevels)     = nsweeps(nlevels:1:-1)
        pf_namelist%nnodes(1:nlevels)      = nnodes(nlevels:1:-1)
        pf_namelist%theta(1:nlevels)       = theta(nlevels:1:-1)
        pf_namelist%directforce(1:nlevels) = directforce(nlevels:1:-1)
        pf_namelist%feval_mode(1:nlevels)  = feval_mode(nlevels:1:-1)
    end subroutine read_in_pf_params

end module pfm_helper
