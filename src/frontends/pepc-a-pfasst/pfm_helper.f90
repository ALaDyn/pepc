module pfm_helper
  use module_debug
  implicit none
  
    type pf_nml_t
        ! Attention: defaults are overridden by settings in readn_inpf_params
        integer :: niter = 1
        integer :: nlevels = 2
        integer :: num_space_instances = 1
        logical :: color_space_div = .true.
        logical :: color_time_div = .false.
        logical :: echo_errors = .true.
        logical :: echo_timings = .true.
        real(kind=8) :: te, res_tol
        integer :: nsteps
        integer, dimension(:), allocatable :: nsweeps, nnodes
    end type pf_nml_t

contains


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!
    !! Init MPI communication, split into TIME and SPACE communicators
    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine init_mpi(pf_nml, MPI_COMM_SPACE, MPI_COMM_TIME)
        use module_pepc_types
        use pf_mod_mpi
        implicit none

        type(pf_nml_t), intent(inout) :: pf_nml
        integer(kind_default), intent(out) :: MPI_COMM_SPACE, MPI_COMM_TIME
        
        integer :: mpi_err, color
        integer(kind_pe) :: mpi_size, mpi_rank, mpi_size_space

        call pepc_status('|--> init_pfasst()')

        ! Global MPI initialization
        call MPI_Init( mpi_err )
        call MPI_Comm_size( MPI_COMM_WORLD, mpi_size, mpi_err )
        call MPI_Comm_rank( MPI_COMM_WORLD, mpi_rank, mpi_err )

        call read_in_pf_params(pf_nml, mpi_rank, MPI_COMM_WORLD)

        if (mod(mpi_size,pf_nml%num_space_instances) .ne. 0) then
            if (mpi_rank == 0) write(*,*) 'Well, this is not going to work, num_space_instances must be a factor of mpi_size:', mpi_size, pf_nml%num_space_instances
            call MPI_ABORT(MPI_COMM_WORLD,1,mpi_err)
        end if

        ! Generate spatial MPI communicator, depending on num_space_instances
        if (pf_nml%color_space_div) then
          color = mpi_rank/(mpi_size/pf_nml%num_space_instances)
        else
          color = mod(mpi_rank,pf_nml%num_space_instances)
        endif

        call MPI_COMM_SPLIT(MPI_COMM_WORLD, color, mpi_rank, MPI_COMM_SPACE, mpi_err)
        call MPI_COMM_SIZE(MPI_COMM_SPACE, mpi_size_space, mpi_err)

        ! Generate temporal MPI communicator
        if (pf_nml%color_space_div) then
          color = mpi_rank/pf_nml%num_space_instances
        else
          color = mod(mpi_rank, mpi_size_space)
        endif

        call MPI_COMM_SPLIT(MPI_COMM_WORLD, color, mpi_rank, MPI_COMM_TIME, mpi_err)       

        if (mpi_rank == 0) write(*,*) 'All right, I can use',mpi_size,'processors and these will be split up into',pf_nml%num_space_instances,'instances with',mpi_size_space,'processors each.'

    end subroutine init_mpi


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!
    !! Fill PFASST object pf (generated earlier), mostly with read-in or derived parameters
    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine fill_pfasst_object(pf, encap, sweeper, pf_nml, wk)
        use pfasst, only : pf_pfasst_t, pf_sweeper_t, PF_WINDOW_BLOCK
        use pfm_encap, only : pf_encap_t, app_data_t
        use pfm_transfer, only : interpolate, restrict
        use iso_c_binding, only : c_loc
        implicit none

        type(pf_pfasst_t), intent(inout) :: pf
        type(pf_nml_t), intent(inout) :: pf_nml
        type(pf_encap_t), target, intent(in) :: encap
        type(app_data_t), pointer, intent(in) :: wk(:)
        type(pf_sweeper_t), target, intent(in) :: sweeper

        integer :: i

        call pepc_status('|--> fill_pfasst_object()')

        do i = 1, pf%nlevels

            pf%levels(i)%nnodes    = pf_nml%nnodes(i)
            pf%levels(i)%nsweeps   = pf_nml%nsweeps(i)
            pf%levels(i)%nvars     = wk(i)%nvar

            pf%levels(i)%interpolate => interpolate
            pf%levels(i)%restrict    => restrict
            pf%levels(i)%encap       => encap
            pf%levels(i)%sweeper     => sweeper

            pf%levels(i)%ctx = c_loc(wk(i))

            pf%levels(i)%Finterp = .true.

        end do

        pf%niters       = pf_nml%niter
        pf%qtype        = 1
        pf%echo_timings = pf_nml%echo_timings .and. (wk(pf%nlevels)%mpi_rank == 0)
        pf%window       = PF_WINDOW_BLOCK
        pf%rel_res_tol  = pf_nml%res_tol
        pf%abs_res_tol  = pf_nml%res_tol

        deallocate(pf_nml%nsweeps, pf_nml%nnodes)

    end subroutine fill_pfasst_object


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!
    !! Simple read-in routine, using the first command line argument
    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine read_in_pf_params(pf_namelist, rank, comm)
        use pf_mod_mpi
        implicit none
        
        type(pf_nml_t), intent(out) :: pf_namelist
        integer, intent(in) :: rank, comm

        integer, parameter :: max_nlevels = 20

        integer :: niter = 1
        integer :: nlevels = max_nlevels
        integer :: num_space_instances = 1
        logical :: echo_errors = .true.
        logical :: echo_timings = .true.
        logical :: color_space_div = .true.
        logical :: color_time_div = .false.
        real(kind=8) :: te = 0.
        real(kind=8) :: res_tol = 0d0
        integer :: nsteps = 0
        integer, dimension(:), allocatable :: nsweeps, nnodes

        namelist /pf_nml/ niter, num_space_instances, nlevels, nnodes, echo_timings, echo_errors, te, nsteps, nsweeps, res_tol, color_space_div, color_time_div

        logical :: available
        character(len=255) :: file_name
        integer :: ierr
        integer, parameter :: para_file_id = 10

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

            allocate(nsweeps(nlevels), nnodes(nlevels)) ! use max_nlevels here only

            call MPI_BCAST( file_name, 255, MPI_CHARACTER, 0, comm, ierr )
            open(para_file_id,file=trim(file_name),action='read')
            rewind(para_file_id)
            read(para_file_id, NML=pf_nml)
            close(para_file_id)

            pf_namelist%nlevels = nlevels
            if (nlevels .gt. max_nlevels) then
                write(*,*) 'Error, nlevels is too large, must be lower than ',max_nlevels
                call MPI_ABORT(MPI_COMM_WORLD, 0, ierr)
            end if
 
            allocate(pf_namelist%nsweeps(nlevels), pf_namelist%nnodes(nlevels)) ! use exact nlevels here
            pf_namelist%nsweeps(:) = nsweeps(1:nlevels)
            pf_namelist%nnodes(:)  = nnodes(1:nlevels)

            deallocate(nsweeps, nnodes)
        else
            allocate(pf_namelist%nsweeps(nlevels), pf_namelist%nnodes(nlevels)) ! use exact nlevels here
            pf_namelist%nsweeps = 5
            pf_namelist%nnodes = 5
        end if

        pf_namelist%niter = niter
        pf_namelist%num_space_instances = num_space_instances
        pf_namelist%echo_timings = echo_timings
        pf_namelist%echo_errors = echo_errors
        pf_namelist%te = te
        pf_namelist%nsteps = nsteps
        pf_namelist%res_tol = res_tol
        pf_namelist%color_space_div = color_space_div
        pf_namelist%color_time_div = color_time_div

    end subroutine read_in_pf_params


end module pfm_helper
