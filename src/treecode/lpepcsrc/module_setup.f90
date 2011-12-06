module module_setup
  implicit none
  private


  public libpepc_setup
  public libpepc_finalize
  public libpepc_get_para_file


contains
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Initializes debug level
    !> Creates and registers user-defined MPI types
    !> Calls initialization for periodic framework
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine libpepc_setup(frontendname, my_rank,n_cpu, db_level_in)
        use treevars
        use treetypes
        use module_branching
        use module_mirror_boxes
        use module_spacefilling
        use module_tree_walk
        use module_tree_domains
        use module_debug, only : pepc_status, debug_level
        use module_calc_force, only : calc_force_init
        implicit none
        include 'mpif.h'
  
        character(*), intent(in) :: frontendname
        integer, intent(out) :: my_rank, n_cpu
        integer, intent(in), optional :: db_level_in

        integer, parameter :: para_file_id = 10
        character(len=255) :: para_file_name
        logical :: para_file_available
        integer :: ierr, provided

        integer, parameter :: MPI_THREAD_LEVEL = MPI_THREAD_FUNNELED ! "The process may be multi-threaded, but the application
                                                                       !  must ensure that only the main thread makes MPI calls."

        namelist /libpepc/ debug_level, np_mult, curve_type, weighted

        call pepc_status('SETUP')

        ! Initialize the MPI system (thread safe version, will fallback automatically if thread safety cannot be guaranteed)
        call MPI_INIT_THREAD(MPI_THREAD_LEVEL, provided, ierr)
        ! Get the id number of the current task
        call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
        ! Get the number of MPI tasks
        call MPI_COMM_size(MPI_COMM_WORLD, n_cpu, ierr)

        if (my_rank == 0) then
          ! verbose startup-output
          write(*,'(a)') "   ____    ____    ____    ____        "
          write(*,'(a)') "  /\  _`\ /\  _`\ /\  _`\ /\  _`\      "
          write(*,'(a)') "  \ \ \L\ \ \ \L\_\ \ \L\ \ \ \/\_\      The Pretty Efficient"
          write(*,'(a)') "   \ \ ,__/\ \  _\L\ \ ,__/\ \ \/_/_           Parallel Coulomb Solver"
          write(*,'(a)') "    \ \ \/  \ \ \L\ \ \ \/  \ \ \L\ \  "
          write(*,'(a)') "     \ \_\   \ \____/\ \_\   \ \____/           p.gibbon@fz-juelich.de"
          write(*,'(a)') "      \/_/    \/___/  \/_/    \/___/   "
          write(*,'(/"Starting PEPC, svn revision [",a,"] with frontend {", a, "} on ", I0, " MPI ranks."//)') &
                         SVNVERSION, frontendname, n_cpu

          if (provided < MPI_THREAD_LEVEL) then
            !inform the user about possible issues concerning MPI thread safety
            write(*,'("Call to MPI_INIT_THREAD failed. Requested/provided level of multithreading:", I2, "/" ,I2)') &
                           MPI_THREAD_LEVEL, provided
            write(*,'(a/)') "Initializing with provided level of multithreading. Usually, this is no problem."
          end if
        endif


        ! copy call parameters to treevars module
        me     = my_rank
        num_pe = n_cpu

        if (present(db_level_in)) then
            debug_level = db_level_in
        else
            debug_level = 0
        endif
        np_mult         = -45
        weighted        =   1

        ! read in parameter file
        call libpepc_get_para_file(para_file_available, para_file_name, my_rank)

        call calc_force_init(para_file_available, para_file_name, my_rank)
        call  tree_walk_init(para_file_available, para_file_name, my_rank)

        if (para_file_available) then

            open(para_file_id,file=para_file_name)

            if(my_rank .eq. 0) write(*,*) "reading parameter file, section libpepc: ", para_file_name
            read(para_file_id,NML=libpepc)

            close(para_file_id)

        else
            if ((my_rank .eq. 0) .and. (debug_level > 0)) write(*,*) "##### using default parameter for libpepc #####"
        end if


        ! create and register mpi types
        call register_lpepc_mpi_types()

        ! initialize mirror boxes
        call calc_neighbour_boxes(ws=1)

        ! initialize data structures in module_branches
        call branches_initialize()

    end subroutine libpepc_setup




    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> checks if the first application argument was set
    !> broadcasts the filename to all mpi processes
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine libpepc_get_para_file(available, file_name, my_rank)

        implicit none
        include 'mpif.h'

        logical,   intent(out)          :: available
        character(len=255), intent(out) :: file_name
        integer,   intent(in)           :: my_rank

        integer :: ierr

        ! rank 0 reads in first command line argument
        available = .false.
        if (my_rank .eq. 0) then
            if( COMMAND_ARGUMENT_COUNT() .ne. 0 ) then
                call GET_COMMAND_ARGUMENT(1, file_name)
                available = .true.
                if(my_rank .eq. 0) write(*,*) "found parameter file: ", file_name
            end if
        end if

        call MPI_BCAST( available, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr )

        ! broadcast file name, read actual inputs from namelist file
        if (available) then
            call MPI_BCAST( file_name, 255, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr )
        end if

    end subroutine libpepc_get_para_file



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Frees internal data structures
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine libpepc_finalize()
        use treevars
        use module_branching
        use module_debug, only : pepc_status
        implicit none
        include 'mpif.h'
        integer :: ierr

        call pepc_status('FINALIZE')

        ! deregister mpi types
        call free_lpepc_mpi_types()

        ! finalize data structures in module_branches
        call branches_finalize()

        call MPI_FINALIZE(ierr)

    end subroutine libpepc_finalize



end module
