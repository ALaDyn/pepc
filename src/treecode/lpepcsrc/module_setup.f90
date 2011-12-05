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
    subroutine libpepc_setup(my_rank,n_cpu, db_level_in)
        use treevars
        use treetypes
        use module_branching
        use module_mirror_boxes
        use module_spacefilling
        use module_tree_walk
        use module_tree_domains
        implicit none
  
        integer, intent(in) :: my_rank, n_cpu
        integer, intent(in), optional :: db_level_in

        integer :: db_level

        integer, parameter :: para_file_id = 10
        integer :: file_start=0
        character(len=255) :: para_file_name
        integer :: para_file_available

        namelist /libpepc/ db_level, np_mult, curve_type, weighted

        call status('SETUP')

        ! copy call parameters to treevars module
        me     = my_rank
        num_pe = n_cpu

        if (present(db_level_in)) then
            db_level = db_level_in
        else
            db_level        =   0
        endif
        np_mult         = -45
        weighted        =   1

        ! read in parameter file
        call libpepc_get_para_file(para_file_available, para_file_name, my_rank)

        if (para_file_available .eq. 1) then

            open(para_file_id,file=para_file_name)

            if(my_rank .eq. 0) write(*,*) "reading parameter file, section libpepc: ", para_file_name
            read(para_file_id,NML=libpepc)

            rewind(para_file_id, iostat=file_start)

            if(my_rank .eq. 0) write(*,*) "reading parameter file, section walk_para: ", para_file_name
            read(para_file_id,NML=walk_para)

            close(para_file_id)

        else
            if(my_rank .eq. 0) write(*,*) "##### using default parameter for libpepc #####"
        end if


        force_debug=.false.
        tree_debug=.false.
        build_debug=.false.
        domain_debug = .false.
        branch_debug=.false.
        walk_summary=.false.
        dump_tree=.false.
        periodic_debug=.false.

        load_file_debug=.false.
        timing_file_debug=.false.


        if (db_level==1) then
            !      domain_debug = .true.
            tree_debug=.true.  ! location information only
	
        else if (db_level==2) then
            tree_debug=.true.  ! location information only
            walk_summary=.true.

        else if (db_level==3) then
            tree_debug=.true.
            force_debug=.false.
            walk_summary=.true.
            domain_debug = .true.
            periodic_debug=.true.

        else if (db_level==4) then
            tree_debug=.true.
            domain_debug = .true.
            build_debug=.true.
            branch_debug=.true.
            props_debug=.true.
            walk_summary=.true.
            force_debug=.true.
            periodic_debug=.true.

        else if (db_level==5) then
            dump_tree=.true.

        else if (db_level==6) then
            tree_debug=.true.
            domain_debug = .true.
            build_debug=.true.
            branch_debug=.true.
            walk_summary=.true.
            force_debug=.true.
            dump_tree=.true.
            periodic_debug=.true.
            load_file_debug=.true.
            timing_file_debug=.true.

        else
          ! all off by default
        endif

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

        integer,   intent(out)          :: available
        character(len=255), intent(out) :: file_name
        integer,   intent(in)           :: my_rank

        integer :: ierr

        ! rank 0 reads in first command line argument
        available = 0
        if (my_rank .eq. 0) then
            if( COMMAND_ARGUMENT_COUNT() .ne. 0 ) then
                call GET_COMMAND_ARGUMENT(1, file_name)
                available = 1
                if(my_rank .eq. 0) write(*,*) "reading parameter file, section pepce: ", file_name
            end if
        end if

        call MPI_BCAST( available, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )

        ! broadcast file name, read actual inputs from namelist file
        if (available .eq. 1) then
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
        implicit none

        call status('FINALIZE')

        ! deregister mpi types
        call free_lpepc_mpi_types()

        ! finalize data structures in module_branches
        call branches_finalize()

    end subroutine libpepc_finalize



end module
