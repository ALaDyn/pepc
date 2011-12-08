!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!> Contains all routines that should be callable from the frontend
!> in most cases
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_pepc
    implicit none
    private

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    public pepc_initialize                !< once per simulation

    public pepc_grow_tree                 !< once per timestep
    public pepc_traverse_tree             !< several times per timestep
    public pepc_statistics                !< once per timestep
    public pepc_restore_particles         !< once or never per timestep
    public pepc_timber_tree               !< once or never per timestep

    public pepc_grow_and_traverse         !< once per timestep, calls pepc_grow_tree, pepc_traverse_tree, pepc_statistics, pepc_restore_particles, pepc_timber_tree

    public pepc_finalize                  !< once per simulation

    public pepc_get_para_file

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  subroutine-implementation  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    contains


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Initializes MPI library and data structures for treecode kernel,
    !> reads several parameters from file, that is given as first parameter
    !> to actual executable, initializes submodules
    !> 
    !> Call this function at program startup before any MPI calls
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine pepc_initialize(frontendname, my_rank,n_cpu, db_level_in)
      use treevars, only : me, num_pe, np_mult
      use treetypes, only : register_lpepc_mpi_types
      use module_walk
      use module_mirror_boxes
      use module_spacefilling
      use module_domains
      use module_debug, only : pepc_status, debug_level
      use module_calc_force, only : calc_force_init
      use module_branching, only : branches_initialize
      implicit none
      include 'mpif.h'
      character(*), intent(in) :: frontendname !< name of the program that uses the treecode (only for output purposes)
      integer, intent(out) :: my_rank !< MPI rank of this instance as returned from MPI
      integer, intent(out) :: n_cpu !< number of MPI ranks as returned from MPI
      integer, intent(in), optional :: db_level_in !< sets debug level for treecode kernel (overrides settings, that may be read from libpepc-section in input file)

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
      np_mult         = -45.
      weighted        =   1

      ! read in parameter file
      call pepc_get_para_file(para_file_available, para_file_name, my_rank)
      call calc_force_init(para_file_available, para_file_name, my_rank)
      call tree_walk_init(para_file_available, para_file_name, my_rank)

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
      call calc_neighbour_boxes()

      ! initialize data structures in module_branches
      call branches_initialize()

    end subroutine



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Finalizes MPI library and reverses all initialization from pepc_initialize
    !> Call this function at program termination after all MPI calls
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine pepc_finalize()
      use module_branching
      use module_debug, only : pepc_status
      use treetypes, only : free_lpepc_mpi_types
      implicit none
      include 'mpif.h'
      integer :: ierr

      call pepc_status('FINALIZE')
      ! deregister mpi types
      call free_lpepc_mpi_types()
      ! finalize data structures in module_branches
      call branches_finalize()

      call MPI_FINALIZE(ierr)
    end subroutine


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> checks if the first application argument was set
    !> broadcasts the filename to all mpi processes
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine pepc_get_para_file(available, file_name, my_rank)

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

    end subroutine pepc_get_para_file


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Builds the tree from the given particles, redistributes particles
    !> to other MPI ranks if necessary (i.e. reallocates particles and changes np_local)
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine pepc_grow_and_traverse(np_local, npart_total, particles, itime, clearresults_before_traversal)
      use treetypes
      use module_libpepc_main
      use module_debug
      implicit none
      integer, intent(inout) :: np_local    !< number of particles on this CPU, i.e. number of particles in particles-array
      integer, intent(in) :: npart_total !< total number of simulation particles (sum over np_local over all MPI ranks)
      type(t_particle), allocatable, intent(inout) :: particles(:) !< input particle data, initializes %x, %data, %work appropriately (and optionally set %label) before calling this function
      integer, intent(in) :: itime !> current timestep (used as filename suffix for statistics output)
      logical, intent(in), optional :: clearresults_before_traversal !< if set to .false., the function does not call particleresults_clear(particles) before traversal

      call pepc_grow_tree(np_local, npart_total, particles)
      call pepc_traverse_tree(np_local, particles, clearresults_before_traversal)
      if (dbg(DBG_STATS)) then
        call pepc_statistics(itime)
      endif
      call pepc_restore_particles(np_local, particles)
      call pepc_timber_tree()

    end subroutine



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Builds the tree from the given particles, redistributes particles
    !> to other MPI ranks if necessary (i.e. reallocates particles and changes np_local)
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine pepc_grow_tree(np_local, npart_total, particles)
      use treetypes
      use module_libpepc_main
      implicit none
      integer, intent(inout) :: np_local    !< number of particles on this CPU, i.e. number of particles in particles-array
      integer, intent(in) :: npart_total !< total number of simulation particles (sum over np_local over all MPI ranks)
      type(t_particle), allocatable, intent(inout) :: particles(:) !< input particle data, initializes %x, %data, %work appropriately (and optionally set %label) before calling this function

      call libpepc_grow_tree(np_local, npart_total, particles)

    end subroutine


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Traverses the complete tree for the given particles, i.e. computes
    !> the field values at their positions. Although missing information
    !> is automatically requested from remote MPI ranks, it is important
    !> that the particle coordinates fit to the local MPI ranks domain
    !> to avoid excessive communication
    !> If field values on some regular grid are needed, they can be
    !> generated using pepc_prepare_local_grid() [TODO: provide this function]
    !> Otherwise, it makes sense to provide the same particles as given/returned
    !> from to pepc_grow_tree()
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine pepc_traverse_tree(nparticles, particles, clearresults_before_traversal)
      use treetypes
      use module_libpepc_main
      implicit none
      integer, intent(in) :: nparticles    !< number of particles on this CPU, i.e. number of particles in particles-array
      type(t_particle), allocatable, intent(inout) :: particles(:) !< input particle data, initializes %x, %data, %work appropriately (and optionally set %label) before calling this function
      logical, intent(in), optional :: clearresults_before_traversal !< if set to .false., the function does not call particleresults_clear(particles) before traversal

      call libpepc_traverse_tree(nparticles, particles, clearresults_before_traversal)

    end subroutine



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Writes detailed statistics on the treecode into stats/stats.ITIME.
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine pepc_statistics(itime)
        use module_debug
        use treevars
        use module_timings
        implicit none
        integer, intent(in) :: itime !< current timestep (used as file suffix)
        character(30) :: cfile

        call timer_start(t_fields_stats)
        call tree_stats(itime)

        if( dbg(DBG_LOADFILE) ) then
            call system("mkdir -p " // "load")
            write(cfile,'("load/load_",i6.6,".dat")') me
            open(60, file=trim(cfile),STATUS='UNKNOWN', POSITION = 'APPEND')
            write(60,'(i5,2f20.10, i12)') itime,interactions_local, mac_evaluations_local,npp
            close(60)
        end if

        call timer_stop(t_fields_stats)

    end subroutine


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Restores the initial particle distribution (before calling pepc_grow_tree() ).
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine pepc_restore_particles(np_local, particles)
      use treetypes
      use module_libpepc_main
      implicit none
      integer, intent(inout) :: np_local    !< number of particles on this CPU, i.e. number of particles in particles-array
      type(t_particle), allocatable, intent(inout) :: particles(:) !< input particle data on local MPI rank - is replaced by original particle data that was given before calling pepc_grow_tree()

      call libpepc_restore_particles(np_local, particles)

    end subroutine



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Frees all tree_specific data fields that were allocated in pepc_groe_tree().
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine pepc_timber_tree()
      use module_timings
      use module_debug, only : pepc_status
      use module_allocation, only : deallocate_tree
      implicit none

      call pepc_status('TIMBER TREE')

     ! deallocate particle and result arrays
      call timer_start(t_deallocate)
      call deallocate_tree()
      call timer_stop(t_deallocate)
      call timer_stop(t_all)

      call pepc_status('TREE HAS FALLEN')

    end subroutine


end module module_pepc
