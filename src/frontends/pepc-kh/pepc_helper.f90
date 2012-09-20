module pepc_helper

   implicit none

   type pepc_nml_t
      integer :: np = 0
      integer :: pdump = 0
      integer :: fdump = 0
      real(kind=8) :: theta = 0.3D0
   end type pepc_nml_t

contains

   subroutine init_pepc(pepc_comm, pepc_nml, MPI_COMM_SPACE)
      use mpi
      use encap
      use module_pepc
      implicit none

      type(pepc_comm_t), intent(out) :: pepc_comm
      type(pepc_nml_t), intent(out) :: pepc_nml
      integer, intent(in) :: MPI_COMM_SPACE

      logical, dimension(1:3) :: mpi_periods
      integer :: mpi_err

      pepc_comm%mpi_comm = MPI_COMM_SPACE

      call MPI_COMM_RANK(pepc_comm%mpi_comm, pepc_comm%mpi_rank, mpi_err)
      call MPI_COMM_SIZE(pepc_comm%mpi_comm, pepc_comm%mpi_size, mpi_err)

      call read_in_pepc_params(pepc_nml, pepc_comm%mpi_rank, pepc_comm%mpi_comm)

      call pepc_initialize("pepc-kh" ,pepc_comm%mpi_rank, pepc_comm%mpi_size, .false., 0, pepc_comm%mpi_comm)

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

      integer :: i
      !double precision, parameter :: pi = 3.1415926535897932384626434D0


      ! Pass MPI stuff to parameters
      pepc_pars%pepc_comm = pepc_comm
      pepc_pars%idump = pepc_nml%idump

      pepc_pars%np = pepc_nml%np
      pepc_pars%npp = ceiling(1.0*pepc_pars(i)%np/pepc_comm%mpi_size) !PRELIMINARY
      pepc_pars%theta = pepc_nml%theta
      pepc_pars%pdump = pepc_nml%pdump
      pepc_pars%fdump = pepc_nml%fdump


      allocate(p(1:pepc_pars%npp))
      call particleresults_clear(p, pepc_pars%npp)

      call special_start(p,pepc_pars)

      call pepc_prepare(3)

   end subroutine pepc_setup


   subroutine special_start(p,pepc_pars)
      use encap
      implicit none

      type(pepc_pars_t), intent(in) :: pepc_pars
      type(t_particle), dimension(pepc_pars%npp), intent(inout) :: p

      ! TODO

   end subroutine


   subroutine read_in_pepc_params(pepc_namelist, mpi_rank, mpi_comm)
      use mpi
      use module_walk
      use module_pepc
      use module_domains, only: weighted
      use treevars, only: interaction_list_length_factor
      use module_walk, only: num_walk_threads
      use module_libpepc_main, only: libpepc_read_parameters
      implicit none

      integer, intent(in) :: mpi_rank, mpi_comm
      type(pepc_nml_t), intent(out) :: pepc_namelist

      logical :: available
      character(len=255) :: file_name
      integer :: ierr
      integer, parameter :: para_file_id = 10

      ! variables for pepc namelist
      integer :: np = 0
      integer :: pdump = 0
      integer :: fdump = 0
      real(kind=8) :: theta = 0.3D0

      namelist /pepc_nml/ np, theta, pdump, fdump

      ! rank 0 reads in first command line argument
      available = .false.
      if (mpi_rank .eq. 0) then
         if( COMMAND_ARGUMENT_COUNT() .ne. 0 ) then
            call GET_COMMAND_ARGUMENT(1, file_name)
            available = .true.
             !if(mpi_rank .eq. 0) write(*,*) "found parameter file: ", file_name
         end if
      end if

      ! broadcast file name, read actual inputs from namelist file (name ist fixed, sorry!)
      call MPI_BCAST( available, 1, MPI_LOGICAL, 0, mpi_comm, ierr )
      if (available) then
         call MPI_BCAST( file_name, 255, MPI_CHARACTER, 0, mpi_comm, ierr )
         open(para_file_id,file=trim(file_name),action='read')
         rewind(para_file_id)
         read(para_file_id, NML=pepc_nml)

         pepc_namelist%np = np
         pepc_namelist%theta = theta
         pepc_namelist%pdump = pdump
         pepc_namelist%fdump = fdump

         interaction_list_length_factor = 4
         weighted = 0
         num_walk_threads = 4

         rewind(para_file_id)
         call tree_walk_read_parameters(para_file_id)
         rewind(para_file_id)
         call libpepc_read_parameters(para_file_id)

         close(para_file_id)

      end if

   end subroutine read_in_pepc_params


end module pepc_helper
