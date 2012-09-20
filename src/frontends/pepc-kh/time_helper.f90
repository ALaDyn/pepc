module time_helper

   implicit none

   type time_nml_t
      real(kind=8) :: te
      integer :: nsteps
   end type time_nml_t

contains


   subroutine setup_time(time_pars,pepc_comm)
      use encap

      type(time_pars_t), intent(out) :: time_pars
      type(pepc_comm_t), intent(in) :: pepc_comm

      type(time_nml_t) :: time_nml

      call read_in_time_params(time_nml,pepc_comm%mpi_rank, pepc_comm%mpi_comm)

      time_pars%te = time_nml%te
      time_pars%nsteps = time_nml%nsteps
      time_pars%dt = time_pars%te/time_pars%nsteps

   end subroutine setup_time

   subroutine push_particles(p, pepc_pars)
      use encap
      implicit none

      type(pepc_pars_t), intent(in) :: pepc_pars
      type(t_particle), dimension(pepc_pars%npp), intent(inout) :: p

      ! TODO

   end subroutine push_particles



   subroutine read_in_time_params(time_namelist, mpi_rank, mpi_comm)
      use mpi
      implicit none

      type(time_nml_t), intent(out) :: time_namelist
      integer, intent(in) :: mpi_rank, mpi_comm

      real(kind=8) :: te = 0.
      integer :: nsteps = 0

      namelist /time_nml/ te, nsteps

      logical :: available
      character(len=255) :: file_name
      integer :: ierr
      integer, parameter :: para_file_id = 10

      ! rank 0 reads in first command line argument
      available = .false.
      if (mpi_rank .eq. 0) then
         if( COMMAND_ARGUMENT_COUNT() .ne. 0 ) then
            call GET_COMMAND_ARGUMENT(1, file_name)
            available = .true.
             !if(mpi_rank .eq. 0) write(*,*) "found parameter file: ", file_name
         end if
      end if

      ! broadcast file name, read actual inputs from namelist file
      call MPI_BCAST( available, 1, MPI_LOGICAL, 0, mpi_comm, ierr )
      if (available) then
         call MPI_BCAST( file_name, 255, MPI_CHARACTER, 0, mpi_comm, ierr )
         open(para_file_id,file=trim(file_name),action='read')
         rewind(para_file_id)
         read(para_file_id, NML=pf_nml)
         close(para_file_id)

         pf_namelist%te = te
         pf_namelist%nsteps = nsteps

      end if

   end subroutine read_in_time_params


end module time_helper
