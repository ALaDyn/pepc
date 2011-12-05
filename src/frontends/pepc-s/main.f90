program pepcs

  use module_pepcs
  implicit none
  include 'mpif.h'

  integer :: nparts, nparts_total
  real*8, allocatable :: xyz(:,:), field(:,:), pot(:), q(:)
  real*8, dimension(3) :: lx, ly, lz
  real*8 :: virial(3,3)
  integer :: lperiod(3), extrcorr

  real*8, parameter :: eps   = 0.05
  real*8, parameter :: theta = 0.30
  integer, parameter :: db_level = 2


  integer :: it, ip

  integer :: my_rank, n_cpu, ierr, provided
  integer, parameter :: MPI_THREAD_LEVEL = MPI_THREAD_FUNNELED ! "The process may be multi-threaded, but the application
                                                                  !  must ensure that only the main thread makes MPI calls."

  ! Initialize the MPI system (thread safe version, will fallback automatically if thread safety cannot be guaranteed)
  call MPI_INIT_THREAD(MPI_THREAD_LEVEL, provided, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
  call MPI_COMM_size(MPI_COMM_WORLD, n_cpu, ierr)

  ! inform the user about possible issues concerning MPI thread safety
  if ((my_rank == 0) .and. (provided < MPI_THREAD_LEVEL)) then
    write(*,'("Call to MPI_INIT_THREAD failed. Requested/provided level of multithreading:", I2, "/" ,I2)') &
         MPI_THREAD_LEVEL, provided
    write(*,*) "Initializing with provided level of multithreading. Stability is possibly not guaranteed."
  end if

  lx = [1, 0, 0]
  ly = [0, 1, 0]
  lz = [0, 0, 1]

  lperiod = [0, 0, 0]
  extrcorr = 0

  do it=1, 10

     write(*,*) "-- doing step ", it

     nparts = (my_rank + 15) * it

     write(*,*) " - number of particles on rank ", my_rank, " is ", nparts

     allocate(xyz(1:3,nparts), field(1:3,nparts), pot(nparts), q(nparts))

     call MPI_ALLREDUCE(nparts, nparts_total, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

     if(my_rank .eq. 0) write(*,*) " - total number of particles on rank ", my_rank, " is ", nparts_total

     do ip=1, nparts
        xyz(1,ip) = ip/(1.0*nparts) + my_rank/(1.0*n_cpu)
        xyz(2,ip) = ip/(1.0*nparts) + my_rank/(1.0*n_cpu) + ip
        xyz(3,ip) = ip/(1.0*nparts) + my_rank/(1.0*n_cpu) + 2
          q(  ip) = 0.13
     end do

     call pepc(nparts, nparts, nparts_total,    &
                 xyz, q, field, pot, virial,    &
                 lx, ly, lz, lperiod, extrcorr, &
                 eps, theta, db_level)

     deallocate(xyz, field, pot, q)

  end do

  call MPI_FINALIZE(ierr)

end program pepcs
