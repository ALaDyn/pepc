program pepcs

  implicit none
  include 'mpif.h'

  integer :: nparts, nparts_total
  real*8, dimension(:), allocatable :: x, y, z, ex, ey, ez, pot, q, m
  real*8, dimension(3) :: lx, ly, lz
  logical :: lperiod(3), extrcorr

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

  lperiod = [.false., .false., .false.]
  extrcorr = .false.

  do it=1, 10

     write(*,*) "-- doing step ", it

     nparts = (my_rank + 15) * it

     write(*,*) " - number of particles on rank ", my_rank, " is ", nparts

     allocate(x(nparts), y(nparts), z(nparts), ex(nparts), ey(nparts), ez(nparts))
     allocate(pot(nparts), q(nparts), m(nparts))

     call MPI_ALLREDUCE(nparts, nparts_total, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

     if(my_rank .eq. 0) write(*,*) " - total number of paricles on rank ", my_rank, " is ", nparts_total

     do ip=1, nparts
        x(ip) = ip/(1.0*nparts) + my_rank/(1.0*n_cpu)
        y(ip) = ip/(1.0*nparts) + my_rank/(1.0*n_cpu) + ip
        z(ip) = ip/(1.0*nparts) + my_rank/(1.0*n_cpu) + 2
        q(ip) = 0.13
        m(ip) = 0.13
     end do

     call pepc(nparts, nparts_total, x, y, z, q, m, ex, ey, ez, pot, lx, ly, lz, lperiod, extrcorr)

     deallocate(x, y, z, ex, ey, ez)
     deallocate(pot, q, m)

  end do

  call MPI_FINALIZE(ierr)

end program pepcs
