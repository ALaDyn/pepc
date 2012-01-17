!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_directsum

      use omp_lib
      implicit none

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


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  private type declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer, parameter :: MPI_TAG_DIRECT_DATA_PACKAGE_SIZE = 47
      integer, parameter :: MPI_TAG_DIRECT_DATA_PACKAGE      = 48


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  subroutine-implementation  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      contains

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> direct computation of coulomb force onto a selection of local particles
        !> due to contributions of all (also remote) other particles
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine directforce(particles, np_local, testidx, ntest, directresults, my_rank, n_cpu, comm)
          use module_pepc_types
          use module_interaction_specific_types
          use module_interaction_specific
          implicit none
          include 'mpif.h'

          type(t_particle), intent(in) :: particles(1:np_local)
          integer, intent(in) :: np_local !< number of local particles
          integer, dimension(:), intent(in) :: testidx !< field with particle indices that direct force has to be computed for
          integer, intent(in) :: ntest !< number of particles in testidx
          type(t_particle_results), dimension(:), allocatable, intent(out) :: directresults !< test results
          integer, intent(in) :: my_rank, n_cpu, comm

          integer :: maxtest !< maximum ntest
          type(t_particle), dimension(:), allocatable :: received, sending
          integer :: nreceived, nsending
          integer :: ierr, req, stat(MPI_STATUS_SIZE), i, j, currank, nextrank, prevrank
          type(t_tree_node_interaction_data), allocatable :: local_nodes(:)
          real*8 :: delta(3)

          call MPI_ALLREDUCE(ntest, maxtest, 1, MPI_INTEGER, MPI_MAX, comm, ierr)
          allocate(received(1:maxtest), sending(1:maxtest))

          ! determine right and left neighbour
          nextrank = modulo(my_rank + 1, n_cpu)
          prevrank = modulo(my_rank - 1 + n_cpu, n_cpu)

          ! insert initial data into input array - these particles will be shipped around later
          nreceived = ntest
          do i=1,ntest
            received(i) = particles(testidx(i))
          end do
           
          call particleresults_clear(received, ntest)

          ! we copy all local particles into a node array to be able to feed them to calc_force_per_interaction
          allocate(local_nodes(np_local))
          do i=1,np_local
            call multipole_from_particle(particles(i)%x, particles(i)%data, local_nodes(i))
          end do

! TODO: loop over vbox-vectors
          ! we will send our data packet to every other mpi rank
          do currank=0,n_cpu-1

            ! calculate force from local particles i onto particles j in received-buffer
            ! loop over all received particles
            do j=1,nreceived
              ! loop over all local particles
              do i=1,np_local
                if ((currank .ne. 0) .or. (testidx(j).ne.i)) then
                  delta = received(j)%x - local_nodes(i)%coc
                  call calc_force_per_interaction(received(j), local_nodes(i), particles(i)%key, delta, dot_product(delta, delta), [0._8, 0._8, 0._8], .true.)
                endif
              end do
            end do

            ! copy particles to process to send-buffer
            nsending = nreceived
            sending(1:nsending) = received(1:nreceived)
            ! send size of current data package to right neighbour, receive size of new data package from left neighbour
            call MPI_ISEND(nsending, 1, MPI_INTEGER, nextrank, MPI_TAG_DIRECT_DATA_PACKAGE_SIZE, comm,  req, ierr)
            call MPI_RECV(nreceived, 1, MPI_INTEGER, prevrank, MPI_TAG_DIRECT_DATA_PACKAGE_SIZE, comm, stat, ierr)
            call MPI_WAIT(req, stat, ierr)
            ! send current data package to right neighbour, receive new data package from left neighbour
            call MPI_ISEND(sending, nsending,  MPI_TYPE_PARTICLE, nextrank, MPI_TAG_DIRECT_DATA_PACKAGE, comm,  req, ierr)
            call MPI_RECV(received, nreceived, MPI_TYPE_PARTICLE, prevrank, MPI_TAG_DIRECT_DATA_PACKAGE, comm, stat, ierr)
            call MPI_WAIT(req, stat, ierr)
          end do

          ! copy results to output array
          allocate(directresults(1:ntest))
          directresults(1:ntest) = received(1:nreceived)%results

          deallocate(received, sending, local_nodes)

        end subroutine

        subroutine omp_init(my_rank)

            use module_walk, only: num_walk_threads
            implicit none

            integer, intent(in) :: my_rank
            integer :: omp_thread_num

            ! Set the number of openmp threads.
            ! Set this only, when compiling with openmp (with !$)
            ! Set number of openmp threads to the same number as pthreads used in the walk
            !$ call omp_set_num_threads(num_walk_threads)

            ! Inform the user that openmp is used, and with how many threads
            !$OMP PARALLEL PRIVATE(omp_thread_num)
            !$ omp_thread_num = OMP_GET_THREAD_NUM()
            !$ if( (my_rank .eq. 0) .and. (omp_thread_num .eq. 0) ) write(*,*) 'Using OpenMP with', OMP_GET_NUM_THREADS(), 'threads.'
            !$OMP END PARALLEL

        end subroutine omp_init

        subroutine omp_clear()

            implicit none
            ! Reset the number of openmp threads to 1.
            ! Set this only, when compiling with openmp (with !$)
            !$ call omp_set_num_threads(1)
        end subroutine omp_clear

end module module_directsum
