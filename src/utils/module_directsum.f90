! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2012 Juelich Supercomputing Centre, 
!                         Forschungszentrum Juelich GmbH,
!                         Germany
! 
! PEPC is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! PEPC is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
! 
! You should have received a copy of the GNU Lesser General Public License
! along with PEPC.  If not, see <http://www.gnu.org/licenses/>.
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_directsum

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
          use omp_lib
          use module_walk, only: num_walk_threads
          use module_timings
          use module_mirror_boxes
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
          integer :: ibox
          type(t_particle) :: latticeparticles(ntest)

          real*8 :: t1
          integer :: omp_thread_num


          call MPI_ALLREDUCE(ntest, maxtest, 1, MPI_INTEGER, MPI_MAX, comm, ierr)
          allocate(received(1:maxtest), sending(1:maxtest))

          call timer_reset(t_direct_force)
          call timer_reset(t_direct_comm)

          ! Set number of openmp threads to the same number as pthreads used in the walk
          !$ call omp_set_num_threads(num_walk_threads)

          ! Inform the user that openmp is used, and with how many threads
          !$OMP PARALLEL PRIVATE(omp_thread_num)
          !$ omp_thread_num = OMP_GET_THREAD_NUM()
          !$ if( (my_rank .eq. 0) .and. (omp_thread_num .eq. 0) ) write(*,*) 'Using OpenMP with', OMP_GET_NUM_THREADS(), 'threads. Adjust by modifying num_walk_threads parameter.'
          !$OMP END PARALLEL

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

          ! we will send our data packet to every other mpi rank
          do currank=0,n_cpu-1

            ! calculate force from local particles i onto particles j in received-buffer
            ! loop over all received particles

            t1 = MPI_WTIME()

            !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(j, i, delta)
            do j=1,nreceived
                do i=1,np_local

                    do ibox = 1,num_neighbour_boxes ! sum over all boxes within ws=1
                      ! if we use our own particles, test for equality; exclude particle itself if we are in central box
                      if ((currank .ne. 0) .or. (ibox < num_neighbour_boxes) .or. (testidx(j).ne.i)) then
                          delta = received(j)%x - (local_nodes(i)%coc - lattice_vect(neighbour_boxes(:,ibox)))

                          if (all(abs(delta) < spatial_interaction_cutoff)) then
                              call calc_force_per_interaction(received(j), local_nodes(i), particles(i)%key, delta, dot_product(delta, delta), [0._8, 0._8, 0._8], .true.)
                          endif

                      endif
                    end do

                end do
            end do
            !$OMP END PARALLEL DO

            call timer_add(t_direct_force,MPI_WTIME()-t1)

            t1 = MPI_WTIME()

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

            call timer_add(t_direct_comm,MPI_WTIME()-t1)

          end do

          ! copy results to output array
          allocate(directresults(1:ntest))
          directresults(1:ntest) = received(1:nreceived)%results

          deallocate(received, sending, local_nodes)

          ! Reset the number of openmp threads to 1.
          !$ call omp_set_num_threads(1)

          !call calc_force_per_particle here: add lattice contribution, compare module_libpepc_main
          call timer_start(t_lattice)
          ! add lattice contribution and other per-particle-forces
          call calc_force_after_grow(particles, np_local)
          latticeparticles(1:ntest)         = particles(testidx)
          latticeparticles(1:ntest)%results = directresults(1:ntest)
          call calc_force_per_particle(latticeparticles, ntest)
          directresults(1:ntest) = latticeparticles(1:ntest)%results
          call timer_stop(t_lattice)

        end subroutine

end module module_directsum
