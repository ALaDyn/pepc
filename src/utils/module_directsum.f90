! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
!
! Copyright (C) 2002-2013 Juelich Supercomputing Centre,
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

module module_directsum

      implicit none

      integer, parameter :: MPI_TAG_DIRECT_DATA_PACKAGE_SIZE = 47
      integer, parameter :: MPI_TAG_DIRECT_DATA_PACKAGE      = 48

      contains

        !>
        !> direct computation of coulomb force onto a selection of local particles
        !> due to contributions of all (also remote) other particles
        !>
        subroutine directforce(particles, testidx, ntest, directresults, comm)
          use module_pepc_types
          use module_interaction_specific_types
          use module_interaction_specific
          use omp_lib
          use treevars, only: num_threads
          use module_timings
          use module_mirror_boxes
          implicit none
          include 'mpif.h'

          type(t_particle), intent(in) :: particles(:)
          integer(kind_particle), dimension(:), intent(in) :: testidx !< field with particle indices that direct force has to be computed for
          integer(kind_particle), intent(in) :: ntest !< number of particles in testidx
          type(t_particle_results), dimension(:), allocatable, intent(out) :: directresults !< test results
          integer, intent(in) :: comm

          integer(kind_particle) :: maxtest !< maximum ntest
          type(t_particle), dimension(:), allocatable :: received, sending
          integer(kind_particle) :: i, j, tile_start, tile_size, nreceived, nsending
          integer :: ierr, stat(MPI_STATUS_SIZE), thread_id, num_threads_
          integer(kind_pe) :: my_rank, n_cpu, currank, nextrank, prevrank
          type(t_tree_node_interaction_data), allocatable :: local_nodes(:)
          real*8, allocatable :: delta(:,:), dist2(:)
          integer :: ibox, id
          type(t_particle) :: latticeparticles(ntest)
          type(t_particle_pack), allocatable :: particle_pack(:)

          real*8 :: t1, vbox(3)
          integer :: omp_thread_num

          call MPI_COMM_RANK(comm, my_rank, ierr)
          call MPI_COMM_SIZE(comm, n_cpu, ierr)

          call MPI_ALLREDUCE(ntest, maxtest, 1, MPI_KIND_PARTICLE, MPI_MAX, comm, ierr)
          allocate(received(1:maxtest), sending(1:maxtest))

          call timer_reset(t_direct_force)
          call timer_reset(t_direct_comm)

          ! Set number of openmp threads to the same number as pthreads used in the walk
          !$ call omp_set_num_threads(num_threads)

          ! Inform the user that openmp is used, and with how many threads
          !$OMP PARALLEL PRIVATE(omp_thread_num)
          !$ omp_thread_num = OMP_GET_THREAD_NUM()
          !;$ if( (my_rank .eq. 0) .and. (omp_thread_num .eq. 0) ) write(*,*) 'Using OpenMP with', OMP_GET_NUM_THREADS(), 'threads. Adjust by modifying num_threads parameter.'
          !$OMP END PARALLEL

          ! determine right and left neighbour
          nextrank = modulo(my_rank + 1_kind_pe, n_cpu)
          prevrank = modulo(my_rank - 1_kind_pe + n_cpu, n_cpu)

          ! insert initial data into input array - these particles will be shipped around later
          nreceived = ntest
          do i=1,ntest
            received(i) = particles(testidx(i))
          end do

          call particleresults_clear(received)

          ! we copy all local particles into a node array to be able to feed them to calc_force_per_interaction
          allocate(local_nodes(size(particles)))

          do i=1,size(particles)
            call multipole_from_particle(particles(i)%x, particles(i)%data, local_nodes(i))
          end do

          ! we will send our data packet to every other mpi rank
          do currank=0_kind_pe,n_cpu-1_kind_pe
            ! calculate force from local particles i onto particles j in received-buffer
            ! loop over all received particles

            t1 = MPI_WTIME()

            !$OMP PARALLEL default(shared) private(vbox, dist2, delta, i, ibox, id, j, tile_start, tile_size, thread_id, num_threads_)
            thread_id = omp_get_thread_num()
            num_threads_ = omp_get_num_threads()

            !$omp master
            allocate(particle_pack(0:num_threads_ - 1))
            !$omp end master
            !$omp barrier

            tile_size = nreceived / num_threads_
            if (thread_id < mod(nreceived, num_threads_)) tile_size = tile_size + 1
            tile_start = 1 + (nreceived / num_threads_ + 1) * min(thread_id, mod(nreceived, num_threads_)) &
              + (nreceived / num_threads_) * max(0, thread_id - mod(nreceived, num_threads_))

            allocate(delta(tile_size,3), dist2(tile_size))
            call pack_particle_list(received(tile_start:tile_start + tile_size - 1), particle_pack(thread_id))

            do i = 1, size(local_nodes)
              do ibox = 1, num_neighbour_boxes ! sum over all boxes within ws=1
                vbox = lattice_vect(neighbour_boxes(:,ibox))
                dist2 = 0.0_8
                do id = 1, 3
                  do j = 1, tile_size
                    delta(j, id) = received(tile_start + j - 1)%x(id) - vbox(id) - local_nodes(i)%coc(id)
                    dist2(j) = dist2(j) + delta(j, id) * delta(j, id)
                  end do
                end do

                call calc_force_per_interaction_with_leaf(delta, dist2, particle_pack(thread_id), local_nodes(i))
              end do
            end do

            deallocate(delta, dist2)
            call unpack_particle_list(particle_pack(thread_id), received(tile_start:tile_start + tile_size - 1))

            !$omp barrier
            !$omp master
            deallocate(particle_pack)
            !$omp end master
            !$OMP END PARALLEL

            call timer_add(t_direct_force,MPI_WTIME()-t1)

            t1 = MPI_WTIME()

            ! copy particles to process to send-buffer
            nsending = nreceived
            sending(1:nsending) = received(1:nreceived)
            ! send size of current data package to right neighbour, receive size of new data package from left neighbour
            call MPI_SENDRECV(nsending, 1, MPI_INTEGER, nextrank, MPI_TAG_DIRECT_DATA_PACKAGE_SIZE, &
              nreceived, 1, MPI_INTEGER, prevrank, MPI_TAG_DIRECT_DATA_PACKAGE_SIZE, &
              comm, stat, ierr)
            ! send current data package to right neighbour, receive new data package from left neighbour
            call MPI_SENDRECV(sending, nsending, MPI_TYPE_PARTICLE, nextrank, MPI_TAG_DIRECT_DATA_PACKAGE, &
              received, nreceived, MPI_TYPE_PARTICLE, prevrank, MPI_TAG_DIRECT_DATA_PACKAGE, &
              comm, stat, ierr)

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
          call calc_force_after_grow(particles)
          latticeparticles(1:ntest)         = particles(testidx)
          latticeparticles(1:ntest)%results = directresults(1:ntest)
          call calc_force_per_particle(latticeparticles(1:ntest))
          directresults(1:ntest) = latticeparticles(1:ntest)%results
          call timer_stop(t_lattice)
        end subroutine
end module module_directsum
