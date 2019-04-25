! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
!
! Copyright (C) 2002-2017 Juelich Supercomputing Centre,
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
      use module_pepc_kinds
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
          use module_pepc_kinds
          use module_interaction_specific_types
          use module_interaction_specific
          use omp_lib
          use treevars, only: num_threads
          use module_timings
          use module_mirror_boxes
          use module_debug
          implicit none
          include 'mpif.h'

          type(t_particle), intent(in) :: particles(:)
          integer(kind_particle), dimension(:), intent(in) :: testidx !< field with particle indices that direct force has to be computed for
          integer(kind_particle), intent(in) :: ntest !< number of particles in testidx
          type(t_particle_results), dimension(:), intent(out) :: directresults !< test results
          integer, intent(in) :: comm

          integer(kind_particle) :: maxtest !< maximum ntest
          type(t_particle), dimension(:), allocatable :: received, sending
          integer(kind_default) :: nreceived, nsending
          integer(kind_particle) :: i, j
          integer :: ierr, stat(MPI_STATUS_SIZE)
          integer(kind_pe) :: my_rank, n_cpu, currank, nextrank, prevrank
          type(t_tree_node_interaction_data), allocatable :: local_nodes(:)
          real(kind_physics) :: delta(3)
          integer :: ibox
          type(t_particle) :: latticeparticles(ntest)

          real*8 :: t1

          call MPI_COMM_RANK(comm, my_rank, ierr)
          call MPI_COMM_SIZE(comm, n_cpu, ierr)

          call MPI_ALLREDUCE(ntest, maxtest, 1, MPI_KIND_PARTICLE, MPI_MAX, comm, ierr)
          allocate(received(1:maxtest), sending(1:maxtest))

          call timer_reset(t_direct_force)
          call timer_reset(t_direct_comm)

          ! Inform the user that openmp is used, and with how many threads
          !$OMP PARALLEL DEFAULT(NONE) SHARED(my_rank) NUM_THREADS(num_threads)
          !$OMP MASTER
          !$ if (my_rank .eq. 0) write(*,*) 'Using OpenMP with', OMP_GET_NUM_THREADS(), 'threads. Adjust by modifying num_threads parameter.'
          !$OMP END MASTER
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

            !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, delta, ibox) SCHEDULE(STATIC) NUM_THREADS(num_threads)
            do j=1,nreceived
                do i=1,size(particles)

                    do ibox = 1,num_neighbour_boxes ! sum over all boxes within ws=1
                      ! if we use our own particles, test for equality; exclude particle itself if we are in central box
                      delta = received(j)%x - lattice_vect(neighbour_boxes(:,ibox)) - local_nodes(i)%coc
                      #ifndef NO_SPATIAL_INTERACTION_CUTOFF
                      if (all(abs(delta) < spatial_interaction_cutoff)) then
                      #endif
                        if (currank == 0) then
                          if ((ibox == num_neighbour_boxes) .and. (testidx(j) == i)) then
                            call calc_force_per_interaction(received(j), local_nodes(i), particles(i)%key, delta, dot_product(delta, delta), lattice_vect(neighbour_boxes(:,ibox)), .false.)
                            cycle
                          end if
                        end if

                        call calc_force_per_interaction(received(j), local_nodes(i), particles(i)%key, delta, dot_product(delta, delta), lattice_vect(neighbour_boxes(:,ibox)), .true.)
                      #ifndef NO_SPATIAL_INTERACTION_CUTOFF
                      endif
                      #endif
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
            call MPI_SENDRECV(nsending, 1, MPI_INTEGER, nextrank, MPI_TAG_DIRECT_DATA_PACKAGE_SIZE, &
              nreceived, 1, MPI_INTEGER, prevrank, MPI_TAG_DIRECT_DATA_PACKAGE_SIZE, &
              comm, stat, ierr)
            ! send current data package to right neighbour, receive new data package from left neighbour
            call MPI_SENDRECV(sending, nsending, MPI_TYPE_particle_vec, nextrank, MPI_TAG_DIRECT_DATA_PACKAGE, &
              received, nreceived, MPI_TYPE_particle_vec, prevrank, MPI_TAG_DIRECT_DATA_PACKAGE, &
              comm, stat, ierr)

            call timer_add(t_direct_comm,MPI_WTIME()-t1)

          end do

          ! copy results to output array
          DEBUG_ASSERT(ntest == nreceived)
          DEBUG_ASSERT(size(directresults, kind = kind_particle) >= ntest)
          directresults(1:ntest) = received(1:nreceived)%results

          deallocate(received, sending, local_nodes)

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
