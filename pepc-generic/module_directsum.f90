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

      type :: direct_particle
         sequence
         real*8, dimension(3) :: r
         real*8, dimension(3) :: field
         real*8 :: potential
      end type direct_particle

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
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine verifydirect(x, y, z, q, ex, ey, ez, pot, np_local, testidx, cf_par, my_rank, n_cpu, comm)
          use treetypes
          implicit none
          include 'mpif.h'

          real*8, dimension(:), intent(in) :: x !< x-coordinates of particles, index 1..np_local
          real*8, dimension(:), intent(in) :: y !< y-coordinates of particles, index 1..np_local
          real*8, dimension(:), intent(in) :: z !< z-coordinates of particles, index 1..np_local
          real*8, dimension(:), intent(in) :: q !< charges of particles, index 1..np_local
          real*8, dimension(:), intent(in) :: ex !< x-coordinates of particles, index 1..np_local
          real*8, dimension(:), intent(in) :: ey !< y-coordinates of particles, index 1..np_local
          real*8, dimension(:), intent(in) :: ez !< z-coordinates of particles, index 1..np_local
          real*8, dimension(:), intent(in) :: pot !< charges of particles, index 1..np_local
          type(calc_force_params), intent(in) :: cf_par !< force calculation parameters
          integer, intent(in) :: np_local !< number of local particles
          integer, dimension(:), intent(in) :: testidx !< field with particle indices that direct force has to be computed for
          integer :: ntest !< number of particles in testidx
          integer, intent(in) :: my_rank, n_cpu, comm

          integer :: i
          type(direct_particle), dimension(:), allocatable :: res !< test results

          ntest = size(testidx)

          if (my_rank ==0) write(*,'("-- DIRECT VERIFICATION --")')

          call directforce(x, y, z, q, np_local, testidx, ntest, cf_par, res, my_rank, n_cpu, comm)

          do i=1,ntest
            associate(p=>testidx(i), re=>res(i))
              write(*,'("[",I6.6,":",I6.6,"]",3(x,F10.4), " | PEPC   ", 4(x,E20.12), /, "[",I6.6,":",I6.6,"]",33x, " | DIRECT ", 4(x,E20.12))') &
                    my_rank, p, x(p), y(p), z(p), ex(p), ey(p), ez(p), pot(p), my_rank, p, re%field, re%potential
            end associate
          end do

          deallocate(res)

        end subroutine

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine directforce(x, y, z, q, np_local, testidx, ntest, cf_par, directresults, my_rank, n_cpu, comm)
          use treetypes
          implicit none
          include 'mpif.h'

          real*8, dimension(:), intent(in) :: x !< x-coordinates of particles, index 1..np_local
          real*8, dimension(:), intent(in) :: y !< y-coordinates of particles, index 1..np_local
          real*8, dimension(:), intent(in) :: z !< z-coordinates of particles, index 1..np_local
          real*8, dimension(:), intent(in) :: q !< charges of particles, index 1..np_local
          type(calc_force_params), intent(in) :: cf_par !< force calculation parameters
          integer, intent(in) :: np_local !< number of local particles
          integer, dimension(:), intent(in) :: testidx !< field with particle indices that direct force has to be computed for
          integer, intent(in) :: ntest !< number of particles in testidx
          type(direct_particle), dimension(:), allocatable, intent(out) :: directresults !< test results
          integer, intent(in) :: my_rank, n_cpu, comm
          real :: eps2

          integer :: maxtest !< maximum ntest
          type(direct_particle), dimension(:), allocatable :: received, sending
          integer :: nreceived, nsending
          integer :: ierr, req, stat(MPI_STATUS_SIZE), i, j, currank, nextrank, prevrank

          call MPI_ALLREDUCE(ntest, maxtest, 1, MPI_INTEGER, MPI_MAX, comm, ierr)
          allocate(received(1:maxtest), sending(1:maxtest))

          ! determine right and left neighbour
          nextrank = modulo(my_rank + 1, n_cpu)
          prevrank = modulo(my_rank - 1 + n_cpu, n_cpu)

          eps2 = cf_par%eps**2

          ! insert initial data into input array
          nreceived = ntest
          do i=1,ntest
            associate (p => testidx(i))
              received(i) = direct_particle([x(p), y(p), z(p)], [0.,   0.,   0.  ], 0.)
            end associate
          end do

          ! we will send our data packet to every other mpi rank
          do currank=0,n_cpu-1

            ! calculate force from local particles i onto particles j in received-buffer
            ! loop over all received particles
            do j=1,nreceived
              ! loop over all local particles
              do i=1,np_local
                if ((currank .ne. 0) .or. (testidx(j).ne.i)) then
                  call calc_direct_coulomb_force(received(j), [x(i), y(i), z(i)], q(i), eps2)
                endif
              end do
            end do

            ! copy particles to process to send-buffer
            nsending = nreceived
            sending(1:nsending) = received(1:nreceived)
            ! send size of current data package to right neighbour, receive size of new data package from left neighbour
            call MPI_ISEND(nsending, 1, MPI_INTEGER, nextrank, 47, comm,  req, ierr)
            call MPI_RECV(nreceived, 1, MPI_INTEGER, prevrank, 47, comm, stat, ierr)
            call MPI_WAIT(req, stat, ierr)
            ! send current data package to right neighbour, receive new data package from left neighbour
            call MPI_ISEND(sending, 7*nsending,  MPI_REAL8, nextrank, 48, comm,  req, ierr)
            call MPI_RECV(received, 7*nreceived, MPI_REAL8, prevrank, 48, comm, stat, ierr)
            call MPI_WAIT(req, stat, ierr)
          end do

          ! copy results to output array
          allocate(directresults(1:ntest))
          directresults(1:ntest) = received(1:nreceived)

          do i=1,ntest
            directresults(i)%field     = directresults(i)%field     * cf_par%force_const
            directresults(i)%potential = directresults(i)%potential * cf_par%force_const
          end do

          deallocate(received, sending)

        end subroutine



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine calc_direct_coulomb_force(p1, r2, q2, eps2)
          implicit none
          real*8, intent(in) :: r2(3), q2
          real, intent(in) :: eps2
          type(direct_particle), intent(inout) :: p1

          real*8 :: delta(3), dr, dr2, dr3

          delta = p1%r - r2
          dr2   = dot_product(delta, delta) + eps2
          dr    = sqrt(dr2)
          dr3   = dr*dr2

          ! electric field
          p1%field     = p1%field     + q2 * delta/dr3
          !electric potential
          p1%potential = p1%potential + q2 * 1./dr

        end subroutine

end module module_directsum