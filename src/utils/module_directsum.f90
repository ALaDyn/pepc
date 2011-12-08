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
        !>
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine verifydirect(x, y, z, q, ex, ey, ez, pot, np_local, testidx, verbosity, my_rank, n_cpu, comm)
          use module_pepc_types
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
          integer, intent(in) :: verbosity !< verbosity level: 0 - only print max. relative deviations, 1 - additionally print all. relative deviations, 2 - additionally print all. calculated forces
          integer, intent(in) :: np_local !< number of local particles
          integer, dimension(:), intent(in) :: testidx !< field with particle indices that direct force has to be computed for
          integer :: ntest !< number of particles in testidx
          integer, intent(in) :: my_rank, n_cpu, comm
          real*8 :: deviation(4), deviation_max(4)
          real*8 :: field_abssum(4), field_average(4)

          integer :: i, ntest_total, ierr
          type(direct_particle), target, dimension(:), allocatable :: res !< test results
          type(direct_particle), pointer :: re
          integer :: p

          ntest = size(testidx)

          if (my_rank ==0) write(*,'("-- DIRECT VERIFICATION --")')

          call directforce(x, y, z, q, np_local, testidx, ntest, res, my_rank, n_cpu, comm)

          deviation     = 0.
          deviation_max = 0.
          field_abssum  = 0.

          do i=1,ntest
            p = testidx(i)
            re=>res(i)
              deviation(1:3)       = abs( re%field - [ex(p), ey(p), ez(p)] )
              field_abssum(1:3)    = field_abssum(1:3)    + abs(re%field)
              deviation(4)         = abs( re%potential - pot(p) )
              field_abssum(4)      = field_abssum(4)      + abs(re%potential)

              deviation_max  = max(deviation, deviation_max)

              if (verbosity > 1) then
                write(*,'("[",I6.6,":",I6.6,"]",3(x,F10.4), " | PEPC    ", 4(x,E20.13))') my_rank, p, x(p), y(p), z(p), ex(p), ey(p), ez(p), pot(p)
                write(*,'("[",I6.6,":",I6.6,"]",33x,        " | DIRECT  ", 4(x,E20.13))') my_rank, p, re%field, re%potential
              endif

              if (verbosity > 0) then
                write(*,'("[",I6.6,":",I6.6,"]",33x,        " | Abs.err ", 4(x,E20.13),"__")') my_rank, p, deviation
              endif
          end do

          call MPI_REDUCE(deviation_max,   deviation,             4, MPI_REAL8,   MPI_MAX, 0, comm, ierr)
          call MPI_REDUCE(field_abssum,    field_average,         4, MPI_REAL8,   MPI_SUM, 0, comm, ierr)
          call MPI_REDUCE(ntest,           ntest_total,           1, MPI_INTEGER, MPI_SUM, 0, comm, ierr)
          field_average         = field_average         / ntest_total

          if ((verbosity > -1) .and. (my_rank == 0)) then
            write(*,'("Maximum absolute deviation (ex, ey, ez, pot): ", 4(2x,E20.12))') deviation
            write(*,'("Average field values       (ex, ey, ez, pot): ", 4(2x,E20.12))') field_average
            write(*,'("Maximum relative deviation (ex, ey, ez, pot): ", 4(2x,F20.3))') deviation / field_average
          endif

          deallocate(res)

        end subroutine

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> direct computation of coulomb force onto a selection of local particles
        !> due to contributions of all (also remote) other particles
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine directforce(x, y, z, q, np_local, testidx, ntest, directresults, my_rank, n_cpu, comm)
          use module_pepc_types
          use module_calc_force, only : eps2
          implicit none
          include 'mpif.h'

          real*8, dimension(:), intent(in) :: x !< x-coordinates of particles, index 1..np_local
          real*8, dimension(:), intent(in) :: y !< y-coordinates of particles, index 1..np_local
          real*8, dimension(:), intent(in) :: z !< z-coordinates of particles, index 1..np_local
          real*8, dimension(:), intent(in) :: q !< charges of particles, index 1..np_local
          integer, intent(in) :: np_local !< number of local particles
          integer, dimension(:), intent(in) :: testidx !< field with particle indices that direct force has to be computed for
          integer, intent(in) :: ntest !< number of particles in testidx
          type(direct_particle), dimension(:), allocatable, intent(out) :: directresults !< test results
          integer, intent(in) :: my_rank, n_cpu, comm

          integer :: maxtest !< maximum ntest
          type(direct_particle), dimension(:), allocatable :: received, sending
          integer :: nreceived, nsending
          integer :: ierr, req, stat(MPI_STATUS_SIZE), i, j, currank, nextrank, prevrank, p

          call MPI_ALLREDUCE(ntest, maxtest, 1, MPI_INTEGER, MPI_MAX, comm, ierr)
          allocate(received(1:maxtest), sending(1:maxtest))

          ! determine right and left neighbour
          nextrank = modulo(my_rank + 1, n_cpu)
          prevrank = modulo(my_rank - 1 + n_cpu, n_cpu)

          ! insert initial data into input array
          nreceived = ntest
          do i=1,ntest
            p = testidx(i)
              received(i) = direct_particle([x(p), y(p), z(p)], [0.,   0.,   0.  ], 0.)
          end do

          ! we will send our data packet to every other mpi rank
          do currank=0,n_cpu-1

            ! calculate force from local particles i onto particles j in received-buffer
            ! loop over all received particles
            do j=1,nreceived
              ! loop over all local particles
              do i=1,np_local
                if ((currank .ne. 0) .or. (testidx(j).ne.i)) then
                  call calc_direct_coulomb_force_3D(received(j), [x(i), y(i), z(i)], q(i), eps2)
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
            call MPI_ISEND(sending, 7*nsending,  MPI_REAL8, nextrank, MPI_TAG_DIRECT_DATA_PACKAGE, comm,  req, ierr)
            call MPI_RECV(received, 7*nreceived, MPI_REAL8, prevrank, MPI_TAG_DIRECT_DATA_PACKAGE, comm, stat, ierr)
            call MPI_WAIT(req, stat, ierr)
          end do

          ! copy results to output array
          allocate(directresults(1:ntest))
          directresults(1:ntest) = received(1:nreceived)

          do i=1,ntest
            directresults(i)%field     = directresults(i)%field
            directresults(i)%potential = directresults(i)%potential
          end do

          deallocate(received, sending)

        end subroutine



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Coulomb force law (3D-case) for direct summation
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine calc_direct_coulomb_force_3D(p1, r2, q2, eps2)
          implicit none
          real*8, intent(in) :: r2(3), q2
          real*8, intent(in) :: eps2
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
