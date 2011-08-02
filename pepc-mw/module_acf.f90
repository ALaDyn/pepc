!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>  Encapsulates ...
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_acf
      implicit none


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      type acf
        private
          integer :: Ntau
          integer :: tau
          integer :: num_pe, my_rank, comm
          real*8, allocatable :: Kt(:)
          integer, allocatable :: ctr(:)
          real*8,  allocatable :: oldvals(:,:)

        contains
          procedure :: initialize => acf_initialize
          procedure :: finalize => acf_finalize
          procedure :: addval => acf_addval
          procedure :: to_file => acf_to_file

      end type acf


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      contains

      subroutine acf_to_file(acf_, filename)
        implicit none
        class(acf) :: acf_
        character(*) :: filename

        if (acf_%my_rank == 0) then
          open(47,file=trim(filename))
          write(47,'(g18.8)') acf_%Kt / acf_%ctr
          close(47)
        endif

      end subroutine


      subroutine acf_initialize(acf_, Nt_, my_rank_, num_pe_, comm_)
        implicit none
        class(acf) :: acf_
        integer, intent(in) :: Nt_
        integer, intent(in) :: num_pe_, my_rank_, comm_

        acf_%Ntau    = Nt_
        acf_%tau     = 0
        acf_%num_pe  = num_pe_
        acf_%my_rank = my_rank_
        acf_%comm    = comm_

        allocate(acf_%Kt(0:acf_%Ntau))
        acf_%Kt = 0.
        allocate(acf_%ctr(0:acf_%Ntau))
        acf_%ctr = 0
        allocate(acf_%oldvals(1:3,1:acf_%Ntau))
      end subroutine


      subroutine acf_finalize(acf_)
        implicit none
        class(acf) :: acf_

        if (allocated(acf_%Kt))      deallocate(acf_%Kt)
        if (allocated(acf_%ctr))     deallocate(acf_%ctr)
        if (allocated(acf_%oldvals)) deallocate(acf_%oldvals)
      end subroutine



      subroutine acf_addval(acf_, val)
        implicit none
        include 'mpif.h'
        class(acf) :: acf_
        real*8, intent(in) :: val(3)
        integer :: s, mystart, myend, ierr
        integer, dimension(:), allocatable :: recvcounts, displs

        allocate(recvcounts(0:acf_%num_pe-1))
        allocate(displs(0:acf_%num_pe-1))

        acf_%tau = acf_%tau + 1

        acf_%oldvals(1:3,acf_%tau) = val

        ! split loop
        !    do s = 0,acf_%tau-1
        ! among pocessors
        mystart = 0
        myend   = acf_%tau / acf_%num_pe
        if (acf_%my_rank < modulo(acf_%tau, acf_%num_pe)) then
          myend = myend + 1
        endif

        recvcounts(acf_%my_rank) = myend
        call MPI_ALLGATHER(MPI_IN_PLACE, 1, MPI_INTEGER, recvcounts, 1, MPI_INTEGER, acf_%comm, ierr)
        displs(0) = 0
        do s=1,acf_%num_pe-1
          displs(s) = displs(s-1) + recvcounts(s-1)
        end do

        mystart = displs(acf_%my_rank)
        myend   = myend + mystart

        do s = mystart,myend-1
          acf_%Kt(s) = acf_%Kt(s) + dot_product(val, acf_%oldvals(1:3,acf_%tau - s))
          acf_%ctr(s) = acf_%ctr(s) + 1
        end do


        call MPI_ALLGATHERV(MPI_IN_PLACE, recvcounts(acf_%my_rank), MPI_REAL8,   acf_%Kt,  recvcounts, displs, MPI_REAL8,   acf_%comm, ierr)
        call MPI_ALLGATHERV(MPI_IN_PLACE, recvcounts(acf_%my_rank), MPI_INTEGER, acf_%ctr, recvcounts, displs, MPI_INTEGER, acf_%comm, ierr)

        deallocate(recvcounts, displs)

      end subroutine


end module module_acf
