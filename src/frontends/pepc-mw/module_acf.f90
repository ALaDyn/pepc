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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>  Provides a class for computation of the auto correlation function
!> of a vectorial quantity
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
          real*8  :: dt
          integer :: tau
          integer :: num_pe, my_rank, comm
          real*8, allocatable :: Kt(:)
          real*8,  allocatable :: oldvals(:,:)
          integer, allocatable :: num_contributions(:)

        contains
          procedure :: initialize => acf_initialize
          procedure :: finalize => acf_finalize
          procedure :: addval => acf_addval
          procedure :: to_file => acf_to_file
          procedure :: from_file => acf_from_file

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

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !> Writes the values of the ACF(delta_t) to a file
      !> first colum:   delta_t in physical time (see parameter dt of acf_initialize() )
      !> second column: ACF(delta_t)
      !> @param filename name of target file
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine acf_to_file(acf_, filename)
        implicit none
        class(acf) :: acf_
        character(*) :: filename
        integer :: i

        if (acf_%my_rank == 0) then
          open(47,file=trim(filename),position='rewind')
          do i=0,acf_%tau-1
            write(47,'(6(g18.8,x),I8)') acf_%dt*i, acf_%Kt(i) / acf_%num_contributions(i), acf_%Kt(i), acf_%oldvals(1:3,i+1), acf_%num_contributions(i)
          end do
          close(47)
        endif

      end subroutine

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !> Reads the values of the ACF(delta_t) from a file
      !> first colum:   delta_t in physical time (see parameter dt of acf_initialize() )
      !> second column: ACF(delta_t)
      !> @param filename name of target file
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine acf_from_file(acf_, filename)
        implicit none
        class(acf) :: acf_
        character(*) :: filename
        integer :: i
        real*8 :: tmp(2)

        acf_%tau = 0

        if (acf_%my_rank == 0) then
          open(47,file=trim(filename),action='read')
          do i=0,acf_%Ntau-1
            acf_%tau = acf_%tau + 1
            read(47,'(6(g18.8,x),I8)', end=300) tmp(1), tmp(2), acf_%Kt(acf_%tau-1), acf_%oldvals(1:3,acf_%tau), acf_%num_contributions(acf_%tau-1)
          end do
          300 close(47)
        endif

      end subroutine


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !> Initializes the ACF class
      !> @param Nt_ total number of timesteps to be processed
      !> @param dt_ physical timestep (just for output purposes)
      !> @param my_rank_ local MPI rank
      !> @param num_pe_ total number of MPI ranks
      !> @param comm_ MPI communicator (e.g. MPI_COMM_WORLD)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine acf_initialize(acf_, Nt_, dt_, my_rank_, num_pe_, comm_)
        implicit none
        class(acf) :: acf_
        integer, intent(in) :: Nt_
        real*8, intent(in) :: dt_
        integer, intent(in) :: num_pe_, my_rank_, comm_

        acf_%Ntau    = Nt_
        acf_%dt      = dt_
        acf_%tau     = 0
        acf_%num_pe  = num_pe_
        acf_%my_rank = my_rank_
        acf_%comm    = comm_

        allocate(acf_%Kt(0:acf_%Ntau-1))
        acf_%Kt = 0.
        allocate(acf_%num_contributions(0:acf_%Ntau-1))
        acf_%num_contributions = 0
        allocate(acf_%oldvals(1:3,1:acf_%Ntau))
      end subroutine


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !> Finalizes the ACF class
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine acf_finalize(acf_)
        implicit none
        class(acf) :: acf_

        if (allocated(acf_%Kt))                deallocate(acf_%Kt)
        if (allocated(acf_%num_contributions)) deallocate(acf_%num_contributions)
        if (allocated(acf_%oldvals))           deallocate(acf_%oldvals)
      end subroutine



      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !> Increments internal timestep counter and includes val as current
      !> value of the quantity under consideration
      !> May only be called once per timestep
      !> @param val value of considered quantity
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
          acf_%Kt(s)                = acf_%Kt(s)  + dot_product(val, acf_%oldvals(1:3,acf_%tau - s))
          acf_%num_contributions(s) = acf_%num_contributions(s) + 1
        end do

        call MPI_ALLGATHERV(MPI_IN_PLACE, recvcounts(acf_%my_rank), MPI_REAL8,   acf_%Kt,                recvcounts, displs, MPI_REAL8,   acf_%comm, ierr)
        call MPI_ALLGATHERV(MPI_IN_PLACE, recvcounts(acf_%my_rank), MPI_INTEGER, acf_%num_contributions, recvcounts, displs, MPI_INTEGER, acf_%comm, ierr)

        deallocate(recvcounts, displs)

      end subroutine


end module module_acf
