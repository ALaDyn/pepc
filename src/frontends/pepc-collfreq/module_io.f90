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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>  Encapsulates any file i/o, opening, closin etc.
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_io
      use treevars
      implicit none
      save
      private

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      character*4, public :: csubme   !< Character string of data subdirectory 'data/peXXXX'

      integer, parameter, public :: file_pepc_out    = 24
      integer, parameter, public :: file_domains_dat = 70
      integer, parameter, public :: file_laser_dat   = 71
      integer, parameter, public :: file_energy_dat  = 75
      integer, parameter, public :: file_parts_info_in = 80

      integer, parameter, public :: file_ipefile     = 20

      integer, parameter, public :: file_stdout = 6

      integer, parameter, public :: file_tempfile_1  = 60
      integer, parameter, public :: file_tempfile_2  = 61
      integer, parameter, public :: file_tempfile_3  = 62

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      public openfiles
      public flushfiles
      public closefiles
      public stamp

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      contains

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!>
		!> Open files for result and debugging output
		!> - Rank 0:
		!>          - tree.out (\e 15): Tree stats
		!>          - pepc.out (\e 24): Physics log
		!>          - domains.dat (\e 70):
		!>          - laser.dat (\e 71): Laser parameters
        !>          - energy.dat (\e 75): energies
		!>          .
		!> - all ranks (only if \code(debug_level>2 .and. idump>0)\endcode) :
		!>          - data/out.1234 (\e 20):
		!>
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		subroutine openfiles

		  use physvars
		  implicit none


		  if (my_rank == 0) then
		     !  master diagnostics output
		     open(file_pepc_out, file='pepc.out') ! Physics log

		     open(file_domains_dat, file='domains.dat')
		     open(file_laser_dat,   file='laser.dat')       ! laser parameters

		     write(*,*) 'debug level: ',debug_level,' idump',idump
		  endif

		end subroutine openfiles


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine flushfiles
      use physvars
      implicit none

      if (my_rank == 0) then
         !  master diagnostics output
         flush(file_pepc_out)

         flush(file_domains_dat)
         flush(file_laser_dat)
      endif

    end subroutine flushfiles


		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!>
		!> Tidy up O/P files
		!>
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		subroutine closefiles
		  use physvars
		  implicit none

		  if (me == 0) then
             close(file_pepc_out)
		     close(file_domains_dat)
		     close(file_laser_dat)
             close(file_energy_dat)
		     close(90)
		  endif
		  !if (db_level > 0) close(file_ipefile)
		  close(file_parts_info_in)  ! initial particle data


		end subroutine closefiles


		!     =========================
		!
		!     Time stamp
		!
		!     =========================

		subroutine stamp(istream,ibegin)
		  implicit none

		  character :: cdate*8, ctime*10, czone*5
		  integer :: ibegin
		  integer :: istream

		     !      call DATE_AND_TIME(cdate,ctime,czone,vals)
		     call DATE_AND_TIME(cdate,ctime,czone)

		     if (ibegin.eq.1) then

		        write(istream,'(//a20,a12/a20,a12/a20,a12//)') 'PEPC run on ' &
		             ,cdate(7:8)//'/'//cdate(5:6)//'/'//cdate(1:4) &
		             ,'Time: ',ctime(1:2)//':'//ctime(3:4)//':'//ctime(5:6),' GMT+',czone

		     else
		        write(istream,'(a,a9)') 'Finished run at time: ',ctime(1:2)//':'//ctime(3:4)//':'//ctime(5:6)
		     endif
		end subroutine stamp

end module module_io
