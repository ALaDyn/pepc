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

! ==============================================
!
!                VIS_CONFIG
!
!  Interactive control of target config
!
! ==============================================

subroutine vis_config

  use physvars
  use treevars
  use utils
  implicit none
  include 'mpif.h'

  integer :: i, p, ierr
  integer :: lvisit_active

  integer :: isteer1 = 0, isteer2 =0, isteer3=0, isteer4=1
  real*8 :: dsteer1, dsteer2, dsteer3, dsteer4

 
  integer, save :: np_beam_dt  ! current # beam particles

  ! First check for VISIT connection

#ifdef VISIT_NBODY
  if (me==0)   call flvisit_nbody2_check_connection(lvisit_active)
  call MPI_BCAST( lvisit_active, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
#endif

  if (lvisit_active==0 )then
     if (me==0) write(*,*) 'VISNB | No connection to visualization'
     return
  endif



#ifdef VISIT_NBODY

  if (me==0) then

     call flvisit_nbody2_check_connection(lvisit_active)

     ! Fetch real-time, user-specified control parameters
     if (lvisit_active /= 0) then 
        call flvisit_nbody2_steering_recv( dsteer1,dsteer2,dsteer3, dsteer4,isteer1,isteer2,isteer3,isteer4)

	ivis = isteer2
	ivis_fields=isteer3
	target_geometry = max(isteer4,0)  ! Target geom 0-8

        stretch: select case(target_geometry)

        case(0,5) ! slab & wedge
           x_plasma = dsteer4 

        case(1,11,2,3,13,4,6,16,26,36) ! spherical/cylindrical
           r_sphere = dsteer4
	   x_plasma = min(dsteer4/2.,dsteer3)

	case(12) ! tube
	   r_sphere = dsteer4
 	   z_plasma = min(dsteer4/2.,dsteer3)

	case default
	   r_sphere = dsteer4
	   x_plasma = dsteer3

        end select stretch

	if (isteer1 == 1) launch = .true.  ! Start simulation if launch box clicked.

	write(*,*) 'VISCO | vis freq = ',ivis
	write(*,*) 'VISCO | geometry = ',isteer4
	write(*,*) 'VISCO | nparts = ',npart_total

     else
        write(*,*) 'VISCO | No connection to visualization'
        return
     endif
     
  endif

#endif

  ! Broadcast beam parameters to all other PEs
  call MPI_BARRIER( MPI_COMM_WORLD, ierr)   ! Synchronize first
  call MPI_BCAST( lvisit_active, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)

  if (lvisit_active /= 0) then
     call MPI_BCAST( ivis, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
     call MPI_BCAST( ivis_fields, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
     call MPI_BCAST( target_geometry, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
     call MPI_BCAST( x_plasma, 1, MPI_REAL, 0, MPI_COMM_WORLD,ierr)
     call MPI_BCAST( r_sphere, 1, MPI_REAL, 0, MPI_COMM_WORLD,ierr)
     call MPI_BCAST( launch, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD,ierr)
  else
     if (me==0) write(*,*) ' No Connection to Visualization'
     return
  endif



end subroutine vis_config
