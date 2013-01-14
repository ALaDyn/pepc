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

! ======================
!
!   VIS_DOMAINS_NBODY
!
!   Send domain data to VISIT-XNBODY for real-time visualisation
!
!
! ======================

subroutine vis_domains_nbody


  use physvars
  use treevars
  use module_spacefilling
  implicit none   


  integer, parameter :: npart_visit_max = 10000  ! Max data points for VIS
  integer, parameter :: ship_max = 10000, attrib_max=22
  real*4, dimension(0:attrib_max-1,npart_visit_max) :: vbuffer
  integer :: lvisit_active, nskip, nproot, ne_buf, ni_buf, npart_buf, nbufe, nbufi
  integer :: i, j, k, ioffset,ixd, iyd, izd, ilev, lcount, wfdatai

  real :: s, simtime, t_display, dummy, xd,yd,zd, dx, dz, dy, epond_max, box_max, epondx, epondy, epondz,phipond, lbox, work_ave
  integer :: nship, ierr, cpuid


  simtime = dt*(itime+itime_start)

  if (beam_config>=3 .and. beam_config<=5) then
     t_display = tlaser*convert_fs
  else
     t_display = simtime
  endif
  work_ave=SUM(work_loads)/num_pe

  if (me==0) then

     ! ship branch nodes to show domains
     do j=1,nbranch_sum
        ilev = level_from_key(branch_key(j))
        ixd = SUM( (/ (2**i*ibits( branch_key(j),3*i,1 ), i=0,ilev-1) /) )
        iyd = SUM( (/ (2**i*ibits( branch_key(j),3*i+1,1 ), i=0,ilev-1) /) )
        izd = SUM( (/ (2**i*ibits( branch_key(j),3*i+2,1 ), i=0,ilev-1) /) )
        lbox = boxsize/2**(ilev)          !  box length
        ! Store attributes for visualizing
        vbuffer(0,j) = t_display
        vbuffer(1,j)= xmin + ixd*lbox + lbox/2 ! corners?
        vbuffer(2,j)= ymin + lbox*(iyd+.5)
        vbuffer(3,j)= zmin + lbox*(izd+.5) 
        vbuffer(4,j) = lbox
        vbuffer(5,j) = branch_owner(j)  ! cpu id of branch node
        cpuid = branch_owner(j)+1
        vbuffer(6,j) = work_loads(cpuid)/work_ave   ! total work load of branch node/cpu
        vbuffer(7,j) = npps(cpuid) ! total # particles on cpu
        vbuffer(8,j) = 0.
        vbuffer(9,j) = 0.
        vbuffer(10,j) = 0.
        vbuffer(11,j) = 0.
        vbuffer(12,j) = 0.
        vbuffer(13,j) = 0.
        vbuffer(14,j) = 0.
        vbuffer(15,j) = 0.
        vbuffer(16,j) = 16    ! Domain type
        vbuffer(17,j) = nbranch_sum  ! Total # branch nodes
        vbuffer(18,j) = 0.
        vbuffer(19,j) = 0.
        vbuffer(20,j) = 0.
        vbuffer(21,j) = 0.

!        write (*,'(7f13.4)') vbuffer(1,j),vbuffer(2,j), vbuffer(4,j), vbuffer(5,j), & 
!             vbuffer(6,j), vbuffer(7,j), vbuffer(17,j)
     end do
! root box
     j=nbranch_sum+1  
        vbuffer(0,j) = t_display
        vbuffer(1,j)= xmin + boxsize/2! corners?
        vbuffer(2,j)= ymin + boxsize/2
        vbuffer(3,j)= zmin + boxsize/2 
        vbuffer(4,j) = boxsize
        vbuffer(5,j) = 0  ! cpu id of branch node
        cpuid = branch_owner(j)+1
        vbuffer(6,j) = 1.   ! total work load of branch node/cpu
        vbuffer(7,j) = npart ! total # particles on cpu
        vbuffer(8,j) = 0.
        vbuffer(9,j) = 0.
        vbuffer(10,j) = 0.
        vbuffer(11,j) = 0.
        vbuffer(12,j) = 0.
        vbuffer(13,j) = 0.
        vbuffer(14,j) = 0.
        vbuffer(15,j) = 0.
        vbuffer(16,j) = 16    ! Domain type
        vbuffer(17,j) = nbranch_sum+1  ! Total # branch nodes
        vbuffer(18,j) = 0.
     
     call flvisit_nbody2_check_connection(lvisit_active)

     call flvisit_nbody2_partstep_send(vbuffer,nbranch_sum,attrib_max)
     


  endif


end subroutine vis_domains_nbody


