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

!  ================================
!
!         DIAGNOSTICS
!
!     Perform diagnostics
!
!
!  ================================

subroutine diagnostics

   use module_physvars
   use module_particle_props
   use module_field_grid
   use module_laser
   use module_diagnostics
   use module_io
   use mpi

   implicit none

   if (mod(itime, ivis_fields) .eq. 0) then
      if (idim .eq. 2) then
         call densities_2d
         call fields_2d
         call dump_fields_2d(itime + itime_start)
      else
         call densities
         call sum_fields
      end if
   end if

   ! Routines for particle and field processing
   ! for VISIT (Online visualisation) and/or Netcdf file

#ifdef VISIT_NBODY
   if (vis_on) then
      !     if ( mod(itime,ivis) ==0 ) call vis_parts
      if (mod(itime, ivis) .eq. 0) then
         call vis_parts_nbody(vcount)
         vcount = vcount + 1
      end if
      if (mod(itime, min(ivis, ivis_fields)) .eq. 0 .and. steering) call vis_control
      if (mod(itime, ivis_fields) .eq. 0) then
         !     call pot_grid
         call vis_fields_nbody(itime + itime_start)
         call vis_vecfields_nbody(itime + itime_start)
      end if

   end if
#endif

   !  if (target_geometry.eq.4) then
   !    do i=1,npp
   !        if (pelabel(i)==60) then
   !        if (itime.eq.0) write(90,'(a)') '! t x ux Ex Ax Axo -dA/dt Bx'
   !       write (90,'(8(1pe15.4))') itime*dt,x(i),ux(i),Ex(i),Ax(i),Axo(i),(Axo(i)-Ax(i))/dt,Bx(i)
   !        endif
   !    enddo
   !  endif

!  if (mod(itime,idump) >= idump-navcycle .or. nt.lt.idump) then
!        call sum_fields    ! Accumulate cycle-averaged fields on grid
!        call sum_fieldave
!  endif

   !  - assume for now that idump > navcycle
   !  - should really use running average to be compatible with online vis.

   if (beam_config .eq. 4 .and. mod(itime, itrack) .eq. 0) call track_nc          ! Gather densities and track critical surface

   if (mod(itime + itime_start, ivis_fields) .eq. 0) then
      if (mod(target_geometry, 10) .eq. 1) call sum_radial(itime + itime_start)  ! Radial moments
      call field_lineout(itime + itime_start) ! lineouts
   end if

   if (itime_start .gt. 0 .and. itime .eq. 0) return  ! Avoid over-writing restart data

   call energy_cons(Ukine, Ukini, Umagnetic, Ubeam)       ! Compute energy balance
   call laser_hist

   if (idump .gt. 0 .and. (mod(itime + itime_start, idump) .eq. 0 .or. itime .eq. nt)) then
      call dump(itime + itime_start)     ! Dump complete set of particle data
      if (idim .eq. 3) then
         call dump_fields(itime + itime_start)  ! Field data
      else
!      call dump_fields_2d(itime+itime_start)
      end if
   end if

end subroutine diagnostics

