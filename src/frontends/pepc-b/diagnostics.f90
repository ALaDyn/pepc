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
  use module_gle

  implicit none
  include 'mpif.h'


 if ( mod(itime,ivis_fields)==0 ) then
  if (idim.eq.2) then
   call densities_2d
   call fields_2d
   call dump_fields_2d(itime+itime_start)
  else
   call densities
   call sum_fields
  endif
 endif 


  ! Routines for particle and field processing
  ! for VISIT (Online visualisation) and/or Netcdf file

#ifdef VISIT_NBODY
  if ( vis_on ) then
     !     if ( mod(itime,ivis) ==0 ) call vis_parts       
     if ( mod(itime,ivis) ==0 ) then
        call vis_parts_nbody(vcount)
        vcount = vcount + 1
     endif
     if ( mod(itime,min(ivis,ivis_fields))==0 .and. steering) call vis_control
     if ( mod(itime,ivis_fields)==0 ) then
        !     call pot_grid
        call vis_fields_nbody(itime+itime_start)
        call vis_vecfields_nbody(itime+itime_start)
     endif

  endif
#endif

  !  if (target_geometry.eq.4) then
  !    do i=1,npp
  !	if (pelabel(i)==60) then
  !	if (itime.eq.0) write(90,'(a)') '! t x ux Ex Ax Axo -dA/dt Bx' 
  !       write (90,'(8(1pe15.4))') itime*dt,x(i),ux(i),Ex(i),Ax(i),Axo(i),(Axo(i)-Ax(i))/dt,Bx(i)
  !	endif
  !    enddo
  !  endif

!  if (mod(itime,idump) >= idump-navcycle .or. nt.lt.idump) then
!	call sum_fields    ! Accumulate cycle-averaged fields on grid
!	call sum_fieldave
!  endif

  !  - assume for now that idump > navcycle
  !  - should really use running average to be compatible with online vis.


  if (beam_config == 4 .and. mod(itime,itrack)==0 ) call track_nc          ! Gather densities and track critical surface 

  if( mod(itime+itime_start,ivis_fields)==0) then
     if (mod(target_geometry,10)==1) call sum_radial(itime+itime_start)  ! Radial moments
     call field_lineout(itime+itime_start) ! lineouts
  endif


  if (itime_start>0 .and. itime==0) return  ! Avoid over-writing restart data

  call energy_cons(Ukine,Ukini,Umagnetic,Ubeam)       ! Compute energy balance
  call laser_hist

  if ( idump>0 .and. (mod(itime+itime_start,idump)==0 .or. itime==nt) ) then
     call dump(itime+itime_start)     ! Dump complete set of particle data
     if (idim.eq.3) then
      call dump_fields(itime+itime_start)  ! Field data
     else
!      call dump_fields_2d(itime+itime_start)
     endif
  endif



end subroutine diagnostics




