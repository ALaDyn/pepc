!  ================================
!
!         DIAGNOSTICS
!
!     Perform diagnostics
!
!  $Revision 1.15$
!
!  ================================


subroutine diagnostics

  use physvars
  use treevars
  use utils

  implicit none
  integer :: i,lvisit_active




! Interface to VISIT (Online visualisation)

  if ( vis_on .and. mod(itime,ivis*2)==0 .and. steering) call beam_control

  if (u_beam>0 .and. beam_config==3) scheme=1  ! Switch off Te control if beam on

  if ( vis_on ) then
     if ( mod(itime,ivis) ==0 ) call vis_parts       
     if ( mod(itime,ivis_fields) ==0 ) call vis_fields       
  endif

  if (initial_config.eq.4) then
    do i=1,npp
	if (pelabel(i)==1) then
  write (90,'(6f15.4)') x(i),y(i),z(i),ux(i),uy(i),uz(i)
	endif
    enddo
  endif

  if (mod(itime,idump) >= idump-navcycle) call sum_fields    ! Accumulate cycle-averaged fields on grid
!  - assume for now that idump > navcycle

  if (beam_config == 4 .and. mod(itime,itrack)==0 ) call track_nc          ! Gather densities and track critical surface 

  if (vis_on .and. mod(itime,ivis_fields)==0 ) then
!     call pot_grid
     call densities
     call vis_fields
  endif

  if (itime_start>0 .and. itime==0) return  ! Avoid over-writing restart data
  call energy_cons       ! Compute energy balance

  if ( dump_tree .and. mod(itime,idump) ==0 ) then
     call diagnose_tree   ! Printed tree info (htable etc)
     call draw_tree2d(xl,yl)     ! Draw PE-trees
     call draw_lists      ! Draw interaction lists
     call draw_domains(itime+itime_start)   ! Domains
  endif

  if ((mod(itime,idump)==0 .or. itime==nt) ) then
     call dump(itime+itime_start)     ! Dump complete set of particle data
     call dump_fields(itime+itime_start)  ! Field data
     if (vis_on)  call vis_fields
  endif



end subroutine diagnostics

