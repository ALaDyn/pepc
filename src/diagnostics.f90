!  ================================
!
!         DIAGNOSTICS
!
!     Perform diagnostics
!
!  ================================


subroutine diagnostics

  use treevars
  use utils

  implicit none
  integer :: lvisit_active

  if ( vis_on .and. mod(itime,ivis)==0 .and. steering) call beam_control


  if ( mod(itime,ivis) ==0 ) then
     if (vis_on) call vis_parts       ! Interface to VISIT
  endif


  if (vis_on .and. mod(itime,ivis_fields)==0 ) then
!     call pot_grid
     call densities
     call vis_fields
  endif

  if (itime_start>0 .and. itime==0) return  ! Avoid over-writing restart data
  call energy_cons       ! Compute energy balance
  if (mod(itime,idens) == 0) call densities    ! Compute electron, ion densities for diagnostics

  if (beam_config == 4 .and. mod(itime,idens)==0 ) call track_nc          ! Gather densities and track critical surface 
  if (perf_anal) return

  if ( dump_tree .and. mod(itime,idump) ==0 ) then
     call diagnose_tree   ! Printed tree info (htable etc)
     call draw_tree2d     ! Draw PE-trees
     call draw_lists      ! Draw interaction lists
     call draw_domains(itime+itime_start)   ! Domains
  endif

  if ((mod(itime,idump)==0 ) ) then
     if (itime.ne.0) call dump(itime+itime_start)     ! Dump complete set of particle data
     call slices(itime+itime_start)  ! 1D lineouts

  endif

end subroutine diagnostics



