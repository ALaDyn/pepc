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
  !VAMPINST include
  INCLUDE 'VTcommon.h'
  INTEGER VTIERR
  !

  implicit none
  integer :: lvisit_active

  !VAMPINST subroutine_start
  CALL VTENTER(IF_diagnostics,VTNOSCL,VTIERR)
  !      write(*,*) 'VT: diagnostics S>',VTIERR,
  !     *    IF_diagnostics,ICLASSH
  !
  if (itime_start>0 .and. itime==0) return  ! Avoid over-writing restart data

  call energy_cons       ! Compute energy balance
!  call densities    ! Compute electron, ion densities for diagnostics

  if (beam_config == 4) call track_nc          ! Gather densities and track critical surface 
  if (perf_anal)  Then

     RETURN
  END IF

  if ( dump_tree .and. mod(itime,idump) ==0 ) then
     call diagnose_tree   ! Printed tree info (htable etc)
     call draw_tree2d     ! Draw PE-trees
     call draw_lists      ! Draw interaction lists
     call draw_domains(itime+itime_start)   ! Domains
  endif

  if (mod(itime,idump)==0 .or. itime==nt ) then
!     call dump(itime+itime_start)     ! Dump complete set of particle data
!     call slices(itime+itime_start)  ! 1D lineouts

  endif

  if ( mod(itime,ivis)==0 .and. steering) call beam_control

  if (vis_on .and. mod(itime,ivis_fields)==0 ) then
     !     call vis_fields
  endif

  if ( mod(itime,ivis) ==0 ) then
     if (vis_on) call vis_parts       ! Interface to VISIT
     !     call visit_dump(itime+itime_start) ! Dump particle data (visit format)
  endif


  !VAMPINST subroutine_end
  CALL VTLEAVE(ICLASSH,VTIERR)
  !      write(*,*) 'VT: diagnostics S<',VTIERR,ICLASSH
  !
end subroutine diagnostics
