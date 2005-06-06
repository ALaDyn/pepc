!  ================================
!
!         DIAGNOSTICS
!
!     Perform diagnostics
!
!  $Revision 1.15$
!
!  ================================


subroutine diagnostics(dump_tree)

  use physvars

  implicit none
  integer :: i,lvisit_active
  logical :: dump_tree




  ! Interface to VISIT (Online visualisation)


  if (u_beam>0 .and. beam_config==3) scheme=1  ! Switch off Te control if beam on

  if ( vis_on ) then
!     if ( mod(itime,ivis) ==0 ) call vis_parts       
     if ( mod(itime,ivis) ==0 ) call vis_parts_nbody       
     if ( mod(itime,ivis_domains) ==0 ) call vis_domains_nbody       
     if ( mod(itime,ivis)==0 .and. steering) call beam_control
     if ( mod(itime,ivis_fields)==0 ) then
        !     call pot_grid
        call densities
        call sum_fields
        call vis_fields
     endif
  endif

!  if (target_geometry.eq.4) then
!    do i=1,npp
!	if (pelabel(i)==60) then
!	if (itime.eq.0) write(90,'(a)') '! t x ux Ex Ax Axo -dA/dt Bx' 
 !       write (90,'(8(1pe15.4))') itime*dt,x(i),ux(i),Ex(i),Ax(i),Axo(i),(Axo(i)-Ax(i))/dt,Bx(i)
!	endif
!    enddo
!  endif

  if (mod(itime,idump) >= idump-navcycle) call sum_fields    ! Accumulate cycle-averaged fields on grid
!  - assume for now that idump > navcycle

  if (beam_config == 4 .and. mod(itime,itrack)==0 ) call track_nc          ! Gather densities and track critical surface 



  if (itime_start>0 .and. itime==0) return  ! Avoid over-writing restart data
  call energy_cons(Ukine,Ukini,Umagnetic,Ubeam)       ! Compute energy balance

  if ( dump_tree .and. mod(itime,iprot) ==0 ) then
     call diagnose_tree   ! Printed tree info (htable etc)
     call draw_tree2d(xl,yl)     ! Draw PE-trees
     call draw_lists      ! Draw interaction lists
     call draw_domains(itime+itime_start)   ! Domains
  endif

  if ((mod(itime+itime_start,idump)==0 .or. itime==nt) ) then
     call dump(itime+itime_start)     ! Dump complete set of particle data
     call dump_fields(itime+itime_start)  ! Field data
     if (vis_on)  call vis_fields
  endif



end subroutine diagnostics





