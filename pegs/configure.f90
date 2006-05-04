! ==============================================
!
!                CONFIGURE
!
!  Perform diagnostics on initial config and
!  construct tree
!
! ==============================================

subroutine configure

  use treevars
  use physvars
  use utils

  implicit none
  integer :: i, ipe, idummy=0, ifile
  real :: t_walk, t_force



!  call diagnostics
 
  if (me==0) then
     do ifile = 6,15,9
        write(ifile,'(//a,i8,a,f10.5)') 'Timestep ',itime,' t=',itime*dt
     end do
  endif

  ! Initial tree construction and force computation
 
  call tree_domains    ! Domain decomposition: allocate particle keys to PEs
  call tree_build      ! Build trees from local particle lists
  call tree_branches   ! Determine and concatenate branch nodes
  call tree_fill       ! Fill in remainder of local tree
  call tree_properties ! Compute multipole moments for local tree
  call forces(1,npp,0.,t_walk,t_force)          ! Calculate initial potentials and forces
  call stars(0.)
  call diagnostics

end subroutine configure


