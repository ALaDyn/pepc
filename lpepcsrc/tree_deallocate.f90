subroutine tree_deallocate(nppm_ori)

  use treevars
  implicit none
  include 'mpif.h'

  integer, intent(in) :: nppm_ori

  nppm = nppm_ori
  if (me==0 .and. tree_debug) write(*,*) 'Deallocating multipole fields'

  deallocate ( htable, all_addr, free_addr, point_free, &
       treekey, branch_key, branch_owner, &
       pebranch, leaf_key, twig_key )

  deallocate ( first_child, node_level )

! multipole moments
  deallocate ( charge, &                    ! charge
       abs_charge, &                ! absolute charge
       xcoc, ycoc, zcoc, size_node, &    ! centre of charge 
       xshift, yshift, zshift, &    ! shift vector
       xdip, ydip, zdip, &          ! dipole moment
       xxquad, yyquad, zzquad, &       ! quadrupole moment
       xyquad, yzquad, zxquad, &
       jx, jy, jz, &      ! current
       magmx, magmy, magmz ) ! magnetic moment 

end subroutine tree_deallocate

  
