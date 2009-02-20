subroutine tree_deallocate(nppm_ori)

  use treevars
  implicit none
  include 'mpif.h'

  integer :: ierr
  integer, intent(in) :: nppm_ori

  nppm = nppm_ori
  if (me==0) write(*,*) 'Deallocating multipole fields'

!  call MPI_BARRIER( MPI_COMM_WORLD, ierr)  ! Wait for everyone to catch up
  
 ! interaction key-, node-lists
  deallocate(nodelist,nterm,intlist)
       
  deallocate ( htable, all_addr, free_addr, point_free, &
       treekey, branch_key, branch_owner, &
       pebranch, leaf_key, twig_key, &
       fetched_owner, fetched_keys, requested_owner, requested_keys )

  deallocate ( first_child, n_children, node_level )

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

  deallocate ( pack_child,get_child )    ! Multipole shipping buffers

end subroutine tree_deallocate

  
