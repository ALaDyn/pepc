subroutine pepc_cleanup(my_rank,n_cpu)
  use treevars
  implicit none
  include 'mpif.h'

  integer, intent(in) :: my_rank  ! MPI cpu rank
  integer, intent(in) :: n_cpu  ! MPI # CPUs

! copy call parameters to treevars module

  me = my_rank
  num_pe = n_cpu
 

  if (me==0) then
     write(*,'(a)') 'LPEPC | De-allocating particle and tree arrays ...'
  endif

  ! particle array deallocation

  deallocate ( x, y, z, ux, uy, uz, & 
       q, m, work, &
       Ex, Ey, Ez, pot, &
       Ax, Ay, Az, &
       Bx, By, Bz,  &
       Axo, Ayo, Azo, &
       pepid, pelabel, pekey )    

  deallocate ( nterm, intlist, nodelist ) ! interaction key-, node-lists

  deallocate ( htable, all_addr, free_addr, point_free, &
       nbranches, igap, &
       treekey, branch_key, branch_owner, &
       pebranch, leaf_key, twig_key, &
       fetched_owner, fetched_keys, requested_owner, requested_keys, &
       nreqs_total, nfetch_total )

  deallocate ( first_child, n_children, node_level )

! multipole moments
  deallocate ( charge, &                    ! charge
       abs_charge, &                ! absolute charge
       xcoc, ycoc, zcoc, &    ! centre of charge 
       xshift, yshift, zshift, &    ! shift vector
       xdip, ydip, zdip, &          ! dipole moment
       xxquad, yyquad, zzquad, &       ! quadrupole moment
       xyquad, yzquad, zxquad, &
       jx, jy, jz, &      ! current
       magmx, magmy, magmz ) ! magnetic moment 

  deallocate ( pack_child,get_child )    ! Multipole shipping buffers

! work balance arrays

  deallocate (work_loads,npps,pivots)

!  if (me==0) then
!     write(*,'(a)') 'LPEPC | De-allocating npps ...'
!  endif

!  deallocate  (npps)  ! Work load & Particle distrib amoung PEs

!  if (me==0) then
!     write(*,'(a)') 'LPEPC | De-allocating pivots ...'
!  endif

!  deallocate (pivots)



end subroutine pepc_cleanup






