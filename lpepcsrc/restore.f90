
!  ================================
!
!         Restore initial particle order
!
!  ================================


subroutine restore(npnew,npold,nppm_ori,indxl,irnkl,islen,irlen,fposts,gposts,&
           pot_tmp,ex_tmp,ey_tmp,ez_tmp,w_tmp,p_pot,p_ex,p_ey,p_ez,p_w)

  use treevars
  implicit none
  include 'mpif.h'

  integer, intent(in) :: npnew,npold,nppm_ori
  integer, intent(in) :: indxl(nppm_ori),irnkl(nppm_ori)
  integer, intent(in) :: islen(num_pe),irlen(num_pe)
  integer, intent(in) :: fposts(num_pe+1),gposts(num_pe+1)
  real*8, intent(in), dimension(npnew) :: ex_tmp,ey_tmp,ez_tmp,pot_tmp,w_tmp
  real*8, intent(out), dimension(npold) :: p_ex, p_ey, p_ez, p_pot,p_w  ! fields and potential to return

  integer :: i, ierr

  type (results) :: ship_parts(npnew), get_parts(npold)

  if (me==0 .and. tree_debug) write(*,'(a)') 'LPEPC | DOMAINS..'

  call MPI_BARRIER( MPI_COMM_WORLD, ierr)  ! Wait for everyone to catch up

  do i=1,npnew
     ship_parts(i) = results( ex_tmp(indxl(i)), ey_tmp(indxl(i)), ez_tmp(indxl(i)), &
                              pot_tmp(indxl(i)), w_tmp(indxl(i)), pelabel(indxl(i)) )
  enddo

  ! perform permute
  call MPI_alltoallv(  ship_parts, islen, fposts, mpi_type_results, &
      get_parts, irlen, gposts, mpi_type_results, &
      MPI_COMM_WORLD,ierr ) 

  do i=1,npold
     p_ex(irnkl(i)) = get_parts(i)%Ex
     p_ey(irnkl(i)) = get_parts(i)%Ey
     p_ez(irnkl(i)) = get_parts(i)%Ez
     p_pot(irnkl(i)) = get_parts(i)%pot
     p_w(irnkl(i)) = get_parts(i)%work
     pelabel(irnkl(i)) = get_parts(i)%label
  enddo

  call MPI_BARRIER( MPI_COMM_WORLD, ierr)  ! Wait for everyone to catch up

end subroutine restore
