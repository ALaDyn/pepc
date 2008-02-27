
!  ================================
!
!         Restore initial particle order
!
!  ================================


subroutine restore_p1(npnew,npold,indxl,irnkl,islen,irlen,fposts,gposts)

  use physvars
  implicit none
  include 'mpif.h'

  integer, intent(in) :: npnew,npold
  integer, intent(in) :: indxl(npold),irnkl(npnew)
  integer, intent(in) :: islen(n_cpu),irlen(n_cpu)
  integer, intent(in) :: fposts(n_cpu+1),gposts(n_cpu+1)

  integer :: i, ierr
  type (particle_p1) :: ship_parts(npold), get_parts(npnew)

  call MPI_BARRIER( MPI_COMM_WORLD, ierr)  ! Wait for everyone to catch up

  do i=1,npold
     ship_parts(i) = particle_p1( x(indxl(i)), y(indxl(i)), z(indxl(i)), &
          ux(indxl(i)), uy(indxl(i)), uz(indxl(i)), &
          q(indxl(i)), m(indxl(i)), work(indxl(i)), pelabel(indxl(i)) )
  enddo

   ! perform permute
  call MPI_alltoallv(  ship_parts, islen, fposts, mpi_type_particle_p1, &
      get_parts, irlen, gposts, mpi_type_particle_p1, &
      MPI_COMM_WORLD,ierr ) 

  do i=1,npnew
     x(irnkl(i)) = get_parts(i)%x
     y(irnkl(i)) = get_parts(i)%y
     z(irnkl(i)) = get_parts(i)%z
     ux(irnkl(i)) = get_parts(i)%ux
     uy(irnkl(i)) = get_parts(i)%uy
     uz(irnkl(i)) = get_parts(i)%uz
     q(irnkl(i)) = get_parts(i)%q
     m(irnkl(i)) = get_parts(i)%m
     work(irnkl(i)) = get_parts(i)%work
     pelabel(irnkl(i)) = get_parts(i)%label
  enddo

  call MPI_BARRIER( MPI_COMM_WORLD, ierr)  ! Wait for everyone to catch up

end subroutine restore_p1
