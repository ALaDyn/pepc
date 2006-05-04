

!  ===================================================================
!
!                              KINENERGY
!
!   Calculate kinetic energy
!
!  ===================================================================

subroutine kinenergy(ekin_dust,ekin_star)
  use treevars
  use physvars

  implicit none
  include 'mpif.h'

  integer :: p,i, ierr
  real :: ekin_dust, ekin_star, sum_dust, sum_star, u2
  real, dimension(nppm) :: uhx, uhy, uhz

  sum_dust = 0.
  sum_star = 0.

  do p=1, npp
    uhx(p) = ux(p)-dt*ax(p)/2.   ! Velocities at previous 1/2-step to synch with P.E.
    uhy(p) = uy(p)-dt*ay(p)/2.
    uhz(p) = uz(p)-dt*az(p)/2.
    u2 = uhx(p)**2+uhy(p)**2+uhz(p)**2
    sum_dust = sum_dust + 0.5*m(p)*u2
  end do

  do i=1,ni
     u2 = ux_star(i)**2+uy_star(i)**2+uz_star(i)**2
     sum_star = sum_star + 0.5*m_star(i)*u2
  end do

! Gather partial sums together for global dust energy
 
  call MPI_ALLREDUCE(sum_dust, ekin_dust, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

  ekin_star = sum_star

end subroutine kinenergy
