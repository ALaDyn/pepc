

!  ===================================================================
!
!                              KINENERGY
!
!   Calculate kinetic energy
!
!  ===================================================================

subroutine kinenergy(ekine,ekini,ebeam)
  use physvars
  use treevars

  implicit none
  include 'mpif.h'

  integer :: p,i,ierr
  real*8 :: ekine, ekini, ebeam, sum_plas_e, sum_plas_i, sum_beam, gamma
  real*8, dimension(nppm) :: uhx, uhy, uhz

  sum_plas_e = 0.
  sum_plas_i = 0.

  sum_beam = 0.


  do p=1, npp
  ! Velocities at previous 1/2-step to synch with P.E.
    uhx(p) = ux(p)-dt*q(p)*Ex(p)/m(p)/2. 
    uhy(p) = uy(p)-dt*q(p)*Ey(p)/m(p)/2.
    uhz(p) = uz(p)-dt*q(p)*Ez(p)/m(p)/2.
    gamma = sqrt(1.0 + uhx(p)**2 + uhy(p)**2 + uhz(p)**2)
    if (pelabel(p) <= ne) then
     !  Sum local plasma electron kinetic energy
      if (scheme.eq.7) then
  ! non-relativistic - u not normalised to c
        sum_plas_e = sum_plas_e + 0.5*m(p)*(uhx(p)**2+uhy(p)**2+uhz(p)**2)
      else
        sum_plas_e = sum_plas_e + m(p)*(gamma - 1.0)
      endif

    else if (pelabel(p) <=ne+ni) then
      if (scheme.eq.7) then
        sum_plas_i = sum_plas_i + 0.5*m(p)*(uhx(p)**2+uhy(p)**2+uhz(p)**2)
      else
        sum_plas_i = sum_plas_i + m(p)*(gamma - 1.0)
      endif

    else
     !  Sum beam energy
      if (scheme.eq.7) then
        sum_beam = sum_beam + 0.5*m(p)*(uhx(p)**2+uhy(p)**2+uhz(p)**2)
      else
        sum_beam = sum_beam + m(p)*(gamma - 1.0)
      endif
    endif
 end do

! Gather partial sums together for global energies
 
  call MPI_ALLREDUCE(sum_plas_e, ekine, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(sum_plas_i, ekini, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(sum_beam, ebeam, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)


end subroutine kinenergy



