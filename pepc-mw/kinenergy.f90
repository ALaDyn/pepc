

!  ===================================================================
!
!                              KINENERGY
!
!   Calculate kinetic energy
!
!  ===================================================================

subroutine kinenergy(ekine,ekini,tempe,tempi)
  use physvars
  use module_units

  implicit none
  include 'mpif.h'

  integer :: p,ierr
  real*8 :: ekine, ekini,tempe,tempi, ebeam, sum_plas_e, sum_plas_i, sum_beam, gamma
  real*8 :: uh(3), uh2
  real*8 :: sum_v2e, sum_v2i, sum_ve(1:3), sum_vi(1:3)

  sum_plas_e = 0.
  sum_plas_i = 0.

  sum_beam = 0.

  sum_v2e=0.0
  sum_v2i=0.0
  sum_ve =0.0
  sum_vi =0.0


  do p=1, np_local
  ! Velocities at previous 1/2-step to synch with P.E.
    uh(1) = ux(p)-dt*q(p)*Ex(p)/m(p)/2.
    uh(2) = uy(p)-dt*q(p)*Ey(p)/m(p)/2.
    uh(3) = uz(p)-dt*q(p)*Ez(p)/m(p)/2.
    uh2   = dot_product(uh,uh)
    gamma = sqrt(1.0 + uh2/unit_c2)

    if (pelabel(p) <= ne) then
     !  Sum local plasma electron kinetic energy
      sum_plas_e = sum_plas_e + m(p)*unit_c2*(gamma - 1.0)
      sum_v2e = sum_v2e + uh2/gamma**2.
      sum_ve  = sum_ve  + uh/gamma

    else if (pelabel(p) <=ne+ni) then
      sum_plas_i = sum_plas_i + m(p)*unit_c2*(gamma - 1.0)
      sum_v2i = sum_v2i + uh2/gamma**2.
      sum_vi  = sum_vi  + uh/gamma

       else
     !  Sum beam energy
      sum_beam = sum_beam + m(p)*(gamma - 1.0)
    endif
 end do

! Gather partial sums together for global energies
 
  call MPI_ALLREDUCE(sum_plas_e, ekine, 1, MPI_REAL8, MPI_SUM, MPI_COMM_PEPC, ierr)
  call MPI_ALLREDUCE(sum_plas_i, ekini, 1, MPI_REAL8, MPI_SUM, MPI_COMM_PEPC, ierr)
  call MPI_ALLREDUCE(sum_beam, ebeam, 1, MPI_REAL8, MPI_SUM, MPI_COMM_PEPC, ierr)

  ! Find global KE sums
  call MPI_ALLREDUCE(MPI_IN_PLACE,sum_v2e, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE,sum_v2i, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE, sum_ve, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE, sum_vi, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

  sum_v2e = sum_v2e / ne
  sum_ve  = sum_ve  / ne
  sum_v2i = sum_v2i / ni
  sum_vi  = sum_vi  / ni

  tempe =  mass_e/(3.*unit_kB)*(sum_v2e - dot_product(sum_ve,sum_ve))
  tempi =  mass_i/(3.*unit_kB)*(sum_v2i - dot_product(sum_vi,sum_vi))

end subroutine kinenergy



