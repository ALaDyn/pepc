module energies
  implicit none
  private

    public energy_cons

  contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !>     Find potential, kinetic energies
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine energy_cons(ekine,ekini)
      use physvars
      use module_fmm_framework
      implicit none

      real*8 :: epot, ekine, ekini, etot

      call potenergy(epot)
      call kinenergy(ekine, ekini)

      etot = epot + ekine + ekini

      if (my_rank == 0) then
    !  if ( my_rank == 0 .and. db_level.ge.1 ) then
    !     do ifile = 6,15,9
    !        write (ifile,'(4(a20,1pe12.5/))') &
    !	     ' P.E. = ',epot, &
    !	     ' Electron K.E. = ',ekine, &
    !             ' Ion K.E. = ',ekini, &
    !	     ' Total: ',etot
    !
    !     end do
         ! Write out to energy.dat file
         open(75,file='energy.dat',STATUS='UNKNOWN', POSITION = 'APPEND')
         if (itime.eq.0)  write(75,'(a)') '! time  Upot(total)  Upot(near field) Upot(far field)  Ukin_e Ukin_i Ukin_e+i Utot '
         write (75,'(f12.5,7(1pe13.4))') trun, epot, potnearfield, potfarfield, ekine, ekini, ekine+ekini, etot
         close(75)
      endif
    end subroutine energy_cons


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !>  Calculate kinetic energy
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine kinenergy(ekine,ekini)
      use physvars

      implicit none
      include 'mpif.h'

      integer :: p,ierr
      real*8 :: ekine, ekini, ebeam, sum_plas_e, sum_plas_i, sum_beam, gamma
      real*8, dimension(nppm) :: uhx, uhy, uhz

      sum_plas_e = 0.
      sum_plas_i = 0.
      sum_beam = 0.

      do p=1, np_local
      ! Velocities at previous 1/2-step to synch with P.E.
        uhx(p) = ux(p)-dt*q(p)*Ex(p)/m(p)/2.
        uhy(p) = uy(p)-dt*q(p)*Ey(p)/m(p)/2.
        uhz(p) = uz(p)-dt*q(p)*Ez(p)/m(p)/2.
        gamma = sqrt(1.0 + uhx(p)**2 + uhy(p)**2 + uhz(p)**2)
        if (pelabel(p) <= ne) then
         !  Sum local plasma electron kinetic energy
          sum_plas_e = sum_plas_e + m(p)*(gamma - 1.0)

        else if (pelabel(p) <=ne+ni) then
          sum_plas_i = sum_plas_i + m(p)*(gamma - 1.0)

           else
         !  Sum beam energy
          sum_beam = sum_beam + m(p)*(gamma - 1.0)
        endif
     end do

      ! Gather partial sums together for global energies
      call MPI_ALLREDUCE(sum_plas_e, ekine, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(sum_plas_i, ekini, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(sum_beam, ebeam, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)


    end subroutine kinenergy



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !>  Calculate E.S. and magnetic potential energies
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine potenergy(epot_total)
      use physvars
      use module_fmm_framework
      implicit none
      include 'mpif.h'

      integer :: p, ierr

      real*8 :: upartial
      real*8, intent(out) :: epot_total
      logical :: pot_debug=.false.

      epot_total = 0.  ! Global potential energy

      !  Sum Potential Energy
      upartial = 0. ! Single PE partial potential energy sum

      do p=1, np_local
         upartial = upartial + 0.5*q(p)*pot(p)

         if (my_rank == 0 .and. pot_debug) then
            write (ifile_cpu,'(a,i5,a,i5,3f10.3,a,f12.4)') &
        'local particle ',p,' label ',pelabel(p),x(p),y(p),q(p),' pot ',pot(p)
         endif
      end do

      if (pot_debug) write (ifile_cpu,'(a,1pe11.4,i2)') 'partial PE sum',upartial,my_rank

      call MPI_ALLREDUCE(upartial, epot_total,1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

      ! this can also be done in fields.f90, but thematically it fits better here :-)
      potfarfield  = potfarfield/2.
      potnearfield = potnearfield/2.

      call MPI_ALLREDUCE(MPI_IN_PLACE, potfarfield,  1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, potnearfield, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

    end subroutine potenergy


end module energies
