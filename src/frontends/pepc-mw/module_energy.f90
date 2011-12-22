!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>  Encapsulates ...
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_energy
    use physvars
    implicit none
    private

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    public energies

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Find potential, kinetic energies
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine energies(ekine,ekini)

        use physvars
        use module_fmm_framework
        use module_laser
        use module_pusher
        use module_io
        use module_units
        use module_param_dump

        implicit none

        real*8 :: epot, ekine, ekini, etot, tempe, tempi
        integer :: ifile
        logical, save :: firstcall = .true.

        call energy_pot(epot)

        call energy_kin(ekine, ekini, tempe, tempi)

        ! rescale total potential and kinetic energy by energy per particle
        ekine        = ekine        / ne
        ekini        = ekini        / ni
        epot         = epot         / npart_total
        potnearfield = potnearfield / npart_total
        potfarfield  = potfarfield  / npart_total

        etot = epot + ekine + ekini

        energy(3,1:np_local) = energy(1,1:np_local) + energy(2,1:np_local)

        if (my_rank == 0) then
            do ifile = 6,15,9
              call PrintEnergies(ifile, epot, ekini, ekine, etot, tempe, tempi)
            end do

            ! Write out to energy.dat file
            if (firstcall) then
              firstcall = .false.
              write(file_energy_dat,'("#",18(1x,a20))') "time", "Upot(total)", "Upot(near field)", "Upot(far field)", &
                                         "Ukin_e" ,"Ukin_i" ,"Ukin_tot", "Ukine(wo drift)" ,"Ukini(wo drift)" ,"Utot" ,&
                                         "Te0", "Te_uncor", "chie", "delta_Te", "Ti0", "Ti_uncor", "chii", "delta_Ti"
            endif

            write(file_energy_dat,'(" ",1(1x,f20.6),17(1x,1pe20.10))') trun*unit_t0_in_fs, epot, potnearfield, potfarfield, &
                                                  ekine, ekini, ekine+ekini, 3./2.*unit_kB*tempe, 3./2.*unit_kB*tempi, etot, &
                                                  Te0, Te_uncor, chie, delta_Te, Ti0, Ti_uncor, chii, delta_Ti
        endif
    end subroutine energies


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Calculate potential energies
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine energy_pot(epot_total)
        use physvars
        use module_fmm_framework
        implicit none
        include 'mpif.h'

        integer :: ierr

        real*8, intent(out) :: epot_total

        energy(1,1:np_local) = 0.5 * particles(1:np_local)%data%q*particles(1:np_local)%results%pot
        epot_total           = sum(energy(1,1:np_local))
        call MPI_ALLREDUCE(MPI_IN_PLACE, epot_total,1, MPI_REAL8, MPI_SUM, MPI_COMM_PEPC, ierr)

        ! this can also be done in fields.f90, but thematically it fits better here :-)
        potfarfield  = potfarfield/2.
        potnearfield = potnearfield/2.

        call MPI_ALLREDUCE(MPI_IN_PLACE, potfarfield,  1, MPI_REAL8, MPI_SUM, MPI_COMM_PEPC, ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE, potnearfield, 1, MPI_REAL8, MPI_SUM, MPI_COMM_PEPC, ierr)

    end subroutine energy_pot



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Calculate kinetic energies
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine energy_kin(ekine,ekini,tempe,tempi)
        use physvars
        use module_units

        implicit none
        include 'mpif.h'

        integer :: p,ierr
        real*8 :: ekine, ekini,tempe,tempi, gamma
        real*8 :: uh(3), uh2
        real*8 :: sum_v2e, sum_v2i, sum_ve(1:3), sum_vi(1:3)
        real*8 :: en

        ekine = 0.0
        ekini = 0.0
        sum_v2e    = 0.0
        sum_v2i    = 0.0
        sum_ve     = 0.0
        sum_vi     = 0.0

        do p=1, np_local
            ! Velocities at previous 1/2-step to synch with P.E.
            uh(1) = particles(p)%data%v(1)-dt*particles(p)%data%q*particles(p)%results%e(1)/particles(p)%data%m/2.
            uh(2) = particles(p)%data%v(2)-dt*particles(p)%data%q*particles(p)%results%e(2)/particles(p)%data%m/2.
            uh(3) = particles(p)%data%v(3)-dt*particles(p)%data%q*particles(p)%results%e(3)/particles(p)%data%m/2.
            uh2   = dot_product(uh,uh)
            gamma = sqrt(1.0 + uh2/unit_c2)

            en          = particles(p)%data%m*unit_c2*(gamma - 1.0)
            energy(2,p) = en

            if (particles(p)%data%q <= 0.) then
                ekine   = ekine   + en
                sum_v2e = sum_v2e + uh2/gamma**2.
                sum_ve  = sum_ve  + uh/gamma
            else
                ekini   = ekini   + en
                sum_v2i = sum_v2i + uh2/gamma**2.
                sum_vi  = sum_vi  + uh/gamma
            endif
        end do

        ! Gather partial sums together for global energies
        call MPI_ALLREDUCE(MPI_IN_PLACE, ekine, 1, MPI_REAL8, MPI_SUM, MPI_COMM_PEPC, ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE, ekini, 1, MPI_REAL8, MPI_SUM, MPI_COMM_PEPC, ierr)
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

    end subroutine energy_kin



end module module_energy
