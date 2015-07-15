! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2015 Juelich Supercomputing Centre, 
!                         Forschungszentrum Juelich GmbH,
!                         Germany
! 
! PEPC is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! PEPC is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
! 
! You should have received a copy of the GNU Lesser General Public License
! along with PEPC.  If not, see <http://www.gnu.org/licenses/>.
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>  Encapsulates ...
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_energy
    use module_pepc_kinds
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
    subroutine energies(particles,ekine,ekini)

        use physvars
        use module_fmm_framework
        use module_laser
        use module_pusher
        use module_io
        use module_units
        use module_param_dump

        implicit none
        type(t_particle), intent(in) :: particles(:)
        real*8 :: epot, ekine, ekini, etot, totalmomentum(3)
        integer :: ifile
        logical, save :: firstcall = .true.

        call energy_pot(particles, epot)
        call energy_kin(particles, ekine, ekini, tempe, tempi, totalmomentum)
        
        vte = sqrt(3*unit_kB*tempe/mass_e)
        vti = sqrt(3*unit_kB*tempi/mass_i)

        ! rescale total potential and kinetic energy by energy per particle
        ekine        = ekine        / ne
        ekini        = ekini        / ni
        epot         = epot         / npart_total
        potnearfield = potnearfield / npart_total
        potfarfield  = potfarfield  / npart_total

        etot = epot + ekine + ekini

        if (my_rank == 0) then
            do ifile = file_stdout,file_pepc_out,file_pepc_out-file_stdout
              call PrintEnergies(ifile, epot, ekini, ekine, etot, tempe, tempi, totalmomentum)
            end do

            ! Write out to energy.dat file
            if (firstcall .and. .not. restart) then
              firstcall = .false.
              open(file_energy_dat, FILE='energy.dat',STATUS='UNKNOWN', POSITION = 'REWIND')
              write(file_energy_dat,'("#",18(1x,a20))') "time", "Upot(total)", "Upot(near field)", "Upot(far field)", &
                                         "Ukin_e" ,"Ukin_i" ,"Ukin_tot", "Ukine(wo drift)" ,"Ukini(wo drift)" ,"Utot" ,&
                                         "Te0", "Te_uncor", "chie", "delta_Te", "Ti0", "Ti_uncor", "chii", "delta_Ti"
            else
              open(file_energy_dat, FILE='energy.dat',STATUS='UNKNOWN', POSITION = 'APPEND')
            endif

            write(file_energy_dat,'(" ",1(1x,f20.6),17(1x,1pe20.10))') trun*unit_t0_in_fs, epot, potnearfield, potfarfield, &
                                                  ekine, ekini, ekine+ekini, 3./2.*unit_kB*tempe, 3./2.*unit_kB*tempi, etot, &
                                                  Te0, Te_uncor, chie, delta_Te, Ti0, Ti_uncor, chii, delta_Ti

            close(file_energy_dat)
        endif
    end subroutine energies


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Calculate potential energies
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine energy_pot(particles, epot_total)
        use physvars
        use module_fmm_framework
        use module_pepc_types
        implicit none
        include 'mpif.h'

        integer :: ierr

        type(t_particle), intent(in) :: particles(:)
        real*8, intent(out) :: epot_total
        integer(kind_particle) :: p

        epot_total = 0.
        
        do p=1,size(particles, kind=kind(p))
          epot_total = epot_total + 0.5 * particles(p)%data%q*particles(p)%results%pot
        end do
        
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
    subroutine energy_kin(particles,ekine,ekini,ltempe,ltempi, totalmomentum)
        use physvars
        use module_units
        use module_pusher
        use module_pepc_types

        implicit none
        include 'mpif.h'

        type(t_particle), intent(in) :: particles(:)
        integer(kind_particle) :: p
        integer(kind_default) :: ierr
        real*8 :: ekine, ekini,ltempe,ltempi, gamma, totalmomentum(3)
        real*8 :: uh(3), uh2, tmp

        integer, parameter :: V2E  =  1
        integer, parameter :: VEX  =  2
        integer, parameter :: VEY  =  3
        integer, parameter :: VEZ  =  4
        integer, parameter :: V2I  =  5
        integer, parameter :: VIX  =  6
        integer, parameter :: VIY  =  7
        integer, parameter :: VIZ  =  8
        integer, parameter :: KINE =  9
        integer, parameter :: KINI = 10

        real*8 :: sums(1:10)
        real*8 :: betae, betai

        sums = 0.0
        totalmomentum = 0.0
        
        betae = chie/(2*chie-1.)
        betai = chii/(2*chii-1.)

        do p=1, size(particles, kind=kind(p))
            ! Velocities at prev. 1/2-step to synch with P.E.
            select case(integrator_scheme)

              case(INTEGRATOR_SCHEME_NVT)
                if (particles(p)%data%q < 0.) then
                    uh(1:3) = betae*(particles(p)%data%v(1:3) - dt*particles(p)%data%q * particles(p)%results%e(1:3) / particles(p)%data%m / 2.)
                elseif (particles(p)%data%q > 0.) then
                    uh(1:3) = betai*(particles(p)%data%v(1:3) - dt*particles(p)%data%q * particles(p)%results%e(1:3) / particles(p)%data%m / 2.)
                endif
                
              case (INTEGRATOR_SCHEME_NVE_IONS_FROZEN)
                ! unconstrained motion for negatively charged particles, frozen positively charged particles
                if (particles(p)%data%q<0.) then
                    uh(1:3) = particles(p)%data%v(1:3) - dt*particles(p)%data%q * particles(p)%results%e(1:3) / particles(p)%data%m / 2.
                else
                    uh(1:3) = 0.
                endif
                
              case default
                uh(1:3) = particles(p)%data%v(1:3) - dt*particles(p)%data%q * particles(p)%results%e(1:3) / particles(p)%data%m / 2.
            end select
            
            uh2   = dot_product(uh,uh)
            gamma = sqrt(1.0 + uh2/unit_c2)

            tmp = particles(p)%data%m*unit_c2*(gamma - 1.0)

            totalmomentum = totalmomentum + particles(p)%data%m * uh

            if (particles(p)%data%q < 0.) then
                ! electrons
                sums(V2E)     = sums(V2E)     + uh2 / (gamma*gamma)
                sums(VEX:VEZ) = sums(VEX:VEZ) + uh  / gamma
                sums(KINE)    = sums(KINE)    + tmp
            else
                ! ions
                sums(V2I)     = sums(V2I)     + uh2 / (gamma*gamma)
                sums(VIX:VIZ) = sums(VIX:VIZ) + uh  / gamma
                sums(KINI)    = sums(KINI)    + tmp
            endif
        end do

        ! Find global KE sums
        call MPI_ALLREDUCE(MPI_IN_PLACE, sums, size(sums), MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

        sums(V2E:VEZ) = sums(V2E:VEZ) / ne
        sums(V2I:VIZ) = sums(V2I:VIZ) / ni

        ekine = sums(KINE)
        ekini = sums(KINI)

        ltempe =  mass_e/(3.*unit_kB)*(sums(V2E) - dot_product(sums(VEX:VEZ),sums(VEX:VEZ)))
        ltempi =  mass_i/(3.*unit_kB)*(sums(V2I) - dot_product(sums(VIX:VIZ),sums(VIX:VIZ)))

    end subroutine energy_kin



end module module_energy
