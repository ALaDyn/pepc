! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
!
! Copyright (C) 2002-2012 Juelich Supercomputing Centre,
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!> Contains types for building the sim domain form walls
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_species
    use variables
    use module_cmdline
    implicit none

    contains

!======================================================================================

    subroutine init_species()
        implicit none

        integer, allocatable :: nfp(:)
        integer, allocatable :: nip(:)
        real(KIND=8), allocatable :: mass(:)
        real(KIND=8), allocatable :: charge(:)
        logical,allocatable :: physical_particle(:)
        character(255),allocatable :: name(:)

        integer :: ns,ns_max
        integer :: rc,ispecies,fid=12

        namelist /species_nml/ ns,nip,nfp,mass,charge,physical_particle,name

        ns_max=100

        allocate(nfp(0:ns_max),stat=rc)
        allocate(nip(0:ns_max),stat=rc)
        allocate(mass(0:ns_max),stat=rc)
        allocate(charge(0:ns_max),stat=rc)
        allocate(physical_particle(0:ns_max),stat=rc)
        allocate(name(0:ns_max),stat=rc)

        ns=0
        nfp=0
        nip=0
        mass=0.
        charge=0.
        physical_particle=.false.
        name=""


        IF(root) write(*,'(a,a)') " == reading parameter file, section species: ", trim(input_file)
        open(fid,file=trim(input_file))
        read(fid,NML=species_nml)
        close(fid)

        IF(ns<=0) THEN
            IF (root) write(*,'(a)') " number of species not set or set to invalid value "
            STOP
        ELSE
            IF (root) write(*,'(a,i3,a)') " == initializing ",ns," species"
        END IF

        nspecies=ns
        allocate(species(0:nspecies-1),stat=rc)

        allocate(tnpps(0:nspecies-1),stat=rc)
        allocate(npps(0:nspecies-1),stat=rc)
        tnpps=0
        npps=0

        DO ispecies=0,nspecies-1
            species(ispecies)%name=trim(name(ispecies))
            species(ispecies)%m=mass(ispecies)
            species(ispecies)%q=charge(ispecies)
            species(ispecies)%indx=ispecies
            species(ispecies)%physical_particle=physical_particle(ispecies)
            species(ispecies)%nfp=nfp(ispecies)
            species(ispecies)%nip=nip(ispecies)
            tnpps(ispecies)=nip(ispecies)
        END DO

        call check_species()

        deallocate(nfp)
        deallocate(mass)
        deallocate(charge)
        deallocate(physical_particle)
        deallocate(name)

    end subroutine init_species


!======================================================================================

    subroutine check_species()
        implicit none

        integer :: ispecies

        IF (species(0)%physical_particle) THEN
            IF (root) write(*,'(a)') " Species 0 have to be wallparticles. physical_particle cannot be .true."
            STOP
        END IF

        IF ((species(0)%nip/=0)) THEN
            IF (root) write(*,'(a)') " Number of initial particles cannot be set for wallparticles."
            IF (root) write(*,'(a)') " The number of wallparticles is set for each wall seperately in the geometry input."
            STOP
        END IF

        DO ispecies=0,nspecies-1
            IF ((species(ispecies)%physical_particle .eqv. .false.) .and. (species(ispecies)%nfp/=0)) THEN
                IF (root) write(*,'(a,i3,a)') " A flux cannot be set for a nonphysical species (species ",ispecies," )."
                STOP
            END IF
        END DO

    end subroutine check_species

!======================================================================================


end module module_species
