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
        use module_geometry
        implicit none

        integer, allocatable :: nfp(:)
        integer, allocatable :: nip(:)
        real(KIND=8), allocatable :: mass(:)
        real(KIND=8), allocatable :: src_t(:)
        real(KIND=8), allocatable :: charge(:)
        logical,allocatable :: physical_particle(:)
        character(255),allocatable :: name(:)
        integer, allocatable :: src_type(:),src_bnd(:)
        real(KIND=8), allocatable :: src_x0(:,:)
        real(KIND=8), allocatable :: src_e1(:,:),src_e2(:,:),src_e3(:,:)
        real(KIND=8), allocatable :: src_v0(:)


        integer :: ns,ns_max
        integer :: rc,ispecies,fid=12

        namelist /species_nml/ src_t,ns,nip,nfp,mass,charge,physical_particle,name,src_type,src_bnd,src_x0,src_e1,src_e2,src_e3,src_v0

        ns=0
        ns_max=1000

                allocate(nfp(0:ns_max),stat=rc)
        allocate(nip(0:ns_max),stat=rc)
        allocate(mass(0:ns_max),stat=rc)
        allocate(src_t(0:ns_max),stat=rc)
        allocate(charge(0:ns_max),stat=rc)
        allocate(physical_particle(0:ns_max),stat=rc)
        allocate(name(0:ns_max),stat=rc)
        allocate(src_type(0:ns_max),stat=rc)
        allocate(src_bnd(0:ns_max),stat=rc)
        allocate(src_x0(0:ns_max,3),stat=rc)
        allocate(src_e1(0:ns_max,3),stat=rc)
        allocate(src_e2(0:ns_max,3),stat=rc)
        allocate(src_e3(0:ns_max,3),stat=rc)
        allocate(src_v0(0:ns_max),stat=rc)


        IF(root) write(*,'(a,a)') " == reading parameter file, section species: ", trim(input_file)
        open(fid,file=trim(input_file))
        read(fid,NML=species_nml)
        rewind(fid)


        deallocate(nfp)
        deallocate(mass)
        deallocate(charge)
        deallocate(physical_particle)
        deallocate(name)
        deallocate(src_type)
        deallocate(src_bnd)
        deallocate(src_x0)
        deallocate(src_e1)
        deallocate(src_e2)
        deallocate(src_e3)
        deallocate(src_v0)


        allocate(nfp(0:ns-1),stat=rc)
        allocate(nip(0:ns-1),stat=rc)
        allocate(mass(0:ns-1),stat=rc)
        allocate(src_t(0:ns-1),stat=rc)
        allocate(charge(0:ns-1),stat=rc)
        allocate(physical_particle(0:ns-1),stat=rc)
        allocate(name(0:ns-1),stat=rc)
        allocate(src_type(0:ns-1),stat=rc)
        allocate(src_bnd(0:ns-1),stat=rc)
        allocate(src_x0(0:ns-1,3),stat=rc)
        allocate(src_e1(0:ns-1,3),stat=rc)
        allocate(src_e2(0:ns-1,3),stat=rc)
        allocate(src_e3(0:ns-1,3),stat=rc)
        allocate(src_v0(0:ns-1),stat=rc)
        nfp=0
        nip=0
        mass=0.
        charge=0.
        physical_particle=.false.
        name=""
        src_t=0.
        src_type=0
        src_bnd=0
        src_v0=0
        src_x0=0
        src_e1=0
        src_e2=0
        src_e3=0

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
            species(ispecies)%src_t=src_t(ispecies)
            species(ispecies)%q=charge(ispecies)
            species(ispecies)%indx=ispecies
            species(ispecies)%physical_particle=physical_particle(ispecies)
            species(ispecies)%nfp=nfp(ispecies)
            species(ispecies)%nip=nip(ispecies)
            tnpps(ispecies)=nip(ispecies)

            IF (species(ispecies)%physical_particle) THEN
                IF (src_type(ispecies)==0) THEN !surface source (whole surface)
                    src_x0(ispecies,:)=0.
                    src_e1(ispecies,:)=0.
                    src_e2(ispecies,:)=0.
                    src_e3(ispecies,:)=0.
                    src_v0(ispecies)=0.
                    IF ((src_bnd(ispecies)<=0).or.(src_bnd(ispecies)>nb)) THEN
                        IF (root) write(*,'(a)') "You have to select one of the boundaries as surface source"
                        STOP
                    END IF
                    IF(boundaries(src_bnd(ispecies))%type==2) THEN
                        IF (root) write(*,'(a)') "Periodic boundary cannot be used as surface source"
                        STOP
                    END IF
                    IF (root) write(*,'(a,i3,a,i3,a,i3)') " == Boundary ",src_bnd(ispecies)," chosen as surface source of type "&
                                                          ,src_type(ispecies)," for species ",ispecies

                ELSE IF (src_type(ispecies)==5) THEN !surface source, exactly like type=0, but with a flux from a drifting Maxwellian
                    src_x0(ispecies,:)=0.
                    src_e1(ispecies,:)=0.
                    src_e2(ispecies,:)=0.
                    src_e3(ispecies,:)=0.
                    IF ((src_bnd(ispecies)<=0).or.(src_bnd(ispecies)>nb)) THEN
                        IF (root) write(*,'(a)') "You have to select one of the boundaries as surface source"
                        STOP
                    END IF
                    IF(boundaries(src_bnd(ispecies))%type==2) THEN
                        IF (root) write(*,'(a)') "Periodic boundary cannot be used as surface source"
                        STOP
                    END IF
                    IF (root) write(*,'(a,i3,a,i3,a,i3)') " == Boundary ",src_bnd(ispecies)," chosen as surface source of type "&
                                                          ,src_type(ispecies)," for species ",ispecies

                ELSE IF (src_type(ispecies)==4) THEN !surface source, circle at bnd_x0 + src_x0(1) * e1 +
                                                     !src_x0(2) * e2 with r = src_x0(3)*|e1|
                                                     !this is a first test and only works if |e1| ~ |e2|
                                                     !and src_x0 reasonable
                    src_e1(ispecies,:)=0.
                    src_e2(ispecies,:)=0.
                    src_e3(ispecies,:)=0.
                    src_v0(ispecies)=0.
                    IF ((src_bnd(ispecies)<=0).or.(src_bnd(ispecies)>nb)) THEN
                        IF (root) write(*,'(a)') "You have to select one of the boundaries as surface source"
                        STOP
                    END IF
                    IF(boundaries(src_bnd(ispecies))%type==2) THEN
                        IF (root) write(*,'(a)') "Periodic boundary cannot be used as surface source"
                        STOP
                    END IF
                    IF (root) write(*,'(a,i3,a,i3,a,i3)') " == Boundary ",src_bnd(ispecies)," chosen as surface source of type "&
                                                          ,src_type(ispecies)," for species ",ispecies

                ELSE IF ((src_type(ispecies)==1).or.(src_type(ispecies)==2).or.(src_type(ispecies)==3)) THEN !Volume Source
                    src_bnd(ispecies)=0
                    src_v0(ispecies)=0.
                    IF (root) write(*,'(a,i2,a,i2,a)') " == Volume source of type ",src_type(ispecies)," for species ",ispecies," set. Parameters:"
                    IF (root) write(*,'(a,3(1pe14.5E3))') " == x0: ",src_x0(ispecies,:)
                    IF (root) write(*,'(a,3(1pe14.5E3))') " == e1: ",src_e1(ispecies,:)
                    IF (root) write(*,'(a,3(1pe14.5E3))') " == e2: ",src_e2(ispecies,:)
                    IF (root) write(*,'(a,3(1pe14.5E3))') " == e3: ",src_e3(ispecies,:)
                ELSE
                    IF (root) write(*,'(a,i3,a)') " Source cannot be set. Type ",src_type(ispecies)," not available."
                    STOP
                END IF
            ELSE
                src_x0(ispecies,:)=0.
                src_e1(ispecies,:)=0.
                src_e2(ispecies,:)=0.
                src_e3(ispecies,:)=0.
                src_type(ispecies)=0
                src_bnd(ispecies)=0
                src_v0(ispecies)=0.
            END IF
            species(ispecies)%src_type=src_type(ispecies)
            species(ispecies)%src_bnd=src_bnd(ispecies)
            species(ispecies)%src_x0=src_x0(ispecies,:)
            species(ispecies)%src_e1=src_e1(ispecies,:)
            species(ispecies)%src_e2=src_e2(ispecies,:)
            species(ispecies)%src_e3=src_e3(ispecies,:)
            species(ispecies)%src_v0=src_v0(ispecies)
        END DO

        allocate(ehit_max(0:nspecies-1),stat=rc)
        allocate(energy_resolved_hits(0:nspecies-1,nb,nbins_energy_resolved_hits+1),stat=rc)
        allocate(angle_resolved_hits(0:nspecies-1,nb,nbins_angle_resolved_hits),stat=rc)
        energy_resolved_hits = 0
        angle_resolved_hits = 0
        ehit_max=0.0_8
        DO ispecies=0,nspecies-1
            IF (species(ispecies)%physical_particle) THEN
                ehit_max(ispecies) = 8. * species(ispecies)%src_t
            END IF
        END DO

        call init_maxw_flux_tables(1000)
        tnpps(0)=count_wallparticles()
        call check_species()

        allocate(probe_start_x(0:nspecies-1),stat=rc)
        allocate(probe_start_y(0:nspecies-1),stat=rc)
        allocate(probe_start_z(0:nspecies-1),stat=rc)
        allocate(probe_end_x(0:nspecies-1),stat=rc)
        allocate(probe_end_y(0:nspecies-1),stat=rc)
        allocate(probe_end_z(0:nspecies-1),stat=rc)

        deallocate(nfp)
        deallocate(mass)
        deallocate(charge)
        deallocate(physical_particle)
        deallocate(name)
        deallocate(src_type)
        deallocate(src_bnd)
        deallocate(src_x0)
        deallocate(src_e1)
        deallocate(src_e2)
        deallocate(src_e3)

    end subroutine init_species


!==================================================================================
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> initialize tables for sampling of random numbers from a drifting
    !> Maxwellian flux
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE init_maxw_flux_tables(N)
        use helper
        implicit none

        integer,intent(in) :: N
        integer :: rc,i,ispecies
        real(KIND=8) :: v0,numer,denomi,v,vth,sqrt2

        allocate(maxw_flux_table_F(0:nspecies-1,N),stat=rc)
        allocate(maxw_flux_table_v(0:nspecies-1,N),stat=rc)
        maxw_flux_table_v = 0._8
        maxw_flux_table_F = 0._8

        sqrt2=sqrt(2._8)

        DO ispecies=0,nspecies-1
            IF (species(ispecies)%src_type == 5) THEN
                vth = sqrt(species(ispecies)%src_t*e/species(ispecies)%m)
                v0 = species(ispecies)%src_v0
                call linspace(0._8,7._8*max(v0,vth),maxw_flux_table_v(ispecies,:))
                DO i=1,N
                    v = maxw_flux_table_v(ispecies,i)
                    numer=exp(-(v0/(sqrt2*vth))**2)-exp(-((v-v0)/(sqrt2*vth))**2)+sqrt(pi)*v0/(sqrt2*vth)*(erf((v-v0)/(sqrt2*vth))+erf(v0/(sqrt2*vth)))
                    denomi = exp(-(v0/(sqrt2*vth))**2) + sqrt(pi)*v0/(sqrt2*vth)*(1+erf(v0/(sqrt2*vth)))
                    maxw_flux_table_F(ispecies,i) = numer / denomi
                END DO
            END IF
        END DO


    end subroutine
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
            IF ((species(ispecies)%physical_particle .eqv. .false.) .and. (species(ispecies)%src_t/=0.)) THEN
                IF (root) write(*,'(a,i3,a)') " Source Temperature cannot be set for a nonphysical species (species ",ispecies," )."
                STOP
            END IF
        END DO

    end subroutine check_species

!======================================================================================


end module module_species
