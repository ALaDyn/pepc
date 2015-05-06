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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!> Contains types for building the sim domain form walls
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_species
    use variables
    use module_cmdline
    use poisson_disc_sampling, ONLY: poisson_sampler
    implicit none

    contains

!======================================================================================

    subroutine init_species()
        use module_geometry
        use module_species_types
        implicit none
        include 'mpif.h'

        integer, allocatable :: nfp(:)
        integer, allocatable :: nip(:)
        real(KIND=8), allocatable :: mass(:)
        real(KIND=8), allocatable :: src_t(:)
        real(KIND=8), allocatable :: charge(:)
        integer,allocatable :: physical_particle(:)
        character(255),allocatable :: name(:)
        integer, allocatable :: src_type_x(:),src_bnd(:),src_type_v(:)
        real(KIND=8), allocatable :: src_x0(:,:)
        real(KIND=8), allocatable :: src_e1(:,:),src_e2(:,:),src_e3(:,:)
        real(KIND=8), allocatable :: src_v0(:)

        integer :: ns,ns_max
        integer :: rc,ispecies,irank,fid=12

        integer, allocatable :: npps_per_rank(:,:), displs(:,:)
        real(KIND=8), allocatable :: starting_positions_x(:),starting_positions_y(:),starting_positions_z(:)
        real(KIND=8), allocatable :: l_starting_positions_x(:),l_starting_positions_y(:),l_starting_positions_z(:)
        real(KIND=8), allocatable :: l_starting_positions(:,:)
        real(KIND=8) :: DomainMin(1:3), DomainMax(1:3), DomainCorners(3,8)
        real(KIND=8) :: eps=1.0e-10
        logical :: source_geometry_ok

        namelist /species_nml/ src_t,ns,nip,nfp,mass,charge,physical_particle,name,src_type_x,src_type_v,src_bnd,src_x0,src_e1,src_e2,src_e3,src_v0

        source_geometry_ok = .true.

        ns=0
        ns_max=1000

        allocate(nfp(0:ns_max),stat=rc)
        allocate(nip(0:ns_max),stat=rc)
        allocate(mass(0:ns_max),stat=rc)
        allocate(src_t(0:ns_max),stat=rc)
        allocate(charge(0:ns_max),stat=rc)
        allocate(physical_particle(0:ns_max),stat=rc)
        allocate(name(0:ns_max),stat=rc)
        allocate(src_type_x(0:ns_max),stat=rc)
        allocate(src_type_v(0:ns_max),stat=rc)
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
        deallocate(src_type_x)
        deallocate(src_type_v)
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
        allocate(src_type_x(0:ns-1),stat=rc)
        allocate(src_type_v(0:ns-1),stat=rc)
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
        physical_particle=0
        name=""
        src_t=0.
        src_type_x=0
        src_type_v=0
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
        allocate(npps_per_rank(0:n_ranks-1,0:nspecies-1))
        allocate(displs(0:n_ranks-1,0:nspecies-1))
        tnpps = 0
        npps = 0
        npps_per_rank(:,:) = 0
        displs(:,:) = 0

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
            IF ((species(ispecies)%physical_particle == 1) .OR. (species(ispecies)%physical_particle == 3)) THEN
                species(ispecies)%moving_particle = .TRUE.
            ELSE
                species(ispecies)%moving_particle = .FALSE.
            END IF

            ! initially set number of local particles
            IF (ispecies > 0) THEN
                IF (species(ispecies)%moving_particle) THEN
                    npps(ispecies) = tnpps(ispecies) / n_ranks
                    if(my_rank.eq.(n_ranks-1)) npps(ispecies) = npps(ispecies) + MOD(tnpps(ispecies), n_ranks)
                ELSE !probes only on root (will be moved to othe ranks in grow_tree anyway)
                    IF (my_rank == 0) THEN
                        npps(ispecies) = tnpps(ispecies)
                    ELSE
                        npps(ispecies) = 0
                    END IF
                END IF
            END IF

            call MPI_ALLGATHER(npps(ispecies), 1, MPI_INTEGER, npps_per_rank(:, ispecies), 1, MPI_INTEGER, MPI_COMM_WORLD,rc)

            IF (species(ispecies)%moving_particle) THEN
                species(ispecies)%v_th = sqrt(species(ispecies)%src_t * e / species(ispecies)%m)
                IF (src_type_x(ispecies) == 1) THEN !surface source (whole surface)
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
                    IF (root) write(*,'(a,i3,a,i3,a,i3)') " == Boundary ",src_bnd(ispecies)," chosen as surface source for species ",ispecies

                ELSE IF (src_type_x(ispecies)==2) THEN !Volume Source, uniformly distributed
                    src_bnd(ispecies)=0
                    IF (root) write(*,'(a,i2,a,i2,a)') " == Volume source for species ",ispecies," set. Parameters:"
                    IF (root) write(*,'(a,3(1pe14.5E3))') " == x0: ",src_x0(ispecies,:)
                    IF (root) write(*,'(a,3(1pe14.5E3))') " == e1: ",src_e1(ispecies,:)
                    IF (root) write(*,'(a,3(1pe14.5E3))') " == e2: ",src_e2(ispecies,:)
                    IF (root) write(*,'(a,3(1pe14.5E3))') " == e3: ",src_e3(ispecies,:)

                ELSE IF (src_type_x(ispecies)==3) THEN !Volume Source, poisson disc sampling
                    src_bnd(ispecies)=0
                    IF (root) write(*,'(a,i2,a,i2,a)') " == Volume source with Poisson disc sampling for species ",ispecies," set. Parameters:"
                    IF (root) write(*,'(a,3(1pe14.5E3))') " == x0: ",src_x0(ispecies,:)
                    IF (root) write(*,'(a,3(1pe14.5E3))') " == e1: ",src_e1(ispecies,:)
                    IF (root) write(*,'(a,3(1pe14.5E3))') " == e2: ",src_e2(ispecies,:)
                    IF (root) write(*,'(a,3(1pe14.5E3))') " == e3: ",src_e3(ispecies,:)

                    IF ( (my_rank == ispecies) .OR. ((my_rank == 0) .AND. (n_ranks-1 < ispecies)) ) THEN
                        allocate(starting_positions_x(species(ispecies)%nip))
                        allocate(starting_positions_y(species(ispecies)%nip))
                        allocate(starting_positions_z(species(ispecies)%nip))
                        IF (real_unequal_zero(dotproduct(src_e1(ispecies,:),src_e2(ispecies,:)), eps)) source_geometry_ok = .false.
                        IF (real_unequal_zero(dotproduct(src_e2(ispecies,:),src_e3(ispecies,:)), eps)) source_geometry_ok = .false.
                        IF (real_unequal_zero(dotproduct(src_e3(ispecies,:),src_e1(ispecies,:)), eps)) source_geometry_ok = .false.

                        IF (source_geometry_ok .eqv. .false.) THEN
                            write(*,*) "Poisson disc sampling routine can only be used cuboid shaped source regions"
                            STOP
                        END IF

                        DomainCorners(:,1) = src_x0(ispecies,:)
                        DomainCorners(:,2) = src_x0(ispecies,:) + src_e1(ispecies,:)
                        DomainCorners(:,3) = src_x0(ispecies,:) + src_e2(ispecies,:)
                        DomainCorners(:,4) = src_x0(ispecies,:) + src_e3(ispecies,:)
                        DomainCorners(:,5) = src_x0(ispecies,:) + src_e1(ispecies,:) + src_e2(ispecies,:)
                        DomainCorners(:,6) = src_x0(ispecies,:) + src_e1(ispecies,:) + src_e3(ispecies,:)
                        DomainCorners(:,7) = src_x0(ispecies,:) + src_e2(ispecies,:) + src_e3(ispecies,:)
                        DomainCorners(:,8) = src_x0(ispecies,:) + src_e1(ispecies,:) + src_e2(ispecies,:) + src_e3(ispecies,:)

                        DomainMin(1) = minval(DomainCorners(1,:))
                        DomainMin(2) = minval(DomainCorners(2,:))
                        DomainMin(3) = minval(DomainCorners(3,:))
                        DomainMax(1) = maxval(DomainCorners(1,:))
                        DomainMax(2) = maxval(DomainCorners(2,:))
                        DomainMax(3) = maxval(DomainCorners(3,:))

                        call poisson_sampler(DomainMin, DomainMax, starting_positions_x, starting_positions_y, starting_positions_z)
                    END IF

                    allocate(l_starting_positions_x(npps(ispecies)))
                    allocate(l_starting_positions_y(npps(ispecies)))
                    allocate(l_starting_positions_z(npps(ispecies)))
                    allocate(l_starting_positions(3,npps(ispecies)))

                    DO irank = 0, n_ranks-1
                        displs(irank, ispecies) = sum(npps_per_rank(0:irank-1, ispecies))
                    END DO

                    IF (n_ranks-1 < ispecies) THEN
                        irank = 0
                    ELSE
                        irank = ispecies
                    END IF



                    call MPI_SCATTERV(starting_positions_x(:), npps_per_rank(:,ispecies), displs(:, ispecies), MPI_REAL8, &
                                      l_starting_positions_x(:), npps(ispecies) , MPI_REAL8, irank, MPI_COMM_WORLD,rc)
                    call MPI_SCATTERV(starting_positions_y(:), npps_per_rank(:,ispecies), displs(:, ispecies), MPI_REAL8, &
                                      l_starting_positions_y(:), npps(ispecies) , MPI_REAL8, irank, MPI_COMM_WORLD,rc)
                    call MPI_SCATTERV(starting_positions_z(:), npps_per_rank(:,ispecies), displs(:, ispecies), MPI_REAL8, &
                                      l_starting_positions_z(:), npps(ispecies) , MPI_REAL8, irank, MPI_COMM_WORLD,rc)


                    l_starting_positions(1,:) = l_starting_positions_x
                    l_starting_positions(2,:) = l_starting_positions_y
                    l_starting_positions(3,:) = l_starting_positions_z
                    allocate(species(ispecies)%starting_positions(3,npps(ispecies)))
                    species(ispecies)%starting_positions = l_starting_positions


                    IF (allocated(starting_positions_x)) deallocate(starting_positions_x)
                    IF (allocated(l_starting_positions_x)) deallocate(l_starting_positions_x)
                    IF (allocated(starting_positions_y)) deallocate(starting_positions_y)
                    IF (allocated(l_starting_positions_y)) deallocate(l_starting_positions_y)
                    IF (allocated(starting_positions_z)) deallocate(starting_positions_z)
                    IF (allocated(l_starting_positions_z)) deallocate(l_starting_positions_z)
                    IF (allocated(l_starting_positions)) deallocate(l_starting_positions)


                ELSE IF (src_type_x(ispecies)==11) THEN !surface source, circle at bnd_x0 + src_x0(1) * e1 +
                                                        !src_x0(2) * e2 with r = src_x0(3)*|e1|
                                                        !this is a first test and only works if |e1| ~ |e2|
                                                        !and src_x0 reasonable
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
                    IF (root) write(*,'(a,i3,a,i3,a,i3)') " == Boundary ",src_bnd(ispecies)," chosen as cylindrical surface source for species ",ispecies

                ELSE IF (src_type_x(ispecies)==12) THEN !Cylindrical Volume Source
                    src_bnd(ispecies)=0
                    src_e2(ispecies,:)=0.
                    src_e3(ispecies,:)=0.
                    IF (root) write(*,'(a,i2,a,i2,a)') " == Cylindrical Volume source of type for species ",ispecies," set. Parameters:"
                    IF (root) write(*,'(a,3(1pe14.5E3))') " == reference point(x0): ",src_x0(ispecies,:)
                    IF (root) write(*,'(a,3(1pe14.5E3))') " == axis(e1): ",src_e1(ispecies,:)
                    IF (root) write(*,'(a,(1pe14.5E3))') " == radius(v0): ",src_v0(ispecies)
                ELSE
                    IF (root) write(*,'(a,i3,a)') " Source cannot be set. Type ",src_type_x(ispecies)," not available."
                    STOP
                END IF

                IF (src_type_v(ispecies) == 1) THEN !Emmert Source along x, Maxwellian along y,z
                    src_v0(ispecies) = 0.
                ELSE IF (src_type_v(ispecies) == 2) THEN !Emmert source along B (along x if B=0)
                    src_v0(ispecies) = 0.
                ELSE IF (src_type_v(ispecies) == 3) THEN !Emmert source along cylinder axis for cylindrical volume source (along x if other src_type_x)
                    src_v0(ispecies) = 0.
                    IF (src_type_x(ispecies) /= 12) src_type_v(ispecies) = 2
                ELSE IF (src_type_v(ispecies) == 4) THEN !Emmert source along surface normal for surface source (along x if other src_type_x)
                    src_v0(ispecies) = 0.
                    IF ( (src_type_x(ispecies) /= 1) .AND. (src_type_x(ispecies) /= 11) )  src_type_v(ispecies) = 2
                ELSE IF (src_type_v(ispecies) == 5) THEN !Maxwellian in all directions (Bissel-Johnson)
                    src_v0(ispecies) = 0.
                ELSE IF (src_type_v(ispecies) == 6) THEN !Maxwellian flux along x
                    src_v0(ispecies) = 0.
                ELSE IF (src_type_v(ispecies) == -6) THEN !Maxwellian flux along -x
                    src_v0(ispecies) = 0.
                ELSE IF (src_type_v(ispecies) == 7) THEN !Maxwellian flux along B (along x if B=0)
                    src_v0(ispecies) = 0.
                ELSE IF (src_type_v(ispecies) == -7) THEN !Maxwellian flux along -B (along -x if B=0)
                    src_v0(ispecies) = 0.
                ELSE IF (src_type_v(ispecies) == 8) THEN !Maxwellian flux along cylinder axis for cylindrical volume source (along x if other src_type_x)
                    src_v0(ispecies) = 0.
                    IF (src_type_x(ispecies) /= 12) src_type_v(ispecies) = 6
                ELSE IF (src_type_v(ispecies) == -8) THEN !Maxwellian flux along -1 * cylinder axis for cylindrical volume source (along -x if other src_type_x)
                    src_v0(ispecies) = 0.
                    IF (src_type_x(ispecies) /= 12) src_type_v(ispecies) = -6
                ELSE IF (src_type_v(ispecies) == 9) THEN !Maxwellian flux along surface normal for surface source (along x if other src_type_x)
                    src_v0(ispecies) = 0.
                    IF ( (src_type_x(ispecies) /= 1) .AND. (src_type_x(ispecies) /= 11) )  src_type_v(ispecies) = 6
                ELSE IF (src_type_v(ispecies) == -9) THEN !Maxwellian flux along -1 * surface normal for surface source (along -x if other src_type_x)
                    src_v0(ispecies) = 0.
                    IF ( (src_type_x(ispecies) /= 1) .AND. (src_type_x(ispecies) /= 11) )  src_type_v(ispecies) = -6
                ELSE IF (src_type_v(ispecies) == 10) THEN !drifting Maxwellian flux along x
                ELSE IF (src_type_v(ispecies) == -10) THEN !drifting Maxwellian flux along -x
                ELSE IF (src_type_v(ispecies) == 11) THEN !drifting Maxwellian flux along B (along x if B=0)
                ELSE IF (src_type_v(ispecies) == -11) THEN !drifting Maxwellian flux along -B (along -x if B=0)
                ELSE IF (src_type_v(ispecies) == 12) THEN !drifting Maxwellian flux along cylinder axis for cylindrical volume source (along x if other src_type_x)
                    IF (src_type_x(ispecies) /= 12) src_type_v(ispecies) = 10
                ELSE IF (src_type_v(ispecies) == -12) THEN !drifting Maxwellian flux along -1 * cylinder axis for cylindrical volume source (along x if other src_type_x)
                    IF (src_type_x(ispecies) /= 12) src_type_v(ispecies) = -10
                ELSE IF (src_type_v(ispecies) == 13) THEN !drifting Maxwellian flux along surface normal for surface source (along x if other src_type_x)
                    IF ( (src_type_x(ispecies) /= 1) .AND. (src_type_x(ispecies) /= 11) )  src_type_v(ispecies) = 10
                ELSE IF (src_type_v(ispecies) == -13) THEN !drifting Maxwellian flux along -1 * surface normal for surface source (along -x if other src_type_x)
                    IF ( (src_type_x(ispecies) /= 1) .AND. (src_type_x(ispecies) /= 11) )  src_type_v(ispecies) = -10
                ELSE
                    IF (root) write(*,'(a,i3,a)') " Source cannot be set. Type ",src_type_v(ispecies)," not available."
                    STOP
                END IF
            ELSE
                species(ispecies)%v_th = 0.0_8
                src_x0(ispecies,:)=0.
                src_e1(ispecies,:)=0.
                src_e2(ispecies,:)=0.
                src_e3(ispecies,:)=0.
                src_type_x(ispecies)=0
                src_type_v(ispecies)=0
                src_bnd(ispecies)=0
                src_v0(ispecies)=0.
            END IF
            species(ispecies)%src_type_x=src_type_x(ispecies)
            species(ispecies)%src_type_v=src_type_v(ispecies)
            species(ispecies)%src_bnd=src_bnd(ispecies)
            species(ispecies)%src_x0=src_x0(ispecies,:)
            species(ispecies)%src_e1=src_e1(ispecies,:)
            species(ispecies)%src_e2=src_e2(ispecies,:)
            species(ispecies)%src_e3=src_e3(ispecies,:)
            species(ispecies)%src_v0=src_v0(ispecies)
        END DO


        IF (bool_energy_resolved_hits) THEN
            allocate(ehit_max(0:nspecies-1),stat=rc)
            allocate(energy_resolved_hits(0:nspecies-1,nb,nbins_energy_resolved_hits+1),stat=rc)
            energy_resolved_hits = 0
            ehit_max=0.0_8
            DO ispecies=0,nspecies-1
                IF (species(ispecies)%moving_particle) THEN
                    ehit_max(ispecies) = ehit_max_in_T * species(ispecies)%src_t
                END IF
            END DO
        END IF
        IF (bool_age_resolved_hits) THEN
            allocate(agehit_max(0:nspecies-1),stat=rc)
            allocate(age_resolved_hits(0:nspecies-1,nb,nbins_age_resolved_hits+1),stat=rc)
            age_resolved_hits = 0
            agehit_max=0.0_8
            DO ispecies=0,nspecies-1
                IF (species(ispecies)%moving_particle) THEN
                    agehit_max(ispecies) = agehit_max_in_t_trav_ion * (dx/2) / species(2)%v_th  !maximum is set to 2 * ion_traversal_time
                END IF
            END DO
        END IF
        IF (bool_angle_resolved_hits) THEN
            allocate(angle_resolved_hits(0:nspecies-1,nb,nbins_angle_resolved_hits),stat=rc)
            angle_resolved_hits = 0
        END IF
        IF (bool_space_resolved_hits) THEN
            allocate(space_resolved_hits(0:nspecies-1,nb,nbins_e1_space_resolved_hits,nbins_e2_space_resolved_hits),stat=rc)
            space_resolved_hits = 0
        END IF


        call init_maxw_flux_tables(1000)
        tnpps(0)=count_wallparticles()
        ! initially set number of local wall particles
        npps(0) = tnpps(0) / n_ranks
        if(my_rank.eq.(n_ranks-1)) npps(0) = npps(0) + MOD(tnpps(0), n_ranks)

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
        deallocate(src_type_x)
        deallocate(src_type_v)
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
            IF ((species(ispecies)%src_type_v == 10) .OR. (species(ispecies)%src_type_v == -10) .OR. &
                (species(ispecies)%src_type_v == 11) .OR. (species(ispecies)%src_type_v == -11) .OR. &
                (species(ispecies)%src_type_v == 12) .OR. (species(ispecies)%src_type_v == -12) .OR. &
                (species(ispecies)%src_type_v == 13) .OR. (species(ispecies)%src_type_v == -13)) THEN
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

    END SUBROUTINE

!======================================================================================

    subroutine check_species()
        implicit none

        integer :: ispecies

        IF (species(0)%physical_particle /= 0) THEN
            IF (root) write(*,'(a)') " Species 0 have to be wallparticles. physical_particle has to be 0"
            STOP
        END IF

        IF ((species(0)%nip/=0)) THEN
            IF (root) write(*,'(a)') " Number of initial particles cannot be set for wallparticles."
            IF (root) write(*,'(a)') " The number of wallparticles is set for each wall seperately in the geometry input."
            STOP
        END IF

        DO ispecies=0,nspecies-1
            IF ((.NOT.(species(ispecies)%moving_particle)) .and. (species(ispecies)%nfp/=0)) THEN
                IF (root) write(*,'(a,i3,a)') " A flux cannot be set for a nonphysical species (species ",ispecies," )."
                STOP
            END IF
            IF ( (.NOT.(species(ispecies)%moving_particle)) .and. ((species(ispecies)%src_t > 0.).or.(species(ispecies)%src_t < 0.)) ) THEN
                IF (root) write(*,'(a,i3,a)') " Source Temperature cannot be set for a nonphysical species (species ",ispecies," )."
                STOP
            END IF
        END DO

    end subroutine check_species

!======================================================================================


end module module_species
