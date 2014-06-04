module particlehandling

    use module_pepc_types
    use helper
    use module_geometry
    use zufall
    use output
    implicit none

    contains

!======================================================================================
    SUBROUTINE count_hits_and_remove_particles(p,hits,reflux)
        implicit none
        include 'mpif.h'

        logical :: hit
        integer :: rp,ib
        integer :: rc
        type(t_particle), allocatable, intent(inout) :: p(:)
        type(t_particle), allocatable :: p_hits_logical_sheath(:,:,:)
        integer, intent(inout) :: hits(0:,:),reflux(0:,:)
        real*8 :: xp(3),xold(3)
        real(KIND=8) v_new(3),mu,sigma,ran,t1(3),t2(3),xi(3)


        IF (ANY(boundaries(1:nb)%type==4)) allocate(p_hits_logical_sheath(0:nspecies-1,nb,sum(npps)),stat=rc)

        !write(bnd_hits_out,*) "Timestep:",step

        rp=sum(npps)
        ib=1
        DO WHILE (rp >= 1)
            ib = 1
            DO WHILE (ib <= nb)
                IF (rp==0) EXIT
                call hit_wall(p(rp),boundaries(ib),hit)

                IF (hit) THEN
                    hits(p(rp)%data%species,ib)=hits(p(rp)%data%species,ib)+1
                    IF (boundaries(ib)%type==3) THEN                                   !Open BC
                        IF (boundaries(ib)%reflux_particles) reflux(p(rp)%data%species,ib)=reflux(p(rp)%data%species,ib)+1
                        npps(p(rp)%data%species) = npps(p(rp)%data%species) - 1
                        p(rp) = p(sum(npps)+1)
                        ib = 0
                        rp = rp - 1
                    ELSE IF (boundaries(ib)%type==4) THEN                              !logical sheath
                        IF (boundaries(ib)%reflux_particles) reflux(p(rp)%data%species,ib)=reflux(p(rp)%data%species,ib)+1
                        p_hits_logical_sheath(p(rp)%data%species,ib,hits(p(rp)%data%species,ib)) = p(rp)
                        npps(p(rp)%data%species) = npps(p(rp)%data%species) - 1
                        p(rp) = p(sum(npps)+1)
                        ib = 0
                        rp = rp - 1
                    ELSE IF (boundaries(ib)%type==5) THEN   !immediate half-Maxwellian refluxing normal to surface, tangential v conserved
                        mu=0.0_8
                        sigma=sqrt(species(p(rp)%data%species)%src_t*e/(p(rp)%data%m/fsup))
                        xold = p(rp)%x-p(rp)%data%v*dt
                        call get_intersect(xold, p(rp)%x, boundaries(ib), xi)
                        p(rp)%data%v = p(rp)%data%v - boundaries(ib)%n*dotproduct(p(rp)%data%v, boundaries(ib)%n)
                        call random_gauss_list(v_new(1:1),mu,sigma)
                        v_new(1) = abs(v_new(1))
                        p(rp)%data%v = p(rp)%data%v + boundaries(ib)%n*v_new(1)
                        ran = rnd_num()
                        p(rp)%x = xi + dt*ran*p(rp)%data%v
                        ib=0
                    ELSE IF (boundaries(ib)%type==6) THEN  !immediate half-Maxwellian refluxing normal to surface, tangential v resampled
                        mu=0.0_8
                        sigma=sqrt(species(p(rp)%data%species)%src_t*e/(p(rp)%data%m/fsup))
                        xold = p(rp)%x-p(rp)%data%v*dt
                        call get_intersect(xold, p(rp)%x, boundaries(ib), xi)
                        call random_gauss_list(v_new(1:3),mu,sigma)
                        v_new(1) = abs(v_new(1))
                        t1=boundaries(ib)%e1/sqrt(dotproduct(boundaries(ib)%e1,boundaries(ib)%e1))
                        t2(1)=boundaries(ib)%n(2)*t1(3) - boundaries(ib)%n(3)*t1(2)
                        t2(2)=boundaries(ib)%n(3)*t1(1) - boundaries(ib)%n(1)*t1(3)
                        t2(3)=boundaries(ib)%n(1)*t1(2) - boundaries(ib)%n(2)*t1(1)
                        t2=t2/sqrt(dotproduct(t2,t2))
                        p(rp)%data%v(1) = boundaries(ib)%n(1)*v_new(1) + t1(1)*v_new(2) + t2(1)*v_new(3)
                        p(rp)%data%v(2) = boundaries(ib)%n(2)*v_new(1) + t1(2)*v_new(2) + t2(2)*v_new(3)
                        p(rp)%data%v(3) = boundaries(ib)%n(3)*v_new(1) + t1(3)*v_new(2) + t2(3)*v_new(3)
                        ran = rnd_num()
                        p(rp)%x = xi + dt*ran*p(rp)%data%v
                        ib=0
                    ELSE IF (boundaries(ib)%type==7) THEN   !immediate Maxwellian flux refluxing normal to surface, tangential v conserved
                        mu=0.0_8
                        sigma=sqrt(species(p(rp)%data%species)%src_t*e/(p(rp)%data%m/fsup))
                        xold = p(rp)%x-p(rp)%data%v*dt
                        call get_intersect(xold, p(rp)%x, boundaries(ib), xi)
                        p(rp)%data%v = p(rp)%data%v - boundaries(ib)%n*dotproduct(p(rp)%data%v, boundaries(ib)%n)
                        call random_gaussian_flux(v_new(1),sigma)
                        p(rp)%data%v = p(rp)%data%v + boundaries(ib)%n*v_new(1)
                        ran = rnd_num()
                        p(rp)%x = xi + dt*ran*p(rp)%data%v
                        ib=0
                    ELSE IF (boundaries(ib)%type==8) THEN  !immediate Maxwellian flux refluxing normal to surface, tangential v resampled
                        mu=0.0_8
                        sigma=sqrt(species(p(rp)%data%species)%src_t*e/(p(rp)%data%m/fsup))
                        xold = p(rp)%x-p(rp)%data%v*dt
                        call get_intersect(xold, p(rp)%x, boundaries(ib), xi)
                        call random_gauss_list(v_new(2:3),mu,sigma)
                        call random_gaussian_flux(v_new(1),sigma)
                        t1=boundaries(ib)%e1/sqrt(dotproduct(boundaries(ib)%e1,boundaries(ib)%e1))
                        t2(1)=boundaries(ib)%n(2)*t1(3) - boundaries(ib)%n(3)*t1(2)
                        t2(2)=boundaries(ib)%n(3)*t1(1) - boundaries(ib)%n(1)*t1(3)
                        t2(3)=boundaries(ib)%n(1)*t1(2) - boundaries(ib)%n(2)*t1(1)
                        t2=t2/sqrt(dotproduct(t2,t2))
                        p(rp)%data%v(1) = boundaries(ib)%n(1)*v_new(1) + t1(1)*v_new(2) + t2(1)*v_new(3)
                        p(rp)%data%v(2) = boundaries(ib)%n(2)*v_new(1) + t1(2)*v_new(2) + t2(2)*v_new(3)
                        p(rp)%data%v(3) = boundaries(ib)%n(3)*v_new(1) + t1(3)*v_new(2) + t2(3)*v_new(3)
                        ran = rnd_num()
                        p(rp)%x = xi + dt*ran*p(rp)%data%v
                        ib=0
                    ELSE IF (boundaries(ib)%type==9) THEN   !immediate drifting Maxwellian flux refluxing normal to surface, tangential v conserved
                        xold = p(rp)%x-p(rp)%data%v*dt
                        call get_intersect(xold, p(rp)%x, boundaries(ib), xi)
                        p(rp)%data%v = p(rp)%data%v - boundaries(ib)%n*dotproduct(p(rp)%data%v, boundaries(ib)%n)
                        call random_drifting_gaussian_flux(v_new(1),maxw_flux_table_F(p(rp)%data%species,:),maxw_flux_table_v(p(rp)%data%species,:))
                        p(rp)%data%v = p(rp)%data%v + boundaries(ib)%n*v_new(1)
                        ran = rnd_num()
                        p(rp)%x = xi + dt*ran*p(rp)%data%v
                        ib=0
                    ELSE IF (boundaries(ib)%type==10) THEN  !immediate drifting Maxwellian flux refluxing normal to surface, tangential v resampled
                        xold = p(rp)%x-p(rp)%data%v*dt
                        call get_intersect(xold, p(rp)%x, boundaries(ib), xi)
                        call random_gauss_list(v_new(2:3),mu,sigma)
                        call random_drifting_gaussian_flux(v_new(1),maxw_flux_table_F(p(rp)%data%species,:),maxw_flux_table_v(p(rp)%data%species,:))
                        t1=boundaries(ib)%e1/sqrt(dotproduct(boundaries(ib)%e1,boundaries(ib)%e1))
                        t2(1)=boundaries(ib)%n(2)*t1(3) - boundaries(ib)%n(3)*t1(2)
                        t2(2)=boundaries(ib)%n(3)*t1(1) - boundaries(ib)%n(1)*t1(3)
                        t2(3)=boundaries(ib)%n(1)*t1(2) - boundaries(ib)%n(2)*t1(1)
                        t2=t2/sqrt(dotproduct(t2,t2))
                        p(rp)%data%v(1) = boundaries(ib)%n(1)*v_new(1) + t1(1)*v_new(2) + t2(1)*v_new(3)
                        p(rp)%data%v(2) = boundaries(ib)%n(2)*v_new(1) + t1(2)*v_new(2) + t2(2)*v_new(3)
                        p(rp)%data%v(3) = boundaries(ib)%n(3)*v_new(1) + t1(3)*v_new(2) + t2(3)*v_new(3)
                        ran = rnd_num()
                        p(rp)%x = xi + dt*ran*p(rp)%data%v
                        ib=0
                    ELSE IF (boundaries(ib)%type==2) THEN                              !Periodic BC
                        p(rp)%x = p(rp)%x + boundaries(ib)%n*boundaries(ib)%dist       !Particle re-enters at opposite boundary
                        ib = 0
                    ELSE IF (boundaries(ib)%type==1) THEN                              !Reflecting BC
                        xp = p(rp)%x - boundaries(ib)%x0
                        p(rp)%x = p(rp)%x - 2.*dotproduct(xp,boundaries(ib)%n)*boundaries(ib)%n
                        p(rp)%data%v = p(rp)%data%v - 2.*dotproduct(p(rp)%data%v,boundaries(ib)%n)*boundaries(ib)%n
                        ib = 0
                    ELSE IF (boundaries(ib)%type==0) THEN                              !Absorbing Wall BC
                        IF (boundaries(ib)%reflux_particles) reflux(p(rp)%data%species,ib)=reflux(p(rp)%data%species,ib)+1
                        npps(p(rp)%data%species) = npps(p(rp)%data%species) - 1
                        p(rp) = p(sum(npps)+1)
                        ib = 0
                        rp = rp - 1
                    END IF
                END IF
                ib = ib+1
            END DO
            rp = rp-1
        END DO

        !write(bnd_hits_out,*)



        IF (ANY(boundaries(1:nb)%type==4)) call treat_logical_sheath_boundaries(p,hits,reflux,p_hits_logical_sheath)
        IF (ANY(boundaries(1:nb)%type==4)) deallocate(p_hits_logical_sheath)

        call reallocate_particles(p,sum(npps), sum(npps))

    END SUBROUTINE

!======================================================================================

    SUBROUTINE rethermalize(p)
        implicit none
        include 'mpif.h'

        integer :: ip
        type(t_particle), allocatable, intent(inout) :: p(:)
        real*8 :: xold(3)
        real(KIND=8) mu,sigma,ran,t1(3),t2(3),v_ran(3),ran2,ran3

        real(KIND=8)       :: n1(3),B_vector(3)
        real(KIND=8)       :: eps=1.0e-10

        IF (retherm /= 0) THEN

            B_vector(1) = Bx
            B_vector(2) = By
            B_vector(3) = Bz

            IF (real_unequal_zero(B,eps)) THEN                                ! With B Field
                n1 = B_vector / sqrt(dotproduct(B_vector,B_vector))           ! normal vector along B
            ELSE                                                              ! Without B Field
                n1(1) = 1.0_8                                                 ! normal vector along x
                n1(2) = 0.0_8
                n1(3) = 0.0_8
            END IF
            IF ((real_unequal_zero(n1(1),eps)) .OR. (real_unequal_zero(n1(2),eps))) THEN
                t1(1) = -n1(2)                                            ! tangential vector
                t1(2) = n1(1)
                t1(3) = 0.0_8
            ELSE
                t1(1) = 1.0_8                                             ! tangential vector
                t1(2) = 0.0_8
                t1(3) = 0.0_8
            END IF
            t1 = t1 / sqrt(dotproduct(t1,t1))
            t2(1) = n1(2)*t1(3) - n1(3)*t1(2)                             ! t2 = n1 x t1 (2nd tangential vector)
            t2(2) = n1(3)*t1(1) - n1(1)*t1(3)
            t2(3) = n1(1)*t1(2) - n1(2)*t1(1)
            t2=t2/sqrt(dotproduct(t2,t2))

            DO ip =1,sum(npps)
                IF (sign(1._8,p(ip)%x(1)) /= sign(1._8,p(ip)%data%v(1))) CYCLE

                xold = p(ip)%x-p(ip)%data%v*dt
                IF (sign(1._8,xold(1)) == sign(1._8,p(ip)%x(1))) CYCLE

                mu=0.0_8
                sigma=sqrt(species(p(ip)%data%species)%src_t*e/(p(ip)%data%m/fsup))
                
                IF (retherm == 1) THEN
                    ! resample vx according to maxwellian flux source, don't touch vy,vz
                    ! resample x = ran*dt*vx, don't touch y,z
                    call random_gaussian_flux(v_ran(1),sigma)
                    p(ip)%data%v(1) = v_ran(1)*sign(1._8,p(ip)%x(1))
                    ran = rnd_num()
                    p(ip)%x(1) = dt*ran*p(ip)%data%v(1)

                ELSE IF (retherm == 2) THEN
                    ! resample vx according to maxwellian flux source, vy,vz according to maxwellian source
                    ! resample x = ran*dt*vx, y,z = ran
                    call random_gauss_list(v_ran(2:3),mu,sigma)
                    call random_gaussian_flux(v_ran(1),sigma)
                    p(ip)%data%v(1) = v_ran(1)*sign(1._8,p(ip)%x(1))
                    p(ip)%data%v(2) = v_ran(2)
                    p(ip)%data%v(3) = v_ran(3)
                    ran = rnd_num()
                    ran2 = rnd_num()
                    ran3 = rnd_num()
                    p(ip)%x(1) = dt*ran*p(ip)%data%v(1)
                    p(ip)%x(2) = ymin + ran2*(ymax-ymin)
                    p(ip)%x(3) = zmin + ran3*(zmax-zmin)


                ELSE IF (retherm == 3) THEN
                    ! resample vpar according to maxwellian flux source, don't touch vperp
                    ! resample x = ran*dt*vx, don't touch y,z
                    call random_gaussian_flux(v_ran(1),sigma)
                    p(ip)%data%v = p(ip)%data%v - (dotproduct(p(ip)%data%v, n1) * n1)                       !this gives v - vpar
                    p(ip)%data%v = p(ip)%data%v + (v_ran(1) * n1 *sign(1._8,p(ip)%x(1)) * sign(1._8,n1(1))) !this adds resampled vpar
                    ran = rnd_num()
                    p(ip)%x(1) = dt*ran*p(ip)%data%v(1)


                ELSE IF (retherm == 4) THEN
                    ! resample vpar according to maxwellian flux source, vperp1 and vperp2 according to Maxwellian source
                    ! resample x = ran*dt*vx, y,z = ran
                    call random_gaussian_flux(v_ran(1),sigma)
                    call random_gauss_list(v_ran(2:3),mu,sigma)
                    v_ran(1) = v_ran(1) * sign(1._8,p(ip)%x(1)) * sign(1._8,n1(1))
                    p(ip)%data%v(1) = n1(1)*v_ran(1) + t1(1)*v_ran(2) + t2(1)*v_ran(3) ! this gives the maxwellian flux
                    p(ip)%data%v(2) = n1(2)*v_ran(1) + t1(2)*v_ran(2) + t2(2)*v_ran(3) ! along n1 and gaussian
                    p(ip)%data%v(3) = n1(3)*v_ran(1) + t1(3)*v_ran(2) + t2(3)*v_ran(3) ! distribution along t1, t2
                    ran = rnd_num()
                    ran2 = rnd_num()
                    ran3 = rnd_num()
                    p(ip)%x(1) = dt*ran*p(ip)%data%v(1)
                    p(ip)%x(2) = ymin + ran2*(ymax-ymin)
                    p(ip)%x(3) = zmin + ran3*(zmax-zmin)


                END IF
            END DO
        END IF
    END SUBROUTINE


!======================================================================================

    SUBROUTINE treat_logical_sheath_boundaries(p,hits,reflux,p_hits_logical_sheath)
        implicit none
        include 'mpif.h'

        type(t_particle), allocatable, intent(inout) :: p(:),p_hits_logical_sheath(:,:,:)
        integer, intent(inout) :: hits(0:,:),reflux(0:,:)

        integer ib,ispecies,i,ip,rc
        integer :: thits(0:n_ranks-1,0:nspecies-1,nb),displs(0:n_ranks-1)
        real*8 :: v_perp_logical_sheath(sum(npps),0:nspecies-1,nb)
        real*8 :: t_v_perp_logical_sheath(sum(tnpps),0:nspecies-1,nb)
        real*8,allocatable :: v_perp_temp(:)
        real*8 :: vcut_el(nb)
        integer :: ihits,ehits


        thits=0
        vcut_el=0
        DO ib = 1,nb
            IF (boundaries(ib)%type==4) THEN
                DO ispecies = 0,nspecies-1
                    IF (species(ispecies)%physical_particle) THEN
                        call MPI_GATHER(hits(ispecies,ib),1,MPI_INTEGER, thits(:,ispecies,ib), 1, MPI_INTEGER, 0,MPI_COMM_WORLD, rc)
                    END IF
                END DO
            END IF
        END DO

        DO ib = 1,nb
            IF (boundaries(ib)%type==4) THEN
                DO ispecies = 0,nspecies-1
                    IF (species(ispecies)%physical_particle) THEN
                        IF (root) THEN
                            displs(0)=0
                            DO i=1,n_ranks-1
                                displs(i)=displs(i-1)+thits(i-1,ispecies,ib)
                            END DO
                        END IF
                        DO ip = 1,hits(ispecies,ib)
                            v_perp_logical_sheath(ip,ispecies,ib)=-dotproduct(p_hits_logical_sheath(ispecies,ib,ip)%data%v &
                                                                              ,boundaries(ib)%n)
                        END DO
                        IF (species(ispecies)%indx==1) THEN !sorting only for electrons
                            call MPI_GATHERV(v_perp_logical_sheath(:,ispecies,ib),hits(ispecies,ib),&
                                             MPI_DOUBLE_PRECISION,t_v_perp_logical_sheath(:,ispecies,ib),&
                                             thits(:,ispecies,ib),displs,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD, rc)
                            IF (root) THEN
                                allocate(v_perp_temp(sum(thits(:,ispecies,ib))),stat=rc)
                                v_perp_temp(:) = t_v_perp_logical_sheath(1:sum(thits(:,ispecies,ib)),ispecies,ib)
                                call QsortC(v_perp_temp)
                                t_v_perp_logical_sheath(1:sum(thits(:,ispecies,ib)),ispecies,ib)=v_perp_temp
                                ehits=sum(thits(:,ispecies,ib))
                                ihits=0
                                DO i=0,nspecies-1
                                    IF (species(i)%physical_particle .and. species(i)%q > 0) THEN
                                        ihits=ihits+sum(thits(:,i,ib))
                                    END IF
                                END DO
                                IF (ehits-ihits>0) THEN
                                    vcut_el(ib)=v_perp_temp(ehits-ihits)
                                ELSE
                                    vcut_el(ib)=0.0_8
                                END IF
                                deallocate(v_perp_temp)
                            END IF
                            call MPI_BCAST(vcut_el(ib),1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,rc)
                        END IF
                    END IF
                END DO
            END IF
        END DO

        DO ib = 1,nb
            IF (boundaries(ib)%type==4) THEN
               call repell_electrons_at_logical_sheath(p,hits,reflux,vcut_el(ib),p_hits_logical_sheath,ib)
            END IF
        END DO


    END SUBROUTINE treat_logical_sheath_boundaries

!======================================================================================

    SUBROUTINE repell_electrons_at_logical_sheath(p,hits,reflux,vcut_el,p_hits_logical_sheath,ib)

        implicit none
        include 'mpif.h'

        type(t_particle), allocatable, intent(inout) :: p(:),p_hits_logical_sheath(:,:,:)
        integer, intent(inout) :: hits(0:,:),reflux(0:,:)
        real(KIND=8),intent(in) :: vcut_el
        integer, intent(in) :: ib
        integer :: ispecies

        integer :: ip
        real(KIND=8) :: v_perp,xp(3)

        ispecies=1

        IF (vcut_el<=0.0_8) return

        DO ip=1,hits(ispecies,ib)
            v_perp=-dotproduct(p_hits_logical_sheath(ispecies,ib,ip)%data%v,boundaries(ib)%n)
            IF (v_perp <= vcut_el) THEN !reflect electrons to get ambipolar flux
                npps(ispecies) = npps(ispecies) + 1
                p(sum(npps))=p_hits_logical_sheath(ispecies,ib,ip)
                xp = p(sum(npps))%x - boundaries(ib)%x0
                p(sum(npps))%x = p(sum(npps))%x - 2.*dotproduct(xp,boundaries(ib)%n)*boundaries(ib)%n
                p(sum(npps))%data%v = p(sum(npps))%data%v - 2.*dotproduct(p(sum(npps))%data%v,boundaries(ib)%n)*boundaries(ib)%n
                hits(ispecies,ib) = hits(ispecies,ib) - 1
                IF (reflux(ispecies,ib)>0) reflux(ispecies,ib) = reflux(ispecies,ib) - 1
            END IF
        END DO

    END SUBROUTINE repell_electrons_at_logical_sheath

!======================================================================================

    SUBROUTINE recycling(p,treflux)

        implicit none
        include 'mpif.h'

        type(t_particle), allocatable, intent(inout) :: p(:)
        integer, intent(in) :: treflux(0:,:)
        integer :: reflux(0:nspecies-1,1:nb)
        integer :: rc,ispecies,ib,new_particles


        reflux=0
        new_particles=0

        DO ispecies=0,nspecies-1
            DO ib=1,nb
                IF (treflux(ispecies,ib)/=0) THEN
                    reflux(ispecies,ib) = treflux(ispecies,ib) / n_ranks
                    next_label = next_label + my_rank*reflux(ispecies,ib)
                    IF (my_rank .eq. (n_ranks-1)) THEN
                        reflux(ispecies,ib) = reflux(ispecies,ib) + MOD(treflux(ispecies,ib), n_ranks)
                    END IF

                    call reflux_particles(p,reflux(ispecies,ib),ispecies)
                    call MPI_BCAST(next_label, 1, MPI_INTEGER, n_ranks-1, MPI_COMM_WORLD, rc)
                END IF

            END DO
        END DO


        DO ispecies=0,nspecies-1
            IF (species(ispecies)%nfp/=0) THEN
                IF (species(ispecies)%physical_particle) THEN
                    new_particles=species(ispecies)%nfp / n_ranks
                    next_label = next_label + my_rank * (new_particles)
                    IF (my_rank .eq. (n_ranks-1)) THEN
                        new_particles = new_particles + MOD((species(ispecies)%nfp), n_ranks)
                    END IF

                    call reflux_particles(p,new_particles,ispecies)
                    call MPI_BCAST(next_label, 1, MPI_INTEGER, n_ranks-1, MPI_COMM_WORLD, rc)
                END IF
            END IF
        END DO

    END SUBROUTINE recycling

!======================================================================================

    SUBROUTINE reflux_particles(p,new_particles,ispecies)

        implicit none

        type(t_particle), allocatable, intent(inout) :: p(:)
        integer, intent(in) :: new_particles,ispecies

        integer rc,ip
        type(t_particle), allocatable                :: p_new(:)


        allocate(p_new(new_particles),stat=rc)

        IF (new_particles==0) THEN
            deallocate(p_new)
            return
        ELSE
            call reallocate_particles(p,sum(npps), sum(npps)+new_particles)

            DO ip=1, new_particles
                p_new(ip)%label       = next_label
                next_label            = next_label + 1
                p_new(ip)%data%q      = species(ispecies)%q*fsup
                p_new(ip)%data%m      = species(ispecies)%m*fsup

                p_new(ip)%results%e   = 0.0_8
                p_new(ip)%results%pot = 0.0_8
                p_new(ip)%work        = 1.0_8
                p_new(ip)%data%species= ispecies

                p_new(ip)%data%B(1)=Bx
                p_new(ip)%data%B(2)=By
                p_new(ip)%data%B(3)=Bz
            END DO
            call source(p_new)

            p(sum(npps)+1:sum(npps)+new_particles)=p_new(:)
            npps(ispecies) = npps(ispecies) + new_particles

            deallocate(p_new)
        END IF

    END SUBROUTINE reflux_particles


!======================================================================================

    SUBROUTINE charge_wall(p,thits)
        implicit none

        type(t_particle), allocatable, intent(inout) :: p(:)
        integer, intent(inout) :: thits(0:,:)
        integer :: ip,ispecies,ib
        real*8  :: dq(nb)
        logical :: hit

        dq=0.0_8

        DO ib=1,nb
            DO ispecies=0,nspecies-1
                IF (species(ispecies)%physical_particle) THEN
                    dq(ib) = dq(ib) + thits(ispecies,ib)*species(ispecies)%q*fsup
                END IF
            END DO
            boundaries(ib)%q_tot = boundaries(ib)%q_tot + dq(ib)
        END DO

        DO ip=1, sum(npps)                                            ! charge wall by adding charge to the wall particles
            IF (p(ip)%data%species/=0) CYCLE
            DO ib=1,nb
                call check_hit(p(ip)%x(1),p(ip)%x(2),p(ip)%x(3),boundaries(ib),hit)
                IF(hit) THEN
                    p(ip)%data%q = boundaries(ib)%q_tot / boundaries(ib)%nwp
                END IF
            END DO
        END DO
    END SUBROUTINE charge_wall


!======================================================================================
    SUBROUTINE set_need_to_reflux()

        implicit none

        need_to_reflux=.false.

        if (MOD(step-last_reflux_step,1)==0)then
            need_to_reflux=.true.
            last_reflux_step=step
        end if

        if(checkp_interval.ne.0) then
            if ((MOD(step,checkp_interval)==0).or.(step==nt+startstep)) then
                need_to_reflux=.true.
                last_reflux_step=step
            end if
        end if
        if(diag_interval.ne.0) then
            if ((MOD(step,diag_interval)==0).or.(step==nt+startstep)) THEN
                need_to_reflux=.true.
                last_reflux_step=step
            end if
        end if

    END SUBROUTINE set_need_to_reflux

!======================================================================================

    SUBROUTINE get_number_of_particles(p)
        implicit none

        type(t_particle), allocatable, intent(inout) :: p(:)
        integer :: ip,ispecies


        npps=0
        DO ip=1,size(p)
           ispecies = p(ip)%data%species
           IF (ispecies >= 0) npps(ispecies) = npps(ispecies) + 1
        END DO





    END SUBROUTINE get_number_of_particles


!======================================================================================
    subroutine hits_on_boundaries(p)
        implicit none
        include 'mpif.h'

        type(t_particle), allocatable, intent(inout) :: p(:)

        integer      :: rc

        integer      :: hits(0:nspecies-1,1:nb),thits(0:nspecies-1,1:nb)
        integer      :: reflux(0:nspecies-1,1:nb),treflux(0:nspecies-1,1:nb)
        integer      :: ispecies,ib


        hits=0
        reflux=0
        thits=0
        treflux=0

        if(root) write(*,'(a)') " == [hits_on_boundaries] count hits and recycle "

        call set_need_to_reflux()
        call count_hits_and_remove_particles(p,hits,reflux)

        call rethermalize(p)

        DO ispecies=0,nspecies-1
            DO ib=1,nb
                call MPI_ALLREDUCE(hits(ispecies,ib), thits(ispecies,ib), 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, rc)
                call MPI_ALLREDUCE(reflux(ispecies,ib), treflux(ispecies,ib), 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, rc)
            END DO
        END DO


        call recycling(p,treflux)

        call charge_wall(p,thits)

        call MPI_ALLREDUCE(npps, tnpps, nspecies, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, rc)

        call set_recycling_output_values(thits,treflux)

    END SUBROUTINE hits_on_boundaries

!======================================================================================
    subroutine check_for_leakage(p)
        implicit none
        include 'mpif.h'

        type(t_particle), allocatable, intent(inout) :: p(:)
        type(t_particle), allocatable :: p1(:)

        integer      :: rc
        integer      :: ip,nlp(0:nspecies-1),tnlp(0:nspecies-1),n

        nlp=0
        tnlp=0

        n=size(p)
        allocate(p1(1),stat=rc)

        DO ip=1,n
            IF (p(ip)%x(1)>=xmin .AND. p(ip)%x(1)<=xmax .AND. p(ip)%x(2)>=ymin .AND. &
                p(ip)%x(2)<=ymax .AND. p(ip)%x(3)>=zmin .AND. p(ip)%x(3)<=zmax)   THEN
                CYCLE
            ELSE
                nlp(p(ip)%data%species) = nlp(p(ip)%data%species) + 1

                IF (species(p(ip)%data%species)%physical_particle) THEN
                    write(*,'(i6,a,i16,i16)')my_rank,": label,species: ", p(ip)%label,p(ip)%data%species
                    write(*,'(i6,a,3(1pe16.7E3))')my_rank,": particle velocity: ", p(ip)%data%v
                    write(*,'(i6,a,3(1pe16.7E3))')my_rank,": particle position: ", p(ip)%x
                    write(*,'(i6,a,3(1pe16.7E3))')my_rank,": dx,dy,dz: ", [dx,dy,dz]
                    write(*,'(i6,a,3(1pe16.7E3))')my_rank,": old particle position: ", p(ip)%x - dt*p(ip)%data%v
                    write(*,'(a)') "Particle refluxed according to chosen source!"
                    p1(1)=p(ip)
                    call source(p1)
                    p(ip)=p1(1)
                END IF
            END IF
        END DO

        IF (nlp(0)>0) call redistribute_wall_particles(p)

        call MPI_ALLREDUCE(nlp, tnlp, nspecies, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, rc)
        IF (root) THEN
            if(sum(tnlp)>0) write(*,*) "Number of leaked particles per species:",tnlp
        END IF

    END SUBROUTINE check_for_leakage


!======================================================================================
    SUBROUTINE init_particles(p)
        implicit none
    
        type(t_particle), intent(inout), allocatable :: p(:)
        integer :: ip,ispecies,i

        ip=0
        DO ispecies=1,nspecies-1
            IF (species(ispecies)%physical_particle) THEN
                DO i=1, npps(ispecies)
                    ip = ip + 1
                    p(ip)%data%B(1)=Bx
                    p(ip)%data%B(2)=By
                    p(ip)%data%B(3)=Bz
                    p(ip)%label       = my_rank * (SUM(tnpps(1:nspecies-1)) / n_ranks) + ip
                    p(ip)%data%q      = species(ispecies)%q*fsup
                    p(ip)%data%m      = species(ispecies)%m*fsup
                    p(ip)%data%species= species(ispecies)%indx
                    p(ip)%results%e   = 0.0_8
                    p(ip)%results%pot = 0.0_8
                    p(ip)%work        = 1.0_8
                END DO
            ELSE
                call init_probes(ispecies,ip,p)
            END IF
        END DO

        ispecies=0
        DO i=1,npps(ispecies)
            ip = ip + 1
            p(ip)%data%B(1)=Bx
            p(ip)%data%B(2)=By
            p(ip)%data%B(3)=Bz
            p(ip)%label       = -(my_rank * (tnpps(0) / n_ranks) + i)
            p(ip)%data%q      = species(ispecies)%q*fsup
            p(ip)%data%m      = species(ispecies)%m*fsup
            p(ip)%data%species= species(ispecies)%indx
            p(ip)%results%e   = 0.0_8
            p(ip)%results%pot = 0.0_8
            p(ip)%work        = 1.0_8
        END DO

        call source(p)
        
        next_label = SUM(tnpps(1:nspecies-1))+1
  
    END SUBROUTINE init_particles

!======================================================================================

    SUBROUTINE init_probes(ispecies,ip,p)
        implicit none

        type(t_particle), intent(inout), allocatable :: p(:)
        integer, intent(in) :: ispecies
        integer, intent(inout) :: ip
        integer :: j,i


        j=0
        DO i=1, npps(ispecies)
            ip = ip + 1
            p(ip)%data%B(1)=Bx
            p(ip)%data%B(2)=By
            p(ip)%data%B(3)=Bz
            p(ip)%label       = my_rank * (SUM(tnpps(1:nspecies-1)) / n_ranks) + ip
            p(ip)%data%q      = species(ispecies)%q*fsup
            p(ip)%data%m      = species(ispecies)%m*fsup
            p(ip)%data%species= species(ispecies)%indx
            p(ip)%results%e   = 0.0_8
            p(ip)%results%pot = 0.0_8
            p(ip)%work        = 1.0_8
            p(ip)%x(1) = probe_start_x(ispecies) + (j+0.5) * (probe_end_x(ispecies) - probe_start_x(ispecies)) / species(ispecies)%nip
            p(ip)%x(2) = probe_start_y(ispecies) + (j+0.5) * (probe_end_y(ispecies) - probe_start_y(ispecies)) / species(ispecies)%nip
            p(ip)%x(3) = probe_start_z(ispecies) + (j+0.5) * (probe_end_z(ispecies) - probe_start_z(ispecies)) / species(ispecies)%nip
            p(ip)%data%v = 0.
            j=j+1
        END DO

    END SUBROUTINE init_probes

!======================================================================================

    SUBROUTINE remove_all_probes(ispecies,p)
        implicit none
        include 'mpif.h'

        type(t_particle), intent(inout), allocatable :: p(:)
        integer, intent(in) :: ispecies
        integer :: i,j

        IF (.not. species(ispecies)%physical_particle) THEN
            j = sum(npps)
            i = 1
            DO WHILE (i <= j)
                IF (i > j) EXIT
                IF (p(i)%data%species  == ispecies) THEN
                    p(i) = p(j)
                    j = j - 1
                    i = i - 1
                    npps(ispecies) = npps(ispecies) - 1
                END IF
                i = i + 1
            END DO
            call reallocate_particles(p,j,j)
            call MPI_ALLREDUCE(npps, tnpps, nspecies, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
        END IF
    END SUBROUTINE remove_all_probes

!======================================================================================

    subroutine redistribute_wall_particles(p)
        implicit none

        type(t_particle), intent(inout) :: p(:)
        real :: ran1,ran2
        integer :: ip,ib
        logical :: hit

        DO ip=1, sum(npps)
            IF (p(ip)%data%species/=0) CYCLE
            DO ib=1,nb
                IF (ANY(boundaries(ib)%wp_labels==p(ip)%label)) THEN
                    call check_hit(p(ip)%x(1),p(ip)%x(2),p(ip)%x(3),boundaries(ib),hit)
                    IF (.not. hit) THEN
                        ran1=rnd_num()
                        ran2=rnd_num()
                        p(ip)%x = boundaries(ib)%x0 + ran1*boundaries(ib)%e1 +ran2*boundaries(ib)%e2
                    END IF
                END IF
            END DO
        END DO

    end subroutine redistribute_wall_particles


!======================================================================================

    SUBROUTINE source(p)
        type(t_particle),intent(inout),allocatable :: p(:)
        real               :: ran,ran1,ran2
        integer            :: n,ip,ib
        real*8             :: mu,sigma
        real(KIND=8)       :: t1(3),t2(3),n1(3),e1(3),vhelp(3),B_vector(3)
        real(KIND=8)       :: eps=1.0e-10

        logical :: init_uniform_gaussian

        init_uniform_gaussian=.false.

        B_vector(1) = Bx
        B_vector(2) = By
        B_vector(3) = Bz

        vhelp=0.0_8
        mu=0.0
        IF (allocated(p)) THEN
            n=size(p)
        ELSE
            n=0
        END IF





        DO ip=1, n
            IF (species(p(ip)%data%species)%physical_particle .eqv. .false.) THEN
                IF (p(ip)%data%species==0) THEN
                    p(ip)%data%v = 0.0_8
                    DO ib=1,nb
                        IF (ANY(boundaries(ib)%wp_labels==p(ip)%label)) THEN
                            ran1=rnd_num()
                            ran2=rnd_num()
                            p(ip)%x = boundaries(ib)%x0 + ran1*boundaries(ib)%e1 +ran2*boundaries(ib)%e2
                        END IF
                    END DO
                END IF
            ELSE
                sigma=sqrt(species(p(ip)%data%species)%src_t*e / (p(ip)%data%m/fsup))
                IF (init_uniform_gaussian) THEN
                    IF (step==0) THEN
                        ran=rnd_num()
                        ran1=rnd_num()
                        ran2=rnd_num()
                        p(ip)%x(1)=ran
                        p(ip)%x(2)=ran1
                        p(ip)%x(3)=ran2

                        p(ip)%x(1)         = p(ip)%x(1)*dx + xmin
                        p(ip)%x(2)         = p(ip)%x(2)*dy + ymin
                        p(ip)%x(3)         = p(ip)%x(3)*dz + zmin

                        call random_gauss_list(p(ip)%data%v(1:3),mu,sigma)
                        CYCLE
                    END IF
                END IF
                IF (species(p(ip)%data%species)%src_type==0) THEN                 !surface source (Maxwellian Flux)
                    call random_gauss_list(vhelp(2:3),mu,sigma)
                    call random_gaussian_flux(vhelp(1),sigma)
                    e1 = boundaries(species(p(ip)%data%species)%src_bnd)%e1
                    n1 = boundaries(species(p(ip)%data%species)%src_bnd)%n        ! n1 = normal vector on src_bnd
                    t1 = e1 / sqrt(dotproduct(e1,e1))                             ! t1 = tangential vector
                    t2(1) = n1(2)*t1(3) - n1(3)*t1(2)                             ! t2 = n1 x t1 (2nd tangential vector)
                    t2(2) = n1(3)*t1(1) - n1(1)*t1(3)
                    t2(3) = n1(1)*t1(2) - n1(2)*t1(1)
                    t2=t2/sqrt(dotproduct(t2,t2))
                    p(ip)%data%v(1) = n1(1)*vhelp(1) + t1(1)*vhelp(2) + t2(1)*vhelp(3) ! this gives the maxwellian flux
                    p(ip)%data%v(2) = n1(2)*vhelp(1) + t1(2)*vhelp(2) + t2(2)*vhelp(3) ! along n1 and gaussian
                    p(ip)%data%v(3) = n1(3)*vhelp(1) + t1(3)*vhelp(2) + t2(3)*vhelp(3) ! distribution along t1, t2
                    ran=rnd_num()
                    ran1=rnd_num()
                    ran2=rnd_num()
                    p(ip)%x = boundaries(species(p(ip)%data%species)%src_bnd)%x0 + &
                              ran1*boundaries(species(p(ip)%data%species)%src_bnd)%e1 + &
                              ran2*boundaries(species(p(ip)%data%species)%src_bnd)%e2
                    p(ip)%x = p(ip)%x + boundaries(species(p(ip)%data%species)%src_bnd)%n * &
                              dotproduct(boundaries(species(p(ip)%data%species)%src_bnd)%n,p(ip)%data%v) * dt * ran
                ELSE IF (species(p(ip)%data%species)%src_type==5) THEN                 !surface source (drifting Maxwellian Flux)
                    call random_gauss_list(vhelp(2:3),mu,sigma)
                    call random_drifting_gaussian_flux(vhelp(1),maxw_flux_table_F(p(ip)%data%species,:),maxw_flux_table_v(p(ip)%data%species,:))
                    e1 = boundaries(species(p(ip)%data%species)%src_bnd)%e1
                    n1 = boundaries(species(p(ip)%data%species)%src_bnd)%n        ! n1 = normal vector on src_bnd
                    t1 = e1 / sqrt(dotproduct(e1,e1))                             ! t1 = tangential vector
                    t2(1) = n1(2)*t1(3) - n1(3)*t1(2)                             ! t2 = n1 x t1 (2nd tangential vector)
                    t2(2) = n1(3)*t1(1) - n1(1)*t1(3)
                    t2(3) = n1(1)*t1(2) - n1(2)*t1(1)
                    t2=t2/sqrt(dotproduct(t2,t2))
                    p(ip)%data%v(1) = n1(1)*vhelp(1) + t1(1)*vhelp(2) + t2(1)*vhelp(3) ! this gives the maxwellian flux
                    p(ip)%data%v(2) = n1(2)*vhelp(1) + t1(2)*vhelp(2) + t2(2)*vhelp(3) ! along n1 and gaussian
                    p(ip)%data%v(3) = n1(3)*vhelp(1) + t1(3)*vhelp(2) + t2(3)*vhelp(3) ! distribution along t1, t2
                    ran=rnd_num()
                    ran1=rnd_num()
                    ran2=rnd_num()
                    p(ip)%x = boundaries(species(p(ip)%data%species)%src_bnd)%x0 + &
                              ran1*boundaries(species(p(ip)%data%species)%src_bnd)%e1 + &
                              ran2*boundaries(species(p(ip)%data%species)%src_bnd)%e2
                    p(ip)%x = p(ip)%x + boundaries(species(p(ip)%data%species)%src_bnd)%n * &
                              dotproduct(boundaries(species(p(ip)%data%species)%src_bnd)%n,p(ip)%data%v) * dt * ran
                ELSE IF (species(p(ip)%data%species)%src_type==4) THEN                 !surface source
                    call random_gauss_list(vhelp(2:3),mu,sigma)
                    call random_gaussian_flux(vhelp(1),sigma)
                    e1 = boundaries(species(p(ip)%data%species)%src_bnd)%e1
                    n1 = boundaries(species(p(ip)%data%species)%src_bnd)%n        ! n1 = normal vector on src_bnd
                    t1 = e1 / sqrt(dotproduct(e1,e1))                             ! t1 = tangential vector
                    t2(1) = n1(2)*t1(3) - n1(3)*t1(2)                             ! t2 = n1 x t1 (2nd tangential vector)
                    t2(2) = n1(3)*t1(1) - n1(1)*t1(3)
                    t2(3) = n1(1)*t1(2) - n1(2)*t1(1)
                    t2=t2/sqrt(dotproduct(t2,t2))
                    p(ip)%data%v(1) = n1(1)*vhelp(1) + t1(1)*vhelp(2) + t2(1)*vhelp(3) ! this gives the maxwellian flux
                    p(ip)%data%v(2) = n1(2)*vhelp(1) + t1(2)*vhelp(2) + t2(2)*vhelp(3) ! along n1 and gaussian
                    p(ip)%data%v(3) = n1(3)*vhelp(1) + t1(3)*vhelp(2) + t2(3)*vhelp(3) ! distribution along t1, t2
                    ran=rnd_num()
                    ran1=rnd_num() * species(p(ip)%data%species)%src_x0(3) !radius
                    ran2=rnd_num() * 2.0_8 * pi              !angle
                    p(ip)%x = boundaries(species(p(ip)%data%species)%src_bnd)%x0 + &
                              species(p(ip)%data%species)%src_x0(1)*boundaries(species(p(ip)%data%species)%src_bnd)%e1 + &
                              species(p(ip)%data%species)%src_x0(2)*boundaries(species(p(ip)%data%species)%src_bnd)%e2
                              !now we are at the center of the plasma column
                    p(ip)%x = p(ip)%x + ran1*sin(ran2)*t1 + ran1*cos(ran2)*t2
                              !now we have a cylindrical plasma
                    p(ip)%x = p(ip)%x + boundaries(species(p(ip)%data%species)%src_bnd)%n * &
                              dotproduct(boundaries(species(p(ip)%data%species)%src_bnd)%n,p(ip)%data%v) * dt * ran
                              !and we moved the particle a little bit along n
                ELSE IF (species(p(ip)%data%species)%src_type==1) THEN            ! Emmert source along x, Maxwellian perpendicular to x
                    call random_gauss_list(p(ip)%data%v(2:3),mu,sigma)
                    call random_gaussian_flux(p(ip)%data%v(1),sigma)
                    ran=rnd_num()
                    IF (ran>0.5) p(ip)%data%v(1)=-p(ip)%data%v(1)
                    ran=rnd_num()
                    ran1=rnd_num()
                    ran2=rnd_num()
                    p(ip)%x = species(p(ip)%data%species)%src_x0 + &
                              ran * species(p(ip)%data%species)%src_e1 + &
                              ran1 * species(p(ip)%data%species)%src_e2 + &
                              ran2 * species(p(ip)%data%species)%src_e3
                ELSE IF (species(p(ip)%data%species)%src_type==3) THEN                ! Emmert source along B, Maxwellian perpendicular to B
                    IF (real_unequal_zero(B,eps)) THEN                                ! With B Field
                        n1 = B_vector / sqrt(dotproduct(B_vector,B_vector))           ! normal vector along B
                    ELSE                                                              ! Without B Field
                        n1(1) = 1.0_8                                                 ! normal vector along x
                        n1(2) = 0.0_8
                        n1(3) = 0.0_8
                    END IF
                    IF ((real_unequal_zero(n1(1),eps)) .OR. (real_unequal_zero(n1(2),eps))) THEN
                        t1(1) = -n1(2)                                            ! tangential vector
                        t1(2) = n1(1)
                        t1(3) = 0.0_8
                    ELSE
                        t1(1) = 1.0_8                                             ! tangential vector
                        t1(2) = 0.0_8
                        t1(3) = 0.0_8
                    END IF
                    t1 = t1 / sqrt(dotproduct(t1,t1))
                    t2(1) = n1(2)*t1(3) - n1(3)*t1(2)                             ! t2 = n1 x t1 (2nd tangential vector)
                    t2(2) = n1(3)*t1(1) - n1(1)*t1(3)
                    t2(3) = n1(1)*t1(2) - n1(2)*t1(1)
                    t2=t2/sqrt(dotproduct(t2,t2))

                    call random_gauss_list(vhelp(2:3),mu,sigma)
                    call random_gaussian_flux(vhelp(1),sigma)
                    ran=rnd_num()
                    IF (ran>0.5) vhelp(1)=-vhelp(1)
                    p(ip)%data%v(1) = n1(1)*vhelp(1) + t1(1)*vhelp(2) + t2(1)*vhelp(3) ! this gives the maxwellian flux
                    p(ip)%data%v(2) = n1(2)*vhelp(1) + t1(2)*vhelp(2) + t2(2)*vhelp(3) ! along n1 and gaussian
                    p(ip)%data%v(3) = n1(3)*vhelp(1) + t1(3)*vhelp(2) + t2(3)*vhelp(3) ! distribution along t1, t2
                    ran=rnd_num()
                    ran1=rnd_num()
                    ran2=rnd_num()
                    p(ip)%x = species(p(ip)%data%species)%src_x0 + &
                              ran * species(p(ip)%data%species)%src_e1 + &
                              ran1 * species(p(ip)%data%species)%src_e2 + &
                              ran2 * species(p(ip)%data%species)%src_e3
                ELSE IF (species(p(ip)%data%species)%src_type==2) THEN            ! Bissel Johnson source along x
                    call random_gauss_list(p(ip)%data%v(1:3),mu,sigma)
                    ran=rnd_num()
                    ran1=rnd_num()
                    ran2=rnd_num()
                    p(ip)%x = species(p(ip)%data%species)%src_x0 + &
                              ran * species(p(ip)%data%species)%src_e1 + &
                              ran1 * species(p(ip)%data%species)%src_e2 + &
                              ran2 * species(p(ip)%data%species)%src_e3
                END IF
            END IF
        END DO

        !If treated in guding centre approximation, electrons just need other startign values in
        !the Boris scheme (see Benjamin Berberichs Diss (pp. 44-56)
        IF (guiding_centre_electrons) THEN
            call transform_electrons_to_gc(p)
        END IF
    END SUBROUTINE


!======================================================================================

    SUBROUTINE transform_electrons_to_gc(p)
        type(t_particle), intent(inout) :: p(:)
        integer            :: n,ip
        real*8             :: vpar(3),b,E_tot(3),E_ext(3),vd(3),omega_c


        E_ext=0.0_8 !Just done for now, because Eext is not introduced yet

        n=size(p)
        DO ip=1, n
            IF (p(ip)%data%species == 1) THEN
                !Electric Field at particle position
                E_tot=E_ext+p(ip)%results%e
                !abs(B)
                b=sqrt(p(ip)%data%b(1)**2+p(ip)%data%b(2)**2+p(ip)%data%b(3))
                !electron gyro frequency
                omega_c=p(ip)%data%q*b/p(ip)%data%m

                IF (b<0.0000001) THEN
                    write(*,*) '##### vanishing B-Field for particle:', ip,sum(npps),my_rank,p(ip)%label
                ELSE
                    !Parallel velocity components
                    vpar(1) = (p(ip)%data%b(1)/b**2) * (p(ip)%data%b(1)*p(ip)%data%v(1) + p(ip)%data%b(2)*p(ip)%data%v(2) + p(ip)%data%b(3)*p(ip)%data%v(3))
                    vpar(2) = (p(ip)%data%b(2)/b**2) * (p(ip)%data%b(1)*p(ip)%data%v(1) + p(ip)%data%b(2)*p(ip)%data%v(2) + p(ip)%data%b(3)*p(ip)%data%v(3))
                    vpar(3) = (p(ip)%data%b(3)/b**2) * (p(ip)%data%b(1)*p(ip)%data%v(1) + p(ip)%data%b(2)*p(ip)%data%v(2) + p(ip)%data%b(3)*p(ip)%data%v(3))
                    !particle velocity without drifts (here only ExB) v - v_ExB
                    vd(1) = p(ip)%data%v(1) - (1/b**2) * (E_tot(2)*p(ip)%data%b(3)-E_tot(3)*p(ip)%data%b(2))
                    vd(2) = p(ip)%data%v(2) - (1/b**2) * (E_tot(3)*p(ip)%data%b(1)-E_tot(1)*p(ip)%data%b(3))
                    vd(3) = p(ip)%data%v(3) - (1/b**2) * (E_tot(1)*p(ip)%data%b(2)-E_tot(2)*p(ip)%data%b(1))
                    !new velocities (vpar + v_ExB)
                    p(ip)%data%v=vpar - vd + p(ip)%data%v
                    !new positions (moved by Lamor raduis to guiding centre position)
                    p(ip)%x(1) = p(ip)%x(1) - (1./(omega_c*b)) * (p(ip)%data%b(2)*vd(3) - p(ip)%data%b(3)*vd(2))
                    p(ip)%x(2) = p(ip)%x(2) - (1./(omega_c*b)) * (p(ip)%data%b(3)*vd(1) - p(ip)%data%b(1)*vd(3))
                    p(ip)%x(3) = p(ip)%x(3) - (1./(omega_c*b)) * (p(ip)%data%b(1)*vd(2) - p(ip)%data%b(2)*vd(1))
                END IF
            END IF
        END DO
    END SUBROUTINE

end module
