module particlehandling

    use module_pepc_types
    use helper
    use module_geometry
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
        real*8 :: xp(3)


        IF (ANY(boundaries(1:nb)%type==4)) allocate(p_hits_logical_sheath(0:nspecies-1,nb,np),stat=rc)

        rp=np
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
                        p(rp) = p(np)
                        np = np - 1
                        ib = 0
                        rp = rp - 1
                    ELSE IF (boundaries(ib)%type==4) THEN                              !logical sheath
                        IF (boundaries(ib)%reflux_particles) reflux(p(rp)%data%species,ib)=reflux(p(rp)%data%species,ib)+1
                        p_hits_logical_sheath(p(rp)%data%species,ib,hits(p(rp)%data%species,ib)) = p(rp)
                        p(rp) = p(np)
                        np = np-1
                        ib = 0
                        rp = rp - 1
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
                        p(rp) = p(np)
                        np = np - 1
                        ib = 0
                        rp = rp - 1
                    END IF
                END IF
                ib = ib+1
            END DO
            rp = rp-1
        END DO

        call reallocate_particles(p,np, np)

        IF (ANY(boundaries(1:nb)%type==4)) call treat_logical_sheath_boundaries(p,hits,reflux,p_hits_logical_sheath)
        IF (ANY(boundaries(1:nb)%type==4)) deallocate(p_hits_logical_sheath)

    END SUBROUTINE

!======================================================================================

    SUBROUTINE treat_logical_sheath_boundaries(p,hits,reflux,p_hits_logical_sheath)
        implicit none
        include 'mpif.h'

        type(t_particle), allocatable, intent(inout) :: p(:),p_hits_logical_sheath(:,:,:)
        integer, intent(inout) :: hits(0:,:),reflux(0:,:)

        integer ib,ispecies,i,ip,rc
        integer :: thits(0:nspecies-1,nb,0:n_ranks-1),displs(0:n_ranks-1)
        real*8 :: v_perp_logical_sheath(0:nspecies-1,nb,np)
        real*8 :: t_v_perp_logical_sheath(0:nspecies-1,nb,tnp)
        real*8,allocatable :: v_perp_temp(:)
        real*8 :: vcut_el(nb)
        integer :: ihits,ehits


        thits=0
        vcut_el=0
        DO ib = 1,nb
            IF (boundaries(ib)%type==4) THEN
                DO ispecies = 0,nspecies-1
                    IF (species(ispecies)%physical_particle) THEN
                        call MPI_GATHER(hits(ispecies,ib),1,MPI_INTEGER, thits(ispecies,ib,:), 1, MPI_INTEGER, 0,MPI_COMM_WORLD, rc)
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
                                displs(i)=displs(i-1)+thits(ispecies,ib,i-1)
                            END DO
                        END IF
                        DO ip = 1,hits(ispecies,ib)
                            v_perp_logical_sheath(ispecies,ib,ip)=-dotproduct(p_hits_logical_sheath(ispecies,ib,ip)%data%v &
                                                                              ,boundaries(ib)%n)
                        END DO
                        IF (species(ispecies)%indx==1) THEN !sorting only for electrons
                            call MPI_GATHERV(v_perp_logical_sheath(ispecies,ib,:),hits(ispecies,ib),&
                                             MPI_DOUBLE_PRECISION,t_v_perp_logical_sheath(ispecies,ib,:),&
                                             thits(ispecies,ib,:),displs,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD, rc)
                            IF (root) THEN
                                allocate(v_perp_temp(sum(thits(ispecies,ib,:))),stat=rc)
                                v_perp_temp=t_v_perp_logical_sheath(ispecies,ib,1:sum(thits(ispecies,ib,:)))
                                call QsortC(v_perp_temp)
                                t_v_perp_logical_sheath(ispecies,ib,1:sum(thits(ispecies,ib,:)))=v_perp_temp
                                ehits=sum(thits(ispecies,ib,:))
                                ihits=0
                                DO i=0,nspecies-1
                                    IF (species(i)%physical_particle .and. species(i)%q > 0) THEN
                                        ihits=ihits+sum(thits(i,ib,:))
                                    END IF
                                END DO
                                vcut_el(ib)=v_perp_temp(ehits-ihits)
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
                DO ispecies = 0,nspecies-1
                    IF (species(ispecies)%physical_particle .and. species(ispecies)%indx==1) THEN
                        call repell_electrons_at_logical_sheath(p,hits,reflux,vcut_el(ib),p_hits_logical_sheath,ib,ispecies)
                    END IF
                END DO
            END IF
        END DO


    END SUBROUTINE treat_logical_sheath_boundaries

!======================================================================================

    SUBROUTINE repell_electrons_at_logical_sheath(p,hits,reflux,vcut_el,p_hits_logical_sheath,ib,ispecies)

        implicit none
        include 'mpif.h'

        type(t_particle), allocatable, intent(inout) :: p(:),p_hits_logical_sheath(:,:,:)
        integer, intent(inout) :: hits(0:,:),reflux(0:,:)
        real(KIND=8),intent(in) :: vcut_el
        integer, intent(in) :: ib,ispecies

        integer :: ip
        real(KIND=8) :: v_perp,xp(3)

        DO ip=1,hits(ispecies,ib)
            v_perp=-dotproduct(p_hits_logical_sheath(ispecies,ib,ip)%data%v,boundaries(ib)%n)
            IF (v_perp <= vcut_el) THEN !reflect electrons to get ambipolar flux
                np=np+1
                p(np)=p_hits_logical_sheath(ispecies,ib,ip)
                xp = p(np)%x - boundaries(ib)%x0
                p(np)%x = p(np)%x - 2.*dotproduct(xp,boundaries(ib)%n)*boundaries(ib)%n
                p(np)%data%v = p(np)%data%v - 2.*dotproduct(p(np)%data%v,boundaries(ib)%n)*boundaries(ib)%n
                hits(ispecies,ib) = hits(ispecies,ib) - 1
                IF (reflux(ispecies,ib)>0) reflux(ispecies,ib) = reflux(ispecies,ib) - 1
            END IF
        END DO

    END SUBROUTINE repell_electrons_at_logical_sheath

!======================================================================================

    SUBROUTINE recycling(p,treflux,treflux_next_ts)

        implicit none
        include 'mpif.h'

        type(t_particle), allocatable, intent(inout) :: p(:)
        integer, intent(in) :: treflux(0:,:)
        integer, intent(inout) :: treflux_next_ts(0:,:)
        integer :: reflux(0:nspecies-1,1:nb)
        integer :: rc,ispecies,ib,new_particles


        reflux=0

        DO ispecies=0,nspecies-1
            DO ib=1,nb
                IF (treflux(ispecies,ib)/=0) THEN
                    reflux(ispecies,ib) = treflux(ispecies,ib) / n_ranks
                    next_label = next_label + my_rank*reflux(ispecies,ib)
                    IF (my_rank .eq. (n_ranks-1)) THEN
                        reflux(ispecies,ib) = reflux(ispecies,ib) + MOD(treflux(ispecies,ib), n_ranks)
                    END IF

                    call reflux_particles(p,reflux(ispecies,ib),ispecies,ib)
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

                    call reflux_particles(p,new_particles,ispecies,0)
                    call MPI_BCAST(next_label, 1, MPI_INTEGER, n_ranks-1, MPI_COMM_WORLD, rc)
                END IF
            END IF
        END DO

    END SUBROUTINE recycling

!======================================================================================

    SUBROUTINE reflux_particles(p,new_particles,ispecies,ib)

        implicit none

        type(t_particle), allocatable, intent(inout) :: p(:)
        integer, intent(in) :: new_particles,ispecies,ib

        integer rc,ip
        type(t_particle), allocatable                :: p_new(:)


        IF (new_particles==0) return

        allocate(p_new(new_particles),stat=rc)
        call reallocate_particles(p,np, np+new_particles)

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
        call source(p_new,quelltyp)

        p(np+1:np+new_particles)=p_new(:)
        np = np + new_particles

        deallocate(p_new)

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

        DO ip=1, np                                            ! charge wall by adding charge to the wall particles
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

    SUBROUTINE get_number_of_particles()
        implicit none

        integer :: ip

        npps=0
        DO ip=1,np
           npps(particles(ip)%data%species) = npps(particles(ip)%data%species) + 1
        END DO

    END SUBROUTINE get_number_of_particles


!======================================================================================
    subroutine hits_on_boundaries(p)
        implicit none
        include 'mpif.h'

        type(t_particle), allocatable, intent(inout) :: p(:)

        integer      :: rc

        integer      :: hits(0:nspecies-1,1:nb),thits(0:nspecies-1,1:nb)
        integer      :: reflux(0:nspecies-1,1:nb),treflux(0:nspecies-1,1:nb),treflux_next_ts(0:nspecies-1,1:nb)
        integer      :: ispecies,ib


        hits=0
        reflux=0
        thits=0
        treflux=0
        treflux_next_ts=0

        if(root) write(*,'(a)') " == [hits_on_boundaries] count hits and recycle "

        call set_need_to_reflux()
        call count_hits_and_remove_particles(p,hits,reflux)


        DO ispecies=0,nspecies-1
            DO ib=1,nb
                call MPI_ALLREDUCE(hits(ispecies,ib), thits(ispecies,ib), 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, rc)
                call MPI_ALLREDUCE(reflux(ispecies,ib), treflux(ispecies,ib), 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, rc)
            END DO
        END DO


        call recycling(p,treflux,treflux_next_ts)

        call charge_wall(p,thits)

        call get_number_of_particles()

        call MPI_ALLREDUCE(npps, tnpps, nspecies, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, rc)
        tnp=SUM(tnpps)

        call set_recycling_output_values(thits,treflux)

    END SUBROUTINE


!======================================================================================
    SUBROUTINE init_particles(p)
        implicit none
    
        type(t_particle), intent(inout) :: p(:)
        integer :: ip,ispecies,i

        ip=0
        DO ispecies=1,nspecies-1
            DO i=1, npps(ispecies)
                ip = ip + 1
                p(ip)%data%B(1)=Bx
                p(ip)%data%B(2)=By
                p(ip)%data%B(3)=Bz
                p(ip)%label       = my_rank * (SUM(tnpps(1:nspecies-1)) / n_ranks) + ip
                p(ip)%data%q      = species(ispecies)%q*fsup    !(-1.0_8 + 2.0_8*MOD(p(ip)%label,2))*e*fsup
                p(ip)%data%m      = species(ispecies)%m*fsup    !me*fsup
                p(ip)%data%species= species(ispecies)%indx      !1
                !if(p(ip)%data%q .gt. 0.0) p(ip)%data%m = mp*fsup
                !if(p(ip)%data%q .gt. 0.0) p(ip)%data%species=2
                p(ip)%results%e   = 0.0_8
                p(ip)%results%pot = 0.0_8
                p(ip)%work        = 1.0_8
            END DO
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

        call source(p,quelltyp)
        
        next_label = SUM(tnpps)+1
  
    END SUBROUTINE init_particles


!======================================================================================
    subroutine init_wall_particles(p)
        implicit none
    
        type(t_particle), intent(inout) :: p(:)
        type(t_boundary) :: wall
        real :: ran1,ran2
        integer :: ip,ib,n

        n=size(p)
        DO ip=1, n
            p(ip)%data%B=0.0_8

            p(ip)%label       = -(my_rank * (tnpps(0) / n_ranks) + ip)
            p(ip)%data%q      = 0.0_8
            p(ip)%data%m      = 10.0_8

            p(ip)%results%e   = 0.0_8
            p(ip)%results%pot = 0.0_8
            p(ip)%work        = 1.0_8
            p(ip)%data%species= 0

            p(ip)%data%v      =0.0_8

            !noch auskommentieren test
            !p(ip)%x(1)        =xmax
            !p(ip)%x(2)        =wall_pos(-p(ip)%label,1)
            !p(ip)%x(3)        =wall_pos(-p(ip)%label,2)
            DO ib=1,nb
                IF (ANY(boundaries(ib)%wp_labels==p(ip)%label)) THEN
                    wall=boundaries(ib)
                    ran1=rnd_num()
                    ran2=rnd_num()
                    p(ip)%x = wall%x0 + ran1*wall%e1 +ran2*wall%e2
                END IF
            END DO
        END DO
 
    end subroutine init_wall_particles

!======================================================================================
    subroutine redistribute_wall_particles(p)
        implicit none

        type(t_particle), intent(inout) :: p(:)
        real :: ran1,ran2
        integer :: ip,ib

        DO ip=1, np
            IF (p(ip)%data%species/=0) CYCLE
            DO ib=1,nb
                IF (ANY(boundaries(ib)%wp_labels==p(ip)%label)) THEN
                    ran1=rnd_num()
                    ran2=rnd_num()
                    p(ip)%x = boundaries(ib)%x0 + ran1*boundaries(ib)%e1 +ran2*boundaries(ib)%e2
                END IF
            END DO
        END DO

    end subroutine redistribute_wall_particles

!======================================================================================

    SUBROUTINE source(p,i)
        type(t_particle), intent(inout) :: p(:)
        integer, intent(in)                           :: i

        quelle: SELECT case(i)
            case(0)
                call source_berberich(p) 
            case(1)
                call source_emmert(p)
            case(2)
                call source_bissel_johnson(p)
            case(3)  !test case with uniform plasma
                call source_uniform_gaussian_plasma(p)
            case default
                !do nothing
        END SELECT quelle

    END SUBROUTINE

!======================================================================================

    SUBROUTINE source_berberich(p)
        type(t_particle), intent(inout) :: p(:)
        real               :: ran,ran1,ran2
        integer            :: n,ip,ib
        real*8             :: mu,sigma


        mu=0.0
        n=size(p) 
        DO ip=1, n
            IF (p(ip)%data%species==0) THEN
                p(ip)%data%v = 0.0_8
                DO ib=1,nb
                    IF (ANY(boundaries(ib)%wp_labels==p(ip)%label)) THEN
                        ran1=rnd_num()
                        ran2=rnd_num()
                        p(ip)%x = boundaries(ib)%x0 + ran1*boundaries(ib)%e1 +ran2*boundaries(ib)%e2
                    END IF
                END DO
            ELSE
                sigma=sqrt(species(p(ip)%data%species)%t_src*e/(p(ip)%data%m/fsup))
                call random_gauss_list(p(ip)%data%v(2:3),mu,sigma)
                call random_gaussian_flux(p(ip)%data%v(1),sigma)
                ran=rnd_num()

                !p(ip)%x(2)=rnd_num()*dy + ymin
                !p(ip)%x(3)=rnd_num()*dz + zmin
                !p(ip)%x(1)=p(ip)%data%v(1)*dt*ran + xmin

                ran1=rnd_num()
                ran2=rnd_num()
                p(ip)%x = (boundaries(src_boundary)%x0 + ran1*boundaries(src_boundary)%e1 + ran2*boundaries(src_boundary)%e2)
                p(ip)%x = p(ip)%x + boundaries(src_boundary)%n*dotproduct(boundaries(src_boundary)%n,p(ip)%data%v)*dt*ran



                !If treated in guding centre approximation, electrons just need other startign values in
                !the Boris scheme (see Benjamin Berberichs Diss (pp. 44-56)
                IF (guiding_centre_electrons) THEN
                    call transform_electrons_to_gc(p)
                END IF
            END IF

        END DO     
    END SUBROUTINE


!======================================================================================

    SUBROUTINE source_emmert(p)
        type(t_particle), intent(inout) :: p(:)
        real               :: ran,ran1,ran2
        integer            :: n,ip,ib
        real*8             :: mu,sigma


        mu=0.0
        n=size(p) 
        DO ip=1, n
            IF (p(ip)%data%species==0) THEN
                p(ip)%data%v = 0.0_8
                DO ib=1,nb
                    IF (ANY(boundaries(ib)%wp_labels==p(ip)%label)) THEN
                        ran1=rnd_num()
                        ran2=rnd_num()
                        p(ip)%x = boundaries(ib)%x0 + ran1*boundaries(ib)%e1 +ran2*boundaries(ib)%e2
                    END IF
                END DO
            ELSE
                sigma=sqrt(species(p(ip)%data%species)%t_src*e / (p(ip)%data%m/fsup))

                IF (p(ip)%data%species==2) THEN    !ions
                    call random_gauss_list(p(ip)%data%v(2:3),mu,sigma)
                    call random_gaussian_flux(p(ip)%data%v(1),sigma)
                    ran=rnd_num()
                    IF (ran>0.5) p(ip)%data%v(1)=-p(ip)%data%v(1)
                ELSE IF (p(ip)%data%species==1) THEN   !electrons
                    call random_gauss_list(p(ip)%data%v(1:3),mu,sigma)
                END IF

                p(ip)%x=x0_src + rnd_num()*e1_src + rnd_num()*e2_src + rnd_num()*e3_src

                !If treated in guding centre approximation, electrons just need other startign values in
                !the Boris scheme (see Benjamin Berberichs Diss (pp. 44-56)
                IF (guiding_centre_electrons) THEN
                    call transform_electrons_to_gc(p)
                END IF
            END IF
        END DO
    END SUBROUTINE

!======================================================================================

    SUBROUTINE source_bissel_johnson(p)
        type(t_particle), intent(inout) :: p(:)
        integer            :: n,ip,ib
        real*8             :: mu,sigma
        real               :: ran1,ran2

        mu=0.0
        n=size(p)
        DO ip=1, n
            IF (p(ip)%data%species==0) THEN
                p(ip)%data%v = 0.0_8
                DO ib=1,nb
                    IF (ANY(boundaries(ib)%wp_labels==p(ip)%label)) THEN
                        ran1=rnd_num()
                        ran2=rnd_num()
                        p(ip)%x = boundaries(ib)%x0 + ran1*boundaries(ib)%e1 +ran2*boundaries(ib)%e2
                    END IF
                END DO
            ELSE
                sigma=sqrt(species(p(ip)%data%species)%t_src*e / (p(ip)%data%m/fsup))
                call random_gauss_list(p(ip)%data%v(1:3),mu,sigma)
                !p(ip)%x(1)=rnd_num()
                !p(ip)%x(2)=rnd_num()
                !p(ip)%x(3)=rnd_num()
                !p(ip)%x(1)         = p(ip)%x(1)*0.5*dx + xmin
                !p(ip)%x(2)         = p(ip)%x(2)*dy + ymin
                !p(ip)%x(3)         = p(ip)%x(3)*dz + zmin

                p(ip)%x=x0_src + rnd_num()*e1_src + rnd_num()*e2_src + rnd_num()*e3_src


                !If treated in guding centre approximation, electrons just need other startign values in
                !the Boris scheme (see Benjamin Berberichs Diss (pp. 44-56)
                IF (guiding_centre_electrons) THEN
                    call transform_electrons_to_gc(p)
                END IF
            END IF
        END DO
    END SUBROUTINE

!======================================================================================

    SUBROUTINE source_uniform_gaussian_plasma(p)
        type(t_particle), intent(inout) :: p(:)
        integer            :: n,ip,ib
        real*8             :: mu,sigma
        real               :: ran1,ran2

        mu=0.0
        n=size(p)
        DO ip=1, n
            IF (p(ip)%data%species==0) THEN
                p(ip)%data%v = 0.0_8
                DO ib=1,nb
                    IF (ANY(boundaries(ib)%wp_labels==p(ip)%label)) THEN
                        ran1=rnd_num()
                        ran2=rnd_num()
                        p(ip)%x = boundaries(ib)%x0 + ran1*boundaries(ib)%e1 +ran2*boundaries(ib)%e2
                    END IF
                END DO
            ELSE
                sigma=sqrt(species(p(ip)%data%species)%t_src*e / (p(ip)%data%m/fsup))
                p(ip)%x=x0_src + rnd_num()*e1_src + rnd_num()*e2_src + rnd_num()*e3_src

                call random_gauss_list(p(ip)%data%v(1:3),mu,sigma)

                !If treated in guding centre approximation, electrons just need other starting values in
                !the Boris scheme (see Benjamin Berberichs Diss (pp. 44-56)
                IF (guiding_centre_electrons) THEN
                    call transform_electrons_to_gc(p)
                END IF
            END IF
        END DO
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
                    write(*,*) '##### vanishing B-Field for particle:', ip,np,my_rank,p(ip)%label
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
