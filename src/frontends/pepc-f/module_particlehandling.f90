module particlehandling

    use module_pepc_types
    use helper
    implicit none

    contains

!==================================== hits_on_boundaries ==========================================================
!hits_on_boundaries identifies the particles that hit the target or that left the simulation domain at 
!x(1)==0 (source)
!those particles are placed in the plasma again, according to the chosen source distribution
!the total charge that hit the target is also determined and the wall is charged
! Christian Salmagne 15.05.2012

    subroutine hits_on_boundaries_with_open_sides(p)
        implicit none
        include 'mpif.h'

        type(t_particle), allocatable, intent(inout) :: p(:)
        type(t_particle), allocatable                :: p_new_e(:)
        type(t_particle), allocatable                :: p_new_i(:)

        integer :: ip, rp, rc, local_ihits,global_ihits,local_ehits,global_ehits
        integer :: local_ihits_side, local_ehits_side, global_ihits_side, global_ehits_side
        integer :: local_ihits_src,global_ihits_src,local_ehits_src,global_ehits_src
        real*8 :: local_q, global_q,local_wall_q,global_wall_q                             !q transfered on target on local proc and total in this timestep, total q on wall at this moment
        integer :: new_e,new_i,new_e_global,new_i_global
        
        local_q=0.0_8
        global_q=0.0_8
        local_wall_q=0.0_8
        global_wall_q=0.0_8
        local_ihits=0
        global_ihits=0
        local_ehits=0
        global_ehits=0
        local_ihits_src=0
        global_ihits_src=0
        local_ehits_src=0
        global_ehits_src=0
        local_ihits_side=0
        global_ihits_side=0
        local_ehits_side=0
        global_ehits_side=0

        rp = 1
        DO WHILE (rp <= np)
            IF( hit_wall(p(rp)) ) THEN                     !hit target
                if (p(rp)%data%q > 0.) local_ihits=local_ihits+1
                if (p(rp)%data%q < 0.) local_ehits=local_ehits+1
                local_q=local_q+p(rp)%data%q
                IF(rp .ne. np) THEN
                    DO
                        IF ( hit_wall(p(np)) .or. hit_src(p(np)).or. hit_side(p(np)) ) THEN !third condition for open walls
                            IF (np==rp) EXIT
                            IF ( hit_wall(p(np)) ) THEN
                                IF (p(np)%data%q > 0.) local_ihits=local_ihits+1  
                                IF (p(np)%data%q < 0.) local_ehits=local_ehits+1
                                local_q=local_q+p(rp)%data%q
                            ELSE IF ( hit_src(p(np)) ) THEN
                                IF (p(np)%data%q > 0.) local_ihits_src=local_ihits_src+1  
                                IF (p(np)%data%q < 0.) local_ehits_src=local_ehits_src+1 
                            ! new part for open walls
                            ELSE IF ( hit_side(p(np)) ) THEN
                                IF (p(np)%data%q > 0.) local_ihits_side=local_ihits_side+1
                                IF (p(np)%data%q < 0.) local_ehits_side=local_ehits_side+1
                            ! end of new part
                            END IF
                            np = np - 1
                        ELSE
                            p(rp) = p(np)
                            EXIT
                        END IF
                    END DO
                END IF

                np = np - 1
            ELSE IF( hit_src(p(rp)) ) THEN               !hit source wall
                if (p(rp)%data%q > 0.) local_ihits_src=local_ihits_src+1
                if (p(rp)%data%q < 0.) local_ehits_src=local_ehits_src+1
                IF(rp .ne. np) THEN
                    DO
                        IF ( hit_wall(p(np)) .or. hit_src(p(np)).or. hit_side(p(np)) ) THEN !third condition for open walls
                            IF (np==rp) EXIT
                            IF ( hit_wall(p(np)) ) THEN
                                IF (p(np)%data%q > 0.) local_ihits=local_ihits+1  
                                IF (p(np)%data%q < 0.) local_ehits=local_ehits+1
                                local_q=local_q+p(rp)%data%q
                            ELSE IF ( hit_src(p(np)) ) THEN
                                IF (p(np)%data%q > 0.) local_ihits_src=local_ihits_src+1  
                                IF (p(np)%data%q < 0.) local_ehits_src=local_ehits_src+1
                            ! new part for open walls
                            ELSE IF ( hit_side(p(np)) ) THEN
                                IF (p(np)%data%q > 0.) local_ihits_side=local_ihits_side+1
                                IF (p(np)%data%q < 0.) local_ehits_side=local_ehits_side+1
                            ! end of new part
                            END IF
                            np = np - 1
                        ELSE
                            p(rp) = p(np)
                            EXIT
                        END IF
                    END DO
                END IF

                np = np - 1

            !open walls, if particles hit side walls, they are gone and will be recreated
            ELSE IF (hit_side(p(rp))) THEN
                if (p(rp)%data%q > 0.) local_ihits_side=local_ihits_side+1
                if (p(rp)%data%q < 0.) local_ehits_side=local_ehits_side+1
                IF(rp .ne. np) THEN
                    DO
                        IF ( hit_wall(p(np)) .or. hit_src(p(np)) .or. hit_side(p(np)) ) THEN
                            IF (np==rp) EXIT
                            IF ( hit_wall(p(np)) ) THEN
                                IF (p(np)%data%q > 0.) local_ihits=local_ihits+1
                                IF (p(np)%data%q < 0.) local_ehits=local_ehits+1
                                local_q=local_q+p(rp)%data%q
                            ELSE IF ( hit_src(p(np)) ) THEN
                                IF (p(np)%data%q > 0.) local_ihits_src=local_ihits_src+1
                                IF (p(np)%data%q < 0.) local_ehits_src=local_ehits_src+1
                            ELSE IF ( hit_side(p(np)) ) THEN
                                IF (p(np)%data%q > 0.) local_ihits_side=local_ihits_side+1
                                IF (p(np)%data%q < 0.) local_ehits_side=local_ehits_side+1
                            END IF
                            np = np - 1
                        ELSE
                            p(rp) = p(np)
                            EXIT
                        END IF
                    END DO
                END IF

                np = np - 1
            END IF

            !Periodic behaviour (seems to cause heating of the electrons if periodic condiations are not used for fields as well)
            !IF(p(rp)%x(2) < ymin) THEN
            !    p(rp)%x(2) = p(rp)%x(2) + dy
            !ELSE IF(p(rp)%x(2) > ymax) THEN
            !    p(rp)%x(2) = p(rp)%x(2) - dy
            !END IF
        
            !IF(p(rp)%x(3) < zmin) THEN
            !    p(rp)%x(3) = p(rp)%x(3) + dz
            !ELSE IF(p(rp)%x(3) > zmax) THEN
            !    p(rp)%x(3) = p(rp)%x(3) - dz
            !END IF
            
            rp = rp + 1
      
        END DO

        call MPI_ALLREDUCE(local_ehits, global_ehits, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(local_ihits, global_ihits, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(local_ehits_src, global_ehits_src, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(local_ihits_src, global_ihits_src, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(local_ehits_side, global_ehits_side, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(local_ihits_side, global_ihits_side, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, rc)
       
        local_q = (local_ihits-local_ehits)*e*fsup
        global_q = (global_ihits-global_ehits)*e*fsup

        new_e_global = global_ehits + global_ehits_src + global_ehits_side
        new_i_global = global_ihits + global_ihits_src + global_ihits_side

        new_e = new_e_global / n_ranks
        new_i = new_i_global / n_ranks
        if (my_rank .eq. (n_ranks-1)) then
            new_e = new_e + MOD(new_e_global, n_ranks)
            new_i = new_i + MOD(new_i_global, n_ranks)
        end if

        next_label = next_label + my_rank * (new_e_global/n_ranks + new_i_global/n_ranks) 

        allocate(p_new_e(new_e),stat=rc)
        allocate(p_new_i(new_i),stat=rc)

        call reallocate_particles(particles,np, np+new_e+new_i)
      
        DO ip=1, new_e   !electrons first
            p_new_e(ip)%label       = next_label
            next_label        = next_label + 1
            p_new_e(ip)%data%q      = -e*fsup
            p_new_e(ip)%data%m      = me*fsup

            p_new_e(ip)%results%e   = 0.0_8
            p_new_e(ip)%results%pot = 0.0_8
            p_new_e(ip)%work        = 1.0_8

            p_new_e(ip)%data%B(1)=Bx
            p_new_e(ip)%data%B(2)=By
            p_new_e(ip)%data%B(3)=Bz
        END DO
        call source(p_new_e,quelltyp)

        DO ip=1, new_i   !ions
            p_new_i(ip)%label       = next_label
            next_label        = next_label + 1  
            p_new_i(ip)%data%q      = e*fsup
            p_new_i(ip)%data%m      = mp*fsup

            p_new_i(ip)%results%e   = 0.0_8
            p_new_i(ip)%results%pot = 0.0_8
            p_new_i(ip)%work        = 1.0_8

            p_new_i(ip)%data%B(1)=Bx
            p_new_i(ip)%data%B(2)=By
            p_new_i(ip)%data%B(3)=Bz
        END DO
        call source(p_new_i,quelltyp)
 

        p(np+1:np+new_e)=p_new_e(:)
        p(np+new_e+1:np+new_e+new_i)=p_new_i(:)
        np = np + new_e + new_i

        call MPI_ALLREDUCE(np, tnp, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_BCAST(next_label, 1, MPI_INTEGER, n_ranks-1, MPI_COMM_WORLD, rc)

        deallocate(p_new_e)
        deallocate(p_new_i)

        DO ip=1, np                                         ! charge wall by adding charge to the wall particles
            IF (p(ip)%label<0) THEN                           ! wall particles have negative labels
                p(ip)%data%q = p(ip)%data%q + global_q/tnwp   
                local_wall_q=local_wall_q+p(ip)%data%q
            END IF
        END DO
        call MPI_ALLREDUCE(local_wall_q, global_wall_q, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        
        if(root) open(unit=99,file='current_on_wall.dat',status='UNKNOWN',position='APPEND')   
        if(root) write(99,*)global_ehits_src,",",global_ihits_src,",",global_ehits,",",global_ihits,",",global_q,",",global_wall_q
        if(root)close(99)
        if(root) write(*,'(a,i12)')    " == [filter] total number of particles            : ", tnp
        if(root) write(*,'(a,i12)')    " == [target_hits] number of ion hits              : ", global_ihits
        if(root) write(*,'(a,i12)')    " == [target_hits] number of electron hits         : ", global_ehits
        if(root) write(*,'(a,i12)')    " == [source_hits] number of ion hits              : ", global_ihits_src
        if(root) write(*,'(a,i12)')    " == [source_hits] number of electron hits         : ", global_ehits_src
        if(root) write(*,'(a,i12)')    " == [side_hits] number of ion hits                : ", global_ihits_side
        if(root) write(*,'(a,i12)')    " == [side_hits] number of electron hits           : ", global_ehits_side

    END SUBROUTINE


!======================================================================================

    subroutine hits_on_boundaries(p)
        implicit none
        include 'mpif.h'

        type(t_particle), allocatable, intent(inout) :: p(:)
        type(t_particle), allocatable                :: p_new_e(:)
        type(t_particle), allocatable                :: p_new_i(:)

        integer :: ip, rp, rc, local_ihits,global_ihits,local_ehits,global_ehits
        integer :: local_ihits_src,global_ihits_src,local_ehits_src,global_ehits_src
        real*8 :: local_q, global_q,local_wall_q,global_wall_q                             !q transfered on target on local proc and total in this timestep, total q on wall at this moment
        integer :: new_e,new_i,new_e_global,new_i_global

        local_q=0.0_8
        global_q=0.0_8
        local_wall_q=0.0_8
        global_wall_q=0.0_8
        local_ihits=0
        global_ihits=0
        local_ehits=0
        global_ehits=0
        local_ihits_src=0
        global_ihits_src=0
        local_ehits_src=0
        global_ehits_src=0


        rp = 1
        DO WHILE (rp <= np)
            IF( hit_wall(p(rp)) ) THEN                     !hit target
                if (p(rp)%data%q > 0.) local_ihits=local_ihits+1
                if (p(rp)%data%q < 0.) local_ehits=local_ehits+1
                local_q=local_q+p(rp)%data%q
                IF(rp .ne. np) THEN
                    DO
                        IF ( hit_wall(p(np)) .or. hit_src(p(np)) ) THEN
                            IF (np==rp) EXIT
                            IF ( hit_wall(p(np)) ) THEN
                                IF (p(np)%data%q > 0.) local_ihits=local_ihits+1
                                IF (p(np)%data%q < 0.) local_ehits=local_ehits+1
                                local_q=local_q+p(rp)%data%q
                            ELSE IF ( hit_src(p(np)) ) THEN
                                IF (p(np)%data%q > 0.) local_ihits_src=local_ihits_src+1
                                IF (p(np)%data%q < 0.) local_ehits_src=local_ehits_src+1
                            END IF
                            np = np - 1
                        ELSE
                            p(rp) = p(np)
                            EXIT
                        END IF
                    END DO
                END IF

                np = np - 1
            ELSE IF( hit_src(p(rp)) ) THEN               !hit source wall
                if (p(rp)%data%q > 0.) local_ihits_src=local_ihits_src+1
                if (p(rp)%data%q < 0.) local_ehits_src=local_ehits_src+1
                IF(rp .ne. np) THEN
                    DO
                        IF ( hit_wall(p(np)) .or. hit_src(p(np))) THEN
                            IF (np==rp) EXIT
                            IF ( hit_wall(p(np)) ) THEN
                                IF (p(np)%data%q > 0.) local_ihits=local_ihits+1
                                IF (p(np)%data%q < 0.) local_ehits=local_ehits+1
                                local_q=local_q+p(rp)%data%q
                            ELSE IF ( hit_src(p(np)) ) THEN
                                IF (p(np)%data%q > 0.) local_ihits_src=local_ihits_src+1
                                IF (p(np)%data%q < 0.) local_ehits_src=local_ehits_src+1
                            END IF
                            np = np - 1
                        ELSE
                            p(rp) = p(np)
                            EXIT
                        END IF
                    END DO
                END IF

                np = np - 1
            END IF

            !Periodic behaviour (seems to cause heating of the electrons if periodic condiations are not used for fields as well)
            IF(p(rp)%x(2) < ymin) THEN
                p(rp)%x(2) = p(rp)%x(2) + dy
            ELSE IF(p(rp)%x(2) > ymax) THEN
                p(rp)%x(2) = p(rp)%x(2) - dy
            END IF

            IF(p(rp)%x(3) < zmin) THEN
                p(rp)%x(3) = p(rp)%x(3) + dz
            ELSE IF(p(rp)%x(3) > zmax) THEN
                p(rp)%x(3) = p(rp)%x(3) - dz
            END IF

            rp = rp + 1

        END DO

        call MPI_ALLREDUCE(local_ehits, global_ehits, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(local_ihits, global_ihits, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(local_ehits_src, global_ehits_src, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(local_ihits_src, global_ihits_src, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, rc)


        local_q = (local_ihits-local_ehits)*e*fsup
        global_q = (global_ihits-global_ehits)*e*fsup

        new_e_global = global_ehits + global_ehits_src
        new_i_global = global_ihits + global_ihits_src

        new_e = new_e_global / n_ranks
        new_i = new_i_global / n_ranks
        if (my_rank .eq. (n_ranks-1)) then
            new_e = new_e + MOD(new_e_global, n_ranks)
            new_i = new_i + MOD(new_i_global, n_ranks)
        end if

        next_label = next_label + my_rank * (new_e_global/n_ranks + new_i_global/n_ranks)

        allocate(p_new_e(new_e),stat=rc)
        allocate(p_new_i(new_i),stat=rc)

        call reallocate_particles(particles,np, np+new_e+new_i)

        DO ip=1, new_e   !electrons first
            p_new_e(ip)%label       = next_label
            next_label        = next_label + 1
            p_new_e(ip)%data%q      = -e*fsup
            p_new_e(ip)%data%m      = me*fsup

            p_new_e(ip)%results%e   = 0.0_8
            p_new_e(ip)%results%pot = 0.0_8
            p_new_e(ip)%work        = 1.0_8

            p_new_e(ip)%data%B(1)=Bx
            p_new_e(ip)%data%B(2)=By
            p_new_e(ip)%data%B(3)=Bz
        END DO
        call source(p_new_e,quelltyp)

        DO ip=1, new_i   !ions
            p_new_i(ip)%label       = next_label
            next_label        = next_label + 1
            p_new_i(ip)%data%q      = e*fsup
            p_new_i(ip)%data%m      = mp*fsup

            p_new_i(ip)%results%e   = 0.0_8
            p_new_i(ip)%results%pot = 0.0_8
            p_new_i(ip)%work        = 1.0_8

            p_new_i(ip)%data%B(1)=Bx
            p_new_i(ip)%data%B(2)=By
            p_new_i(ip)%data%B(3)=Bz
        END DO
        call source(p_new_i,quelltyp)


        p(np+1:np+new_e)=p_new_e(:)
        p(np+new_e+1:np+new_e+new_i)=p_new_i(:)
        np = np + new_e + new_i

        call MPI_ALLREDUCE(np, tnp, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_BCAST(next_label, 1, MPI_INTEGER, n_ranks-1, MPI_COMM_WORLD, rc)

        deallocate(p_new_e)
        deallocate(p_new_i)

        DO ip=1, np                                         ! charge wall by adding charge to the wall particles
            IF (p(ip)%label<0) THEN                           ! wall particles have negative labels
                p(ip)%data%q = p(ip)%data%q + global_q/tnwp
                local_wall_q=local_wall_q+p(ip)%data%q
            END IF
        END DO
        call MPI_ALLREDUCE(local_wall_q, global_wall_q, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)

        if(root) open(unit=99,file='current_on_wall.dat',status='UNKNOWN',position='APPEND')
        if(root) write(99,*)global_ehits_src,",",global_ihits_src,",",global_ehits,",",global_ihits,",",global_q,",",global_wall_q
        if(root)close(99)
        if(root) write(*,'(a,i12)')    " == [filter] total number of particles            : ", tnp
        if(root) write(*,'(a,i12)')    " == [target_hits] number of ion hits              : ", global_ihits
        if(root) write(*,'(a,i12)')    " == [target_hits] number of electron hits         : ", global_ehits
        if(root) write(*,'(a,i12)')    " == [source_hits] number of ion hits              : ", global_ihits_src
        if(root) write(*,'(a,i12)')    " == [source_hits] number of electron hits         : ", global_ehits_src

    END SUBROUTINE


!======================================================================================
    SUBROUTINE init_particles(p)
        implicit none
    
        type(t_particle), allocatable, intent(inout) :: p(:)
        integer :: ip


        DO ip=1, npp
            p(ip)%data%B(1)=Bx
            p(ip)%data%B(2)=By
            p(ip)%data%B(3)=Bz

            p(ip)%label       = my_rank * (tnpp / n_ranks) + ip
            p(ip)%data%q      = (-1.0_8 + 2.0_8*MOD(p(ip)%label,2))*e*fsup
            p(ip)%data%m      = me*fsup
            if(p(ip)%data%q .gt. 0.0) p(ip)%data%m = mp*fsup
            p(ip)%results%e   = 0.0_8
            p(ip)%results%pot = 0.0_8
            p(ip)%work        = 1.0_8

        END DO

        call source(p,quelltyp)
        
        next_label = tnpp+1
  
    END SUBROUTINE init_particles


!======================================================================================
    subroutine init_wall_particles(p)
        implicit none
    
        type(t_particle), allocatable, intent(inout) :: p(:)
        integer :: ip

        DO ip=1, nwp
            p(ip)%data%B=0.0_8

            p(ip)%label       = -(my_rank * (tnwp / n_ranks) + ip)
            p(ip)%data%q      = 0.0_8
            p(ip)%data%m      = 10.0_8

            p(ip)%results%e   = 0.0_8
            p(ip)%results%pot = 0.0_8
            p(ip)%work        = 1.0_8

            p(ip)%data%v      =0.0_8
            p(ip)%x(1)        =xmax
            p(ip)%x(2)        =wall_pos(-p(ip)%label,1)
            p(ip)%x(3)        =wall_pos(-p(ip)%label,2)

        END DO
 
    end subroutine init_wall_particles

!======================================================================================

    SUBROUTINE source(p,i)
        type(t_particle), allocatable, intent(inout) :: p(:)
        integer, intent(in)                          :: i

        quelle: SELECT case(i)
            case(0)
                call source_berberich(p) 
                !write(*,*)"BERBERICH"
            case(1)
                call source_emmert(p)
                !write(*,*)"EMMERT"
            case default
                !do nothing
        END SELECT quelle

    END SUBROUTINE

!======================================================================================

    SUBROUTINE source_berberich(p)
        type(t_particle), allocatable, intent(inout) :: p(:)
        real               :: ran
        integer            :: n,ip
        real*8             :: mu,sigma
        real*8             :: vpar(3),b,E_tot(3),E_ext(3),vd(3),omega_c

        !Just done now, because Eext is not introduced yet
        E_ext=0.0_8
        mu=0.0
        n=size(p) 
        DO ip=1, n
            sigma=sqrt(ti_ev*e/(p(ip)%data%m/fsup))
            p(ip)%x(2)=rnd_num()
            p(ip)%x(3)=rnd_num()
            !call random_number(p(ip)%x(2:3))
            p(ip)%x(2)         = p(ip)%x(2)*dy + ymin
            p(ip)%x(3)         = p(ip)%x(3)*dz + zmin

            call random_gauss_list(p(ip)%data%v(2:3),mu,sigma) 
            call random_gaussian_flux(p(ip)%data%v(1),sigma)
            ran=rnd_num()
            !call random_number(ran)
            p(ip)%x(1)=p(ip)%data%v(1)*dt*ran + xmin

            !If treated in guding centre approximation, electrons just need other startign values in
            !the Boris scheme (see Benjamin Berberichs Diss (pp. 44-56)
            IF (guiding_centre_electrons) THEN
                IF ((p(ip)%label > 0) .AND. (p(ip)%data%q < 0.)) THEN
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
            END IF

        END DO     
    END SUBROUTINE


!======================================================================================

    SUBROUTINE source_emmert(p)
        type(t_particle), allocatable, intent(inout) :: p(:)
        real               :: ran
        integer            :: n,ip
        real*8             :: mu,sigma
        real*8             :: vpar(3),b,E_tot(3),E_ext(3),vd(3),omega_c

        !Just done now, because Eext is not introduced yet
        E_ext=0.0_8

        mu=0.0
        n=size(p) 
        DO ip=1, n
            sigma=sqrt(ti_ev*e / (p(ip)%data%m/fsup))
            p(ip)%x(1)=rnd_num()
            p(ip)%x(2)=rnd_num()
            p(ip)%x(3)=rnd_num()
            !call random_number(p(ip)%x)
            p(ip)%x(1)         = p(ip)%x(1)*0.5*dx + xmin 
            p(ip)%x(2)         = p(ip)%x(2)*dy + ymin
            p(ip)%x(3)         = p(ip)%x(3)*dz + zmin

            IF (p(ip)%data%q>0) THEN
                call random_gauss_list(p(ip)%data%v(2:3),mu,sigma) 
                call random_gaussian_flux(p(ip)%data%v(1),sigma)
                ran=rnd_num()
                !call random_number(ran)
                IF (ran>0.5) p(ip)%data%v(1)=-p(ip)%data%v(1)
            ELSE
                call random_gauss_list(p(ip)%data%v(1:3),mu,sigma) 
            END IF

            !If treated in guding centre approximation, electrons just need other startign values in
            !the Boris scheme (see Benjamin Berberichs Diss (pp. 44-56)
            IF (guiding_centre_electrons) THEN
                IF ((p(ip)%label > 0) .AND. (p(ip)%data%q < 0.)) THEN
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
            END IF

        END DO     
    END SUBROUTINE

end module
