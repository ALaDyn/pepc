module particlehandling

    use module_pepc_types
    use helper
    implicit none

    contains

!==================================== hits_on_boundaries ===============================================================
!hits_on_boundaries identifies the particles that hit the target or that left the simulation domain at x(1)==0 (source)
!those particles are placed in the plasma again, according to the chosen source distribution
!the total charge that hit the target is also determined and the wall is charge
! Christian Salmagne 15.05.2012

    subroutine hits_on_boundaries(p)
        implicit none
        include 'mpif.h'

        type(t_particle), allocatable, intent(inout) :: p(:)
        type(t_particle), allocatable                :: p_new_e(:)
        type(t_particle), allocatable                :: p_new_i(:)

        integer :: ip, rp, rc, local_ihits,global_ihits,local_ehits,global_ehits
        integer :: local_ihits_src,global_ihits_src,local_ehits_src,global_ehits_src
        real*8 :: local_q, global_q,local_wall_q,global_wall_q                             !q transfered on target on local proc and total in this timestep, total q on wall at this moment
        integer :: new_e,new_i
        
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
        DO ip=1, np
            IF( p(rp)%x(1) > dx ) THEN                     !hit target
                if (p(rp)%data%q > 0 .and. p(rp)%label>0) local_ihits=local_ihits+1
                if (p(rp)%data%q < 0 .and. p(rp)%label>0) local_ehits=local_ehits+1
                local_q=local_q+p(rp)%data%q
                IF(rp .ne. np) p(rp) = p(np)
                np = np - 1
            ELSE IF( p(rp)%x(1) < 0.0 ) THEN               !hit source wall
                if (p(rp)%data%q > 0 .and. p(rp)%label>0) then
                    local_ihits_src=local_ihits_src+1
                end if
                if (p(rp)%data%q < 0 .and. p(rp)%label>0) then
                    local_ehits_src=local_ehits_src+1

                end if

                IF(rp .ne. np) p(rp) = p(np)
                np = np - 1
            END IF

            IF(p(rp)%x(2) < 0.0_8) THEN
                p(rp)%x(2) = p(rp)%x(2) + dy
            ELSE IF(p(rp)%x(2) > dy) THEN
                p(rp)%x(2) = p(rp)%x(2) - dy
            END IF
        
            IF(p(rp)%x(3) < 0.0_8) THEN
                p(rp)%x(3) = p(rp)%x(3) + dz
            ELSE IF(p(rp)%x(3) > dz) THEN
                p(rp)%x(3) = p(rp)%x(3) - dz
            END IF
            rp = rp + 1
      
        END DO

        call MPI_ALLREDUCE(local_ehits, global_ehits, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(local_ihits, global_ihits, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(local_ehits_src, global_ehits_src, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(local_ihits_src, global_ihits_src, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, rc)
       
        local_q = (local_ihits-local_ehits)*e
        global_q = (global_ihits-global_ehits)*e

        new_e=local_ehits+local_ehits_src
        new_i=local_ihits+local_ihits_src

        allocate(p_new_e(new_e),stat=rc)
        allocate(p_new_i(new_i),stat=rc)

        call reallocate_particles(particles,np, np+new_e+new_i)
      
        DO ip=1, new_e   !electrons first
            p_new_e(ip)%label       = next_label
            next_label        = next_label + 1
            p_new_e(ip)%data%q      = -e
            p_new_e(ip)%data%m      = me

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
            p_new_i(ip)%data%q      = e
            p_new_i(ip)%data%m      = mp

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

        deallocate(p_new_e)
        deallocate(p_new_i)

        
        DO ip=1, np                                         ! charge wall by adding charge to the wall particles
            IF (p(ip)%label<0) THEN                           ! wall particles have negative labels
                p(ip)%data%q = p(ip)%data%q + global_q/tnwp   
                local_wall_q=local_wall_q+p(ip)%data%q
            END IF
        END DO

        call MPI_ALLREDUCE(local_wall_q, global_wall_q, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)

        if(root) write(99,*)global_ehits,",",global_ihits,",",global_q,",",global_wall_q
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
            p(ip)%work        = 0

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
            p(ip)%work        = 0

            p(ip)%data%v      =0.0_8
            p(ip)%x(1)        =dx
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

        mu=0.0
        n=size(p) 
        DO ip=1, n
            sigma=sqrt(ti_ev*e/p(ip)%data%m)
            call random_number(p(ip)%x(2:3))
            p(ip)%x(2)         = p(ip)%x(2)*dy
            p(ip)%x(3)         = p(ip)%x(3)*dz

            call random_gauss_list(p(ip)%data%v(2:3),mu,sigma) 
            call random_gaussian_flux(p(ip)%data%v(1),sigma)
            call random_number(ran)
            p(ip)%x(1)=p(ip)%data%v(1)*dt*ran

        END DO     
    END SUBROUTINE


!======================================================================================

    SUBROUTINE source_emmert(p)
        type(t_particle), allocatable, intent(inout) :: p(:)
        real               :: ran
        integer            :: n,ip
        real*8             :: mu,sigma

        mu=0.0
        n=size(p) 
        DO ip=1, n
            sigma=sqrt(ti_ev*e/p(ip)%data%m)
            call random_number(p(ip)%x)
            p(ip)%x(1)         = p(ip)%x(1)*0.5*dx
            p(ip)%x(2)         = p(ip)%x(2)*dy
            p(ip)%x(3)         = p(ip)%x(3)*dz

            call random_gauss_list(p(ip)%data%v(2:3),mu,sigma) 
            call random_gaussian_flux(p(ip)%data%v(1),sigma)
            call random_number(ran)
            IF (ran>0.5) p(ip)%data%v(1)=-p(ip)%data%v(1)

        END DO     
    END SUBROUTINE

end module
