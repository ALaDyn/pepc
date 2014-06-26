!!!!!!!!!!!!!!!!!!!!
!! Diagnose module
!!
!! Enthaelt Methoden fuer die Diagnose des Plasmas
!!!!!!!!!!!!!!!!!!!!

!erster test

MODULE diagnostics


    use variables
    use helper

    implicit none

    CONTAINS
!===============================================================================

    function get_epot(p,ispecies)
        implicit none
        include 'mpif.h'

        type(t_particle),  intent(in) :: p(:)
        integer, intent(in) :: ispecies
        real(KIND=8) :: get_epot

        integer :: ip,rc
        real(KIND=8) :: esum,tesum

        esum=0.0
        tesum=0.0

        DO ip=1,size(p)
            IF(p(ip)%data%species==ispecies) THEN
                esum=esum+0.5*p(ip)%results%pot*fc*species(ispecies)%q/e
            END IF
        END DO
        esum=esum/tnpps(ispecies)

        call MPI_ALLREDUCE(esum, tesum, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, rc)

        get_epot=tesum
        return

    end function get_epot


!===============================================================================


    subroutine bin_v(ispecies,vsum_bin,v2sum_bin,n_bin)
        implicit none

        integer, intent(in) :: ispecies
        integer, intent(out):: n_bin(:,:,:)
        real(KIND=8), intent(out) :: vsum_bin(:,:,:,:)
        real(KIND=8), intent(out) :: v2sum_bin(:,:,:,:)

        integer :: ip,n
        integer                 :: ix,iy,iz

        real(KIND=8)            :: n1(3), t1(3), t2(3), B_vector(3)
        real(KIND=8)            :: eps=1.0e-10

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

        dx = (xmax-xmin)/diag_bins_x
        dy = (ymax-ymin)/diag_bins_y
        dz = (zmax-zmin)/diag_bins_z

        n = size(particles)
        n_bin = 0
        vsum_bin = 0.0_8
        v2sum_bin = 0.0_8

        DO ip=1, n
            IF (particles(ip)%data%species == ispecies) THEN
                ix = int((particles(ip)%x(1) - xmin) / dx) + 1
                iy = int((particles(ip)%x(2) - ymin) / dy) + 1
                iz = int((particles(ip)%x(3) - zmin) / dz) + 1

                n_bin(ix,iy,iz) = n_bin(ix,iy,iz) + 1
                vsum_bin(1,ix,iy,iz) = vsum_bin(1,ix,iy,iz) + particles(ip)%data%v(1)
                vsum_bin(2,ix,iy,iz) = vsum_bin(2,ix,iy,iz) + particles(ip)%data%v(2)
                vsum_bin(3,ix,iy,iz) = vsum_bin(3,ix,iy,iz) + particles(ip)%data%v(3)
                vsum_bin(4,ix,iy,iz) = vsum_bin(4,ix,iy,iz) + dotproduct(particles(ip)%data%v,n1)
                vsum_bin(5,ix,iy,iz) = vsum_bin(5,ix,iy,iz) + dotproduct(particles(ip)%data%v,t1)
                vsum_bin(6,ix,iy,iz) = vsum_bin(6,ix,iy,iz) + dotproduct(particles(ip)%data%v,t2)

                v2sum_bin(1,ix,iy,iz) = v2sum_bin(1,ix,iy,iz) + particles(ip)%data%v(1)**2                                                !vx  vx
                v2sum_bin(2,ix,iy,iz) = v2sum_bin(2,ix,iy,iz) + particles(ip)%data%v(2)**2                                                !vy  vy
                v2sum_bin(3,ix,iy,iz) = v2sum_bin(3,ix,iy,iz) + particles(ip)%data%v(3)**2                                                !vz  vz
                v2sum_bin(4,ix,iy,iz) = v2sum_bin(4,ix,iy,iz) + particles(ip)%data%v(1)*particles(ip)%data%v(2)                           !vx  vy
                v2sum_bin(5,ix,iy,iz) = v2sum_bin(5,ix,iy,iz) + particles(ip)%data%v(2)*particles(ip)%data%v(3)                           !vy  vz
                v2sum_bin(6,ix,iy,iz) = v2sum_bin(6,ix,iy,iz) + particles(ip)%data%v(3)*particles(ip)%data%v(1)                           !vz  vx

                v2sum_bin(7,ix,iy,iz) = v2sum_bin(7,ix,iy,iz) + dotproduct(particles(ip)%data%v,n1)**2                                    !v||  v||
                v2sum_bin(8,ix,iy,iz) = v2sum_bin(8,ix,iy,iz) + dotproduct(particles(ip)%data%v,t1)**2                                    !v_|_1  v_|_1
                v2sum_bin(9,ix,iy,iz) = v2sum_bin(9,ix,iy,iz) + dotproduct(particles(ip)%data%v,t2)**2                                    !v_|_2  v_|_2
                v2sum_bin(10,ix,iy,iz) = v2sum_bin(10,ix,iy,iz) + dotproduct(particles(ip)%data%v,n1)*dotproduct(particles(ip)%data%v,t1) !v||  v_|_1
                v2sum_bin(11,ix,iy,iz) = v2sum_bin(11,ix,iy,iz) + dotproduct(particles(ip)%data%v,t1)*dotproduct(particles(ip)%data%v,t2) !v_|_1  v_|_2
                v2sum_bin(12,ix,iy,iz) = v2sum_bin(12,ix,iy,iz) + dotproduct(particles(ip)%data%v,t2)*dotproduct(particles(ip)%data%v,n1) !v_|_2  v||
            END IF
        END DO
    end subroutine

!===============================================================================


!    subroutine bin_v(ispecies,npoints,x_bin,vsum_bin,v2sum_bin,n_bin)
!        implicit none

!        integer, intent(in) :: ispecies,npoints
!        real(KIND=8), intent(out) :: x_bin(:)
!        integer, intent(out):: n_bin(:)
!        real(KIND=8), intent(out) :: vsum_bin(:,:)
!        real(KIND=8), intent(out) :: v2sum_bin(:,:)

!        integer :: ip,n
!        real(KIND=8)            :: dx,x0
!        integer                 :: i

!        real(KIND=8)            :: n1(3), t1(3), t2(3), B_vector(3)
!        real(KIND=8)            :: eps=1.0e-10

!        B_vector(1) = Bx
!        B_vector(2) = By
!        B_vector(3) = Bz
!        IF (real_unequal_zero(B,eps)) THEN                                ! With B Field
!            n1 = B_vector / sqrt(dotproduct(B_vector,B_vector))           ! normal vector along B
!        ELSE                                                              ! Without B Field
!            n1(1) = 1.0_8                                                 ! normal vector along x
!            n1(2) = 0.0_8
!            n1(3) = 0.0_8
!        END IF
!        IF ((real_unequal_zero(n1(1),eps)) .OR. (real_unequal_zero(n1(2),eps))) THEN
!            t1(1) = -n1(2)                                            ! tangential vector
!            t1(2) = n1(1)
!            t1(3) = 0.0_8
!        ELSE
!            t1(1) = 1.0_8                                             ! tangential vector
!            t1(2) = 0.0_8
!            t1(3) = 0.0_8
!        END IF
!        t1 = t1 / sqrt(dotproduct(t1,t1))
!        t2(1) = n1(2)*t1(3) - n1(3)*t1(2)                             ! t2 = n1 x t1 (2nd tangential vector)
!        t2(2) = n1(3)*t1(1) - n1(1)*t1(3)
!        t2(3) = n1(1)*t1(2) - n1(2)*t1(1)
!        t2=t2/sqrt(dotproduct(t2,t2))

!        dx = (xmax-xmin)/npoints
!        x0 = xmin + dx/2.
!        x_bin = x0
!        DO i=1,npoints
!            x_bin(i) = x_bin(i) + (i-1)*dx
!        END DO

!        n = size(particles)
!        n_bin = 0
!        vsum_bin = 0.0_8
!        v2sum_bin = 0.0_8

!        DO i=1, npoints
!            DO ip=1, n
!                IF (particles(ip)%data%species == ispecies) THEN
!                    IF ((particles(ip)%x(1) < x_bin(i)+dx/2.) .AND. (particles(ip)%x(1) >= x_bin(i)-dx/2.)) THEN
!                        n_bin(i) = n_bin(i) + 1
!                        vsum_bin(1,i) = vsum_bin(1,i) + particles(ip)%data%v(1)
!                        vsum_bin(2,i) = vsum_bin(2,i) + particles(ip)%data%v(2)
!                        vsum_bin(3,i) = vsum_bin(3,i) + particles(ip)%data%v(3)
!                        vsum_bin(4,i) = vsum_bin(4,i) + dotproduct(particles(ip)%data%v,n1)
!                        vsum_bin(5,i) = vsum_bin(5,i) + dotproduct(particles(ip)%data%v,t1)
!                        vsum_bin(6,i) = vsum_bin(6,i) + dotproduct(particles(ip)%data%v,t2)
!                        v2sum_bin(1,i) = v2sum_bin(1,i) + particles(ip)%data%v(1)**2
!                        v2sum_bin(2,i) = v2sum_bin(2,i) + particles(ip)%data%v(2)**2
!                        v2sum_bin(3,i) = v2sum_bin(3,i) + particles(ip)%data%v(3)**2
!                        v2sum_bin(4,i) = v2sum_bin(4,i) + dotproduct(particles(ip)%data%v,n1)**2
!                        v2sum_bin(5,i) = v2sum_bin(5,i) + dotproduct(particles(ip)%data%v,t1)**2
!                        v2sum_bin(6,i) = v2sum_bin(6,i) + dotproduct(particles(ip)%data%v,t2)**2
!                    END IF
!                END IF
!            END DO
!        END DO
!    end subroutine

!===============================================================================

    subroutine get_probe_particles(probes,ispecies)
        implicit none
        include 'mpif.h'

        integer, intent(in)        :: ispecies
        integer                    :: rc,ip,n,i

        type(t_particle),  intent(inout),allocatable :: probes(:)


        n=size(particles)

        allocate(probes(npps(ispecies)),stat=rc)

        i=1
        DO ip=1,n
            IF (particles(ip)%data%species == ispecies) THEN
                probes(i) = particles(ip)
                i = i+1
            END IF
        END DO


    end subroutine get_probe_particles

!===============================================================================

    function get_avg_wallpotential(p,ib)
        use module_geometry
        implicit none
        include 'mpif.h'

        type(t_particle),  intent(in) :: p(:)
        integer, intent(in) :: ib
        integer :: tn_wallparticles
        integer :: n_particles
        integer :: ip, rc
        logical :: hit
        real(KIND=8) :: x_hit_rel(2)

        real(KIND=8) :: get_avg_wallpotential
        real(KIND=8) :: potsum

        tn_wallparticles = tnpps(0)
        potsum=0.
        n_particles = size(p)


        DO ip=1,n_particles
            IF(p(ip)%data%species/=0) CYCLE
            call check_hit(p(ip)%x(1),p(ip)%x(2),p(ip)%x(3),boundaries(ib),hit,x_hit_rel)
            IF (hit) potsum=potsum+p(ip)%results%pot*fc
        END DO

        potsum = potsum/tn_wallparticles

        call MPI_ALLREDUCE(potsum, get_avg_wallpotential ,1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, rc)

    end function get_avg_wallpotential

!===============================================================================

    function get_v2_mean(p,ispecies)
        implicit none
        include 'mpif.h'

        type(t_particle),  intent(in) :: p(:)
        integer, intent(in) :: ispecies
        real(KIND=8),dimension(3) :: get_v2_mean

        integer :: ip,rc
        real(KIND=8) :: v2sum(3),tv2sum(3)

        v2sum=0.0
        tv2sum=0.0

        DO ip=1,size(p)
            IF(p(ip)%data%species==ispecies) THEN
                v2sum=v2sum+p(ip)%data%v*p(ip)%data%v
            END IF
        END DO
        v2sum=v2sum/tnpps(ispecies)

        call MPI_ALLREDUCE(v2sum, tv2sum, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, rc)

        get_v2_mean=tv2sum
        return

    end function get_v2_mean

!===============================================================================

    function get_v_mean(p,ispecies)
        implicit none
        include 'mpif.h'

        type(t_particle),  intent(in) :: p(:)
        integer, intent(in) :: ispecies
        real(KIND=8),dimension(3) :: get_v_mean

        integer :: ip,rc
        real(KIND=8) :: vsum(3),tvsum(3)

        vsum=0.0
        tvsum=0.0

        DO ip=1,size(p)
            IF(p(ip)%data%species==ispecies) THEN
                vsum=vsum+p(ip)%data%v
            END IF
        END DO
        vsum=vsum/tnpps(ispecies)

        call MPI_ALLREDUCE(vsum, tvsum, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, rc)

        get_v_mean=tvsum
        return

    end function get_v_mean


!===============================================================================

    subroutine write_velocities(p)

        implicit none 

        type(t_particle), allocatable, intent(in) :: p(:)
        integer                                     :: n,ip

        n=size(p)

        open(20, file='test.ion', action='write')
        open(21, file='test.electron', action='write')
        DO ip=1,n
            IF (p(ip)%label>0) THEN                  !only particles  
                IF (p(ip)%data%q>0) THEN             !ions
                    write(20,*)p(ip)%label,p(ip)%x,p(ip)%data%v,p(ip)%data%q,p(ip)%data%m
                ELSE
                    write(20,*)p(ip)%label,p(ip)%x,p(ip)%data%v,p(ip)%data%q,p(ip)%data%m
                END IF     
            END IF
        END DO 
        close(20)
        close(21)

    end subroutine 




END MODULE
