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

    real(KIND=8), allocatable :: initial_velocities(:,:)

    CONTAINS
!===============================================================================

    subroutine store_initial_velocities(p)
        implicit none
        include 'mpif.h'

        type(t_particle),  intent(in) :: p(:)

        integer :: ip, ierr
        real(KIND=8), allocatable :: initial_velocities_local(:,:)

        if (.not. allocated(initial_velocities)) then
            allocate(initial_velocities(3, sum(tnpps)))
        end if
        allocate(initial_velocities_local(3, sum(tnpps)))

        initial_velocities = 0.0_8
        initial_velocities_local = 0.0_8

        do ip=1, sum(npps)
            initial_velocities_local(:, p(ip)%label) = p(ip)%data%v
        end do

        call MPI_ALLREDUCE(initial_velocities_local, initial_velocities, 3*sum(tnpps), MPI_REAL8, MPI_SUM, &
                           MPI_COMM_WORLD, ierr)

    end subroutine store_initial_velocities
!===============================================================================

    subroutine hockney_diag(p, avg_1, avg_2, avg_3, avg_4, avg_5)
        implicit none
        include 'mpif.h'

        type(t_particle),  intent(in) :: p(:)
        real(KIND=8), intent(inout) :: avg_1(:), avg_2(:), avg_3(:), avg_4(:), avg_5(:)
        !avg_1 = <beta(t)^2>^0.5
        !avg_2 = <vperp(t)^2>^0.5
        !avg_3 = <vpar(t)>
        !avg_4 = <h(t)>
        !avg_5 = <h(t)^2>^0.5
        !< x > = 1/N * sum(1..N) (x)
        !h = 0.5 * m * (v(t)^2 - v(0)^2)

        integer :: ip, ierr, ispecies
        real(KIND=8) :: v0_p(3), h_p, vpar_p, vperp_p, v_p(3), e0_p(3)
        real(KIND=8) :: beta_p !deflection angle; angle between v0 and v

        real(KIND=8) :: avg_beta2(nspecies-1), avg_vperp2(nspecies-1), avg_vpar(nspecies-1),&
                        avg_h(nspecies-1), avg_h2(nspecies-1)

        avg_beta2 = 0.
        avg_vperp2 = 0.
        avg_vpar = 0.
        avg_h = 0.
        avg_h2 = 0.

        avg_1 = 0.
        avg_2 = 0.
        avg_3 = 0.
        avg_4 = 0.
        avg_5 = 0.


        do ip=1, sum(npps)
            ispecies = p(ip)%data%species
            if (.not. species(ispecies)%physical_particle) cycle
            v0_p = initial_velocities(:, p(ip)%label)
            e0_p = v0_p / norm(v0_p)
            v_p = p(ip)%data%v

            h_p = 0.5 * species(ispecies)%m/e * (sum(v_p**2) - sum(v0_p**2))
            vpar_p = dotproduct(v_p, e0_p)
            vperp_p = sqrt(sum(v_p**2) - vpar_p**2)
            beta_p = acos( dotproduct(v0_p, v_p) / (norm(v0_p) * norm(v_p)) ) !range of acos func: 0..pi

            avg_beta2(ispecies) = avg_beta2(ispecies) + beta_p**2
            avg_vperp2(ispecies) = avg_vperp2(ispecies) + vperp_p**2
            avg_vpar(ispecies) = avg_vpar(ispecies) + vpar_p
            avg_h(ispecies) = avg_h(ispecies) + h_p
            avg_h2(ispecies) = avg_h2(ispecies) + h_p**2
        end do

        call MPI_REDUCE(avg_beta2, avg_1, nspecies-1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(avg_vperp2, avg_2, nspecies-1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(avg_vpar, avg_3, nspecies-1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(avg_h, avg_4, nspecies-1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(avg_h2, avg_5, nspecies-1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

        avg_1 = sqrt(avg_1(:) / tnpps(1:))
        avg_2 = sqrt(avg_2(:) / tnpps(1:))
        avg_3 = avg_3(:) / tnpps(1:)
        avg_4 = avg_4(:) / tnpps(1:)
        avg_5 = sqrt(avg_5 / tnpps(1:))


    end subroutine hockney_diag

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
! This routine bins v and v^2 in cylinder coordinates.
! The cylinder axis is at y,z = 0, xmin < x <xmax
! diag_bins_x bins for x: xmin < x < xmax
! diag_bins_y bins for r: 0 < r < ymax
! diag_bins_z bins for phi: 0 < phi < 2Pi

    subroutine fill_data_bins_cylindrical(ispecies,nbs,dbs)
        implicit none

        integer, intent(in)         :: ispecies
        integer, intent(inout)      :: nbs(:,:,:)
        real(KIND=8), intent(inout) :: dbs(:,:,:,:)
        real(KIND=8)                :: cylinder_coords(3)
        real(KIND=8)                :: delx,delr,deltheta

        integer :: ip,n
        integer                 :: ix,ir,itheta

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

        delx = (xmax-xmin)/diag_bins_x
        delr = ymax/diag_bins_y
        deltheta = 2*pi/diag_bins_z

        n = size(particles)

        DO ip=1, n
            IF (particles(ip)%data%species == ispecies) THEN
                cylinder_coords(1) = particles(ip)%x(1)
                cylinder_coords(2) = sqrt(particles(ip)%x(2)**2 + particles(ip)%x(3)**2)
                IF (real_unequal_zero(cylinder_coords(2),eps)) THEN
                    IF (particles(ip)%x(3) >= 0._8) THEN
                        cylinder_coords(3) = acos(particles(ip)%x(2) / cylinder_coords(2))
                    ELSE
                        cylinder_coords(3) = 2*pi - acos(particles(ip)%x(2) / cylinder_coords(2))
                    END IF
                ELSE
                    cylinder_coords(3) = 0.
                END IF

                ix = int((cylinder_coords(1) - xmin) / delx) + 1
                ir = int(cylinder_coords(2) / delr) + 1
                itheta = int(cylinder_coords(3) / deltheta) + 1

                !number of particles
                nbs(ix,ir,itheta) = nbs(ix,ir,itheta) + 1

                !vx, vy, vz
                dbs(1,ix,ir,itheta) = dbs(1,ix,ir,itheta) + particles(ip)%data%v(1)
                dbs(2,ix,ir,itheta) = dbs(2,ix,ir,itheta) + particles(ip)%data%v(2)
                dbs(3,ix,ir,itheta) = dbs(3,ix,ir,itheta) + particles(ip)%data%v(3)

                !vpar, vperp1, vperp2
                dbs(4,ix,ir,itheta) = dbs(4,ix,ir,itheta) + dotproduct(particles(ip)%data%v,n1)
                dbs(5,ix,ir,itheta) = dbs(5,ix,ir,itheta) + dotproduct(particles(ip)%data%v,t1)
                dbs(6,ix,ir,itheta) = dbs(6,ix,ir,itheta) + dotproduct(particles(ip)%data%v,t2)

                !2nd moments: vxvx vyvy vzvz vxvy vyvz vzvx
                dbs(7,ix,ir,itheta) = dbs(7,ix,ir,itheta) + particles(ip)%data%v(1)**2                                                !vx  vx
                dbs(8,ix,ir,itheta) = dbs(8,ix,ir,itheta) + particles(ip)%data%v(2)**2                                                !vy  vy
                dbs(9,ix,ir,itheta) = dbs(9,ix,ir,itheta) + particles(ip)%data%v(3)**2                                                !vz  vz
                dbs(10,ix,ir,itheta) = dbs(10,ix,ir,itheta) + particles(ip)%data%v(1)*particles(ip)%data%v(2)                         !vx  vy
                dbs(11,ix,ir,itheta) = dbs(11,ix,ir,itheta) + particles(ip)%data%v(2)*particles(ip)%data%v(3)                         !vy  vz
                dbs(12,ix,ir,itheta) = dbs(12,ix,ir,itheta) + particles(ip)%data%v(3)*particles(ip)%data%v(1)                         !vz  vx

                !2nd moments: vparvpar, vperp1vperp1, vperp2vperp2, vparvperp1, vperp1vperp2, vperp2,vpar
                dbs(13,ix,ir,itheta) = dbs(13,ix,ir,itheta) + dotproduct(particles(ip)%data%v,n1)**2                                  !v||  v||
                dbs(14,ix,ir,itheta) = dbs(14,ix,ir,itheta) + dotproduct(particles(ip)%data%v,t1)**2                                  !v_|_1  v_|_1
                dbs(15,ix,ir,itheta) = dbs(15,ix,ir,itheta) + dotproduct(particles(ip)%data%v,t2)**2                                  !v_|_2  v_|_2
                dbs(16,ix,ir,itheta) = dbs(16,ix,ir,itheta) + dotproduct(particles(ip)%data%v,n1)*dotproduct(particles(ip)%data%v,t1) !v||  v_|_1
                dbs(17,ix,ir,itheta) = dbs(17,ix,ir,itheta) + dotproduct(particles(ip)%data%v,t1)*dotproduct(particles(ip)%data%v,t2) !v_|_1  v_|_2
                dbs(18,ix,ir,itheta) = dbs(18,ix,ir,itheta) + dotproduct(particles(ip)%data%v,t2)*dotproduct(particles(ip)%data%v,n1) !v_|_2  v||

                !phi
                dbs(19,ix,ir,itheta) = dbs(19,ix,ir,itheta) + particles(ip)%results%pot*fc

                !Ex Ey Ez
                dbs(20:22,ix,ir,itheta) = dbs(20:22,ix,ir,itheta) + particles(ip)%results%E*fc

                !Epar, Eperp1, Eperp2
                dbs(23,ix,ir,itheta) = dbs(23,ix,ir,itheta) + dotproduct(particles(ip)%results%E,n1)*fc
                dbs(24,ix,ir,itheta) = dbs(24,ix,ir,itheta) + dotproduct(particles(ip)%results%E,t1)*fc
                dbs(25,ix,ir,itheta) = dbs(25,ix,ir,itheta) + dotproduct(particles(ip)%results%E,t2)*fc

                !2nd moment: phiphi
                dbs(26,ix,ir,itheta) = dbs(26,ix,ir,itheta) + (particles(ip)%results%pot*fc)**2

                !2nd moments: ExEx, EyEy, EzEz, ExEy, EyEz, EzEx
                dbs(27:29,ix,ir,itheta) = dbs(27:29,ix,ir,itheta) + (particles(ip)%results%E*fc) * (particles(ip)%results%E*fc)
                dbs(30,ix,ir,itheta) = dbs(30,ix,ir,itheta) + particles(ip)%results%E(1)*fc * particles(ip)%results%E(2)*fc
                dbs(31,ix,ir,itheta) = dbs(31,ix,ir,itheta) + particles(ip)%results%E(2)*fc * particles(ip)%results%E(3)*fc
                dbs(32,ix,ir,itheta) = dbs(32,ix,ir,itheta) + particles(ip)%results%E(3)*fc * particles(ip)%results%E(1)*fc

                !2nd moments: EparEpar, Eperp1Eperp1, Eperp2Eperp2, EparEperp1, Eperp1Eperp2, Eperp2Epar
                dbs(33,ix,ir,itheta) = dbs(33,ix,ir,itheta) + dotproduct(particles(ip)%results%E,n1)**2
                dbs(34,ix,ir,itheta) = dbs(34,ix,ir,itheta) + dotproduct(particles(ip)%results%E,t1)**2
                dbs(35,ix,ir,itheta) = dbs(35,ix,ir,itheta) + dotproduct(particles(ip)%results%E,t2)**2
                dbs(36,ix,ir,itheta) = dbs(36,ix,ir,itheta) + dotproduct(particles(ip)%results%E,n1)*dotproduct(particles(ip)%results%E,t1)
                dbs(37,ix,ir,itheta) = dbs(37,ix,ir,itheta) + dotproduct(particles(ip)%results%E,t1)*dotproduct(particles(ip)%results%E,t2)
                dbs(38,ix,ir,itheta) = dbs(38,ix,ir,itheta) + dotproduct(particles(ip)%results%E,t2)*dotproduct(particles(ip)%results%E,n1)

                !age
                dbs(39,ix,ir,itheta) = dbs(39,ix,ir,itheta) + particles(ip)%data%age

            END IF
        END DO
    end subroutine

!===============================================================================


    subroutine fill_data_bins(ispecies,nbs,dbs)
        implicit none

        integer, intent(in)         :: ispecies
        integer, intent(inout)      :: nbs(:,:,:)
        real(KIND=8), intent(inout) :: dbs(:,:,:,:)

        integer :: ip,n
        integer                 :: ix,iy,iz
        real(KIND=8)            :: delx,dely,delz

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

        delx = (xmax-xmin)/diag_bins_x
        dely = (ymax-ymin)/diag_bins_y
        delz = (zmax-zmin)/diag_bins_z

        n = size(particles)

        DO ip=1, n
            IF (particles(ip)%data%species == ispecies) THEN
                ix = int((particles(ip)%x(1) - xmin) / delx) + 1
                iy = int((particles(ip)%x(2) - ymin) / dely) + 1
                iz = int((particles(ip)%x(3) - zmin) / delz) + 1

                nbs(ix,iy,iz) = nbs(ix,iy,iz) + 1

                !vx, vy, vz
                dbs(1,ix,iy,iz) = dbs(1,ix,iy,iz) + particles(ip)%data%v(1)
                dbs(2,ix,iy,iz) = dbs(2,ix,iy,iz) + particles(ip)%data%v(2)
                dbs(3,ix,iy,iz) = dbs(3,ix,iy,iz) + particles(ip)%data%v(3)

                !vpar, vperp1, vperp2
                dbs(4,ix,iy,iz) = dbs(4,ix,iy,iz) + dotproduct(particles(ip)%data%v,n1)
                dbs(5,ix,iy,iz) = dbs(5,ix,iy,iz) + dotproduct(particles(ip)%data%v,t1)
                dbs(6,ix,iy,iz) = dbs(6,ix,iy,iz) + dotproduct(particles(ip)%data%v,t2)

                !2nd moments: vxvx vyvy vzvz vxvy vyvz vzvx
                dbs(7,ix,iy,iz) = dbs(7,ix,iy,iz) + particles(ip)%data%v(1)**2                                                !vx  vx
                dbs(8,ix,iy,iz) = dbs(8,ix,iy,iz) + particles(ip)%data%v(2)**2                                                !vy  vy
                dbs(9,ix,iy,iz) = dbs(9,ix,iy,iz) + particles(ip)%data%v(3)**2                                                !vz  vz
                dbs(10,ix,iy,iz) = dbs(10,ix,iy,iz) + particles(ip)%data%v(1)*particles(ip)%data%v(2)                         !vx  vy
                dbs(11,ix,iy,iz) = dbs(11,ix,iy,iz) + particles(ip)%data%v(2)*particles(ip)%data%v(3)                         !vy  vz
                dbs(12,ix,iy,iz) = dbs(12,ix,iy,iz) + particles(ip)%data%v(3)*particles(ip)%data%v(1)                         !vz  vx

                !2nd moments: vparvpar, vperp1vperp1, vperp2vperp2, vparvperp1, vperp1vperp2, vperp2,vpar
                dbs(13,ix,iy,iz) = dbs(13,ix,iy,iz) + dotproduct(particles(ip)%data%v,n1)**2                                  !v||  v||
                dbs(14,ix,iy,iz) = dbs(14,ix,iy,iz) + dotproduct(particles(ip)%data%v,t1)**2                                  !v_|_1  v_|_1
                dbs(15,ix,iy,iz) = dbs(15,ix,iy,iz) + dotproduct(particles(ip)%data%v,t2)**2                                  !v_|_2  v_|_2
                dbs(16,ix,iy,iz) = dbs(16,ix,iy,iz) + dotproduct(particles(ip)%data%v,n1)*dotproduct(particles(ip)%data%v,t1) !v||  v_|_1
                dbs(17,ix,iy,iz) = dbs(17,ix,iy,iz) + dotproduct(particles(ip)%data%v,t1)*dotproduct(particles(ip)%data%v,t2) !v_|_1  v_|_2
                dbs(18,ix,iy,iz) = dbs(18,ix,iy,iz) + dotproduct(particles(ip)%data%v,t2)*dotproduct(particles(ip)%data%v,n1) !v_|_2  v||

                !phi
                dbs(19,ix,iy,iz) = dbs(19,ix,iy,iz) + particles(ip)%results%pot*fc

                !Ex Ey Ez
                dbs(20:22,ix,iy,iz) = dbs(20:22,ix,iy,iz) + particles(ip)%results%E*fc

                !Epar, Eperp1, Eperp2
                dbs(23,ix,iy,iz) = dbs(23,ix,iy,iz) + dotproduct(particles(ip)%results%E,n1)*fc
                dbs(24,ix,iy,iz) = dbs(24,ix,iy,iz) + dotproduct(particles(ip)%results%E,t1)*fc
                dbs(25,ix,iy,iz) = dbs(25,ix,iy,iz) + dotproduct(particles(ip)%results%E,t2)*fc

                !2nd moment: phiphi
                dbs(26,ix,iy,iz) = dbs(26,ix,iy,iz) + (particles(ip)%results%pot*fc)**2

                !2nd moments: ExEx, EyEy, EzEz, ExEy, EyEz, EzEx
                dbs(27:29,ix,iy,iz) = dbs(27:29,ix,iy,iz) + (particles(ip)%results%E*fc) * (particles(ip)%results%E*fc)
                dbs(30,ix,iy,iz) = dbs(30,ix,iy,iz) + particles(ip)%results%E(1)*fc * particles(ip)%results%E(2)*fc
                dbs(31,ix,iy,iz) = dbs(31,ix,iy,iz) + particles(ip)%results%E(2)*fc * particles(ip)%results%E(3)*fc
                dbs(32,ix,iy,iz) = dbs(32,ix,iy,iz) + particles(ip)%results%E(3)*fc * particles(ip)%results%E(1)*fc

                !2nd moments: EparEpar, Eperp1Eperp1, Eperp2Eperp2, EparEperp1, Eperp1Eperp2, Eperp2Epar
                dbs(33,ix,iy,iz) = dbs(33,ix,iy,iz) + dotproduct(particles(ip)%results%E,n1)**2
                dbs(34,ix,iy,iz) = dbs(34,ix,iy,iz) + dotproduct(particles(ip)%results%E,t1)**2
                dbs(35,ix,iy,iz) = dbs(35,ix,iy,iz) + dotproduct(particles(ip)%results%E,t2)**2
                dbs(36,ix,iy,iz) = dbs(36,ix,iy,iz) + dotproduct(particles(ip)%results%E,n1)*dotproduct(particles(ip)%results%E,t1)
                dbs(37,ix,iy,iz) = dbs(37,ix,iy,iz) + dotproduct(particles(ip)%results%E,t1)*dotproduct(particles(ip)%results%E,t2)
                dbs(38,ix,iy,iz) = dbs(38,ix,iy,iz) + dotproduct(particles(ip)%results%E,t2)*dotproduct(particles(ip)%results%E,n1)

                !age
                dbs(39,ix,iy,iz) = dbs(39,ix,iy,iz) + particles(ip)%data%age
            END IF
        END DO
    end subroutine


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
