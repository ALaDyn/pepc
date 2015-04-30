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

    subroutine hockney_diag(p, avg_1, avg_2, avg_3, avg_4, avg_5, avg_6, avg_7, &
                            avg_8, avg_9, avg_10, avg_11, avg_12, avg_13, avg_14, avg_15, avg_16, avg_fields)
        implicit none
        include 'mpif.h'

        type(t_particle),  intent(in) :: p(:)
        real(KIND=8), intent(inout) :: avg_1(:), avg_2(:), avg_3(:), avg_4(:), avg_5(:), avg_6(:), avg_7(:)
        real(KIND=8), intent(inout) :: avg_8(:), avg_9(:), avg_10(:), avg_11(:), avg_12(:), avg_13(:), avg_14(:)
        real(KIND=8), intent(inout) :: avg_15(:), avg_16(:), avg_fields(:,:)
        !avg_1 = <beta(t)>
        !avg_2 = <|beta(t)|>
        !avg_3 = <beta(t)^2>^0.5
        !avg_4 = <h(t)>
        !avg_5 = <|h(t)|>
        !avg_6 = <h(t)^2>^0.5
        !avg_7 = <vpar(t)>
        !avg_8 = <|vpar(t)|>
        !avg_9 = <vpar(t)^2>^0.5
        !avg_10 = <|vperp(t)|>
        !avg_11 = <vperp(t)^2>^0.5
        !avg_12 = <|v(t)|>
        !avg_13 = <v(t)^2>^0.5
        !avg_14 = <Ekin_par(t)> = 1/2 * m * <vpar(t)^2>
        !avg_15 = <Ekin_perp(t)> = 1/2 * m * <vperp(t)^2>
        !avg_16 = <Ekin(t)> = <Ekin_par(t)> + <Ekin_perp(t)>
        !avg_fields(1,:) = <Ex(t)>
        !avg_fields(2,:) = <Ey(t)>
        !avg_fields(3,:) = <Ez(t)>
        !avg_fields(4,:) = <ExEx(t)>
        !avg_fields(5,:) = <EyEy(t)>
        !avg_fields(6,:) = <EzEz(t)>
        !avg_fields(7,:) = <ExEy(t)>
        !avg_fields(8,:) = <EyEz(t)>
        !avg_fields(9,:) = <EzEx(t)>
        !avg_fields(10,:) = <|E)(t|>
        !avg_fields(11,:) = <Phi(t)>
        !avg_fields(12,:) = <PhiPhi(t)>

        !with:
        !< x > = 1/N * sum(1..N) (x)
        !h = 0.5 * m * (v(t)^2 - v(0)^2)

        integer :: ip, ierr, ispecies, indx
        real(KIND=8) :: v0_p(3), h_p, vpar_p, vperp_p, v_p(3), e0_p(3)
        real(KIND=8) :: beta_p !deflection angle; angle between v0 and v

        real(KIND=8) :: avg_beta(nspecies-1), avg_absbeta(nspecies-1), avg_beta2(nspecies-1)
        real(KIND=8) :: avg_h(nspecies-1), avg_absh(nspecies-1), avg_h2(nspecies-1)
        real(KIND=8) :: avg_vpar(nspecies-1), avg_absvpar(nspecies-1), avg_vpar2(nspecies-1)
        real(KIND=8) :: avg_absvperp(nspecies-1), avg_vperp2(nspecies-1)
        real(KIND=8) :: avg_absv(nspecies-1)
        real(KIND=8) :: avg_fields_l(12, nspecies-1)

        avg_beta = 0.
        avg_absbeta = 0.
        avg_beta2 = 0.
        avg_h = 0.
        avg_absh = 0.
        avg_h2 = 0.
        avg_vpar = 0.
        avg_absvpar = 0.
        avg_vpar2 = 0.
        avg_absvperp = 0.
        avg_vperp2 = 0.
        avg_absv = 0.

        avg_fields_l = 0.

        avg_1 = 0.
        avg_2 = 0.
        avg_3 = 0.
        avg_4 = 0.
        avg_5 = 0.
        avg_6 = 0.
        avg_7 = 0.
        avg_8 = 0.
        avg_9 = 0.
        avg_10 = 0.
        avg_11 = 0.
        avg_12 = 0.
        avg_13 = 0.
        avg_14 = 0.
        avg_15 = 0.
        avg_16 = 0.
        avg_fields = 0.

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

            avg_beta(ispecies) = avg_beta(ispecies) + beta_p
            avg_absbeta(ispecies) = avg_absbeta(ispecies) + abs(beta_p)
            avg_beta2(ispecies) = avg_beta2(ispecies) + beta_p**2

            avg_h(ispecies) = avg_h(ispecies) + h_p
            avg_absh(ispecies) = avg_absh(ispecies) + abs(h_p)
            avg_h2(ispecies) = avg_h2(ispecies) + h_p**2

            avg_vpar(ispecies) = avg_vpar(ispecies) + vpar_p
            avg_absvpar(ispecies) = avg_absvpar(ispecies) + abs(vpar_p)
            avg_vpar2(ispecies) = avg_vpar2(ispecies) + vpar_p**2

            avg_absvperp(ispecies) = avg_absvperp(ispecies) + vperp_p
            avg_vperp2(ispecies) = avg_vperp2(ispecies) + vperp_p**2

            avg_absv(ispecies) = avg_absv(ispecies) + norm(v_p)

            avg_fields_l(1,ispecies) = avg_fields_l(1,ispecies) + p(ip)%results%E(1) * fc
            avg_fields_l(2,ispecies) = avg_fields_l(2,ispecies) + p(ip)%results%E(2) * fc
            avg_fields_l(3,ispecies) = avg_fields_l(3,ispecies) + p(ip)%results%E(3) * fc
            avg_fields_l(4,ispecies) = avg_fields_l(4,ispecies) + p(ip)%results%E(1)*p(ip)%results%E(1) * fc**2
            avg_fields_l(5,ispecies) = avg_fields_l(5,ispecies) + p(ip)%results%E(2)*p(ip)%results%E(2) * fc**2
            avg_fields_l(6,ispecies) = avg_fields_l(6,ispecies) + p(ip)%results%E(3)*p(ip)%results%E(3) * fc**2
            avg_fields_l(7,ispecies) = avg_fields_l(7,ispecies) + p(ip)%results%E(1)*p(ip)%results%E(2) * fc**2
            avg_fields_l(8,ispecies) = avg_fields_l(8,ispecies) + p(ip)%results%E(2)*p(ip)%results%E(3) * fc**2
            avg_fields_l(9,ispecies) = avg_fields_l(9,ispecies) + p(ip)%results%E(3)*p(ip)%results%E(1) * fc**2
            avg_fields_l(10,ispecies) = avg_fields_l(10,ispecies) + norm(p(ip)%results%E) * fc
            avg_fields_l(11,ispecies) = avg_fields_l(11,ispecies) + p(ip)%results%pot * fc
            avg_fields_l(12,ispecies) = avg_fields_l(12,ispecies) + p(ip)%results%pot**2 * fc**2
        end do

        call MPI_REDUCE(avg_beta, avg_1, nspecies-1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(avg_absbeta, avg_2, nspecies-1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(avg_beta2, avg_3, nspecies-1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(avg_h, avg_4, nspecies-1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(avg_absh, avg_5, nspecies-1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(avg_h2, avg_6, nspecies-1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(avg_vpar, avg_7, nspecies-1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(avg_absvpar, avg_8, nspecies-1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(avg_vpar2, avg_9, nspecies-1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(avg_absvperp, avg_10, nspecies-1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(avg_vperp2, avg_11, nspecies-1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(avg_absv, avg_12, nspecies-1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(avg_fields_l, avg_fields, (nspecies-1)*12, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)


        avg_14 = 0.5 * species(1:)%m/e * avg_9(:) / tnpps(1:)
        avg_15 = 0.5 * species(1:)%m/e * avg_11(:) / tnpps(1:)
        avg_16 = avg_14 + avg_15
        avg_13 = sqrt( (avg_9(:) + avg_11(:)) / tnpps(1:) )


        avg_1 = avg_1(:) / tnpps(1:)
        avg_2 = avg_2(:) / tnpps(1:)
        avg_3 = sqrt(avg_3(:) / tnpps(1:))
        avg_4 = avg_4(:) / tnpps(1:)
        avg_5 = avg_5(:) / tnpps(1:)
        avg_6 = sqrt(avg_6(:) / tnpps(1:))
        avg_7 = avg_7(:) / tnpps(1:)
        avg_8 = avg_8(:) / tnpps(1:)
        avg_9 = sqrt(avg_9(:) / tnpps(1:))
        avg_10 = avg_10(:) / tnpps(1:)
        avg_11 = sqrt(avg_11(:) / tnpps(1:))
        avg_12 = avg_12(:) / tnpps(1:)
        do indx = 1, 11
            avg_fields(indx,:) = avg_fields(indx, :) / tnpps(1:)
        end do

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
                dbs(33,ix,ir,itheta) = dbs(33,ix,ir,itheta) + dotproduct(particles(ip)%results%E,n1)**2 * fc**2
                dbs(34,ix,ir,itheta) = dbs(34,ix,ir,itheta) + dotproduct(particles(ip)%results%E,t1)**2 * fc**2
                dbs(35,ix,ir,itheta) = dbs(35,ix,ir,itheta) + dotproduct(particles(ip)%results%E,t2)**2 * fc**2
                dbs(36,ix,ir,itheta) = dbs(36,ix,ir,itheta) + dotproduct(particles(ip)%results%E,n1)*dotproduct(particles(ip)%results%E,t1) * fc**2
                dbs(37,ix,ir,itheta) = dbs(37,ix,ir,itheta) + dotproduct(particles(ip)%results%E,t1)*dotproduct(particles(ip)%results%E,t2) * fc**2
                dbs(38,ix,ir,itheta) = dbs(38,ix,ir,itheta) + dotproduct(particles(ip)%results%E,t2)*dotproduct(particles(ip)%results%E,n1) * fc**2

                !age
                dbs(39,ix,ir,itheta) = dbs(39,ix,ir,itheta) + particles(ip)%data%age

            END IF
        END DO
    end subroutine

!===============================================================================
    subroutine fill_velocity_bins(ispecies, nbs, dbs)
        implicit none

        integer, intent(in)         :: ispecies
        integer, intent(inout)      :: nbs(:,:)
        real(KIND=8), intent(inout) :: dbs(:,:,:)

        integer :: ip, n
        integer :: ivx, iv2
        real(KIND=8) :: vth, cellsizevx, cellsizev2
        real(KIND=8) :: v(3), v2
        real(KIND=8) :: q,m,E(3)

        n = size(particles)
        vth = species(ispecies)%v_th*sqrt(2.)  !note: vth = sqrt(2T/m); Unterscheidet sich von sonst verwendeten v_th = sqrt(T/m) um sqrt(2)
        cellsizevx = (v_grid_max*2)/diag_bins_vx
        cellsizev2 = (v_grid_max*3)/diag_bins_v2

        DO ip=1, n
            IF (particles(ip)%data%species == ispecies) THEN
                q = particles(ip)%data%q
                m = particles(ip)%data%m
                E = particles(ip)%results%e * fc
                v = particles(ip)%data%v
                v2 = v(2)*v(2) + v(3)*v(3)

                ivx = int((v(1)/vth - (-v_grid_max)) / cellsizevx) + 1
                iv2 = int((v2/vth**2) / cellsizev2) + 1
                ivx = min( max(0, ivx), diag_bins_vx+1 )
                iv2 = min( iv2, diag_bins_v2+1 )

                nbs(ivx,iv2) = nbs(ivx,iv2) + 1

                dbs(1,ivx,iv2) = dbs(1,ivx,iv2) + q/m*E(1)                                                !d(v_x)/dt
                dbs(2,ivx,iv2) = dbs(2,ivx,iv2) + 2*q/m*v(1)*E(1) + dt*q**2/m**2*E(1)**2                  !d(v_x**2)/dt
                dbs(3,ivx,iv2) = dbs(3,ivx,iv2) + 2*q/m*dotproduct(v,E) + dt*q**2/m**2*dotproduct(E,E)    !d(v**2)/dt
            END IF
        END DO


    end subroutine fill_velocity_bins

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
                dbs(33,ix,iy,iz) = dbs(33,ix,iy,iz) + dotproduct(particles(ip)%results%E,n1)**2 * fc**2
                dbs(34,ix,iy,iz) = dbs(34,ix,iy,iz) + dotproduct(particles(ip)%results%E,t1)**2 * fc**2
                dbs(35,ix,iy,iz) = dbs(35,ix,iy,iz) + dotproduct(particles(ip)%results%E,t2)**2 * fc**2
                dbs(36,ix,iy,iz) = dbs(36,ix,iy,iz) + dotproduct(particles(ip)%results%E,n1)*dotproduct(particles(ip)%results%E,t1) * fc**2
                dbs(37,ix,iy,iz) = dbs(37,ix,iy,iz) + dotproduct(particles(ip)%results%E,t1)*dotproduct(particles(ip)%results%E,t2) * fc**2
                dbs(38,ix,iy,iz) = dbs(38,ix,iy,iz) + dotproduct(particles(ip)%results%E,t2)*dotproduct(particles(ip)%results%E,n1) * fc**2

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
