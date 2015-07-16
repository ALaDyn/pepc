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
                            avg_8, avg_9, avg_10, avg_11, avg_12, avg_13, avg_14, avg_15, avg_16, &
                            avg_17, avg_18, avg_19, avg_20, avg_21, avg_22, avg_23, avg_24, avg_25, &
                            avg_26, avg_27, avg_28, avg_29, avg_30, avg_31, avg_32, avg_33, avg_34, &
                            avg_35, avg_36, avg_37, avg_38, avg_39, avg_40, avg_41, avg_42, &
                            avg_fields)
        implicit none
        include 'mpif.h'

        type(t_particle),  intent(in) :: p(:)
        real(KIND=8), intent(inout) :: avg_1(:), avg_2(:), avg_3(:), avg_4(:), avg_5(:), avg_6(:), avg_7(:)
        real(KIND=8), intent(inout) :: avg_8(:), avg_9(:), avg_10(:), avg_11(:), avg_12(:), avg_13(:), avg_14(:)
        real(KIND=8), intent(inout) :: avg_15(:), avg_16(:), avg_17(:), avg_18(:), avg_19(:), avg_20(:),avg_21(:)
        real(KIND=8), intent(inout) :: avg_22(:), avg_23(:), avg_24(:), avg_25(:), avg_26(:), avg_27(:),avg_28(:)
        real(KIND=8), intent(inout) :: avg_29(:), avg_30(:), avg_31(:), avg_32(:), avg_33(:), avg_34(:),avg_35(:)
        real(KIND=8), intent(inout) :: avg_36(:), avg_37(:), avg_38(:), avg_39(:), avg_40(:), avg_41(:),avg_42(:)
        real(KIND=8), intent(inout) :: avg_fields(:,:)
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
        !avg_17 = <|v(t)|^3>
        !avg_18 = <|v(t)|^4>
        !avg_19 = <|v(t)|^5>
        !avg_20 = <|v(t)|^6>
        !avg_21 = <vx(t)>
        !avg_22 = <vx(t)^2>
        !avg_23 = <vx(t)^3>
        !avg_24 = <vx(t)^4>
        !avg_25 = <vx(t)^5>
        !avg_26 = <vx(t)^6>
        !avg_27 = <|vx(t)|>
        !avg_28 = <|vx(t)|^3>
        !avg_29 = <|vx(t)|^5>
        !avg_30 = <vpar(t)^3>
        !avg_31 = <vpar(t)^4>
        !avg_32 = <vpar(t)^5>
        !avg_33 = <vpar(t)^6>
        !avg_34 = <|vpar(t)|^3>
        !avg_35 = <|vpar(t)|^5>
        !avg_36 = <v(t)^2 v(0)^2>
        !avg_37 = <vx(t) vx(0)>
        !avg_38 = <|vx(t)| |vx(0)|>
        !avg_39 = <vx(t)^2 vx(0)^2>
        !avg_40 = <vpar(t) vpar(0)>
        !avg_41 = <|vpar(t)| |vpar(0)|>
        !avg_42 = <vpar(t)^2 vpar(0)^2>

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
        real(KIND=8) :: avg_absv(nspecies-1),avg_absv3(nspecies-1),avg_absv4(nspecies-1)
        real(KIND=8) :: avg_absv5(nspecies-1),avg_absv6(nspecies-1)
        real(KIND=8) :: avg_vt2v02(nspecies-1)
        real(KIND=8) :: avg_fields_l(12, nspecies-1)
        real(KIND=8) :: avg_vx1(nspecies-1), avg_vx2(nspecies-1), avg_vx3(nspecies-1), &
                        avg_vx4(nspecies-1), avg_vx5(nspecies-1), avg_vx6(nspecies-1)
        real(KIND=8) :: avg_absvx1(nspecies-1), avg_absvx3(nspecies-1), avg_absvx5(nspecies-1)
        real(KIND=8) :: avg_vpar3(nspecies-1), &
                        avg_vpar4(nspecies-1), avg_vpar5(nspecies-1), avg_vpar6(nspecies-1)
        real(KIND=8) :: avg_absvpar3(nspecies-1), avg_absvpar5(nspecies-1)
        real(KIND=8) :: avg_vxt1vx01(nspecies-1), avg_absvxt1absvx01(nspecies-1), avg_vxt2vx02(nspecies-1)
        real(KIND=8) :: avg_vpart1vpar01(nspecies-1), avg_absvpart1absvpar01(nspecies-1), avg_vpart2vpar02(nspecies-1)



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
        avg_absv3 = 0.
        avg_absv4 = 0.
        avg_absv5 = 0.
        avg_absv6 = 0.
        avg_vx1 = 0.
        avg_vx2 = 0.
        avg_vx3 = 0.
        avg_vx4 = 0.
        avg_vx5 = 0.
        avg_vx6 = 0.
        avg_absvx1 = 0.
        avg_absvx3 = 0.
        avg_absvx5 = 0.
        avg_vpar3 = 0.
        avg_vpar4 = 0.
        avg_vpar5 = 0.
        avg_vpar6 = 0.
        avg_absvpar3 = 0.
        avg_absvpar5 = 0.
        avg_vxt1vx01 = 0.
        avg_absvxt1absvx01 = 0.
        avg_vxt2vx02 = 0.
        avg_vpart1vpar01 = 0.
        avg_absvpart1absvpar01 = 0.
        avg_vpart2vpar02 = 0.

        avg_vt2v02 = 0.

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
        avg_17 = 0.
        avg_18 = 0.
        avg_19 = 0.
        avg_20 = 0.
        avg_21 = 0.
        avg_22 = 0.
        avg_23 = 0.
        avg_24 = 0.
        avg_25 = 0.
        avg_26 = 0.
        avg_27 = 0.
        avg_28 = 0.
        avg_29 = 0.
        avg_30 = 0.
        avg_31 = 0.
        avg_32 = 0.
        avg_33 = 0.
        avg_34 = 0.
        avg_35 = 0.
        avg_36 = 0.
        avg_37 = 0.
        avg_38 = 0.
        avg_39 = 0.
        avg_40 = 0.
        avg_41 = 0.
        avg_42 = 0.
        avg_fields = 0.

        do ip=1, sum(npps)
            ispecies = p(ip)%data%species
            IF (.NOT.(species(p(ip)%data%species)%moving_particle)) cycle

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
            avg_absv3(ispecies) = avg_absv3(ispecies) + norm(v_p)**3
            avg_absv4(ispecies) = avg_absv4(ispecies) + norm(v_p)**4
            avg_absv5(ispecies) = avg_absv5(ispecies) + norm(v_p)**5
            avg_absv6(ispecies) = avg_absv6(ispecies) + norm(v_p)**6

            avg_vx1(ispecies) = avg_vx1(ispecies) + v_p(1)
            avg_vx2(ispecies) = avg_vx2(ispecies) + v_p(1)**2
            avg_vx3(ispecies) = avg_vx3(ispecies) + v_p(1)**3
            avg_vx4(ispecies) = avg_vx4(ispecies) + v_p(1)**4
            avg_vx5(ispecies) = avg_vx5(ispecies) + v_p(1)**5
            avg_vx6(ispecies) = avg_vx6(ispecies) + v_p(1)**6
            avg_absvx1(ispecies) = avg_absvx1(ispecies) + abs(v_p(1))
            avg_absvx3(ispecies) = avg_absvx3(ispecies) + abs(v_p(1))**3
            avg_absvx5(ispecies) = avg_absvx5(ispecies) + abs(v_p(1))**5

            avg_vpar3(ispecies) = avg_vpar3(ispecies) + vpar_p**3
            avg_vpar4(ispecies) = avg_vpar4(ispecies) + vpar_p**4
            avg_vpar5(ispecies) = avg_vpar5(ispecies) + vpar_p**5
            avg_vpar6(ispecies) = avg_vpar6(ispecies) + vpar_p**6
            avg_absvpar3(ispecies) = avg_absvpar3(ispecies) + abs(vpar_p)**3
            avg_absvpar5(ispecies) = avg_absvpar5(ispecies) + abs(vpar_p)**5

            avg_vxt1vx01(ispecies) = avg_vxt1vx01(ispecies) + v_p(1)*v0_p(1)
            avg_absvxt1absvx01(ispecies) = avg_absvxt1absvx01(ispecies) + abs(v_p(1))*abs(v0_p(1))
            avg_vxt2vx02(ispecies) = avg_vxt2vx02(ispecies) + v_p(1)**2*v0_p(1)**2
            avg_vpart1vpar01(ispecies) = avg_vpart1vpar01(ispecies) + vpar_p*norm(v0_p)
            avg_absvpart1absvpar01(ispecies) = avg_absvpart1absvpar01(ispecies) + abs(vpar_p)*norm(v0_p)
            avg_vpart2vpar02(ispecies) = avg_vpart2vpar02(ispecies) + vpar_p**2*norm(v0_p)**2
            avg_vt2v02(ispecies) = avg_vt2v02(ispecies) + dotproduct(v_p,v_p)*dotproduct(v0_p,v0_p)

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
        call MPI_REDUCE(avg_absv3, avg_17, nspecies-1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(avg_absv4, avg_18, nspecies-1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(avg_absv5, avg_19, nspecies-1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(avg_absv6, avg_20, nspecies-1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(avg_vx1, avg_21, nspecies-1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(avg_vx2, avg_22, nspecies-1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(avg_vx3, avg_23, nspecies-1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(avg_vx4, avg_24, nspecies-1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(avg_vx5, avg_25, nspecies-1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(avg_vx6, avg_26, nspecies-1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(avg_absvx1, avg_27, nspecies-1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(avg_absvx3, avg_28, nspecies-1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(avg_absvx5, avg_29, nspecies-1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(avg_vpar3, avg_30, nspecies-1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(avg_vpar4, avg_31, nspecies-1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(avg_vpar5, avg_32, nspecies-1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(avg_vpar6, avg_33, nspecies-1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(avg_absvpar3, avg_34, nspecies-1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(avg_absvpar5, avg_35, nspecies-1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(avg_vt2v02, avg_36, nspecies-1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(avg_vxt1vx01, avg_37, nspecies-1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(avg_absvxt1absvx01, avg_38, nspecies-1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(avg_vxt2vx02, avg_39, nspecies-1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(avg_vpart1vpar01, avg_40, nspecies-1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(avg_absvpart1absvpar01, avg_41, nspecies-1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(avg_vpart2vpar02, avg_42, nspecies-1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

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
        avg_17 = avg_17(:) / tnpps(1:)
        avg_18 = avg_18(:) / tnpps(1:)
        avg_19 = avg_19(:) / tnpps(1:)
        avg_20 = avg_20(:) / tnpps(1:)
        avg_21 = avg_21(:) / tnpps(1:)
        avg_22 = avg_22(:) / tnpps(1:)
        avg_23 = avg_23(:) / tnpps(1:)
        avg_24 = avg_24(:) / tnpps(1:)
        avg_25 = avg_25(:) / tnpps(1:)
        avg_26 = avg_26(:) / tnpps(1:)
        avg_27 = avg_27(:) / tnpps(1:)
        avg_28 = avg_28(:) / tnpps(1:)
        avg_29 = avg_29(:) / tnpps(1:)
        avg_30 = avg_30(:) / tnpps(1:)
        avg_31 = avg_31(:) / tnpps(1:)
        avg_32 = avg_32(:) / tnpps(1:)
        avg_33 = avg_33(:) / tnpps(1:)
        avg_34 = avg_34(:) / tnpps(1:)
        avg_35 = avg_35(:) / tnpps(1:)
        avg_36 = avg_36(:) / tnpps(1:)
        avg_37 = avg_37(:) / tnpps(1:)
        avg_38 = avg_38(:) / tnpps(1:)
        avg_39 = avg_39(:) / tnpps(1:)
        avg_40 = avg_40(:) / tnpps(1:)
        avg_41 = avg_41(:) / tnpps(1:)
        avg_42 = avg_42(:) / tnpps(1:)
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
        integer, intent(inout)      :: nbs(:,0:)
        real(KIND=8), intent(inout) :: dbs(:,0:)

        integer :: ip, n
        integer :: ivpar
        real(KIND=8) :: cellsizevpar
        real(KIND=8) :: v(3), v0_p(3), e0_p(3), h_p, vpar_p, vperp_p, beta_p, vpar0_p, vth
        real(KIND=8) :: q,m
        n = size(particles)
        cellsizevpar = (v_grid_max*2)/diag_bins_vpar
        vth = species(ispecies)%v_th

        DO ip=1, n
            IF (particles(ip)%data%species == ispecies) THEN
                q = species(ispecies)%q
                m = species(ispecies)%m
                v = particles(ip)%data%v

                v0_p = initial_velocities(:, particles(ip)%label)
                e0_p = v0_p / norm(v0_p)
                h_p = 0.5 * species(ispecies)%m/e * (sum(v**2) - sum(v0_p**2))
                vpar_p = dotproduct(v, e0_p)
                vpar0_p = dotproduct(v0_p, e0_p)
                vperp_p = sqrt(sum(v**2) - vpar_p**2)
                beta_p = acos( dotproduct(v0_p, v) / (norm(v0_p) * norm(v)) ) !range of acos func: 0..pi

                !binning for the pdf is based on current parallel velocity vpar_p
                ivpar = int((vpar_p/vth - (-v_grid_max)) / cellsizevpar) + 1
                ivpar = min( max(0, ivpar), diag_bins_vpar+1 )
                nbs(1,ivpar) = nbs(1,ivpar) + 1

                !binning for velocity resolved Hockney diag is based on the initial parallel velocity vpar0_p
                ivpar = int((vpar0_p/vth - (-v_grid_max)) / cellsizevpar) + 1
                ivpar = min( max(0, ivpar), diag_bins_vpar+1 )
                nbs(2,ivpar) = nbs(2,ivpar) + 1
                dbs(1,ivpar) = dbs(1,ivpar) + vpar_p
                dbs(2,ivpar) = dbs(2,ivpar) + vpar_p**2
                dbs(3,ivpar) = dbs(3,ivpar) + vperp_p**2
                dbs(4,ivpar) = dbs(4,ivpar) + beta_p**2
                dbs(5,ivpar) = dbs(5,ivpar) + h_p
                dbs(6,ivpar) = dbs(6,ivpar) + h_p**2

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
