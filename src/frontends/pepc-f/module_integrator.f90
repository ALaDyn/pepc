!!!!!!!!!!!!!!!!!!!!
!! integrator module
!!!!!!!!!!!!!!!!!!!!

MODULE integrator

    USE module_pepc_types
    use helper
    IMPLICIT NONE

    CONTAINS


!=======================================================================================================
!takes array of type t_particle, calculates new velocities and positions using the boris algorithm
    SUBROUTINE boris_nonrel(p, deltat)

        IMPLICIT NONE
  
        type(t_particle), allocatable, intent(inout) :: p(:)
        real(KIND=8), intent(in) :: deltat
        real*8                :: beta,tt, sx, sy, sz, tz, ty, tx
        real*8                :: uxm,uym,uzm, uxp,uyp,uzp, uxd,uyd,uzd
        integer               :: ip


        if(root) write(*,'(a)') " == [boris_nonrel] calculate velocities "

        DO ip=1, size(p)
            IF (species(p(ip)%data%species)%physical_particle) THEN ! only move physical particals
                beta=p(ip)%data%q / p(ip)%data%m *deltat*0.5     ! charge/mass constant needed in Boris-algorithm
                !   first half-accn <-> first part of Boris-algorithm
                uxm = p(ip)%data%v(1) + beta * p(ip)%results%e(1) * fc
                uym = p(ip)%data%v(2) + beta * p(ip)%results%e(2) * fc
                uzm = p(ip)%data%v(3) + beta * p(ip)%results%e(3) * fc

                !   rotation
                tx = beta * p(ip)%data%B(1)
                ty = beta * p(ip)%data%B(2)
                tz = beta * p(ip)%data%B(3)
                tt = 1.0 + tx**2 + ty**2 + tz**2

                sx = 2.0 * tx/tt
                sy = 2.0 * ty/tt
                sz = 2.0 * tz/tt

                uxd = uxm + uym*tz - uzm*ty
                uyd = uym + uzm*tx - uxm*tz
                uzd = uzm + uxm*ty - uym*tx

                uxp = uxm + uyd*sz - uzd*sy
                uyp = uym + uzd*sx - uxd*sz
                uzp = uzm + uxd*sy - uyd*sx

                !   second half-accn  <-> second part of Boris-algorithm
                p(ip)%data%v(1) =  uxp + beta * p(ip)%results%e(1) * fc
                p(ip)%data%v(2) =  uyp + beta * p(ip)%results%e(2) * fc
                p(ip)%data%v(3) =  uzp + beta * p(ip)%results%e(3) * fc

            END IF
        END DO

    END SUBROUTINE boris_nonrel 

!=======================================================================================================
!takes array of type t_particle, calculates new velocities and positions, standard from pepc-mini
    SUBROUTINE standard_integrator(p, deltat)

        IMPLICIT NONE
  
        type(t_particle), allocatable, intent(inout) :: p(:)
        real(KIND=8), intent(in) :: deltat
        integer               :: ip

        if(root) write(*,'(a)') " == [standard_integrator] calculate velocities "


        DO ip=1, size(p)
            IF (species(p(ip)%data%species)%physical_particle) THEN ! only move physical particals
                p(ip)%data%v = p(ip)%data%v + deltat * p(ip)%data%q / p(ip)%data%m * p(ip)%results%e * fc
            END IF
        END DO

    END SUBROUTINE standard_integrator 


!======================================================================================================    
!takes array of type t_particle and calculates new positions
    SUBROUTINE push_particles(p, deltat)
        USE module_mirror_boxes
        IMPLICIT NONE
    
        type(t_particle), allocatable, intent(inout) :: p(:)
        real(KIND=8), intent(in) :: deltat
        integer :: ip

        if(root) write(*,'(a)') " == [pusher] push particles "

        DO ip=1, size(p)
            IF (species(p(ip)%data%species)%physical_particle) THEN ! only move physical particals
                p(ip)%x      = p(ip)%x      + deltat * p(ip)%data%v
            END IF
        END DO
    END SUBROUTINE push_particles

END MODULE
