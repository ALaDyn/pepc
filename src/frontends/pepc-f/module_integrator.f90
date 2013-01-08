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
    SUBROUTINE boris_nonrel(p)

        IMPLICIT NONE
  
        type(t_particle), allocatable, intent(inout) :: p(:)
        real*8                :: fact
        real*8                :: beta,tt, sx, sy, sz, tz, ty, tx
        real*8                :: uxm,uym,uzm, uxp,uyp,uzp, uxd,uyd,uzd
        integer               :: ip


        if(root) write(*,'(a)') " == [boris_nonrel] calculate velocities "

        IF (step==0) THEN
            fact=dt/2.0
        ELSE
            fact=dt
        END IF

        DO ip=1, np
            IF (p(ip)%label>0) THEN                            ! wall particles have negative labels
                beta=p(ip)%data%q / p(ip)%data%m *fact*0.5     ! charge/mass constant needed in Boris-algorithm
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

        call push_particles(p)

    END SUBROUTINE boris_nonrel 

!=======================================================================================================
!takes array of type t_particle, calculates new velocities and positions, standard from pepc-mini
    SUBROUTINE standard_integrator(p)

        IMPLICIT NONE
  
        type(t_particle), allocatable, intent(inout) :: p(:)
        real*8                :: fact
        integer               :: ip

        fact=dt

        DO ip=1, np
            IF (p(ip)%label>0) THEN
                p(ip)%data%v = p(ip)%data%v + fact * p(ip)%data%q / p(ip)%data%m * p(ip)%results%e * fc
            END IF
        END DO

        call push_particles(p)

    END SUBROUTINE standard_integrator 


!======================================================================================================    
!takes array of type t_particle and calculates new positions
    SUBROUTINE push_particles(p)
        USE module_mirror_boxes
        IMPLICIT NONE
    
        type(t_particle), allocatable, intent(inout) :: p(:)
        integer :: ip
        real*8  :: fact

        fact = dt
        if(root) write(*,'(a)') " == [pusher] push particles "

        DO ip=1, np
            IF (p(ip)%label>0) THEN
                p(ip)%x      = p(ip)%x      + fact * p(ip)%data%v    
            END IF
        END DO
    END SUBROUTINE push_particles

END MODULE
