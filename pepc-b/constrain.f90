subroutine constrain
    use physvars
    use treevars
    use utils

    implicit none

    real*8, dimension(1:3)		:: r_new, r_old, r_test, v
    real*8 				:: c_status		! particle in or out?
    integer 			        :: p, face_nr

    ! bisections
    real*8, dimension(1:3)		:: r_d, temp	        ! crossing point of particle
    real*8				:: p_new, p_old, diff

    ! reflection plane
    real*8, dimension(1:3)		:: n

    ! new direction of v
    real*8, dimension(1:3)                :: c

    ! file id
    integer                             :: c_file = 76, nr_out
    logical :: constrain_debug=.false.

    if (constrain_debug) write(*, *) 'in constrain'

    ! walk through all particles
    nr_out = 0
    do p = 1, npp
        r_new = (/x(p), y(p), z(p)/)
        v = (/ux(p), uy(p), uz(p)/)

        do face_nr = 1, number_faces ! walk through all faces
            call face(r_new, c_status, face_nr, & 
              target_geometry, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre )

            if (c_status .ge. 0.) then
                nr_out = nr_out + 1

                ! get a good r_old
                r_test = r_new - dt * v
                call face(r_test, c_status, face_nr, & 
                  target_geometry, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre )
                if (c_status .gt. 0.) then
                    call cutvector(r_test, face_nr, n, r_old)
                    temp = r_new - r_old
                    v = sqrt(dot_product(v, v)) * temp / sqrt(dot_product(temp, temp))
                else
                    r_old = r_test
                end if
        
                ! bisection for the particle
                p_new = 1.
                p_old = 0.
                r_test = r_old + p_new * (r_new - r_old)
                call face(r_test, c_status, face_nr, &
                  target_geometry, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre )
                do while (abs(c_status) .gt. constrain_proof)
                    diff = abs(p_new - p_old) / 2.
                    p_old = p_new
                    if (c_status .gt. 0.) then
                        p_new = p_new - diff
                    else
                        p_new = p_new + diff
                    end if
                    r_test = r_old + p_new * (r_new - r_old)
                    call face(r_test, c_status, face_nr, &
                       target_geometry, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre )
                end do
!                if (constrain_debug) write(*, *) 'Bisection for particle done.'
                call cutvector(r_test, face_nr, n, r_d)

                ! Reflect
                v = v - 2. * dot_product(v, n) * n
                if (constrain_debug) write(*, *) 'Reflecting .. new v: v     = ', v(1:3)
                c = v / sqrt(dot_product(v, v))
                r_new = r_d + sqrt(dot_product(r_d - r_new, r_d - r_new)) * c
                x(p) = r_new(1)
                y(p) = r_new(2)
                z(p) = r_new(3)
                ux(p) = v(1)
                uy(p) = v(2)
                uz(p) = v(3)
            end if
        end do
    end do
 if (constrain_debug)    write(*, *) 'Number of reflections: ', nr_out
 
end subroutine constrain
