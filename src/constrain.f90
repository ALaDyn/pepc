subroutine constrain
    use treevars
    use utils

    implicit none

    real, dimension(1:3)		:: r_new, r_old, r_test, v
    real 				:: c_status		! particle in or out?
    integer 			        :: p, face_nr

    ! bisections
    real, dimension(1:3)		:: r_d, temp	        ! crossing point of particle
    real				:: p_new, p_old, diff

    ! reflection plane
    real, dimension(1:3)		:: n

    ! new direction of v
    real, dimension(1:3)                :: c

    ! file id
    integer                             :: c_file = 76, nr_out

    write(c_file, *) 'in constrain'

    ! walk through all particles
    nr_out = 0
    do p = 1, npp
        r_new = (/x(p), y(p), z(p)/)
        v = (/ux(p), uy(p), uz(p)/)

        do face_nr = 1, number_faces ! walk through all faces
            call face(r_new, c_status, face_nr)
            if (c_status .ge. 0.) then
                nr_out = nr_out + 1

                ! get a good r_old
                r_test = r_new - dt * v
                call face(r_test, c_status, face_nr) 
                if (c_status .gt. 0.) then
                    call cutvector(r_test, face_nr, n, r_old)
                    temp = r_new - r_old
                    v = sqrt(dot_product(v, v)) * temp / sqrt(dot_product(temp, temp))
                else
                    r_old = r_test
                end if
        
                ! besection for the particle
                p_new = 1.
                p_old = 0.
                r_test = r_old + p_new * (r_new - r_old)
                call face(r_test, c_status, face_nr)
                do while (abs(c_status) .gt. constrain_proof)
                    diff = abs(p_new - p_old) / 2.
                    p_old = p_new
                    if (c_status .gt. 0.) then
                        p_new = p_new - diff
                    else
                        p_new = p_new + diff
                    end if
                    r_test = r_old + p_new * (r_new - r_old)
                    call face(r_test, c_status, face_nr)
                end do
                write(c_file, *) 'Bisection for particle done.'
                call cutvector(r_test, face_nr, n, r_d)

                ! Reflect
                write(c_file, *) 'Reflecting particle.'
                v = v - 2. * dot_product(v, n) * n
                write(c_file, *) 'New v: v     = ', v(1:3)
                c = v / sqrt(dot_product(v, v))
                r_new = r_d + sqrt(dot_product(r_d - r_new, r_d - r_new)) * c
                x(p) = r_new(1)
                y(p) = r_new(2)
                z(p) = r_new(3)
                ux(p) = v(1)
                uy(p) = v(2)
                uz(p) = v(3)
                write(c_file, *) 'Reflection done.'
                write(c_file, *) 'Bisection for particle nr ', p, ' done.'
                write(c_file, *) ''
            end if
        end do
    end do
    write(c_file, *) 'Number of reflections: ', nr_out
 
end subroutine constrain
