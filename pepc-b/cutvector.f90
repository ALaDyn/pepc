subroutine cutvector(r_in, face_nr, n, r_out)
    !
    ! cuts a vector on the face face_nr
    ! and optionally puts out thr normal vector
    !
    use physvars

    implicit none
    
    real*8, dimension(1:3), intent(in)  :: r_in
    real*8, dimension(1:3), intent(out) :: n, r_out
    real*8, dimension(1:3)              :: r_affin, temp
    integer, intent(in)               :: face_nr
    real*8, dimension(1:3, 1:3)	      :: T
    real*8                              :: c_status, phi, psi, gamma

    r_affin = r_in - plasma_centre
    select case (target_geometry) 
    case (0) ! slab
        select case (face_nr)
        case (1, 2) ! x-y-planes
            n = (/0., 0., 1./)
            if (face_nr == 2) n = -n
            r_affin = .5 * z_plasma * r_affin / abs(dot_product(r_affin, n))
        case (3, 4) ! x-z-planes
            n = (/0., 1., 0./)
            if (face_nr == 4) n = -n
            r_affin = .5 * y_plasma * r_affin / abs(dot_product(r_affin, n))
        case (5, 6) ! y-z-planes
            n = (/1., 0., 0./)
            if (face_nr == 6) n = -n
            r_affin = .5 * x_plasma * r_affin / abs(dot_product(r_affin, n))
        end select
        r_out = r_affin + plasma_centre
        
    case (1) ! sphere
        n = - r_affin / sqrt(dot_product(r_affin, r_affin))
        r_out = r_sphere * r_affin / abs(dot_product(r_affin, n)) + plasma_centre
        
    case (2) ! disc
        select case (face_nr)
        case (1) ! zylinder
            n = (/0.d0, -r_affin(2), -r_affin(3)/)
            n = n / sqrt(dot_product(n, n))
            r_affin = r_sphere * r_affin / abs(dot_product(r_affin, n))
        case (2, 3) ! y-z-planes
            n = (/1., 0., 0./)
            if (face_nr == 3) n = -n
            r_affin = .5 * x_plasma * r_affin / abs(dot_product(r_affin, n))
        end select
        r_out = r_affin + plasma_centre
        
    case (3) ! wire
        select case (face_nr)
        case (1) ! zylinder
            n = (/-r_affin(1), -r_affin(2), 0.d0/)
            n = n / sqrt(dot_product(n, n))
            r_affin = r_sphere * r_affin / abs(dot_product(r_affin, n))
        case (2, 3) ! x-y-planes
            n = (/0., 0., 1./)
            if (face_nr == 3) n = -n
            r_affin = .5 * z_plasma * r_affin / abs(dot_product(r_affin, n))
        end select
        r_out = r_affin + plasma_centre
        
    case (4) ! ellipsoid
        T = reshape((/1., 0., 0., 0., x_plasma / y_plasma, 0., 0., 0., x_plasma / z_plasma/), (/3, 3/))
        temp = matmul(T, r_affin)
        r_out = r_sphere * x_plasma * r_affin / sqrt(dot_product(temp, temp))
        temp = r_sphere * x_plasma * temp / sqrt(dot_product(temp, temp))
        phi = atan(temp(1) / sqrt(dot_product(temp, temp) - temp(1)**2))
        psi = atan(temp(1) * sqrt(dot_product(r_out, r_out) - temp(1)**2) / (dot_product(temp, temp) - temp(1)**2))
        phi = phi - psi
        call rotate(phi, temp, n)
        n = -n / sqrt(dot_product(n, n))
        r_out = r_out + plasma_centre
        
    case (5) ! prism
        gamma = atan(y_plasma / (2 * x_plasma))
        select case (face_nr)
        case (1, 2) ! x-y-planes
            n = (/0., 0.,-1./)
            if (face_nr == 2) n = -n
            r_affin = .5 * z_plasma * r_affin / abs(dot_product(r_affin, n))
        case (3) ! y-z-plane
            n = (/1., 0., 0./)
            r_affin = .5 * x_plasma * r_affin / abs(dot_product(r_affin, n))
        case (4, 5) ! the other planes
            n = (/-sin(gamma), +cos(gamma),0.d0/)
            r_affin(2) = tan(gamma) * r_affin(1) - y_plasma / 4.
            if (face_nr == 5) then
                n(2) = -n(2)
                r_affin(2) = -r_affin(2)
            end if
        end select
        r_out = r_affin + plasma_centre
        
    case (6) ! semisphere
        temp = r_affin - (/-r_sphere / 2.d0, 0.d0, 0.d0/)
        select case (face_nr)
        case (1) ! semisphere
            n = -temp
            n = n / sqrt(dot_product(n, n))
            r_affin = r_sphere * temp / abs(dot_product(temp, n)) + (/-r_sphere / 2.d0, 0.d0, 0.d0/)
        case (2) ! plane
            n = (/-1., 0., 0./)
            r_affin = .5 * r_sphere * r_affin / abs(dot_product(r_affin, n))
        end select
        r_out = r_affin + plasma_centre
        
    case (7) ! hollow ball
        select case (face_nr)
        case (1) ! outer sphere
            n = -r_affin
            n = n / sqrt(dot_product(n, n))
            r_affin = r_sphere * r_affin / abs(dot_product(r_affin, n))
        case (2) ! inner sphere
            n = r_affin
            n = n / sqrt(dot_product(n, n))
            r_affin = (r_sphere - x_plasma) * r_affin / abs(dot_product(r_affin, n))
        end select
        r_out = r_affin + plasma_centre
        
    case (8) ! hollow semisphere
        temp = r_affin - (/-r_sphere / 2.d0, 0.d0, 0.d0/)
        select case (face_nr)
        case (1) ! outer semisphere
            n = -temp
            n = n / sqrt(dot_product(n, n))
            r_affin = r_sphere * temp / abs(dot_product(temp, n)) + (/-r_sphere / 2., 0., 0./)
        case (2) ! inner semisphere
            n = temp
            n = n / sqrt(dot_product(n, n))
            r_affin = (r_sphere - x_plasma) * temp / abs(dot_product(temp, n)) + (/-r_sphere / 2., 0., 0./)
        case (3) ! y-z-plane
            n = (/-1., 0., 0./)
            r_affin = .5 * r_sphere * r_affin / abs(dot_product(r_affin, n))
        end select
        r_out = r_affin + plasma_centre
        
    end select
    
contains

    !
    ! Rotates a vector around the y-axis
    ! with an angle phi (needed for the 
    ! ellipsoid)
    !
    subroutine rotate(phi, v_in, v_out)
        real*8				:: phi
        real*8, dimension(1:3)		:: v_in, v_out
        real*8				:: c, s
        real*8, dimension(1:3, 1:3)	:: R

        c = cos(phi)
        s = sin(phi)

        R = reshape((/1.d0, 0.d0, 0.d0, 0.d0, c, s, 0.d0, -s, c/), (/3, 3/))
        v_out = matmul(R, v_in)
    end subroutine rotate

end subroutine cutvector
