!
!   FACE
!
!  Function to set up bounded geometry targets in conjunction with 
!  randion. 
!  M. Hammes, U. Wuppertal
!
!  $Revision$
!
! Generalised description of convex faces described 
! by a function f(x, y, z) = 0. 
! Saves f(r) in c_status
! and means c_status > 0.: particle is out,
! c_status < 0.: particle is in.
!
subroutine face(r, c_status, face_nr)
    use treevars
    use utils

    implicit none 

    real, dimension(1:3)	    :: r, r_diff
    real			    :: c_status
    integer                         :: face_nr

    ! vectors and other values for the slab
    real, dimension(1:3)            :: normal_vector, offset_vector
    real                            :: gamma

    ! matrix for the ellipsoid
    real, dimension(1:3, 1:3)       :: A

    select case (initial_config)
    case (0) ! slab
        select case (face_nr)
        case(1) ! x-y-plane, negative z
            normal_vector = (/0., 0., +1./)
            offset_vector = (/0., 0., -z_plasma / 2./) + plasma_centre
        case(2) ! x-y-plane, positive z
            normal_vector = (/0., 0., -1./)
            offset_vector = (/0., 0., +z_plasma / 2./) + plasma_centre
        case(3) ! x-z-plane, negative y
            normal_vector = (/0., +1., 0./)
            offset_vector = (/0., -y_plasma / 2., 0./) + plasma_centre
        case(4) ! x-z-plane, positive y
            normal_vector = (/0., -1., 0./)
            offset_vector = (/0., +y_plasma / 2., 0./) + plasma_centre
        case(5) ! y-z-plane, negative x
            normal_vector = (/+1., 0., 0./)
            offset_vector = (/-x_plasma / 2., 0., 0./) + plasma_centre
        case(6) ! y-z-plane, positive x
            normal_vector = (/-1., 0., 0./)
            offset_vector = (/+x_plasma / 2., 0., 0./) + plasma_centre
        end select

        c_status = dot_product(normal_vector, offset_vector - r)
    case(1) ! sphere
        r_diff = r - plasma_centre
        c_status = dot_product(r_diff, r_diff) - r_sphere**2
    case(2) ! disc
        select case (face_nr)
        case(1) ! the tube in x-direction
            r_diff = r - plasma_centre
            c_status = r_diff(2)**2 + r_diff(3)**2 - r_sphere**2
        case(2) ! y-z-plane, negative x
            normal_vector = (/+1., 0., 0./)
            offset_vector = (/-x_plasma / 2., 0., 0./) + plasma_centre
            c_status = dot_product(normal_vector, offset_vector - r)
        case(3) ! y-z-plane, positive x
            normal_vector = (/-1., 0., 0./)
            offset_vector = (/+x_plasma / 2., 0., 0./) + plasma_centre
            c_status = dot_product(normal_vector, offset_vector - r)
        end select
    case(3)  ! wire
        select case (face_nr)
        case(1) ! the tube in z-direction
            r_diff = r - plasma_centre
            c_status = r_diff(1)**2 + r_diff(2)**2 - r_sphere**2
        case(2) ! x-y-plane, negative z
            normal_vector = (/0., 0., +1./)
            offset_vector = (/0., 0., -z_plasma / 2./) + plasma_centre
            c_status = dot_product(normal_vector, offset_vector - r)
        case(3) ! x-y-plane, positive z
            normal_vector = (/0., 0., -1./)
            offset_vector = (/0., 0., +z_plasma / 2./) + plasma_centre
            c_status = dot_product(normal_vector, offset_vector - r)
        end select
    case (4) ! ellipsoid
        A = reshape((/1 / x_plasma, 0., 0., 0., 1 / y_plasma, 0., 0., 0., 1 / z_plasma/), (/3, 3/))
        r_diff = matmul(A, r - plasma_centre)
        c_status = dot_product(r_diff, r_diff) - r_sphere**2
    case (5) ! prism in x-y plane
        gamma = atan(y_plasma / (2 * x_plasma))
        select case (face_nr)
        case(1) ! x-y-plane, positive z
            normal_vector = (/0., 0., -1./)
            offset_vector = (/0., 0., +z_plasma / 2./) + plasma_centre
        case(2) ! x-y-plane, negative z
            normal_vector = (/0., 0., +1./)
            offset_vector = (/0., 0., -z_plasma / 2./) + plasma_centre
        case(3) ! y-z-plane, negative x
            normal_vector = (/+1., 0., 0./)
            offset_vector = (/-x_plasma / 2., 0., 0./) + plasma_centre
        case(4) ! negative y
            normal_vector = (/-sin(gamma), +cos(gamma), 0./)
            offset_vector = (/-x_plasma / 2., -y_plasma / 2., 0./) + plasma_centre
        case(5) ! positive y
            normal_vector = (/-sin(gamma), -cos(gamma), 0./)
            offset_vector = (/-x_plasma / 2., +y_plasma / 2.,0./) + plasma_centre
        end select
        c_status = dot_product(normal_vector, offset_vector - r)
    case (6) ! semisphere
        select case (face_nr)
        case(1) ! semisphere
            r_diff = r + (/r_sphere / 2., 0., 0./) - plasma_centre
            c_status = dot_product(r_diff, r_diff) - r_sphere**2
        case(2) ! y-z-plane
            normal_vector = (/+1., 0., 0./)
            offset_vector = (/-r_sphere / 2., 0., 0./) + plasma_centre
            c_status = dot_product(normal_vector, offset_vector - r)
        end select
    case (7) ! hollow sphere
        r_diff = r - plasma_centre
        select case (face_nr)
        case(1) ! outer_sphere
            c_status = dot_product(r_diff, r_diff) - r_sphere**2
        case(2) ! inner_sphere
            c_status = (r_sphere - x_plasma)**2 - dot_product(r_diff, r_diff)
        end select
    case (8) ! hollow semisphere
        select case (face_nr)
        case(1) ! outer semisphere
            r_diff = r + (/r_sphere / 2., 0., 0./) - plasma_centre
            c_status = dot_product(r_diff, r_diff) - r_sphere**2
        case(2) ! inner semisphere
            r_diff = r + (/r_sphere / 2., 0., 0./) - plasma_centre
            c_status = (r_sphere - x_plasma)**2 - dot_product(r_diff, r_diff)
        case(3) ! y-z-plane
            normal_vector = (/+1., 0., 0./)
            offset_vector = (/-r_sphere / 2., 0., 0./) + plasma_centre
            c_status = dot_product(normal_vector, offset_vector - r)
        end select
    end select

end subroutine face
