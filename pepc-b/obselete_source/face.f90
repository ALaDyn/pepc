!
!   FACE
!
!  Function to set up bounded geometry targets in conjunction with 
!  PLASMA_START. 
!  M. Hammes, U. Wuppertal
!
!  $Revision: 608 $
!
! Generalised description of convex faces described 
! by a function f(x, y, z) = 0. 
! Saves f(r) in c_status
! and means c_status > 0.: particle is out,
! c_status < 0.: particle is in.

subroutine face(r, c_status, face_nr, &
    target_geometry, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre )

    implicit none 

    real*8, intent(in), dimension(1:3)	    :: r
    real*8, intent(out)			    :: c_status
    integer, intent(in)                     :: face_nr
    integer, intent(in)	                    :: target_geometry
    real, intent(in) :: x_plasma, y_plasma, z_plasma, r_sphere
    real, intent(in) :: plasma_centre(3)

    ! vectors and other values for the slab
    real*8, dimension(1:3)	    :: r_diff
    real*8, dimension(1:3)            :: normal_vector, offset_vector
    real*8                            :: gamma

    ! matrix for the ellipsoid
    real, dimension(1:3, 1:3)       :: A

    select case (target_geometry)
    case (0) ! slab
        select case (face_nr)
        case(1) ! x-y-plane, negative z
            normal_vector = (/ 0., 0., +1. /)
            offset_vector = (/ 0., 0., -z_plasma/2. /) + plasma_centre
        case(2) ! x-y-plane, positive z
            normal_vector = (/ 0., 0., -1./)
            offset_vector = (/ 0., 0., +z_plasma/2. /) + plasma_centre
        case(3) ! x-z-plane, negative y
            normal_vector = (/ 0., +1., 0. /)
            offset_vector = (/ 0., -y_plasma/2., 0. /) + plasma_centre
        case(4) ! x-z-plane, positive y
            normal_vector = (/ 0., -1., 0. /)
            offset_vector = (/ 0., +y_plasma/2., 0. /) + plasma_centre
        case(5) ! y-z-plane, negative x
            normal_vector = (/ +1., 0., 0. /)
            offset_vector = (/ -x_plasma/2., 0., 0. /) + plasma_centre
        case(6) ! y-z-plane, positive x
            normal_vector = (/ -1., 0., 0. /)
            offset_vector = (/ +x_plasma/2., 0., 0. /) + plasma_centre
        end select

        c_status = dot_product(normal_vector, offset_vector - r)

    case(1) ! sphere
        r_diff = r - plasma_centre
        c_status = dot_product(r_diff, r_diff) - r_sphere**2

    case(11) ! hollow sphere
        r_diff = r - plasma_centre
        select case (face_nr)
        case(1) ! outer_sphere
            c_status = dot_product(r_diff, r_diff) - r_sphere**2
        case(2) ! inner_sphere
            c_status = (r_sphere - x_plasma)**2 - dot_product(r_diff, r_diff)
        end select

    case(2) ! disc
        select case (face_nr)
        case(1) ! the tube in x-direction
            r_diff = r - plasma_centre
            c_status = r_diff(2)**2 + r_diff(3)**2 - r_sphere**2
        case(2) ! y-z-plane, negative x
            normal_vector = (/ +1., 0., 0. /)
            offset_vector = (/ -x_plasma/2., 0., 0. /) + plasma_centre
            c_status = dot_product(normal_vector, offset_vector - r)
        case(3) ! y-z-plane, positive x
            normal_vector = (/ -1., 0., 0. /)
            offset_vector = (/ +x_plasma/2., 0., 0. /) + plasma_centre
            c_status = dot_product(normal_vector, offset_vector - r)
        end select

    case(12)  ! hollow disc/tube, thickness z_plasma
        select case (face_nr)
        case(1) ! the tube in x-direction
            r_diff = r - plasma_centre
            c_status = r_diff(2)**2 + r_diff(3)**2 - r_sphere**2
        case(2) ! the inner tube in z-direction
            r_diff = r - plasma_centre
            c_status =  (r_sphere-z_plasma)**2 - (r_diff(2)**2 + r_diff(3)**2)
        case(3) ! y-z-plane, negative x
            normal_vector = (/ +1., 0., 0. /)
            offset_vector = (/ -x_plasma/2., 0., 0. /) + plasma_centre
            c_status = dot_product(normal_vector, offset_vector - r)
        case(4) ! y-z-plane, positive x
            normal_vector = (/ -1., 0., 0. /)
            offset_vector = (/ +x_plasma/2., 0., 0. /) + plasma_centre
            c_status = dot_product(normal_vector, offset_vector - r)
        end select

    case(3)  ! wire
        select case (face_nr)
        case(1) ! the tube in z-direction
            r_diff = r - plasma_centre
            c_status = r_diff(1)**2 + r_diff(2)**2 - r_sphere**2
        case(2) ! x-y-plane, negative z
            normal_vector = (/ 0., 0., +1. /)
            offset_vector = (/ 0., 0., -z_plasma/2. /) + plasma_centre
            c_status = dot_product(normal_vector, offset_vector - r)
        case(3) ! x-y-plane, positive z
            normal_vector = (/ 0., 0., -1. /)
            offset_vector = (/ 0., 0., +z_plasma/2. /) + plasma_centre
            c_status = dot_product(normal_vector, offset_vector - r)
        end select

    case(13)  ! hollow wire
        select case (face_nr)
        case(1) ! the tube in z-direction
            r_diff = r - plasma_centre
            c_status = r_diff(1)**2 + r_diff(2)**2 - r_sphere**2
        case(2) ! the inner tube in z-direction
            r_diff = r - plasma_centre
            c_status =  (r_sphere-x_plasma)**2 - (r_diff(1)**2 + r_diff(2)**2)
        case(3) ! x-y-plane, negative z
            normal_vector = (/ 0., 0., +1. /)
            offset_vector = (/ 0., 0., -z_plasma/2. /) + plasma_centre
            c_status = dot_product(normal_vector, offset_vector - r)
        case(4) ! x-y-plane, positive z
            normal_vector = (/ 0., 0., -1. /)
            offset_vector = (/ 0., 0., +z_plasma/2. /) + plasma_centre
            c_status = dot_product(normal_vector, offset_vector - r)
        end select

    case(4) ! ellipsoid
        A = reshape((/ 1./x_plasma, 0., 0., 0., 1./y_plasma, 0., 0., 0., 1./z_plasma /), (/3, 3/))
        r_diff = matmul(A, r - plasma_centre)
        c_status = dot_product(r_diff, r_diff) - r_sphere**2

    case(5) ! prism in x-y plane
        gamma = atan(y_plasma / (2 * x_plasma))
        select case (face_nr)
        case(1) ! x-y-plane, positive z
            normal_vector = (/ 0., 0., -1. /)
            offset_vector = (/ 0., 0., +z_plasma/2. /) + plasma_centre
        case(2) ! x-y-plane, negative z
            normal_vector = (/ 0., 0., +1. /)
            offset_vector = (/ 0., 0., -z_plasma/2. /) + plasma_centre
        case(3) ! y-z-plane, negative x
            normal_vector = (/ +1., 0., 0. /)
            offset_vector = (/ -x_plasma/2., 0., 0. /) + plasma_centre
        case(4) ! negative y
            normal_vector = (/ -sin(gamma), +cos(gamma), 0.d0 /)
            offset_vector = (/ -x_plasma/2., -y_plasma/2., 0.0 /) + plasma_centre
        case(5) ! positive y
            normal_vector = (/ -sin(gamma), -cos(gamma), 0.d0 /)
            offset_vector = (/ -x_plasma/2., +y_plasma/2.,0. /) + plasma_centre
        end select
        c_status = dot_product(normal_vector, offset_vector - r)


    case(6) ! hemisphere: flat side facing laser
        select case (face_nr)
        case(1) ! hemisphere
            r_diff = r + (/ r_sphere/2., 0., 0. /) - plasma_centre
            c_status = dot_product(r_diff, r_diff) - r_sphere**2
        case(2) ! y-z-plane
            normal_vector = (/ +1., 0., 0. /)
            offset_vector = (/ -r_sphere/2., 0., 0. /) + plasma_centre
            c_status = dot_product(normal_vector, offset_vector - r)
        end select

    case(16) ! hollow hemisphere
        select case (face_nr)
        case(1) ! outer hemisphere
            r_diff = r + (/ r_sphere/2., 0., 0. /) - plasma_centre
            c_status = dot_product(r_diff, r_diff) - r_sphere**2
        case(2) ! inner hemisphere
            r_diff = r + (/ r_sphere/2., 0., 0. /) - plasma_centre
            c_status = (r_sphere - x_plasma)**2 - dot_product(r_diff, r_diff)
        case(3) ! y-z-plane
            normal_vector = (/ 1., 0., 0. /)
            offset_vector = (/-r_sphere/2., 0., 0./) + plasma_centre
            c_status = dot_product(normal_vector, offset_vector - r)
        end select

    case(26) ! flipped hemisphere: curved side facing laser
        select case (face_nr)
        case(1) ! hemisphere
            r_diff = r - (/ r_sphere/2., 0., 0. /) - plasma_centre
            c_status = dot_product(r_diff, r_diff) - r_sphere**2
        case(2) ! y-z-plane
            normal_vector = (/-1., 0., 0./)
            offset_vector = (/ +r_sphere/2., 0., 0./) + plasma_centre
            c_status = dot_product(normal_vector, offset_vector - r)
        end select


    case(36) ! flipped hollow hemisphere: curved surface facing laser
        select case (face_nr)
        case(1) ! outer hemisphere
            r_diff = r - (/r_sphere/2., 0., 0./) - plasma_centre
            c_status = dot_product(r_diff, r_diff) - r_sphere**2
        case(2) ! inner hemisphere
            r_diff = r - (/r_sphere/2., 0., 0./) - plasma_centre
            c_status = (r_sphere - x_plasma)**2 - dot_product(r_diff, r_diff)
        case(3) ! y-z-plane
            normal_vector = (/-1., 0., 0./)
            offset_vector = (/+r_sphere/2., 0., 0./) + plasma_centre
            c_status = dot_product(normal_vector, offset_vector - r)
        end select

    case default  ! revert to slab
        select case (face_nr)
        case(1) ! x-y-plane, negative z
            normal_vector = (/ 0., 0., +1. /)
            offset_vector = (/ 0., 0., -z_plasma/2. /) + plasma_centre
        case(2) ! x-y-plane, positive z
            normal_vector = (/ 0., 0., -1./)
            offset_vector = (/ 0., 0., +z_plasma/2. /) + plasma_centre
        case(3) ! x-z-plane, negative y
            normal_vector = (/ 0., +1., 0./)
            offset_vector = (/ 0., -y_plasma/2., 0. /) + plasma_centre
        case(4) ! x-z-plane, positive y
            normal_vector = (/ 0., -1., 0./)
            offset_vector = (/ 0., +y_plasma/2., 0. /) + plasma_centre
        case(5) ! y-z-plane, negative x
            normal_vector = (/ +1., 0., 0. /)
            offset_vector = (/ -x_plasma/2., 0., 0. /) + plasma_centre
        case(6) ! y-z-plane, positive x
            normal_vector = (/ -1., 0., 0. /)
            offset_vector = (/ +x_plasma/2., 0., 0. /) + plasma_centre
        end select

        c_status = dot_product(normal_vector, offset_vector - r)

    end select

end subroutine face
