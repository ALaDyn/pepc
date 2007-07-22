! ==============================================
!
!                RANDION
!
!   $Revision: 1.5
!
!  Sets up 3D random distribution of particles
! 
! ==============================================

subroutine randion

    use physvars
    use treevars
    use utils

    implicit none

    integer              :: i, j, face_nr
    integer              :: idum, iseed1, iseed2, iseed3, i1, n1,p, k, nsphere, nx, ny
    real*8               :: xt, yt, zt, radius, dpx, s
    real*8 :: c_status
    real, dimension(1:3) :: r_temp

    iseed1 = -11 - me      ! Select seed depending on PE
    iseed2 = -1001 - me      ! Select seed depending on PE
    iseed3 = -2111 - me      ! Select seed depending on PE
    if (debug_level>1) write (ipefile, '(a,3i8)') 'Seeds: ', iseed1, iseed2, iseed3
    
    !  Initialise particles according to target geometry
    !       0 = random slab
    !       1 = random sphere
    !       2 = random disc
    !       3 = random wire
    !       4 = ellipsoid
    !       5 = wedge
    !       6 = hemisphere
    !       7 = hollow sphere
    !       8 = hollow hemisphere

    nsphere = npp

    p = 0
    do while (p < nsphere)
        geometry: select case (target_geometry)
        case(0, 5) ! slab or wedge 
            xt = .5 * x_plasma * (2 * rano(iseed1) - 1.) + plasma_center(1)
            yt = .5 * y_plasma * (2 * rano(iseed1) - 1.) + plasma_center(2)         
            zt = .5 * z_plasma * (2 * rano(iseed1) - 1.) + plasma_center(3)
 
        case(1, 7) ! sphere and hollow sphere
            xt = r_sphere * (2 * rano(iseed1) - 1.) + plasma_center(1)
            yt = r_sphere * (2 * rano(iseed1) - 1.) + plasma_center(2)        
            zt = r_sphere * (2 * rano(iseed1) - 1.) + plasma_center(3)
            
        case(2) ! disc
            xt = .5 * x_plasma * (2 * rano(iseed1) - 1.) + plasma_center(1)
            yt = r_sphere * (2 * rano(iseed1) - 1.) + plasma_center(2)         
            zt = r_sphere * (2 * rano(iseed1) - 1.) + plasma_center(3)
            
        case(3) ! wire
            xt = r_sphere * (2 * rano(iseed1) - 1.) + plasma_center(1)         
            yt = r_sphere * (2 * rano(iseed1) - 1.) + plasma_center(2)
            zt = .5 * z_plasma * (2 * rano(iseed1) - 1.) + plasma_center(3)
            
        case(4) ! ellipsoid
            xt = r_sphere * (2 * rano(iseed1) - 1.) * x_plasma + plasma_center(1)
            yt = r_sphere * (2 * rano(iseed1) - 1.) * y_plasma + plasma_center(2)
            zt = r_sphere * (2 * rano(iseed1) - 1.) * z_plasma + plasma_center(3)

        case(6, 8) ! hemisphere and hollow hemisphere
            xt = .5 * r_sphere * (2 * rano(iseed1) - 1.) + plasma_center(1)
            yt = r_sphere * (2 * rano(iseed1) - 1.) + plasma_center(2)
            zt = r_sphere * (2 * rano(iseed1) - 1.) + plasma_center(3)

        end select geometry
            
        ! check if the new particle is in
        ! and if yes then add it
        do face_nr = 1, number_faces
            call face((/xt, yt, zt/), c_status, face_nr)
            if (c_status .ge. 0.) exit
        end do
        if (c_status .lt. 0) then
            p = p + 1
            z(p) = zt 
            y(p) = yt
            x(p) = xt
        end if
    end do
        
    ! scramble to remove correlations
    iseed3 = -17 - 4 * me
    n1 = nsphere
        
    !  exclude odd one out
    if (mod(nsphere,2) .ne. 0) then
        n1 = n1 - 1
    end if

    do i = 1, n1
        k = n1 * rano(iseed3) + 1
        xt = x(i)
        yt = y(i)
        zt = z(i)
        x(i) = x(k)
        y(i) = y(k)
        z(i) = z(k)
        x(k) = xt
        y(k) = yt
        z(k) = zt
    end do

    q(1:nep)             = qe        ! plasma electrons
    q(nep + 1:npp)       = qi        ! plasma ions (need Z* here)
    m(1:nep)             = mass_e    ! electron mass
    m(nep + 1:npp)       = mass_i    ! ion mass
    pepid(1:npp)         = me        ! processor ID
    pelabel(1:nep)       = me * nep + (/(i, i = 1, nep)/)      ! Electron labels
    pelabel(nep + 1:npp) = ne + me * nip + (/(i, i = 1, nip)/) ! Ion labels
! zero fields
  Ex(1:npp) = 0
  Ey(1:npp) = 0
  Ez(1:npp) = 0
  Bx(1:npp) = 0
  By(1:npp) = 0
  Bz(1:npp) = 0
  Ax(1:npp) = 0
  Ay(1:npp) = 0
  Az(1:npp) = 0
  Axo(1:npp) = 0
  Ayo(1:npp) = 0
  Azo(1:npp) = 0
  pot(1:npp) = 0
  work(1:npp) = 1.   ! set work load balanced initially


end subroutine randion

