! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2012 Juelich Supercomputing Centre, 
!                         Forschungszentrum Juelich GmbH,
!                         Germany
! 
! PEPC is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! PEPC is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
! 
! You should have received a copy of the GNU Lesser General Public License
! along with PEPC.  If not, see <http://www.gnu.org/licenses/>.
!

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
    use module_velocity_setup

    implicit none

    integer              :: i, face_nr
    integer              :: iseed1, iseed2, iseed3, n1,p, k, nsphere
    real*8                :: xt, yt, zt
    real                  ::c_status

    iseed1 = -11 - my_rank      ! Select seed depending on PE
    iseed2 = -1001 - my_rank     ! Select seed depending on PE
    iseed3 = -2111 - my_rank      ! Select seed depending on PE
    write (*, '(a,3i8)') 'Seeds: ', iseed1, iseed2, iseed3
    
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

    nsphere = np_local

!     open(90,file='distr0.dat')
!     open(91,file='distr1.dat')
!     open(92,file='distr2.dat')
!     open(93,file='distr3.dat')


    p = 0
    do while (p < nsphere)
        geometry: select case (target_geometry)
        case(0, 5) ! slab or wedge          
            xt = .5 * x_plasma * (2 * rano(iseed1) - 1.) + plasma_centre(1)
            yt = .5 * y_plasma * (2 * rano(iseed1) - 1.) + plasma_centre(2)         
            zt = .5 * z_plasma * (2 * rano(iseed1) - 1.) + plasma_centre(3)
 
        case(1, 7) ! sphere and hollow sphere
            xt = r_sphere * (2 * rano(iseed1) - 1.) + plasma_centre(1)
            yt = r_sphere * (2 * rano(iseed1) - 1.) + plasma_centre(2)        
            zt = r_sphere * (2 * rano(iseed1) - 1.) + plasma_centre(3)
            
        case(2) ! disc
            xt = .5 * x_plasma * (2 * rano(iseed1) - 1.) + plasma_centre(1)
            yt = r_sphere * (2 * rano(iseed1) - 1.) + plasma_centre(2)         
            zt = r_sphere * (2 * rano(iseed1) - 1.) + plasma_centre(3)
            
        case(3) ! wire
            xt = r_sphere * (2 * rano(iseed1) - 1.) + plasma_centre(1)         
            yt = r_sphere * (2 * rano(iseed1) - 1.) + plasma_centre(2)
            zt = .5 * z_plasma * (2 * rano(iseed1) - 1.) + plasma_centre(3)
            
        case(4) ! ellipsoid
            xt = r_sphere * (2 * rano(iseed1) - 1.) * x_plasma + plasma_centre(1)
            yt = r_sphere * (2 * rano(iseed1) - 1.) * y_plasma + plasma_centre(2)
            zt = r_sphere * (2 * rano(iseed1) - 1.) * z_plasma + plasma_centre(3)

        case(6, 8) ! hemisphere and hollow hemisphere
            xt = .5 * r_sphere * (2 * rano(iseed1) - 1.) + plasma_centre(1)
            yt = r_sphere * (2 * rano(iseed1) - 1.) + plasma_centre(2)
            zt = r_sphere * (2 * rano(iseed1) - 1.) + plasma_centre(3)

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
    iseed3 = -17 - 4 * my_rank
    n1 = nsphere
        
    !  exclude odd one out
    if (mod(nsphere,2) .ne. 0) then
        n1 = n1 - 1
    end if

    do i = 1, n1
        k = int(n1 * rano(iseed3) + 1)
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
    q(nep + 1:np_local)       = qi        ! plasma ions (need Z* here)
    m(1:nep)             = mass_e    ! electron mass
    m(nep + 1:np_local)       = mass_i    ! ion mass
    pelabel(1:nep)       = my_rank * nep + (/(i, i = 1, nep)/)      ! Electron labels
    pelabel(nep + 1:np_local) = ne + my_rank * nip + (/(i, i = 1, nip)/) ! Ion labels
! zero fields
  Ex(1:np_local) = 0.
  Ey(1:np_local) = 0.
  Ez(1:np_local) = 0.

  pot(1:np_local) = 0.
  work(1:np_local) = 1.

end subroutine randion



