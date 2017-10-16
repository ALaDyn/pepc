! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2017 Juelich Supercomputing Centre, 
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Module for creating Mackay icosahedra
!> see
!> http://de.wikipedia.org/wiki/Ikosaeder
!> Acta Crystallographia 15, 916 (1962)
!>
!> original C-source:
!> http://www.pas.rochester.edu/~wangyt/algorithms/ih/index.html
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_icosahedron

  implicit none
  private

  ! particle types for get_particle
  integer, parameter, public :: type_particle_undef  =-1 
  integer, parameter, public :: type_particle_center = 0 
  integer, parameter, public :: type_particle_corner = 1 
  integer, parameter, public :: type_particle_edge   = 2 
  integer, parameter, public :: type_particle_face   = 3 
  ! half of golden ratio
  real, parameter :: HT = ( sqrt(5.0) + 1.0 ) / 4.0 
  ! corners
  real, parameter :: b(3,12) = reshape([ HT,  0.0,  0.5, &
                                         HT,  0.0, -0.5, &
                                        0.5,   HT,  0.0, &
                                       -0.5,   HT,  0.0, &
                                        0.0,  0.5,   HT, &
                                        0.0, -0.5,   HT, &
                                        0.5,  -HT,  0.0, &
                                        0.0,  0.5,  -HT, &
                                        -HT,  0.0,  0.5, &
                                        0.0, -0.5,  -HT, &
                                        -HT,  0.0, -0.5, &
                                       -0.5,  -HT,  0.0 ], shape(b))
  ! edges
  integer, parameter :: e(2,30) = reshape([ 0, 1, 0, 2, 0, 4, 0, 5, 0, 6, &
                                           10, 3,10, 7,10, 8,10, 9,10,11, &
                                            1, 2, 1, 6, 1, 7, 1, 9, &
                                            8, 3, 8, 4, 8, 5, 8,11, &
                                            2, 3, 2, 4, 2, 7, &
                                           11, 5,11, 6,11, 9, &
                                            6, 5, 6, 9, &
                                            3, 4, 3, 7, &
                                            7, 9, &
                                            5, 4], shape(e))
  ! facets
  integer, parameter :: f(3,20) = reshape([ 0, 1, 2, &
                                            0, 2, 4, &
                                            0, 4, 5, &
                                            0, 5, 6, &
                                            0, 1, 6, &
                                           10, 3, 7, &
                                           10, 3, 8, &
                                           10, 8,11, &
                                           10, 9,11, &
                                           10, 7, 9, &
                                            1, 2, 7, &
                                            1, 7, 9, &
                                            1, 6, 9, &
                                            8, 5,11, &
                                            8, 4, 5, &
                                            8, 3, 4, &
                                            2, 3, 7, &
                                            2, 3, 4, &
                                           11, 5, 6, &
                                           11, 6, 9], shape(f))

  public num_shell_particles
  public num_total_particles
  public get_particle
  public get_nextlower_particles

  contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !> Returns the number of particles in a certain shell, first shell is shell=0
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function num_shell_particles(shell)
      implicit none
      integer :: num_shell_particles
      integer, intent(in) :: shell
      integer :: i

        if    ( shell <  0 ) then
          num_shell_particles =  -1
        elseif( shell == 0 ) then
          num_shell_particles =   1
        elseif( shell == 1 ) then
          num_shell_particles =  12
        elseif( shell == 2 ) then
          num_shell_particles =  42
        else
          num_shell_particles = 0
          num_shell_particles = num_shell_particles + 12            ! corner particles
          num_shell_particles = num_shell_particles + (shell-1)*30  ! edge particles

          do i=1,shell-2
            num_shell_particles = num_shell_particles + i*20 ! side particles
          end do
        endif

    end function num_shell_particles


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !> Returns the total number of particles in an icosahedron with a certain number of shells, first shell is shell=0
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function num_total_particles(shell)
      implicit none
      integer :: num_total_particles
      integer, intent(in) :: shell
      integer :: c
      ! this could also be done by summing over all contributions from num_shell_particles()
      c = shell + 1
      num_total_particles = (10*c*c*c - 15*c*c + 11*c -3)/3

    end function num_total_particles

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !> Returns the number of icosahedron particles right below the given particle number
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function get_nextlower_particles(n)
      implicit none
      integer :: get_nextlower_particles
      integer, intent(in) :: n
      integer :: c
      ! this could also be done by summing over all contributions from num_shell_particles()
      c = 0

      do while (n>num_total_particles(c+1))
        c = c+1
      end do

      get_nextlower_particles = num_total_particles(c)

    end function get_nextlower_particles



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !> Returns the nth (zero-based) particle of a filled icosahedron, number of layers is determined automatically
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function get_particle(n, currlayer, t)
      implicit none
      integer, intent(in) :: n
      integer, intent(out) :: t !< type of particle: corner, edge or face (see above)
      integer, intent(out) :: currlayer
      real, dimension(3) :: get_particle

      if (n < 0) then
        write(*,*) "Requested negative particle index in icosahedron, possibly datatype overflow??"
        stop
      endif

      currlayer = 0
      ! find current layer
      do while ((n > num_total_particles(currlayer)-1) .and. (currlayer < 863))
        currlayer = currlayer + 1
      end do

      if (currlayer == 863) then
        write(*,*) "Requested icosahedron size is to large (more than 2 billion particles in total)"
        stop
      endif

      get_particle = get_particle_from_layer(n - max(num_total_particles(currlayer - 1), 0), currlayer, t)

    end function get_particle

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !> Returns the nth (zero-based) particle of the given icosahedron-shell
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function get_particle_from_layer(n, currlayer, t)
      implicit none
      integer, intent(in) :: n
      integer, intent(out) :: t !< type of particle: corner, edge or face (see above)
      integer, intent(in) :: currlayer
      real, dimension(3) :: get_particle_from_layer
      real, dimension(3) :: e1, e2, e3, v1, v2
      integer :: i, j, k, c

      t = type_particle_undef

      if ((currlayer == 0) .and. (n == 0)) then
        ! only central particle for innermost shell
        get_particle_from_layer = [0, 0, 0]
        t = type_particle_center
      elseif ((currlayer > 0) .and. (n >= 0) .and. (n<12)) then
        ! corner particles
        get_particle_from_layer = b(:,n+1) * currlayer
        t = type_particle_corner
      elseif ((currlayer > 1) .and. (n >= 12) .and. (n<12+30*(currlayer-1))) then
        ! edge particles
        i  = (n-12)  /  (currlayer - 1)
        j  = mod((n-12),(currlayer - 1)) + 1
        e1 = b(:, e(1,i + 1) + 1 ) * currlayer
        e2 = b(:, e(2,i + 1) + 1 ) * currlayer
        get_particle_from_layer = e1 + (e2-e1) * j / real(currlayer)
        t = type_particle_edge
      elseif ((currlayer > 2) .and. (n >= 12+30*(currlayer-1))) then
        ! face particles
        get_particle_from_layer = [0,0,0]
        c = 12+30*(currlayer-1)-1

        do i=0,19
          do j=1,currlayer-2
            do k=1,j
              c = c+1
              if (c==n) then
                e1 = b(:, f(1,i+1)+1 ) * currlayer
                e2 = b(:, f(2,i+1)+1 ) * currlayer
                e3 = b(:, f(3,i+1)+1 ) * currlayer
                v1 = e1 + (e2-e1) * (j+1) / real(currlayer)
                v2 = e1 + (e3-e1) * (j+1) / real(currlayer)
                get_particle_from_layer = v1 + (v2-v1) * k / real(j+1);
                t = type_particle_face
              endif
            end do
          end do
        end do
      else
        write(*,*) "get_particle_from_layer: invalid combination of layer =", currlayer, " and particle =", n
        stop
      endif

    end function get_particle_from_layer

end module module_icosahedron
