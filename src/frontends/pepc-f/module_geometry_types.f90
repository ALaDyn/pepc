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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!> Contains types for building the sim domain form walls
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_geometry_types
      implicit none

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public type declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      !> Data structure for storing interaction-specific particle data
      type t_boundary
         real(KIND=8) :: x0(3)                !< reference point
         real(KIND=8) :: e1(3)                !< edge1: vector to adjacent corner
         real(KIND=8) :: e2(3)                !< edge2: vector to other adjacent corner; adjacent and opposite corner are computed
         real(KIND=8) :: n(3)                 !< surface normal
         real(KIND=8) :: A = 0._8             !< surface Area
         logical :: rectangle=.false.         !< TRUE if e1 and e2 are perpendicular
         integer :: type                !< surface type
                                        !<  0 absorbing wall with uniform charge distribution realized with wall particles
                                        !<  1 absorbing wall with uniform charge distribution realized with analytic field and potential
                                        !<  2 open boundary
                                        !<  3 logical sheath (not available at the moment)
                                        !<  4 reflecting boundary
                                        !<  5 immediate half-Maxwellian refluxing normal to surface, tangential v conserved
                                        !<  6 immediate half-Maxwellian refluxing normal to surface, tangential v resampled
                                        !<  7 immediate Maxwellian flux refluxing normal to surface, tangential v conserved
                                        !<  8 immediate Maxwellian flux refluxing normal to surface, tangential v resampled
                                        !<  9 immediate drifting Maxwellian flux refluxing normal to surface, tangential v conserved
                                        !< 10 immediate drifting Maxwellian flux refluxing normal to surface, tangential v resampled
                                        !< 11 periodic boundary (an according boundary on the other side of the system has to be set)
                                        !< -1 virtual boundary (used for diagnostic purposes)
         logical :: reflux_particles    !< if true, particles hitting the boundary will be refluxed according to chosen source
         logical :: accumulate_charge   !> if true, incident charge is accumulated (and can be create external fields)
         integer :: indx                !< index (should be the same as in the boundary array)
         integer :: opp_bnd=0           !< opposite boundary if periodic bc's
         real(KIND=8) :: dist=0.              !< distance to opposite boundary (only set if type=2)
         integer :: nwp=0, nwpe1=0, nwpe2=0               !< number of wall particles (total, along e1 and e2) (can only be set if type=0)
         integer, allocatable :: wp_labels(:)             !< wp labels for this wall
         integer :: wp_label_max=0, wp_label_min=0        !< min and max value of wp_label. Is used during initialization to find out
                                                          !< whether a particle is part of th boundary and in the wp source function to
                                                          !< find the particle position in wppe1 and wppe2
         real(kind=8), allocatable :: wppe1(:), wppe2(:)  !< wall particle positions along e1 and e2 (relative values from 0 to 1)
         real*8 :: q_tot=0.             !< total caharge of wall (can only be set if type=0)

      end type t_boundary


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  subroutine-implementation  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      contains

end module module_geometry_types
