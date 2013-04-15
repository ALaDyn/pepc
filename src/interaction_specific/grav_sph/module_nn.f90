! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2013 Juelich Supercomputing Centre, 
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

module module_nn
  implicit none
  private

  public :: nn_prepare_particleresults

  contains

  subroutine nn_prepare_particleresults(t, particles, nparticles)
    use module_tree, only: t_tree, tree_lookup_node_critical, tree_node_get_parent
    use treevars, only: nlev
    use module_pepc_types, only: t_particle, t_tree_node
    use module_interaction_specific_types, only: num_neighbour_particles
    implicit none
    type(t_tree), intent(in) :: t
    type(t_particle), intent(inout) :: particles(nparticles)
    integer, intent(in) :: nparticles

    type(t_tree_node), pointer :: node, parent
    integer :: i
    real*8, dimension(:), allocatable :: boxdiag2

    allocate(boxdiag2(0:nlev))
    boxdiag2(0) = 3.*dot_product(t%bounding_box%boxsize,t%bounding_box%boxsize)
    do i=1,nlev
       boxdiag2(i) =  boxdiag2(i-1)/4.
    end do


    ! for each particle, we traverse the tree upwards, until the current twig
    ! contains more leaves than number of necessary neighbours - as a first guess for the
    ! search radius, we use its diameter
    do i=1,nparticles

       call tree_lookup_node_critical(t, particles(i)%key_leaf, node, 'nn_prepare_particleresults()')

       do
          if (node%leaves >= num_neighbour_particles) then
            ! this twig contains enough particles --> we use its diameter as search radius
            particles(i)%results%maxdist2 = boxdiag2(node%level)
            particles(i)%results%neighbour_keys(1:num_neighbour_particles) = node%key

            exit ! from this loop
          endif

          if (.not. tree_node_get_parent(t, node, parent)) then
             write(*,*) particles(i)
             exit
          endif
          node => parent
       end do

       particles(i)%results%dist2(1:num_neighbour_particles) = particles(i)%results%maxdist2
       particles(i)%results%dist_vector(:,1:num_neighbour_particles) = -13._8 

    end do
  end subroutine nn_prepare_particleresults
end module module_nn
