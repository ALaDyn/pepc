! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
!
! Copyright (C) 2002-2014 Juelich Supercomputing Centre,
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

!>
!> Encapsulates anything that is directly involved in force calculation
!> and multipole expansion manipulation
!> i.e. shifting along the tree, computing forces between particles and cluster, etc.
!>
module module_interaction_specific
  use module_pepc_types
  use module_interaction_specific_types
  implicit none
  save
  private

  integer, public :: force_law    = 5      !< 5=NN-list "interaction"
  integer, public :: mac_select   = 3      !< selector for multipole acceptance criterion
  ! 0 = barnes-Hut-MAC
  ! 1 = Bmax-MAC
  ! 2 = N^2 MAC
  ! 3 = asymmetric NN-MAC
  ! 4 = symmetric NN-MAC

  real*8, public  :: theta2       = 0.36  !< square of multipole opening angle
  real*8, public  :: eps2         = 0.0    !< square of short-distance cutoff parameter for plummer potential (0.0 corresponds to classical Coulomb)


  namelist /calc_force_nearestneighbour/ force_law

  ! currently, all public functions in module_interaction_specific are obligatory
  public multipole_from_particle
  public shift_multipoles_up
  public results_add
  public calc_force_per_interaction_with_leaf
  public calc_force_per_interaction_with_twig
  public calc_force_per_particle
  public mac
  public particleresults_clear
  public calc_force_read_parameters
  public calc_force_write_parameters
  public calc_force_finalize
  public calc_force_prepare
  public calc_force_after_grow
  public get_number_of_interactions_per_particle
  public pack_particle_list
  public unpack_particle_list

contains

  !>
  !> Computes multipole properties of a single particle
  !>
  subroutine multipole_from_particle(particle_pos, particle, multipole)
    implicit none
    real*8, intent(in) :: particle_pos(3)
    type(t_particle_data), intent(in) :: particle
    type(t_tree_node_interaction_data), intent(out) :: multipole

    ! use velocity (v) at same time step as coordinate, not v_minus_half
    multipole = t_tree_node_interaction_data(particle_pos, & ! coc
         particle%particle_id, &                             ! unique particle id for neighbour-test
         particle%q, &                                       ! charge (here mass)
         [0., 0., 0.], &                                     ! dipole moment
         [0., 0., 0.], &                                     ! diagonal quadrupole moments
         0., 0., 0., &                                       ! other quadrupole moments
         0., &                                               ! bmax
         particle%v, &                                       ! particle velocity
         particle%temperature, &                             ! particle temperature
         0., &                                               ! sph density
         0. )                                                ! sph smoothing-length
  end subroutine multipole_from_particle


  !>
  !> Accumulates multipole properties of child nodes to parent node
  !>
  subroutine shift_multipoles_up(parent, children)
    implicit none
    type(t_tree_node_interaction_data), intent(out) :: parent
    type(t_tree_node_interaction_data), intent(in) :: children(:)

    integer :: nchild, j

    real*8 :: shift(1:3)

    nchild = size(children)

    parent%charge     = SUM( children(1:nchild)%charge )

    parent%coc = [0._8, 0._8, 0._8]
    parent%particle_id = 0_8


    if (parent%charge .ne. 0.) then
       ! use center-of-charge because we may divide by charge
       do j=1,nchild
          parent%coc(1:3) = parent%coc(1:3) + ( children(j)%coc(1:3) * children(j)%charge )
       end do

       parent%coc(1:3) = parent%coc(1:3) / parent%charge
    else
       ! use geometric center
       do j=1,nchild
          parent%coc(1:3) = parent%coc(1:3) +   children(j)%coc(1:3)
       end do

       parent%coc(1:3) = parent%coc(1:3) / nchild
    endif

    ! multipole properties
    parent%dip    = [0., 0., 0.]
    parent%quad   = [0., 0., 0.]
    parent%xyquad = 0.
    parent%yzquad = 0.
    parent%zxquad = 0.

    do j=1,nchild
       shift(1:3) = parent%coc(1:3) - children(j)%coc

       ! dipole moment
       parent%dip = parent%dip + children(j)%dip - children(j)%charge*shift(1:3)

       ! quadrupole moment
       parent%quad(1:3) = parent%quad(1:3) + children(j)%quad(1:3) - 2*children(j)%dip(1:3)*shift(1:3) + children(j)%charge*shift(1:3)**2

       parent%xyquad = parent%xyquad + children(j)%xyquad - children(j)%dip(1)*shift(2) - children(j)%dip(2)*shift(1) + children(j)%charge*shift(1)*shift(2)
       parent%yzquad = parent%yzquad + children(j)%yzquad - children(j)%dip(2)*shift(3) - children(j)%dip(3)*shift(2) + children(j)%charge*shift(2)*shift(3)
       parent%zxquad = parent%zxquad + children(j)%zxquad - children(j)%dip(3)*shift(1) - children(j)%dip(1)*shift(3) + children(j)%charge*shift(3)*shift(1)
    end do

    parent%bmax = maxval(sqrt((parent%coc(1)-children(1:nchild)%coc(1))**2+(parent%coc(2)-children(1:nchild)%coc(2))**2+(parent%coc(3)-children(1:nchild)%coc(3))**2) + children(1:nchild)%bmax)

    ! set velocity, temperature and rho for tree node to a nonsense value which may be recognised as nonsense when one tries to use them
    parent%v = [-13._8,-13._8,-13._8]
    parent%temperature = -13._8
    parent%rho = -13._8

    ! this is stupid in the construction of the tree, because all h are zero. But it is necessary to rebuild this multipole property later for the symmetric neighbour search
    parent%h = maxval(children(1:nchild)%h)
  end subroutine shift_multipoles_up


  !>
  !> adds res2 to res1
  !>
  subroutine results_add(res1, res2)
    implicit none
    type(t_particle_results), intent(inout) :: res1
    type(t_particle_results), intent(in) :: res2

    if( (force_law .eq. 2) .or. (force_law .eq. 3) ) then ! for gravity
       res1%e    = res1%e    + res2%e
       res1%pot  = res1%pot  + res2%pot

    else if( force_law .eq. 5) then
       ! TODONN: ist das wirklich alles?
       res1 = res2
    else
       write(*,*) "!!! value of force_law is not allowed in results_add:", force_law
    end if
  end subroutine results_add


  !>
  !> reads interaction-specific parameters from file
  !>
  subroutine calc_force_read_parameters(filehandle)
    use module_debug, only: pepc_status
    implicit none
    integer, intent(in) :: filehandle

    call pepc_status("READ PARAMETERS, section calc_force_nearestneighbour")
    read(filehandle, NML=calc_force_nearestneighbour)
  end subroutine


  !>
  !> writes interaction-specific parameters to file
  !>
  subroutine calc_force_write_parameters(filehandle)
    use module_debug, only: pepc_status
    implicit none
    integer, intent(in) :: filehandle

    write(filehandle, NML=calc_force_nearestneighbour)
  end subroutine


  !>
  !> computes derived parameters for calc force module
  !>
  subroutine calc_force_prepare()
    implicit none
    ! nothing to do here
    !        call fmm_framework_init(my_rank, wellsep=mirror_box_layers)
    ! TODO: what ist that good for?
    ! currently no periodic gravity
  end subroutine calc_force_prepare


  !>
  !> initializes static variables of calc force module that depend
  !> on particle data and might be reused on subsequent traversals
  !>
  subroutine calc_force_after_grow(particles)
    use module_pepc_types
    implicit none
    type(t_particle), dimension(:), intent(in) :: particles

    ! nothing to be done here for now
  end subroutine


  !>
  !> subroutine must return the estimated number of iteractions per particle
  !> for the current mac and/or parameters and the supplied total number of particles
  !>
  subroutine get_number_of_interactions_per_particle(npart_total, nintmax)
    implicit none
    integer(kind_particle), intent(in) :: npart_total !< total number of particles
    integer(kind_node), intent(out) :: nintmax !< maximum number of interactions per particle

    real*8 :: invnintmax !< inverse of nintmax to avoid division by zero for theta == 0.0

    real*8, parameter :: theta2 = 0.3**2 !TODO: this function must be adapted to the actual MAC

    ! Estimate of interaction list length - Hernquist expression
    ! applies for BH-MAC
    invnintmax = max(theta2 / (35._8*log(1._8*npart_total)) , 1._8/npart_total)
    nintmax    = int(1._8/invnintmax)
  end subroutine get_number_of_interactions_per_particle


  !>
  !> finalizes the calc force module at end of simulation
  !>
  subroutine calc_force_finalize()
    implicit none
    ! nothing to do here
  end subroutine calc_force_finalize


  !>
  !> generic Multipole Acceptance Criterion
  !>
  function mac(particle, node, dist2, boxlength2)
    implicit none

    logical :: mac
    type(t_tree_node_interaction_data), intent(in) :: node
    type(t_particle), intent(in) :: particle
    real*8, intent(in) :: dist2
    real*8, intent(in) :: boxlength2

    select case (mac_select)
    case(0)
       ! Barnes-Hut-MAC
       ! mac = (theta2 * dist2 > boxlength2)
    case(1)
       ! Bmax-MAC
       mac = (theta2 * dist2 > min(node%bmax**2,3.0*boxlength2)) !TODO: Can we put the min into bmax itself? And **2?
    case(2)
       ! N^2 code
       mac = .false.
    case(3)
       ! NN-MAC: we may "interact" with the node if it is further away than maxdist2 --> this leads to the node *not* being put onto the NN-list (strange, i know)
       ! first line: original formulation, last line: after transition to formulation with only one square root
       ! mac = sqrt(dist2) - sqrt(3.* boxlength2)  >  sqrt(results%maxdist2)                                ! + sqrt(3.*boxlength2)
       !     = sqrt(dist2)                         >  sqrt(results%maxdist2) + sqrt(3.*boxlength2)          ! ^2
       !     =      dist2                          > (sqrt(results%maxdist2) + sqrt(3.*boxlength2))**2
       !     =      dist2                          > results%maxdist2 + 2.*sqrt( 3.*results%maxdist2*boxlength2) + 3.*boxlength2
       mac =      dist2                          > particle%results%maxdist2 +    sqrt(12.*particle%results%maxdist2*boxlength2) + 3.*boxlength2
       ! TODO NN: this estimation should be evaluated without (!!) any square roots for performance reasons (which does not seem to be trivial)
    case(4)
       ! dist - node_diagonal > maxdist2                             .OR. dist - node_diagonal > 2*h2
       ! sqrt(dist2) - sqrt(3.*boxlength2) > sqrt(results%maxdist2)  .OR. sqrt(dist2) - sqrt(3.*boxlength2) > tree_nodes(node)%h*2.  ! + sqrt(3. * boxlength2)
       ! sqrt(dist2) > sqrt(results%maxdist2) + sqrt(3.*boxlength2)  .OR. sqrt(dist2) > tree_nodes(node)%h*2. + sqrt(3.*boxlength2)  ! ^2
       ! dist2 > results%maxdist2 + 2.*sqrt(results%maxdist2*3.*boxlength2) + 3.*boxlength2 .OR. dist2 > 4.*tree_nodes(node)%h**2 + 2.*sqrt(2.*tree_nodes(node)%h*3.*boxlength2) + 3. *boxlength2
       !
       mac = ( (dist2 > particle%results%maxdist2 + sqrt(12.*particle%results%maxdist2*boxlength2) + 3.*boxlength2) .or. &
            (   dist2 > 4.*node%h**2 + sqrt(24.*node%h *boxlength2) + 3.*boxlength2) )

    case default
       write(*,*) "!!! value of mac_select is not allowed in mac:", mac_select
    end select
  end function mac


  !>
  !> clears result in t_particle datatype - usually, this function does not need to be touched
  !> due to dependency on module_pepc_types and(!) on module_interaction_specific, the
  !> function cannot reside in module_interaction_specific that may not include module_pepc_types
  !>
  subroutine particleresults_clear(particles)
    use module_pepc_types, only: t_particle
    implicit none
    type(t_particle), intent(inout) :: particles(:)

    integer(kind_particle) :: i,j

    do i=1,size(particles, kind=kind(i))
       particles(i)%results%h         = 0._8
       particles(i)%results%rho       = 0._8
       particles(i)%results%sph_force = [0._8, 0._8, 0._8]
       particles(i)%results%e         = [0._8, 0._8, 0._8]
       particles(i)%results%pot       = 0._8

       particles(i)%results%maxdist2             = huge(0._8)
       do j=1,size(particles(i)%results%neighbour_nodes)
         particles(i)%results%neighbour_nodes(j)%p => null()
       end do
       particles(i)%results%maxidx               = 1
    end do
  end subroutine particleresults_clear


  !>
  !> Force calculation wrapper.
  !> This function is thought for pre- and postprocessing of
  !> calculated fields, and for being able to call several
  !> (different) force calculation routines
  !>
  subroutine calc_force_per_interaction_with_leaf(delta, dist2, particle_pack, node_data)
    use module_pepc_types
    use module_debug
    use module_coulomb_kernels
    use treevars
    implicit none

    real*8, intent(in) :: delta(:,:)
    real*8, intent(in) :: dist2(:)
    type(t_particle_pack), intent(inout) :: particle_pack
    type(t_tree_node_interaction_data), target, intent(in) :: node_data
    integer(kind_particle) :: np
    
    np = size(dist2,kind=kind_particle)

    select case (force_law)
    case (2)  !  compute 2D-Coulomb fields and potential of particle p from its interaction list
       call calc_force_coulomb_2D_direct(np, delta, particle_pack%ex, particle_pack%ey, particle_pack%pot, node_data, eps2)
    case (3)  !  compute 3D-Coulomb fields and potential of particle p from its interaction list
       call calc_force_coulomb_3D_direct(np, delta, dist2, particle_pack%ex, particle_pack%ey, particle_pack%ez, particle_pack%pot, node_data, eps2)
    case (4)  ! LJ potential for quiet start
       call calc_force_LJ(np, delta, dist2, particle_pack%ex, particle_pack%ey, particle_pack%ez, node_data, eps2)
    case (5)
        call update_nn_list(particle_pack, node_data, delta, dist2)
    case default
       DEBUG_ERROR(*, "value of force_law is not allowed in calc_force_per_interaction:", force_law)
    end select
  end subroutine


  !>
  !> Force calculation wrapper.
  !> This function is thought for pre- and postprocessing of
  !> calculated fields, and for being able to call several
  !> (different) force calculation routines
  !>
  subroutine calc_force_per_interaction_with_twig(delta, dist2, particle_pack, node_data)
    use module_pepc_types
    use module_debug
    use module_coulomb_kernels
    use treevars
    implicit none

    real*8, intent(in) :: delta(:,:)
    real*8, intent(in) :: dist2(:)
    type(t_particle_pack), intent(inout) :: particle_pack
    type(t_tree_node_interaction_data), target, intent(in) :: node_data
    integer(kind_particle) :: np

    np = size(dist2,kind=kind_particle)

    select case (force_law)
    case (2)  !  compute 2D-Coulomb fields and potential of particle p from its interaction list
       call calc_force_coulomb_2D(np, delta, particle_pack%ex, particle_pack%ey, particle_pack%pot, node_data, eps2)
    case (3)  !  compute 3D-Coulomb fields and potential of particle p from its interaction list
       call calc_force_coulomb_3D(np, delta, dist2, particle_pack%ex, particle_pack%ey, particle_pack%ez, particle_pack%pot, node_data, eps2)
    case (4)  ! LJ potential for quiet start
       call calc_force_LJ(np, delta, dist2, particle_pack%ex, particle_pack%ey, particle_pack%ez, node_data, eps2)
    case (5)
       call update_nn_list(particle_pack, node_data, delta, dist2)
    case default
      DEBUG_ERROR(*, "value of force_law is not allowed in calc_force_per_interaction:", force_law)
    end select
  end subroutine


  !>
  !> Force calculation wrapper for contributions that only have
  !> to be added once per particle
  !>
  subroutine calc_force_per_particle(particles)
    implicit none

    type(t_particle), intent(inout) :: particles(:)

    ! currently nothing to do here
  end subroutine calc_force_per_particle


  ! >
  ! > If particle closer than the furthest of the current n closest particles
  ! > put in on the lists instead of the formerly furthest particle and update
  ! > the distance to the furthest particle.
  ! >
  subroutine update_nn_list(particle_pack, node_data, delta, dist2)
    use module_pepc_types
    use treevars
    use module_debug
    implicit none
    include 'mpif.h'

    real*8, intent(in) :: delta(:,:)
    real*8, intent(in) :: dist2(:)
    type(t_particle_pack), intent(inout) :: particle_pack
    type(t_tree_node_interaction_data), target, intent(in) :: node_data

    integer :: tmp(1), p


    select case (mac_select)
    case(3)

       do p=1,size(dist2)
         if (dist2(p) < particle_pack%maxdist2(p)) then
           ! add node to NN_list
           particle_pack%neighbour_nodes(p, particle_pack%maxidx(p))%p => node_data
           particle_pack%dist2(p, particle_pack%maxidx(p))           = dist2(p)
           particle_pack%dist_vector(p, :,particle_pack%maxidx(p))   = delta(p,:)
           tmp = maxloc(particle_pack%dist2(p, 1:num_neighbour_particles)) ! this is really ugly, but maxloc returns a 1-by-1 vector instead of the expected scalar
           particle_pack%maxidx(p)   = tmp(1)
           particle_pack%maxdist2(p) = particle_pack%dist2(p, particle_pack%maxidx(p))
         else
           ! node is further away than farest particle in nn-list --> can be ignored
         endif
       end do

    case(4)

       do p=1,size(dist2)
         if ( (dist2(p) < particle_pack%maxdist2(p)) .or. (dist2(p) < 4.*node_data%h*node_data%h)) then
           ! add node to NN_list
           particle_pack%neighbour_nodes(p, particle_pack%maxidx(p))%p => node_data
           particle_pack%dist2(p, particle_pack%maxidx(p))           = dist2(p)
           particle_pack%dist_vector(p, :,particle_pack%maxidx(p))   = delta(p,:)
           
           particle_pack%maxidx(p) = particle_pack%maxidx(p) + 1
           if( particle_pack%maxidx(p) > max_neighbour_particles ) then
             DEBUG_ERROR(*, "Number of neighbours found in symmetric neighbour search bigger than max_neighbour_particles. Increase this value, recompile and try again.")
           end if
         else
           ! node is further away than farest particle in nn-list --> can be ignored
         endif
       end do

    case default
       DEBUG_ERROR(*, "value of mac_select not allowed in update_nn_list:", mac_select)

    end select
  end subroutine update_nn_list

  subroutine pack_particle_list(particles, packed)
    use module_pepc_types, only: t_particle, kind_particle
    implicit none

    type(t_particle), intent(in) :: particles(:)
    type(t_particle_pack), intent(inout) :: packed

    integer(kind_particle) :: ip, np

    np = size(particles, kind = kind_particle)

    allocate(packed%maxdist2(np), &
             packed%maxidx(np), &
             packed%neighbour_nodes(np, max_neighbour_particles), &
             packed%dist2(np, max_neighbour_particles), &
             packed%dist_vector(np, 3, max_neighbour_particles), &
             packed%q(np), &
             packed%ex(np), &
             packed%ey(np), &
             packed%ez(np), &
             packed%pot(np) &
             )

    do ip = 1, np
      packed%maxdist2(ip) = particles(ip)%results%maxdist2
      packed%maxidx(ip)   = particles(ip)%results%maxidx
      packed%neighbour_nodes(ip,:) = particles(ip)%results%neighbour_nodes
      packed%dist2(ip,:)  = particles(ip)%results%dist2
      packed%dist_vector(ip,:,:)   = particles(ip)%results%dist_vector(:,:)
      packed%q(ip) = particles(ip)%data%q
      packed%ex(ip)  = particles(ip)%results%e(1)
      packed%ey(ip)  = particles(ip)%results%e(2)
      packed%ez(ip)  = particles(ip)%results%e(3)
      packed%pot(ip) = particles(ip)%results%pot
    end do
  end subroutine pack_particle_list

  subroutine unpack_particle_list(packed, particles)
    use module_pepc_types, only: t_particle, kind_particle
    use module_debug
    implicit none

    type(t_particle_pack), intent(inout) :: packed
    type(t_particle), intent(inout) :: particles(:)

    integer(kind_particle) :: ip, np

    np = size(particles, kind = kind_particle)
    DEBUG_ASSERT(np == size(packed%maxdist2, kind = kind_particle))

    do ip = 1, np
      particles(ip)%results%maxdist2 = packed%maxdist2(ip)
      particles(ip)%results%maxidx = packed%maxidx(ip)
      particles(ip)%results%neighbour_nodes = packed%neighbour_nodes(ip,:)
      particles(ip)%results%dist2 = packed%dist2(ip,:)
      particles(ip)%results%dist_vector(:,:) = packed%dist_vector(ip,:,:)
      particles(ip)%data%q = packed%q(ip)
      particles(ip)%results%e(1) = packed%ex(ip)
      particles(ip)%results%e(2) = packed%ey(ip)
      particles(ip)%results%e(3) = packed%ez(ip)
      particles(ip)%results%pot = packed%pot(ip)
    end do

    deallocate(packed%maxdist2, &
               packed%maxidx, &
               packed%neighbour_nodes, &
               packed%dist2, &
               packed%dist_vector, &
               packed%q, &
               packed%ex, &
               packed%ey, &
               packed%ez, &
               packed%pot &
               )
  end subroutine unpack_particle_list
end module module_interaction_specific
