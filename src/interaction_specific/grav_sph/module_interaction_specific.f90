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

!>
!> Encapsulates anything that is directly involved in force calculation
!> and multipole expansion manipulation
!> i.e. shifting along the tree, computing forces between particles and cluster, etc.
!>
module module_interaction_specific
  use module_pepc_kinds
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
  public calc_force_per_interaction_with_self
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

    integer(kind_particle) :: i

    do i=1,size(particles, kind=kind(i))
       particles(i)%results%h         = 0._8
       particles(i)%results%rho       = 0._8
       particles(i)%results%sph_force = [0._8, 0._8, 0._8]
       particles(i)%results%e         = [0._8, 0._8, 0._8]
       particles(i)%results%pot       = 0._8

       particles(i)%results%maxdist2           = huge(0._8)
       particles(i)%results%neighbour_nodes(:) = -1 ! FIXME: this should be NODE_INVALID
       particles(i)%results%maxidx             = 1
    end do
  end subroutine particleresults_clear


  !>
  !> Force calculation wrapper.
  !> This function is thought for pre- and postprocessing of
  !> calculated fields, and for being able to call several
  !> (different) force calculation routines
  !>
  subroutine calc_force_per_interaction_with_self(particle, node, node_idx, delta, dist2, vbox)
    use module_pepc_types
    use treevars
    implicit none

    type(t_tree_node_interaction_data), intent(in) :: node
    integer(kind_node), intent(in) :: node_idx
    type(t_particle), intent(inout) :: particle
    real*8, intent(in) :: vbox(3), delta(3), dist2
  end subroutine


  !>
  !> Force calculation wrapper.
  !> This function is thought for pre- and postprocessing of
  !> calculated fields, and for being able to call several
  !> (different) force calculation routines
  !>
  subroutine calc_force_per_interaction_with_leaf(particle, node, node_idx, delta, dist2, vbox)
    use module_pepc_types
    use treevars
    implicit none

    type(t_tree_node_interaction_data), intent(in) :: node
    integer(kind_node), intent(in) :: node_idx
    type(t_particle), intent(inout) :: particle
    real*8, intent(in) :: vbox(3), delta(3), dist2

    real*8 :: exyz(3), phic

    select case (force_law)
    case (2)  !  compute 2D-Coulomb fields and potential of particle p from its interaction list
       call calc_force_coulomb_2D_direct(node, delta(1:2), dot_product(delta(1:2), delta(1:2)), exyz(1), exyz(2),phic)
       exyz(3) = 0.

       particle%results%e         = particle%results%e    + exyz
       particle%results%pot       = particle%results%pot  + phic
    case (3)  !  compute 3D-Coulomb fields and potential of particle p from its interaction list
       call calc_force_coulomb_3D_direct(node, delta, dist2, exyz(1), exyz(2), exyz(3), phic)

       particle%results%e         = particle%results%e    + exyz
       particle%results%pot       = particle%results%pot  + phic
    case (4)  ! LJ potential for quiet start
       call calc_force_LJ(node, delta, dist2, exyz(1), exyz(2), exyz(3), phic)
       exyz(3) = 0.

       particle%results%e         = particle%results%e    + exyz
       particle%results%pot       = particle%results%pot  + phic
    case (5)
       call update_nn_list(particle, node, node_idx, delta, dist2)
    case default
       write(*,*) "value of force_law is not allowed in calc_force_per_interaction:", force_law
    end select
  end subroutine


  !>
  !> Force calculation wrapper.
  !> This function is thought for pre- and postprocessing of
  !> calculated fields, and for being able to call several
  !> (different) force calculation routines
  !>
  subroutine calc_force_per_interaction_with_twig(particle, node, node_idx, delta, dist2, vbox)
    use module_pepc_types
    use treevars
    implicit none

    type(t_tree_node_interaction_data), intent(in) :: node
    integer(kind_node), intent(in) :: node_idx
    type(t_particle), intent(inout) :: particle
    real*8, intent(in) :: vbox(3), delta(3), dist2

    real*8 :: exyz(3), phic

    select case (force_law)
    case (2)  !  compute 2D-Coulomb fields and potential of particle p from its interaction list
       call calc_force_coulomb_2D(node, delta(1:2), dot_product(delta(1:2), delta(1:2)), exyz(1), exyz(2),phic)
       exyz(3) = 0.

       particle%results%e         = particle%results%e    + exyz
       particle%results%pot       = particle%results%pot  + phic
    case (3)  !  compute 3D-Coulomb fields and potential of particle p from its interaction list
       call calc_force_coulomb_3D(node, delta, dist2, exyz(1), exyz(2), exyz(3), phic)

       particle%results%e         = particle%results%e    + exyz
       particle%results%pot       = particle%results%pot  + phic
    case (4)  ! LJ potential for quiet start
       call calc_force_LJ(node, delta, dist2, exyz(1), exyz(2), exyz(3), phic)
       exyz(3) = 0.

       particle%results%e         = particle%results%e    + exyz
       particle%results%pot       = particle%results%pot  + phic
    case (5)
       call update_nn_list(particle, node, node_idx, delta, dist2)
    case default
       write(*,*) "value of force_law is not allowed in calc_force_per_interaction:", force_law
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
  subroutine update_nn_list(particle, node, node_idx, d, dist2)
    use module_interaction_specific_types, only: &
         max_neighbour_particles

    use module_pepc_types
    use treevars
    use mpi
    implicit none

    integer*8, intent(in) :: node_idx !< node index of particle to interact with
    type(t_tree_node_interaction_data), intent(in) :: node
    real*8, intent(in) :: d(3), dist2 !< separation vector and magnitude**2 precomputed in walk_single_particle
    type(t_particle), intent(inout) :: particle

    integer :: ierr, tmp(1)

    select case (mac_select)
    case(3)

       if (dist2 < particle%results%maxdist2) then
          ! add node to NN_list
          particle%results%neighbour_nodes(particle%results%maxidx) = node_idx
          particle%results%dist2(particle%results%maxidx)           = dist2
          particle%results%dist_vector(:,particle%results%maxidx)   = d
          tmp                       = maxloc(particle%results%dist2(1:num_neighbour_particles)) ! this is really ugly, but maxloc returns a 1-by-1 vector instead of the expected scalar
          particle%results%maxidx   = tmp(1)
          particle%results%maxdist2 = particle%results%dist2(particle%results%maxidx)
       else
          ! node is further away than furthest particle in nn-list --> can be ignored
       endif

    case(4)

       if ( (dist2 < particle%results%maxdist2) .or. ( dist2 < 4.*node%h*node%h ) ) then
          ! add node to NN_list
          particle%results%neighbour_nodes(particle%results%maxidx) = node_idx
          particle%results%dist2(particle%results%maxidx)           = dist2
          particle%results%dist_vector(:,particle%results%maxidx)   = d
          particle%results%maxidx   = particle%results%maxidx + 1
          if( particle%results%maxidx > max_neighbour_particles ) then
             write(*,*) "Number of neighbours found in symmetric neighbour search bigger than max_neighbour_particles. Increase this value, recompile and try again."
             call MPI_ABORT(MPI_COMM_WORLD, ierr)
          end if

       else
          ! node to far away --> can be ignored
       endif

    case default
       write(*,*) "value of mac_select not allowed in update_nn_list:", mac_select

    end select
  end subroutine update_nn_list


        !>
        !> Calculates 3D Coulomb interaction of particle p with tree node inode
        !> that is shifted by the lattice vector vbox
        !> results are returned in eps, sumfx, sumfy, sumfz, sumphi
        !>
        subroutine calc_force_coulomb_3D(t, d, dist2, sumfx, sumfy, sumfz, sumphi)
          use module_pepc_types
          use treevars
          use mpi
          implicit none

          type(t_tree_node_interaction_data), intent(in) :: t !< index of particle to interact with
          real*8, intent(in) :: d(3), dist2 !< separation vector and magnitude**2 precomputed in walk_single_particle
          real*8, intent(out) ::  sumfx,sumfy,sumfz,sumphi

          real*8 :: rd,dx,dy,dz,r,dx2,dy2,dz2
          real*8 :: dx3,dy3,dz3,rd3,rd5,rd7,fd1,fd2,fd3,fd4,fd5,fd6

             sumfx  = 0.
             sumfy  = 0.
             sumfz  = 0.
             sumphi = 0.

             !  preprocess distances
             dx = d(1)
             dy = d(2)
             dz = d(3)


             r = sqrt(dist2+eps2)
             rd = 1./r
             rd3 = rd**3
             rd5 = rd**5
             rd7 = rd**7

             dx2 = dx**2
             dy2 = dy**2
             dz2 = dz**2
             dx3 = dx**3
             dy3 = dy**3
             dz3 = dz**3

             fd1 = 3.*dx2*rd5 - rd3
             fd2 = 3.*dy2*rd5 - rd3
             fd3 = 3.*dz2*rd5 - rd3
             fd4 = 3.*dx*dy*rd5
             fd5 = 3.*dy*dz*rd5
             fd6 = 3.*dx*dz*rd5

             ! potential

             sumphi = sumphi + t%charge*rd    &                           !  monopole term
                                        !
                  + (dx*t%dip(1) + dy*t%dip(2) + dz*t%dip(3))*rd3  &    !  dipole
                                        !     Dx             Dy            Dz
                  + 0.5*fd1*t%quad(1) + 0.5*fd2*t%quad(2) + 0.5*fd3*t%quad(3)  &  !  quadrupole
                                        !           Qxx                 Qyy                 Qzz
                  + fd4*t%xyquad + fd5*t%yzquad + fd6*t%zxquad
             !   Qxy            Qyz             Qzx

             !  forces

             sumfx = sumfx + t%charge*dx*rd3 &      ! monopole term
                                        !
                  + fd1*t%dip(1) + fd4*t%dip(2) + fd6*t%dip(3)   &   !  dipole term
                                        !
                  + (15.*dx3*rd7 - 9.*dx*rd5 )*0.5*t%quad(1) &     !
                  + ( 15.*dy*dx2*rd7 - 3.*dy*rd5 )*t%xyquad &     !
                  + ( 15.*dz*dx2*rd7 - 3.*dz*rd5 )*t%zxquad &     !   quadrupole term
                  + ( 15*dx*dy*dz*rd7 )*t%yzquad &                !
                  + ( 15.*dx*dy2*rd7 - 3.*dx*rd5 )*0.5*t%quad(2) & !
                  + ( 15.*dx*dz2*rd7 - 3.*dx*rd5 )*0.5*t%quad(3)   !

             sumfy = sumfy + t%charge*dy*rd3 &
                  + fd2*t%dip(2) + fd4*t%dip(1) + fd5*t%dip(3)  &
                  + ( 15.*dy3*rd7 - 9.*dy*rd5 )*0.5*t%quad(2) &
                  + ( 15.*dx*dy2*rd7 - 3.*dx*rd5 )*t%xyquad &
                  + ( 15.*dz*dy2*rd7 - 3.*dz*rd5 )*t%yzquad &
                  + ( 15.*dx*dy*dz*rd7 )*t%zxquad &
                  + ( 15.*dy*dx2*rd7 - 3.*dy*rd5 )*0.5*t%quad(1) &
                  + ( 15.*dy*dz2*rd7 - 3.*dy*rd5 )*0.5*t%quad(3)

             sumfz = sumfz + t%charge*dz*rd3 &
                  + fd3*t%dip(3) + fd5*t%dip(2) + fd6*t%dip(1)  &
                  + ( 15.*dz3*rd7 - 9.*dz*rd5 )*0.5*t%quad(3) &
                  + ( 15.*dx*dz2*rd7 - 3.*dx*rd5 )*t%zxquad &
                  + ( 15.*dy*dz2*rd7 - 3.*dy*rd5 )*t%yzquad &
                  + ( 15.*dx*dy*dz*rd7 )*t%xyquad &
                  + ( 15.*dz*dy2*rd7 - 3.*dz*rd5 )*0.5*t%quad(2) &
                  + ( 15.*dz*dx2*rd7 - 3.*dz*rd5 )*0.5*t%quad(1)
        end subroutine calc_force_coulomb_3D


        !>
        !> Calculates 2D Coulomb interaction of particle p with tree node inode
        !> that is shifted by the lattice vector vbox
        !> results are returned in eps, sumfx, sumfy, sumphi
        !> Unregularized force law is:
        !>   Phi = -2q log R
        !>   Ex = -dPhi/dx = 2 q x/R^2 etc
        !>
        subroutine calc_force_coulomb_2D(t, d, dist2, sumfx, sumfy, sumphi)
          use module_pepc_types
          use treevars
          use mpi
          implicit none

          type(t_tree_node_interaction_data), intent(in) :: t !< index of particle to interact with
          real*8, intent(in) :: d(2), dist2 !< separation vector and magnitude**2 precomputed in walk_single_particle
          real*8, intent(out) ::  sumfx,sumfy,sumphi

          real*8 :: dx,dy,d2,rd2,rd4,rd6,dx2,dy2,dx3,dy3

          sumfx  = 0.
          sumfy  = 0.
          sumphi = 0.

          !  preprocess distances and reciprocals
          dx = d(1)
          dy = d(2)

          d2  = dist2+eps2
          rd2 = 1./d2
          rd4 = rd2**2
          rd6 = rd2**3
          dx2 = dx**2
          dy2 = dy**2
          dx3 = dx**3
          dy3 = dy**3

          sumphi = sumphi - 0.5*t%charge*log(d2)    &                           !  monopole term
               !
               + (dx*t%dip(1) + dy*t%dip(2) )*rd2  &    !  dipole
               !
               + 0.5*t%quad(1)*(dx2*rd4 - rd2) + 0.5*t%quad(2)*(dy2*rd4 - rd2) + t%xyquad*dx*dy*rd4  !  quadrupole

          sumfx = sumfx + t%charge*dx*rd2  &   ! monopole
               !
               + t%dip(1)*(2*dx2*rd4 - rd2) + t%dip(2)*2*dx*dy*rd4  &  ! dipole
               !
               + 0.5*t%quad(1)*(8*dx3*rd6 - 6*dx*rd4) &                    ! quadrupole
               + 0.5*t%quad(2)*(8*dx*dy**2*rd6 - 2*dx*rd4) &
               +     t%xyquad*(8*dx2*dy*rd6 - 2*dy*rd4)

          sumfy = sumfy + t%charge*dy*rd2  &   ! monopole
               !
               + t%dip(2)*(2*dy2*rd4 - rd2) + t%dip(1)*2*dx*dy*rd4  &  ! dipole
               !
               + 0.5*t%quad(2)*(8*dy3*rd6 - 6*dy*rd4) &                    ! quadrupole
               + 0.5*t%quad(1)*(8*dy*dx**2*rd6 - 2*dy*rd4) &
               +     t%xyquad*(8*dy2*dx*rd6 - 2*dx*rd4)

        end subroutine calc_force_coulomb_2D


        !>
        !> CALC_FORCE_LJ
        !>
        !> Calculates 3D Lennard-Jones interaction of particle p with tree node inode
        !> shifted by the lattice vector vbox
        !> results are returned sumfx, sumfy, sumfz, sumphi
        !>
        subroutine calc_force_LJ(t, d, dist2, sumfx, sumfy, sumfz, sumphi)
          use module_pepc_types
          use treevars
          use mpi
          implicit none

          type(t_tree_node_interaction_data), intent(in) :: t !< index of particle to interact with
          real*8, intent(in) :: d(3), dist2 !< separation vector and magnitude**2 precomputed in walk_single_particle
          real*8, intent(out) ::  sumfx,sumfy,sumfz,sumphi
          real*8 :: dx,dy,dz,r2
          real*8 :: flj, epsc2, plj, aii2, aii2_r2, r

          sumfx  = 0.
          sumfy  = 0.
          sumfz  = 0.
          sumphi = 0.

          !  preprocess distances
          dx  = d(1)
          dy  = d(2)
          dz  = d(3)
          r2 = dist2

          !    epsc should be > a_ii to get evenly spaced ions
          aii2  = eps2
          epsc2 = 0.8*aii2
          plj   = 0.

          ! Force is repulsive up to and just beyond aii
          if (r2 > epsc2) then
              aii2_r2 = aii2/r2
          else
              aii2_r2 = aii2/epsc2
          endif

          flj = 2.*(aii2_r2)**4 - 1.*(aii2_r2  )**2

          ! potential
          sumphi = sumphi + plj

          !  forces
          r     = sqrt(r2)
          sumfx = sumfx + dx/r*flj
          sumfy = sumfy + dy/r*flj
          !       sumfz = sumfz + dz/r*flj
          sumfz=0.
      end subroutine calc_force_LJ


      !>
      !> Calculates 3D Coulomb interaction of particle p with particle inode
      !> that is shifted by the lattice vector vbox
      !> results are returned in eps, sumfx, sumfy, sumfz, sumphi
      !>
      subroutine calc_force_coulomb_3D_direct(t, d, dist2, sumfx, sumfy, sumfz, sumphi)
          use module_pepc_types
          use treevars
          use mpi
          implicit none

          type(t_tree_node_interaction_data), intent(in) :: t !< index of particle to interact with
          real*8, intent(in) :: d(3), dist2 !< separation vector and magnitude**2 precomputed in walk_single_particle
          real*8, intent(out) ::  sumfx,sumfy,sumfz,sumphi

          real*8 :: rd,dx,dy,dz,r,charge, rd3

          !  preprocess distances
          dx = d(1)
          dy = d(2)
          dz = d(3)

          r = sqrt(dist2+eps2)
          rd = 1./r
          rd3 = rd**3

          charge = t%charge

          ! potential
          sumphi = charge*rd

          !  forces

          sumfx = charge*dx*rd3

          sumfy = charge*dy*rd3

          sumfz = charge*dz*rd3
      end subroutine calc_force_coulomb_3D_direct


      !>
      !> Calculates 2D Coulomb interaction of particle p with tree node inode
      !> that is shifted by the lattice vector vbox
      !> results are returned in eps, sumfx, sumfy, sumphi
      !> Unregularized force law is:
      !>   Phi = -2q log R
      !>   Ex = -dPhi/dx = 2 q x/R^2 etc
      !>
      subroutine calc_force_coulomb_2D_direct(t, d, dist2, sumfx, sumfy, sumphi)
          use module_pepc_types
          use treevars
          use mpi
          implicit none

          type(t_tree_node_interaction_data), intent(in) :: t !< index of particle to interact with
          real*8, intent(in) :: d(2), dist2 !< separation vector and magnitude**2 precomputed in walk_single_particle
          real*8, intent(out) ::  sumfx,sumfy,sumphi

          real*8 :: dx,dy,d2,rd2,charge

          !  preprocess distances and reciprocals
          dx = d(1)
          dy = d(2)

          d2  = dist2+eps2
          rd2 = 1./d2


          charge = t%charge

          sumphi = - 0.5*charge*log(d2)

          sumfx = charge*dx*rd2

          sumfy = charge*dy*rd2
      end subroutine calc_force_coulomb_2D_direct
end module module_interaction_specific
