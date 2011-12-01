!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>  Encapsulates anything that is directly involved in force calculation
!>
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module module_calc_force
     use treetypes
     implicit none
     save
     private

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real*8, parameter :: WORKLOAD_PENALTY_MAC  = 1._8 !< TODO: currently unused
      real*8, parameter :: WORKLOAD_PENALTY_INTERACTION = 30._8


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      public calc_force_per_interaction
      public calc_force_per_particle
      public mac
      public particleresults_clear

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

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !>
      !> generic Multipole Acceptance Criterion
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      function mac(particle, node, cf_par, dist2, boxlength2)
        implicit none

        logical :: mac
        integer, intent(in) :: node
        type(t_particle), intent(in) :: particle
        type(t_calc_force_params), intent(in) :: cf_par
        real*8, intent(in) :: dist2
        real*8, intent(in) :: boxlength2

        select case (cf_par%mac)
            case (0)
              ! Barnes-Hut-MAC
              mac = (cf_par%theta2 * dist2 > boxlength2)
            case (1)
              ! NN-MAC: we may "interact" with the node if it is further away than maxdist2 --> this leads to the node *not* being put onto the NN-list (strange, i know)
              ! first line: original formulation, last line: after transition to formulation with only one square root
              ! mac = sqrt(dist2) - sqrt(3.* boxlength2)  >  sqrt(results%maxdist2)                                ! + sqrt(3.*boxlength2)
              !     = sqrt(dist2)                         >  sqrt(results%maxdist2) + sqrt(3.*boxlength2)          ! ^2
              !     =      dist2                          > (sqrt(results%maxdist2) + sqrt(3.*boxlength2))**2
              !     =      dist2                          > results%maxdist2 + 2.*sqrt( 3.*results%maxdist2*boxlength2) + 3.*boxlength2
                mac =      dist2                          > particle%results%maxdist2 +    sqrt(12.*particle%results%maxdist2*boxlength2) + 3.*boxlength2
              ! TODO NN: this estimation should be evaluated without (!!) any square roots for performance reasons (which does not seem to be trivial)
            case default
              ! N^2 code
              mac = .false.
        end select

      end function

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !>
      !> clears result in t_particle datatype - usually, this function does not need to be touched
      !> due to dependency on treetypes and(!) on module_interaction_specific, the
      !> function cannot reside in module_interaction_specific that may not include treetypes
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      elemental subroutine particleresults_clear(particle)
        implicit none
        type(t_particle), intent(inout) :: particle

        call results_clear(particle%results)

      end subroutine

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Force calculation wrapper.
        !> This function is thought for pre- and postprocessing of
        !> calculated fields, and for being able to call several
        !> (different) force calculation routines
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine calc_force_per_interaction(particle, inode, delta, dist2, vbox, cf_par)
          use module_interaction_specific
          use treevars
          implicit none

          integer, intent(in) :: inode
          type(t_particle), intent(inout) :: particle
          real*8, intent(in) :: vbox(3), delta(3), dist2

          type(t_calc_force_params), intent(in) :: cf_par

          select case (cf_par%force_law)
            case (5)
                call update_nn_list(particle, inode, delta, dist2, cf_par)
                particle%work = particle%work + WORKLOAD_PENALTY_INTERACTION
          end select


        end subroutine calc_force_per_interaction


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Force calculation wrapper for contributions that only have
        !> to be added once per particle
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine calc_force_per_particle(particles, nparticles, cf_par)
          use module_interaction_specific
          use treevars, only : me
          implicit none

          integer, intent(in) :: nparticles
          type(t_particle), intent(inout) :: particles(:)
          type(t_calc_force_params), intent(in) :: cf_par

          ! currently nothing to do here

        end subroutine calc_force_per_particle


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !>
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine update_nn_list(particle, inode, d, dist2, cf_par)
          use treetypes
          use treevars
          implicit none
          include 'mpif.h'

          integer, intent(in) :: inode !< index of particle to interact with
          real*8, intent(in) :: d(3), dist2 !< separation vector and magnitude**2 precomputed in walk_single_particle
          type(t_calc_force_params), intent(in) :: cf_par !< Force parameters - see module_treetypes
          type(t_particle), intent(inout) :: particle

          integer :: ierr, tmp(1)

          if (dist2 < particle%results%maxdist2) then
            ! add node to NN_list
            particle%results%neighbour_nodes(particle%results%maxidx) = inode
            particle%results%dist2(particle%results%maxidx)           = dist2
            tmp                       = maxloc(particle%results%dist2(:)) ! this is really ugly, but maxloc returns a 1-by-1 vector instead of the expected scalar
            particle%results%maxidx   = tmp(1)
            particle%results%maxdist2 = particle%results%dist2(particle%results%maxidx)
          else
            ! node is further away than farest particle in nn-list --> can be ignored
          endif

        end subroutine update_nn_list

        
  end module module_calc_force
