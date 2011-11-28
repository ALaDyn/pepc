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
        !> Force calculation wrapper.
        !> This function is thought for pre- and postprocessing of
        !> calculated fields, and for being able to call several
        !> (different) force calculation routines
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine calc_force_per_interaction(particle, res, inode, delta, dist2, vbox, cf_par)
          use module_interaction_specific
          use treevars
          implicit none

          integer, intent(in) :: inode
          type(t_particle_data), intent(in) :: particle
          type(t_particle_results), intent(inout) :: res
          real*8, intent(in) :: vbox(3), delta(3), dist2

          type(t_calc_force_params), intent(in) :: cf_par

          select case (cf_par%force_law)
            case (5)
                call update_nn_list(inode, delta, dist2, cf_par, res)
          end select

          res%work = res%work + WORKLOAD_PENALTY_INTERACTION

        end subroutine calc_force_per_interaction


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Force calculation wrapper for contributions that only have
        !> to be added once per particle
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine calc_force_per_particle(particles, nparticles, res, cf_par)
          use module_interaction_specific
          use treevars, only : me
          use module_fmm_framework
          implicit none

          integer, intent(in) :: nparticles
          type(t_particle), intent(in) :: particles(:)
          type(t_calc_force_params), intent(in) :: cf_par
          type(t_particle_results), intent(inout) :: res(:)
          real*8 :: ex_lattice, ey_lattice, ez_lattice, phi_lattice
          integer :: p

          potfarfield  = 0.
          potnearfield = 0.

        end subroutine calc_force_per_particle


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !>
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine update_nn_list(inode, d, dist2, cf_par, res)
          use treetypes
          use treevars
          use module_interaction_specific, only : t_particle_results
          implicit none
          include 'mpif.h'

          integer, intent(in) :: inode !< index of particle to interact with
          real*8, intent(in) :: d(3), dist2 !< separation vector and magnitude**2 precomputed in walk_single_particle
          type(t_calc_force_params), intent(in) :: cf_par !< Force parameters - see module_treetypes
          type(t_particle_results), intent(inout) :: res

          integer :: ierr, tmp(1)

          if (dist2 < res%maxdist2) then
            ! add node to NN_list
            res%neighbour_keys(res%maxidx) = inode !TODONN replace with key
            res%dist2(res%maxidx)          = dist2
            tmp                            = maxloc(res%dist2(:)) ! this is really ugly, but maxloc returns a 1-by-1 vector instead of the expected scalar
            res%maxidx                     = tmp(1)
            res%maxdist2                   = res%dist2(res%maxidx)
            res%work                       = res%work + WORKLOAD_PENALTY_INTERACTION
          else
!            ! node is further away than farest particle in nn-list --> this should have been avoided by the MAC
!            write(*,*) "update_nn_list(): node is further away than farest particle in nn-list --> this should have been avoided by the MAC"
!            write(*,*) inode, d, dist2
!            write(*,*) res
!            call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
          endif

        end subroutine update_nn_list
        
  end module module_calc_force
