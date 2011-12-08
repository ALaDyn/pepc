!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>  Encapsulates anything that is directly involved in force calculation
!>
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module module_calc_force
     use module_pepc_types
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

      integer, public :: force_law    = 5      !< 5=NN-list "interaction"
      integer, public :: mac_select   = 1      !< selector for multipole acceptance criterion, 1 = NN-MAC

      namelist /calc_force_nearestneighbour/ force_law, mac_select

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      public calc_force_per_interaction
      public calc_force_per_particle
      public mac
      public particleresults_clear
      public calc_force_init
      public calc_force_finalize

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
      !> initializes interaction specific parameters, redas them from file
      !> if optional argument para_file_name is given
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calc_force_init(para_file_available, para_file_name, my_rank)
        implicit none
        logical, intent(in) :: para_file_available
        character(*), intent(in) :: para_file_name
        integer, intent(in) :: my_rank
        integer, parameter :: para_file_id = 47

        if (para_file_available) then
            open(para_file_id,file=para_file_name)

            if(my_rank .eq. 0) write(*,*) "reading parameter file, section calc_force_nearestneighbour: ", para_file_name
            read(para_file_id,NML=calc_force_nearestneighbour)

            close(para_file_id)
        endif

      end subroutine

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !>
      !> finalizes the calc force module at end of simulation
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calc_force_finalize()
      end subroutine calc_force_finalize


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !>
      !> generic Multipole Acceptance Criterion
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      function mac(particle, node, dist2, boxlength2)
        implicit none

        logical :: mac
        integer, intent(in) :: node
        type(t_particle), intent(in) :: particle
        real*8, intent(in) :: dist2
        real*8, intent(in) :: boxlength2

        select case (mac_select)
            case (0)
              ! Barnes-Hut-MAC
              ! mac = (theta2 * dist2 > boxlength2)
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
      !> due to dependency on module_pepc_types and(!) on module_interaction_specific, the
      !> function cannot reside in module_interaction_specific that may not include module_pepc_types
      !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine particleresults_clear(particles, nparticles)
        use module_pepc_types
        use module_htable
        use treevars
        use module_spacefilling
        implicit none
        type(t_particle), intent(inout) :: particles(nparticles)
        integer, intent(in) :: nparticles

        integer*8 :: key
        integer :: addr, i
        real*8, dimension(:), allocatable :: boxdiag2

        allocate(boxdiag2(0:nlev))
        boxdiag2(0) = (sqrt(3.)*boxsize)**2
        do i=1,nlev
           boxdiag2(i) =  boxdiag2(i-1)/4.
        end do


        ! for each particle, we traverse the tree upwards, until the current twig
        ! contains more leaves than number of necessary neighbours - as a first guess for the
        ! search radius, we use its diameter
        do i=1,nparticles
            key = particles(i)%key

            particles(i)%results%maxdist2 = huge(0._8)
            particles(i)%results%neighbour_nodes(:) = 0
            particles(i)%results%maxidx             = 1

            do while (key .ne. 0)
              if (testaddr(key, addr)) then
                if (htable(addr)%leaves > num_neighbour_particles) then
                  ! this twig contains enough particles --> we use its diameter as search radius
                  particles(i)%results%maxdist2 = boxdiag2(level_from_key(key))
                  particles(i)%results%neighbour_nodes(:) = htable(addr)%node

                  exit ! from this loop
                endif
              endif

              key = ishft(key, -3)
              if (key.eq.0) then
               write(*,*) particles(i)
              endif
            end do

            particles(i)%results%dist2(:) = particles(i)%results%maxdist2
            particles(i)%results%dist_vector(:,:) = -13._8 
        end do


      end subroutine

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Force calculation wrapper.
        !> This function is thought for pre- and postprocessing of
        !> calculated fields, and for being able to call several
        !> (different) force calculation routines
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine calc_force_per_interaction(particle, inode, delta, dist2, vbox)
          use module_interaction_specific
          use treevars
          implicit none

          integer, intent(in) :: inode
          type(t_particle), intent(inout) :: particle
          real*8, intent(in) :: vbox(3), delta(3), dist2

          select case (force_law)
            case (5)
                call update_nn_list(particle, inode, delta, dist2)
                particle%work = particle%work + WORKLOAD_PENALTY_INTERACTION
          end select


        end subroutine calc_force_per_interaction


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Force calculation wrapper for contributions that only have
        !> to be added once per particle
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine calc_force_per_particle(particles, nparticles)
          use module_interaction_specific
          use treevars, only : me
          implicit none

          integer, intent(in) :: nparticles
          type(t_particle), intent(inout) :: particles(:)

          ! currently nothing to do here

        end subroutine calc_force_per_particle


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !>
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine update_nn_list(particle, inode, d, dist2)
          use module_pepc_types
          use treevars
          implicit none
          include 'mpif.h'

          integer, intent(in) :: inode !< index of particle to interact with
          real*8, intent(in) :: d(3), dist2 !< separation vector and magnitude**2 precomputed in walk_single_particle
          type(t_particle), intent(inout) :: particle

          integer :: ierr, tmp(1)

          if (dist2 < particle%results%maxdist2) then
            ! add node to NN_list
            particle%results%neighbour_nodes(particle%results%maxidx) = inode
            particle%results%dist2(particle%results%maxidx)           = dist2
            particle%results%dist_vector(:,particle%results%maxidx) = d
            tmp                       = maxloc(particle%results%dist2(:)) ! this is really ugly, but maxloc returns a 1-by-1 vector instead of the expected scalar
            particle%results%maxidx   = tmp(1)
            particle%results%maxdist2 = particle%results%dist2(particle%results%maxidx)
          else
            ! node is further away than farest particle in nn-list --> can be ignored
          endif

        end subroutine update_nn_list

        
  end module module_calc_force
