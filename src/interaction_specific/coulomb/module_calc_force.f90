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

      integer, public :: force_law    = 3      !< 3 = 3D-Coulomb, 2 = 2D-Coulomb
      integer, public :: mac_select   = 0      !< selector for multipole acceptance criterion, mac_select==0: Barnes-Hut
      logical, public :: include_far_field_if_periodic = .true. !< if set to false, the far-field contribution to periodic boundaries is ignored (aka 'minimum-image-mode')
      real*8, public  :: theta2       = 0.6**2.  !< square of multipole opening angle
      real*8, public  :: eps2         = 0.0    !< square of short-distance cutoff parameter for plummer potential (0.0 corresponds to classical Coulomb)

      namelist /calc_force_coulomb/ force_law, mac_select, include_far_field_if_periodic, theta2, eps2


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! currently, all public functions in module_calc_force are obligatory
      public calc_force_per_interaction
      public calc_force_per_particle
      public mac
      public particleresults_clear
      public calc_force_init

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

            if(my_rank .eq. 0) write(*,*) "reading parameter file, section calc_force_coulomb: ", para_file_name
            read(para_file_id,NML=calc_force_coulomb)

            close(para_file_id)
        endif

      end subroutine


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
              mac = (theta2 * dist2 > boxlength2)
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
      subroutine particleresults_clear(particles, nparticles)
        use treetypes
        implicit none
        type(t_particle), intent(inout) :: particles(nparticles)
        integer, intent(in) :: nparticles

        particles(1:nparticles)%results = EMPTY_PARTICLE_RESULTS

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


          real*8 :: exyz(3), phic

          select case (force_law)
            case (2)  !  compute 2D-Coulomb fields and potential of particle p from its interaction list
                call calc_force_coulomb_2D(inode, delta(1:2), dot_product(delta(1:2), delta(1:2)), exyz(1), exyz(2),phic)
                exyz(3) = 0.

            case (3)  !  compute 3D-Coulomb fields and potential of particle p from its interaction list
                call calc_force_coulomb_3D(inode, delta, dist2, exyz(1), exyz(2), exyz(3), phic)

            case (4)  ! LJ potential for quiet start
                call calc_force_LJ(inode, delta, dist2, exyz(1), exyz(2), exyz(3), phic)
                exyz(3) = 0.

            case default
              exyz = 0.
              phic = 0.
          end select

          particle%results%e         = particle%results%e    + exyz
          particle%results%pot       = particle%results%pot  + phic
          particle%work              = particle%work         + WORKLOAD_PENALTY_INTERACTION

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
          use module_fmm_framework
          use module_mirror_boxes
          implicit none

          integer, intent(in) :: nparticles
          type(t_particle), intent(inout) :: particles(:)
          real*8 :: ex_lattice, ey_lattice, ez_lattice, phi_lattice
          integer :: p

          ! calculate spherical multipole expansion of central box
          if (include_far_field_if_periodic) call fmm_framework_timestep(particles)

          potfarfield  = 0.
          potnearfield = 0.

          if ((do_periodic) .and. (include_far_field_if_periodic)) then

             if ((me==0) .and. (force_law .ne. 3)) write(*,*) "Warning: far-field lattice contribution is currently only supported for force_law==3"

             do p=1,nparticles
                call fmm_sum_lattice_force(particles(p), ex_lattice, ey_lattice, ez_lattice, phi_lattice) !TODO: use coordinates from particles

                potfarfield  = potfarfield  + phi_lattice * particles(p)%data%q
                potnearfield = potnearfield + particles(p)%results%pot  * particles(p)%data%q

                particles(p)%results%e     = particles(p)%results%e     + [ex_lattice, ey_lattice, ez_lattice]
                particles(p)%results%pot   = particles(p)%results%pot   +  phi_lattice
             end do

          end if

        end subroutine calc_force_per_particle


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Calculates 3D Coulomb interaction of particle p with tree node inode
        !> that is shifted by the lattice vector vbox
        !> results are returned in eps, sumfx, sumfy, sumfz, sumphi
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine calc_force_coulomb_3D(inode, d, dist2, sumfx, sumfy, sumfz, sumphi)
          use treetypes
          use treevars
          implicit none

          include 'mpif.h'

          integer, intent(in) :: inode !< index of particle to interact with
          real*8, intent(in) :: d(3), dist2 !< separation vector and magnitude**2 precomputed in walk_single_particle
          real*8, intent(out) ::  sumfx,sumfy,sumfz,sumphi

          real*8 :: rd,dx,dy,dz,r,dx2,dy2,dz2
          real*8 :: dx3,dy3,dz3,rd3,rd5,rd7,fd1,fd2,fd3,fd4,fd5,fd6
          type(t_tree_node_interaction_data), pointer :: t

             t=>tree_nodes(inode)

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


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Calculates 2D Coulomb interaction of particle p with tree node inode
        !> that is shifted by the lattice vector vbox
        !> results are returned in eps, sumfx, sumfy, sumphi 
        !> Unregularized force law is: 
        !>   Phi = -2q log R 
        !>   Ex = -dPhi/dx = 2 q x/R^2 etc 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine calc_force_coulomb_2D(inode, d, dist2, sumfx, sumfy, sumphi)
          use module_interaction_specific
          use treevars
          implicit none

          include 'mpif.h'

          integer, intent(in) :: inode !< index of particle to interact with
          real*8, intent(in) :: d(2), dist2 !< separation vector and magnitude**2 precomputed in walk_single_particle
          real*8, intent(out) ::  sumfx,sumfy,sumphi

          real*8 :: dx,dy,d2,rd2,rd4,rd6,dx2,dy2,dx3,dy3
          type(t_tree_node_interaction_data), pointer :: t

          sumfx  = 0.
          sumfy  = 0.
          sumphi = 0.

          t=>tree_nodes(inode)

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
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> CALC_FORCE_LJ
        !>
        !> Calculates 3D Lennard-Jones interaction of particle p with tree node inode
        !> shifted by the lattice vector vbox
        !> results are returned sumfx, sumfy, sumfz, sumphi
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine calc_force_LJ(inode, d, dist2, sumfx, sumfy, sumfz, sumphi)
          use treetypes
          use treevars
          implicit none

          include 'mpif.h'

          integer, intent(in) :: inode !< index of particle to interact with
          real*8, intent(in) :: d(3), dist2 !< separation vector and magnitude**2 precomputed in walk_single_particle
          real*8, intent(out) ::  sumfx,sumfy,sumfz,sumphi
          real*8 :: dx,dy,dz,r2
          real*8 :: flj, epsc2, plj, aii2, aii2_r2, r

          type(t_tree_node_interaction_data), pointer :: t

          t=>tree_nodes(inode)

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
          !    	  sumfz = sumfz + dz/r*flj
          sumfz=0.

      end subroutine calc_force_LJ


  end module module_calc_force
