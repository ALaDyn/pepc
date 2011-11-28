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
      function mac(cf_par, dist2, boxlength2, results)
        implicit none

        logical :: mac
        type(t_calc_force_params), intent(in) :: cf_par
        real*8, intent(in) :: dist2
        real*8, intent(in) :: boxlength2
        type(t_particle_results), intent(in) :: results

        select case (cf_par%mac)
            case (0)
              ! Barnes-Hut-MAC
              mac = (cf_par%theta2 * dist2 > boxlength2)
            case default
              ! N^2 code
              mac = .false.
        end select

      end function

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
          !> Force law struct has following content (defined in module treetypes)
          !> These need to be included/defined in call to fields from frontend
          !>    real    :: eps
          !>    real    :: force_const
          !>    integer :: force_law   0= no interaction (default); 2=2D Coulomb; 3=3D Coulomb
          type(t_calc_force_params), intent(in) :: cf_par

          real*8 :: exyz(3), phic

          select case (cf_par%force_law)
            case (2)  !  compute 2D-Coulomb fields and potential of particle p from its interaction list
                call calc_force_coulomb_2D(inode, delta(1:2), dot_product(delta(1:2), delta(1:2)), cf_par, exyz(1), exyz(2),phic)
                exyz(3) = 0.

            case (3)  !  compute 3D-Coulomb fields and potential of particle p from its interaction list
                call calc_force_coulomb_3D(inode, delta, dist2, cf_par, exyz(1), exyz(2), exyz(3), phic)

            case (4)  ! LJ potential for quiet start
                call calc_force_LJ(inode, delta, dist2, cf_par, exyz(1), exyz(2), exyz(3), phic)
                exyz(3) = 0.

            case default
              exyz = 0.
              phic = 0.
          end select

          res%e    = res%e    + cf_par%force_const * exyz
          res%pot  = res%pot  + cf_par%force_const * phic
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

          if ((do_periodic) .and. (cf_par%include_far_field_if_periodic)) then

             if ((me==0) .and. (cf_par%force_law .ne. 3)) write(*,*) "Warning: far-field lattice contribution is currently only supported for force_law==3"

             do p=1,nparticles
                call fmm_sum_lattice_force(particles(p), ex_lattice, ey_lattice, ez_lattice, phi_lattice) !TODO: use coordinates from particles

                potfarfield  = potfarfield  + phi_lattice * particles(p)%data%q
                potnearfield = potnearfield + res(p)%pot  * particles(p)%data%q

                res(p)%e     = res(p)%e     + cf_par%force_const * [ex_lattice, ey_lattice, ez_lattice]
                res(p)%pot   = res(p)%pot   + cf_par%force_const * phi_lattice
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
        subroutine calc_force_coulomb_3D(inode, d, dist2, cf_par, sumfx, sumfy, sumfz, sumphi)
          use treetypes
          use treevars
          implicit none

          include 'mpif.h'

          integer, intent(in) :: inode !< index of particle to interact with
          real*8, intent(in) :: d(3), dist2 !< separation vector and magnitude**2 precomputed in walk_single_particle
          type(t_calc_force_params), intent(in) :: cf_par !< Force parameters - see module_treetypes
          real*8, intent(out) ::  sumfx,sumfy,sumfz,sumphi

          real*8 :: rd,dx,dy,dz,r,dx2,dy2,dz2
          real*8 :: dx3,dy3,dz3,rd3,rd5,rd7,fd1,fd2,fd3,fd4,fd5,fd6
          type(t_multipole_data), pointer :: t

             t=>tree_nodes(inode)

             sumfx  = 0.
             sumfy  = 0.
             sumfz  = 0.
             sumphi = 0.

             !  preprocess distances
             dx = d(1)
             dy = d(2)
             dz = d(3)


             r = sqrt(dist2+cf_par%eps**2)
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
        subroutine calc_force_coulomb_2D(inode, d, dist2, cf_par, sumfx, sumfy, sumphi)
          use module_interaction_specific
          use treevars
          implicit none

          include 'mpif.h'

          integer, intent(in) :: inode !< index of particle to interact with
          real*8, intent(in) :: d(2), dist2 !< separation vector and magnitude**2 precomputed in walk_single_particle
          type(t_calc_force_params), intent(in) :: cf_par
          real*8, intent(out) ::  sumfx,sumfy,sumphi

          real*8 :: dx,dy,d2,rd2,rd4,rd6,dx2,dy2,dx3,dy3
          real :: eps2
          type(t_multipole_data), pointer :: t

          eps2   = cf_par%eps**2
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
        subroutine calc_force_LJ(inode, d, dist2, cf_par, sumfx, sumfy, sumfz, sumphi)
          use treetypes
          use treevars
          implicit none

          include 'mpif.h'

          integer, intent(in) :: inode !< index of particle to interact with
          real*8, intent(in) :: d(3), dist2 !< separation vector and magnitude**2 precomputed in walk_single_particle
          type(t_calc_force_params), intent(in) :: cf_par
          real*8, intent(out) ::  sumfx,sumfy,sumfz,sumphi
          real*8 :: dx,dy,dz,r
          real*8 :: flj, epsc, plj, aii

          type(t_multipole_data), pointer :: t

          t=>tree_nodes(inode)

          sumfx  = 0.
          sumfy  = 0.
          sumfz  = 0.
          sumphi = 0.

          !  preprocess distances
          dx = d(1)
          dy = d(2)
          dz = d(3)
          r = sqrt(dist2)

          !    epsc should be > a_ii to get evenly spaced ions
          aii = cf_par%eps
          epsc = 0.8*aii
          plj =0.

          ! Force is repulsive up to and just beyond aii
          if (r > epsc) then
              flj =2.*aii**8/r**8 - 1.*aii**4/r**4
          else
              flj = 2.*aii**8/epsc**8 - 1.*aii**4/epsc**4
          endif


          ! potential
          sumphi = sumphi + plj

          !  forces

          sumfx = sumfx + dx/r*flj
          sumfy = sumfy + dy/r*flj
          !    	  sumfz = sumfz + dz/r*flj
          sumfz=0.

      end subroutine calc_force_LJ
        
  end module module_calc_force
