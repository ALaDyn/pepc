!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>  Encapsulates anything that is directly involved in force calculation
!>
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module module_calc_force
     implicit none
     save
     private

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real*8, allocatable, public :: ex_tmp(:),ey_tmp(:),ez_tmp(:),pot_tmp(:),w_tmp(:)


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
	subroutine calc_force_per_interaction(p, inode, delta, dist2, vbox, cf_par)
          use treetypes
	  use treevars
	  implicit none

	  integer, intent(in) :: p, inode
	  real*8, intent(in) :: vbox(3), delta(3), dist2
	  !> Force law struct has following content (defined in module treetypes) 
	  !> These need to be included/defined in call to fields from frontend 
	  !>    real    :: eps 
	  !>    real    :: force_const 
	  !>    integer :: force_law   0= no interaction (default); 2=2D Coulomb; 3=3D Coulomb
          type(calc_force_params), intent(in) :: cf_par

	  real*8 :: exc, eyc, ezc, phic

          select case (cf_par%force_law)
            case (2)  !  compute 2D-Coulomb fields and potential of particle p from its interaction list
                call calc_force_coulomb_2D(p, inode, vbox, cf_par, exc, eyc, phic)
                ezc = 0.

            case (3)  !  compute 3D-Coulomb fields and potential of particle p from its interaction list
                call calc_force_coulomb_3D(inode, delta, dist2, cf_par, exc, eyc, ezc, phic)

		    case default
		      exc  = 0.
		      eyc  = 0.
		      ezc  = 0.
		      phic = 0.
		  end select

		  pot_tmp(p) = pot_tmp(p)+ cf_par%force_const * phic
		  ex_tmp(p)  = ex_tmp(p) + cf_par%force_const * exc
		  ey_tmp(p)  = ey_tmp(p) + cf_par%force_const * eyc
		  ez_tmp(p)  = ez_tmp(p) + cf_par%force_const * ezc


		end subroutine calc_force_per_interaction


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Force calculation wrapper for contributions that only have
        !> to be added once per particle
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine calc_force_per_particle(cf_par)
          use treetypes
          use treevars
          use module_fmm_framework
          implicit none

          type(calc_force_params), intent(in) :: cf_par
          real*8 :: ex_lattice, ey_lattice, ez_lattice, phi_lattice
          integer :: p

          potfarfield  = 0.
          potnearfield = 0.

          if (do_periodic) then

             do p=1,npp
                call fmm_sum_lattice_force(p, ex_lattice, ey_lattice, ez_lattice, phi_lattice)
                
                ! write(*,*) p,          "|", q(p), "|", pot_tmp(p), "|", phi_lattice
                ! write(*,*) "            |", x(p), "|", ex_tmp(p),  "|", ex_lattice
                ! write(*,*) "            |", y(p), "|", ey_tmp(p),  "|", ey_lattice
                ! write(*,*) "            |", z(p), "|", ez_tmp(p),  "|", ez_lattice

                potfarfield  = potfarfield  + phi_lattice * particles(p)%q
                potnearfield = potnearfield + pot_tmp(p) * particles(p)%q
                pot_tmp(p) = pot_tmp(p) + cf_par%force_const * phi_lattice
                ex_tmp(p)  = ex_tmp(p)  + cf_par%force_const * ex_lattice
                ey_tmp(p)  = ey_tmp(p)  + cf_par%force_const * ey_lattice
                ez_tmp(p)  = ez_tmp(p)  + cf_par%force_const * ez_lattice
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
          real*8, intent(in) :: d(3), dist2
          type(calc_force_params), intent(in) :: cf_par
          real*8, intent(out) ::  sumfx,sumfy,sumfz,sumphi

          real*8 :: rd,dx,dy,dz,r,dx2,dy2,dz2
          real*8 :: dx3,dy3,dz3,rd3,rd5,rd7,fd1,fd2,fd3,fd4,fd5,fd6
          type(multipole), pointer :: t

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
                  + (dx*t%xdip + dy*t%ydip + dz*t%zdip)*rd3  &    !  dipole
                                        !     Dx             Dy            Dz
                  + 0.5*fd1*t%xxquad + 0.5*fd2*t%yyquad + 0.5*fd3*t%zzquad  &  !  quadrupole
                                        !           Qxx                 Qyy                 Qzz
                  + fd4*t%xyquad + fd5*t%yzquad + fd6*t%zxquad
             !   Qxy            Qyz             Qzx

             !  forces

             sumfx = sumfx + t%charge*dx*rd3 &      ! monopole term
                                        !
                  + fd1*t%xdip + fd4*t%ydip + fd6*t%zdip   &   !  dipole term
                                        !
                  + (15.*dx3*rd7 - 9.*dx*rd5 )*0.5*t%xxquad &     !
                  + ( 15.*dy*dx2*rd7 - 3.*dy*rd5 )*t%xyquad &     !
                  + ( 15.*dz*dx2*rd7 - 3.*dz*rd5 )*t%zxquad &     !   quadrupole term
                  + ( 15*dx*dy*dz*rd7 )*t%yzquad &                !
                  + ( 15.*dx*dy2*rd7 - 3.*dx*rd5 )*0.5*t%yyquad & !
                  + ( 15.*dx*dz2*rd7 - 3.*dx*rd5 )*0.5*t%zzquad   !

             sumfy = sumfy + t%charge*dy*rd3 &
                  + fd2*t%ydip + fd4*t%xdip + fd5*t%zdip  &
                  + ( 15.*dy3*rd7 - 9.*dy*rd5 )*0.5*t%yyquad &
                  + ( 15.*dx*dy2*rd7 - 3.*dx*rd5 )*t%xyquad &
                  + ( 15.*dz*dy2*rd7 - 3.*dz*rd5 )*t%yzquad &
                  + ( 15.*dx*dy*dz*rd7 )*t%zxquad &
                  + ( 15.*dy*dx2*rd7 - 3.*dy*rd5 )*0.5*t%xxquad &
                  + ( 15.*dy*dz2*rd7 - 3.*dy*rd5 )*0.5*t%zzquad

             sumfz = sumfz + t%charge*dz*rd3 &
                  + fd3*t%zdip + fd5*t%ydip + fd6*t%xdip  &
                  + ( 15.*dz3*rd7 - 9.*dz*rd5 )*0.5*t%zzquad &
                  + ( 15.*dx*dz2*rd7 - 3.*dx*rd5 )*t%zxquad &
                  + ( 15.*dy*dz2*rd7 - 3.*dy*rd5 )*t%yzquad &
                  + ( 15.*dx*dy*dz*rd7 )*t%xyquad &
                  + ( 15.*dz*dy2*rd7 - 3.*dz*rd5 )*0.5*t%yyquad &
                  + ( 15.*dz*dx2*rd7 - 3.*dz*rd5 )*0.5*t%xxquad

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
        subroutine calc_force_coulomb_2D(p, inode, vbox, cf_par, sumfx, sumfy, sumphi)
          use treetypes
          use treevars
          implicit none

          include 'mpif.h'

          integer, intent(in) :: p  !< particle label
          integer, intent(in) :: inode !< index of particle to interact with
          real*8, intent(in) :: vbox(3) !< vector to neighbour box that is currently processed
          type(calc_force_params), intent(in) :: cf_par
          real*8, intent(out) ::  sumfx,sumfy,sumphi

          real*8 :: dx,dy,d2,rd2,rd4,rd6,dx2,dy2,dx3,dy3
          real :: eps2
          type(multipole), pointer :: t

          eps2   = cf_par%eps**2
          sumfx  = 0.
          sumfy  = 0.
          sumphi = 0.

          t=>tree_nodes(inode)

          !  preprocess distances and reciprocals
          dx = particles(p)%x(1) - ( t%xcoc + vbox(1) )
          dy = particles(p)%x(2) - ( t%ycoc + vbox(2) )

          d2  = dx**2+dy**2+eps2 
          rd2 = 1./d2 
          rd4 = rd2**2 
          rd6 = rd2**3 
          dx2 = dx**2 
          dy2 = dy**2 
          dx3 = dx**3 
          dy3 = dy**3 
	  
          sumphi = sumphi - 0.5*t%charge*log(d2)    &                           !  monopole term 
               ! 
               + (dx*t%xdip + dy*t%ydip )*rd2  &    !  dipole 
               !                               
               + 0.5*t%xxquad*(dx2*rd4 - rd2) + 0.5*t%yyquad*(dy2*rd4 - rd2) + t%xyquad*dx*dy*rd4  !  quadrupole 
          
          sumfx = sumfx + t%charge*dx*rd2  &   ! monopole 
               ! 
               + t%xdip*(2*dx2*rd4 - rd2) + t%ydip*2*dx*dy*rd4  &  ! dipole 
               ! 
               + 0.5*t%xxquad*(8*dx3*rd6 - 6*dx*rd4) &                    ! quadrupole 
               + 0.5*t%yyquad*(8*dx*dy**2*rd6 - 2*dx*rd4) & 
               +     t%xyquad*(8*dx2*dy*rd6 - 2*dy*rd4) 
          
          sumfy = sumfy + t%charge*dy*rd2  &   ! monopole 
               ! 
               + t%ydip*(2*dy2*rd4 - rd2) + t%xdip*2*dx*dy*rd4  &  ! dipole 
               ! 
               + 0.5*t%yyquad*(8*dy3*rd6 - 6*dy*rd4) &                    ! quadrupole 
               + 0.5*t%xxquad*(8*dy*dx**2*rd6 - 2*dy*rd4) & 
               +     t%xyquad*(8*dy2*dx*rd6 - 2*dx*rd4) 

        end subroutine calc_force_coulomb_2D
        
        
        
end module module_calc_force
