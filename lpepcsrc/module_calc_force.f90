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
		subroutine calc_force_per_interaction(p, inode, vbox, eps, force_const)
		  use treevars
		  implicit none

		  integer, intent(in) :: p, inode
		  real*8, intent(in) :: vbox(3)
		  real, intent(in) :: eps
		  real, intent(in) :: force_const

		  real*8 :: exc, eyc, ezc, phic

		  !  compute Coulomb fields and potential of particle p from its interaction list
		  call calc_force_coulomb(p, inode, vbox, eps, exc, eyc, ezc, phic)

		  pot_tmp(p) = pot_tmp(p)+ force_const * phic
		  ex_tmp(p)  = ex_tmp(p) + force_const * exc
		  ey_tmp(p)  = ey_tmp(p) + force_const * eyc
		  ez_tmp(p)  = ez_tmp(p) + force_const * ezc


		end subroutine calc_force_per_interaction


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Force calculation wrapper for contributions that only have
        !> to be added once per particle
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine calc_force_per_particle(eps, force_const)
          use treevars
          use module_fmm_framework
          implicit none

          real, intent(in) :: eps
          real, intent(in) :: force_const
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

			    potfarfield  = potfarfield  + phi_lattice * q(p)
			    potnearfield = potnearfield + pot_tmp(p) * q(p)
			    pot_tmp(p) = pot_tmp(p) + force_const * phi_lattice
			    ex_tmp(p)  = ex_tmp(p)  + force_const * ex_lattice
			    ey_tmp(p)  = ey_tmp(p)  + force_const * ey_lattice
			    ez_tmp(p)  = ez_tmp(p)  + force_const * ez_lattice
			  end do

          end if

        end subroutine calc_force_per_particle


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Calculates Coulomb interaction of particle p with tree node inode
        !> that is shifted by the lattice vector vbox
        !> results are returned in eps, sumfx, sumfy, sumfz, sumphi
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		subroutine calc_force_coulomb(p, inode, vbox, eps, sumfx, sumfy, sumfz, sumphi)
		  use treevars
		  implicit none

		  include 'mpif.h'

		  integer, intent(in) :: p  !< particle label
		  integer, intent(in) :: inode !< index of particle to interact with
		  real*8, intent(in) :: vbox(3) !< vector to neighbour box that is currently processed
		  real, intent(in) :: eps ! smoothing parameter
		  real*8, intent(out) ::  sumfx,sumfy,sumfz,sumphi

		  real*8 :: rd,dx,dy,dz,d,dx2,dy2,dz2
		  real*8 :: dx3,dy3,dz3,rd3,rd5,rd7,fd1,fd2,fd3,fd4,fd5,fd6
		  real :: eps2

		  eps2=eps**2
		  sumfx = 0.
		  sumfy = 0.
		  sumfz = 0.
		  sumphi = 0.

		     !  preprocess distances
		     dx = x(p) - ( xcoc(inode) + vbox(1) )
		     dy = y(p) - ( ycoc(inode) + vbox(2) )
		     dz = z(p) - ( zcoc(inode) + vbox(3) )


		     d = sqrt(dx**2+dy**2+dz**2+eps2)
		     rd = 1./d
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

		     sumphi = sumphi + charge(inode)*rd    &                           !  monopole term
		                                !
		          + (dx*xdip(inode) + dy*ydip(inode) + dz*zdip(inode))*rd3  &    !  dipole
		                                !     Dx             Dy            Dz
		          + 0.5*fd1*xxquad(inode) + 0.5*fd2*yyquad(inode) + 0.5*fd3*zzquad(inode)  &  !  quadrupole
		                                !           Qxx                 Qyy                 Qzz
		          + fd4*xyquad(inode) + fd5*yzquad(inode) + fd6*zxquad(inode)
		     !   Qxy            Qyz             Qzx

		     !  forces

		     sumfx = sumfx + charge(inode)*dx*rd3 &      ! monopole term
		                                !
		          + fd1*xdip(inode) + fd4*ydip(inode) + fd6*zdip(inode)   &   !  dipole term
		                                !
		          + (15.*dx3*rd7 - 9.*dx*rd5 )*0.5*xxquad(inode) &     !
		          + ( 15.*dy*dx2*rd7 - 3.*dy*rd5 )*xyquad(inode) &     !
		          + ( 15.*dz*dx2*rd7 - 3.*dz*rd5 )*zxquad(inode) &     !   quadrupole term
		          + ( 15*dx*dy*dz*rd7 )*yzquad(inode) &                !
		          + ( 15.*dx*dy2*rd7 - 3.*dx*rd5 )*0.5*yyquad(inode) & !
		          + ( 15.*dx*dz2*rd7 - 3.*dx*rd5 )*0.5*zzquad(inode)   !

		     sumfy = sumfy + charge(inode)*dy*rd3 &
		          + fd2*ydip(inode) + fd4*xdip(inode) + fd5*zdip(inode)  &
		          + ( 15.*dy3*rd7 - 9.*dy*rd5 )*0.5*yyquad(inode) &
		          + ( 15.*dx*dy2*rd7 - 3.*dx*rd5 )*xyquad(inode) &
		          + ( 15.*dz*dy2*rd7 - 3.*dz*rd5 )*yzquad(inode) &
		          + ( 15.*dx*dy*dz*rd7 )*zxquad(inode) &
		          + ( 15.*dy*dx2*rd7 - 3.*dy*rd5 )*0.5*xxquad(inode) &
		          + ( 15.*dy*dz2*rd7 - 3.*dy*rd5 )*0.5*zzquad(inode)

		     sumfz = sumfz + charge(inode)*dz*rd3 &
		          + fd3*zdip(inode) + fd5*ydip(inode) + fd6*xdip(inode)  &
		          + ( 15.*dz3*rd7 - 9.*dz*rd5 )*0.5*zzquad(inode) &
		          + ( 15.*dx*dz2*rd7 - 3.*dx*rd5 )*zxquad(inode) &
		          + ( 15.*dy*dz2*rd7 - 3.*dy*rd5 )*yzquad(inode) &
		          + ( 15.*dx*dy*dz*rd7 )*xyquad(inode) &
		          + ( 15.*dz*dy2*rd7 - 3.*dz*rd5 )*0.5*yyquad(inode) &
		          + ( 15.*dz*dx2*rd7 - 3.*dz*rd5 )*0.5*xxquad(inode)

		end subroutine calc_force_coulomb



end module module_calc_force
