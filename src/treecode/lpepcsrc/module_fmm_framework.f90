!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!> Encapsulates calculation of the lattice contribution by means
!> of the FMM-approach to the lattice
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_fmm_framework
      use module_debug
      implicit none
      save
      private

      !! TODO: set dipole- and low-order stuff in MLattice to zero
      !! TODO: add formal dipole correction as surface term



      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! lattice basis vectors
      real*8, public :: t_lattice_1(3) = [1, 0, 0] !< 1st vector of lattice basis
      real*8, public :: t_lattice_2(3) = [0, 1, 0] !< 2nd vector of lattice basis
      real*8, public :: t_lattice_3(3) = [0, 0, 1] !< 3rd vector of lattice basis
      real*8, public :: Lattice(3,3) !< holds the lattice transformation matrix
      real*8, public :: LatticeInv(3,3) !< holds the inverse lattice transformation matrix
      real*8, public :: LatticeCenter(3) !< holds the central point of the lattice box
      logical, public :: periodicity(3) = [.false., .false., .false.]  !< boolean switches for determining periodicity directions
      !> variables that should not be written to
      integer, public :: num_neighbour_boxes = 1
      integer, dimension(:,:), allocatable, public :: neighbour_boxes !dimensions in 3D case (3,27)
      !> is set to OR(periodicity), to be able to exit from all procedures as fast as possible
      !> if nothing is to be done
      logical, public :: do_periodic
      !> far- and near-field contribution to potential energy (has to be calculated in fields.p90)
      real*8, public :: potfarfield, potnearfield
      !> whether to do dipole correction or not, see [J.Chem.Phys. 107, 10131, eq. (19,20)]
      logical, public :: do_extrinsic_correction = .false.
      !> set to .true. to activate debug mode (table output etc.)
      logical, public :: periodic_debug =.false.

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      public fmm_framework_init
      public fmm_framework_finalize
      public fmm_framework_timestep
      public fmm_sum_lattice_force
      public lattice_vect
      public fmm_framework_param_dump
      public constrain_periodic

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! general stuff
      integer :: myrank
      ! FMM-PARAMETERS
      integer, parameter :: Lmax = 8
      integer, parameter :: LmaxL = 4
      integer, parameter :: MaxIter = 16
      integer :: ws = 1
      ! FMM-VARIABLES
      complex*16 :: mu_cent(1:Lmax*(Lmax+1)/2+Lmax+1)
      ! externally calculated FMM variables
      complex*16 :: omega_tilde(1:Lmax*(Lmax+1)/2+Lmax+1)
      complex*16 :: MLattice(1:Lmax*(Lmax+1)/2+Lmax+1)
      !> variables for extrinsic to intrinsic correction
      real*8 :: box_dipole(3) = 0.
      real*8 :: quad_trace    = 0.

      ! Variables for lattice coefficient calculation
      complex*16 :: Mstar(1:Lmax*(Lmax+1)/2+Lmax+1)
      complex*16 :: Lstar(1:Lmax*(Lmax+1)/2+Lmax+1)
      ! number of boxes to include into each direction
      integer :: periodicity_switches(3)
      ! true for lattices with axes parallel to cartesian system
      logical :: simplelattice

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  subroutine-implementation  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      contains

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Module Initialization, should be called on program startup
        !> after setting up all geometric parameters etc.
        !> @param[in] mpi_rank MPI rank of process for controlling debug output
        !> @param[in] wellsep well-separation criterion parmater
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine fmm_framework_init(mpi_rank, wellsep)
          use module_math_tools
          implicit none
          integer, intent(in) :: mpi_rank
          integer, intent(in) :: wellsep

          myrank = mpi_rank
          ws     = wellsep

          call init_movement_constraint

          LatticeCenter = 0.5*(t_lattice_1 + t_lattice_2 + t_lattice_3)

          do_periodic = periodicity(1) .or. periodicity(2) .or. periodicity(3)

          call calc_neighbour_boxes

           ! anything above has to be done in any case
          if (do_periodic) then
            call calc_lattice_coefficients
          end if

        end subroutine fmm_framework_init


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Clean up the allocated memory
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine fmm_framework_finalize
          implicit none
          deallocate(neighbour_boxes)

        end subroutine fmm_framework_finalize


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Prepares the list of neighbour boxes within ws
        !> stores their number in num_neighbour_boxes and their logical
        !> indices/coordinates in neighbour_boxes
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine calc_neighbour_boxes
        implicit none
        integer :: i,j,k,idx

          do i = 1,3
            if (periodicity(i)) then
               periodicity_switches(i) = ws
             else
               periodicity_switches(i) = 0
             end if
          end do

          idx = 0

          num_neighbour_boxes = product(2*periodicity_switches+1)
          allocate(neighbour_boxes(3,num_neighbour_boxes))

          ! the boxes of the shells that surround the central box
          ! are ordered in a way, that a box and its central-symmetric counterpart are grouped together
          do i = -periodicity_switches(1),-1
            do j = -periodicity_switches(2),periodicity_switches(2)
              do k = -periodicity_switches(3),periodicity_switches(3)
                  idx = idx + 1
                  neighbour_boxes(:,idx) = [ i,  j ,  k]
                  idx = idx + 1
                  neighbour_boxes(:,idx) = [-i, -j , -k]
              end do
            end do
          end do

          i = 0
          do j = -periodicity_switches(2),-1
            do k = -periodicity_switches(3),periodicity_switches(3)
                idx = idx + 1
                neighbour_boxes(:,idx) = [ i,  j ,  k]
                idx = idx + 1
                neighbour_boxes(:,idx) = [-i, -j , -k]
            end do
          end do

          i = 0
          j = 0
          do k = -periodicity_switches(3),-1
            idx = idx + 1
            neighbour_boxes(:,idx) = [ i,  j ,  k]
            idx = idx + 1
            neighbour_boxes(:,idx) = [-i, -j , -k]
          end do

          neighbour_boxes(:,idx+1) = [0, 0, 0] ! center box is put to the back of the boxlist for easier iteration


        end subroutine calc_neighbour_boxes

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Initializes transformation matrices between cartesian system and lattice basis
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine init_movement_constraint
          use module_math_tools
          implicit none

          integer :: i,j

          Lattice(1,:) = t_lattice_1
          Lattice(2,:) = t_lattice_2
          Lattice(3,:) = t_lattice_3
          LatticeInv = Inverse3(Lattice)

          ! simplify the movement constraint if the lattice is really simple
          simplelattice = .true.
          do i = 1,3
            do j = 1,3
              if (i.ne.j) then
                simplelattice = simplelattice .and. (Lattice(i,j) == 0)
              else
                simplelattice = simplelattice .and. (Lattice(i,j) >  0)
              end if
            end do
          end do


        end subroutine init_movement_constraint


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Refreshes Multipole information and Taylor coefficients,
        !> has to be called every timestep
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine fmm_framework_timestep(particles)
          use treetypes
          implicit none
          type(t_particle), dimension(:), intent(in) :: particles
          if (.not. do_periodic) return
          call calc_omega_tilde(particles)
          call calc_mu_cent
          call calc_extrinsic_correction(particles)
        end subroutine fmm_framework_timestep


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Calculates the lattice coefficients for computing mu_cent
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine calc_lattice_coefficients
          use module_debug
          implicit none
          include 'mpif.h'

          integer :: l, m
          integer :: iter

          write (debug, *) 'LATTICE COEFFICIENTS: Starting calculation'
          call print_debug(.true., 2, [6, 15], myrank.eq.0 .and. [.true., .true.])

          Mstar    = 0
          Lstar    = 0
          MLattice = 0

          ! pretabulation of necessary values
          do l = 0,Lmax
            do m = 0,l
              Mstar( tblinv(l, m) ) = MstarFunc(l, m)
              Lstar( tblinv(l, m) ) = LstarFunc(l, m)
            end do
          end do

          if ((myrank == 0) .and. periodic_debug) then
            call WriteTableToFile('Mstar.tab', Mstar)
            call WriteTableToFile('Lstar.tab', Lstar)
          end if

          ! zeroth step of iteration
          MLattice = Lstar

          do iter = 1,MaxIter

            MLattice = M2M( UL( MLattice ) , MStar ) + Lstar

            !DEBUG
            !write(*,*) "----------- After Iteration ", iter
            !write(*,*) "MLattice = ", MLattice
          end do

          ! MLattice(1:tblinv(3,3))=0

          if ((myrank == 0) .and. periodic_debug) then
            call WriteTableToFile('MLattice.tab', MLattice)
          end if

        end subroutine calc_lattice_coefficients


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Calculates the overall multipole expansion of the whole
        !> central box
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine calc_omega_tilde(particles)
          use treevars
          use treetypes
          use module_debug
          implicit none
          include 'mpif.h'

          type(t_particle), dimension(:), intent(in) :: particles

          integer :: ll, mm, p
          integer :: ierr
          complex*16 :: tmp

          omega_tilde = 0

          ! calculate multipole contributions of all local particles
          do ll=0,LmaxL
            do mm=0,ll
              tmp = 0

              do p=1,npp
                tmp = tmp + particles(p)%data%q * O(ll, mm, particles(p)%x - LatticeCenter)
              end do

              omega_tilde( tblinv(ll, mm) ) = tmp
            end do
          end do

          ! sum multipole contributions from all processors
          call MPI_ALLREDUCE(MPI_IN_PLACE, omega_tilde, Lmax*(Lmax+1)/2+Lmax+1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)

          if ((myrank == 0) .and. periodic_debug) then
            call WriteTableToFile('omega_tilde.tab', omega_tilde)
          end if

          if (real(omega_tilde( tblinv(0, 0))) > 1.E-14) then

            write (debug, *) 'WARNIG: The central box is not charge-neutral. Switching off calculation of lattice contribution. omega_tilde( tblinv(0, 0))=', omega_tilde( tblinv(0, 0))
            call print_debug(.true., 2, [ipefile, 6], [.true., myrank.eq.0])

            do_periodic         = .false.
            num_neighbour_boxes = 1
          end if

        end subroutine calc_omega_tilde


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Calculates the (cartesian) overall dipole moment
        !> \f$\frac{4\pi}{3}\sum_p q(p){\vec r}_p\f$ and the
        !> trace of the quadrupole matrix
        !> \f$\frac{2\pi}{3}\sum_p q(p){\vec r}_p\cdot{\vec r}_p\f$
        !> for performing the extrinsic-to-intrinsic correction
        !> (see [J.Chem.Phys. 107, 10131, eqn.(19,20)] for details
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine calc_extrinsic_correction(particles)
          use treevars
          use module_debug
          use treetypes
          implicit none
          include 'mpif.h'

          type(t_particle), intent(in), dimension(:) :: particles(:)

          real, parameter :: pi=3.141592654
          real*8 :: box_dipole(3) = 0.
          real*8 :: quad_trace    = 0.
          real*8 :: r(3)

          integer :: p
          integer :: ierr

          box_dipole = 0.
          quad_trace = 0.

          if (do_extrinsic_correction) then

	          ! calculate multipole contributions of all local particles
	          do p=1,npp
	            r = particles(p)%x - LatticeCenter

	            box_dipole = box_dipole + particles(p)%data%q * r
	            quad_trace = quad_trace + particles(p)%data%q * dot_product(r, r)
	          end do

	          ! sum contributions from all processors
	          call MPI_ALLREDUCE(MPI_IN_PLACE, box_dipole, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
	          call MPI_ALLREDUCE(MPI_IN_PLACE, quad_trace, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

	          box_dipole = 4.*pi/3. * box_dipole
	          quad_trace = 2.*pi/3. * quad_trace

	      end if

        end subroutine calc_extrinsic_correction


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Calculates the lattice contribution with respect to the
        !> centre of the original box
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine calc_mu_cent
          implicit none

          integer :: j, k, lambda, nu
          complex*16 :: tmp

          mu_cent = 0;

          ! contribution of all outer lattice cells, with regards to the centre of the original box
          do j = 0,LmaxL
            do k = 0,j

              tmp = 0

              do lambda = 0,LmaxL
                do nu = -lambda,lambda
                  tmp = tmp + (-1.)**lambda*tbl(MLattice,j+lambda,k+nu) * tbl(omega_tilde, lambda, nu)
                end do
              end do

              mu_cent( tblinv(j, k) ) = tmp
            end do
          end do

          if ((myrank == 0) .and. periodic_debug) then
            call WriteTableToFile('mu_cent.tab', mu_cent)
          end if

        end subroutine calc_mu_cent


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Calculates the chargeless moments of the multipole expansion
        !> for a certain particle
        !> this is formally identical to #Mvec, however #O also treats negative m
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		complex*16 function O(l, m, r)
		  use module_math_tools
          implicit none
          integer, intent(in) :: l, m
          real*8, intent(in) :: r(3)
          integer :: mm

          real*8 :: s(3)
          complex*16, parameter :: i = (0,1)

          if ((l<0) .or. (m<-l) .or. (m>l)) then
              O = 0 ! these are zero per definition
          else
	          call cartesian_to_spherical(r, s)
	          mm = abs(m)

	          O = div_by_fac(LegendreP(l, mm, s(2)) * Exp(-i*mm*s(3) ), l+mm)

	          if (s(1) > 0) O = O * s(1)**l
	          if (m    < 0) O = (-1)**mm * conjg(O)
          endif

		end function O


		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!>
		!> Calculates force onto individual particles that results
		!> from mirror boxes beyond the near field region,
		!> i.e. the lattice contribution
		!>
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		subroutine fmm_sum_lattice_force(particle, ex_lattice, ey_lattice, ez_lattice, phi_lattice)
		  use treevars
		  implicit none

		  type(t_particle), intent(in) :: particle
		  real*8, intent(out) ::  ex_lattice, ey_lattice, ez_lattice, phi_lattice

		  complex*16 :: mu_shift(0:1,0:1), tmp
		  integer :: ll, mm, j, k
		  real*8 :: r(3)

		  if (.not. do_periodic) then
            ex_lattice = 0
            ey_lattice = 0
            ez_lattice = 0
            phi_lattice = 0
		    return
		  end if

		 ! shift mu_cent to the position of our particle
		  mu_shift =  0
		  r        = particle%x - LatticeCenter

		  do ll = 0,1
		    do mm = 0,ll

		      tmp = 0

		      do j = ll,LmaxL
		        do k=-j,j
		           tmp = tmp + O(j-ll, k-mm, r) * tbl( mu_cent, j, k )
		        end do
		      end do

		      mu_shift(ll, mm) = tmp

		    end do
		  end do

                         ! E = -grad(Phi)
          if (do_extrinsic_correction) then    ! extrinsic correction
		    ex_lattice  = -real(mu_shift(1,1))  + box_dipole(1)
		    ey_lattice  = -aimag(mu_shift(1,1)) + box_dipole(2)
            ez_lattice  = -real(mu_shift(1,0))  + box_dipole(3)
		    phi_lattice =  real(mu_shift(0,0))  - dot_product(r, box_dipole) + quad_trace
		  else
            ex_lattice  = -real(mu_shift(1,1))
            ey_lattice  = -aimag(mu_shift(1,1))
            ez_lattice  = -real(mu_shift(1,0))
            phi_lattice =  real(mu_shift(0,0))
          endif
		end subroutine fmm_sum_lattice_force



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Calculates vector with respect to lattice base vectors
        !> @param[in] ijk lattice indices
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function lattice_vect(ijk)
          implicit none
          integer, intent(in) :: ijk(3)
          real*8 :: lattice_vect(3)

          lattice_vect = ijk(1)*t_lattice_1 + ijk(2)*t_lattice_2 + ijk(3)*t_lattice_3

        end function lattice_vect


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Calculate Spherical Multipole coefficients \f$\mathcal{M}\f$ for a cartesian vector
        !> @param[in] l multipole order
        !> @param[in] m
        !> @param[in] r cartesian vector [x, y, z]
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        complex*16 function Mvec(l,m,r)
          use module_math_tools
          implicit none
          integer, intent(in) :: l, m
          real*8, intent(in) :: r(3)

          real*8 :: s(3)
          complex*16, parameter :: i = (0,1)

          call cartesian_to_spherical(r, s)

          Mvec = div_by_fac(LegendreP(l, m, s(2)) * Exp(-i*m*s(3) ), l+m)

          if (s(1) > 0) Mvec = Mvec * s(1)**l

          ! DEBUG (for comparison with Mathematica worksheet)
          !write(*, '("Mvec[", I2.2, ", ", I2.2, ", {", F8.5,",",F8.5,",",F8.5, "}],  ", D20.10, D20.10)') &
          !                l, m, r(1), r(2), r(3), real(Mvec), aimag(Mvec)
          !write(*, '("M[", I2.2, ", ", I2.2, ", ", F8.5,",",F8.5,",",F8.5, "],  ", D20.10, D20.10)') &
          !                l, m, s(1), s(2), s(3), real(Mvec), aimag(Mvec)

        end function Mvec


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Calculate Spherical Taylor coefficients \f$\mathcal{L}\f$ for a cartesian vector
        !> @param[in] l multipole order
        !> @param[in] m
        !> @param[in] r cartesian vector [x, y, z]
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        complex*16 function Lvec(l,m,r)
          use module_math_tools
          implicit none
          integer, intent(in) :: l, m
          real*8, intent(in) :: r(3)
          real*8 :: s(3)
          complex*16, parameter :: i = (0,1)

          call cartesian_to_spherical(r, s)

          Lvec = mult_by_fac(LegendreP(l, m, s(2)) / s(1)**(l+1) * Exp( i*m*s(3) ), l-m)

        end function Lvec


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Table access function for giving arbitrary l and m
        !>
        !> The function cares for picking the right indices and
        !> respects symmetry
        !>
        !> @param[in] l multipole order
        !> @param[in] m
        !> @param[in] A table
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        complex*16 function tbl(A,l,m)
          implicit none
          integer, intent(in) :: l, m
          complex*16, intent(in) :: A(1:Lmax*(Lmax+1)/2+Lmax+1)

          if (l>Lmax) then
            tbl = 0
          else
            tbl = A( tblinv(l, abs(m)) )

            if (m<0) tbl = (-1)**m * conjg(tbl)
          end if

        end function tbl


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Table access function for giving arbitrary l and m
        !>
        !> Calculates the flat index, where an entry for l and m
        !> has to be stored in a one-dimensional array
        !>
        !> @param[in] l
        !> @param[in] m
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        integer function tblinv(l,m)
          implicit none
          integer, intent(in) :: l, m

          tblinv = l*(l+1)/2 + 1 + m

        end function tblinv


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> M2M-Operator (denoted with \f$\otimes\f$ )
        !>
        !> @param[in] l multipole order
        !> @param[in] m
        !> @param[in] A table
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function M2M(L, M)
          implicit none
          complex*16, intent(in) :: L(1:Lmax*(Lmax+1)/2+Lmax+1)
          complex*16, intent(in) :: M(1:Lmax*(Lmax+1)/2+Lmax+1)
          complex*16, dimension(1:Lmax*(Lmax+1)/2+Lmax+1) :: M2M

          complex*16 :: t
          integer :: ll, mm, j, k

          ! DEBUG
          !write(*,*) "L = ", L
          !write(*,*) "M = ", M

          do ll = 0,Lmax
            do mm = 0,ll

              t = 0

              do j = 0,Lmax
                do k = -j,j
                  t = t + tbl(L,ll+j, mm+k) * tbl(M,j,k)
                end do
              end do

              M2M( tblinv(ll, mm) ) = t
            end do
          end do

          ! DEBUG
          !write(*,*) "M2M = ", M2M

        end function M2M


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Scaling Operator \f$\mathcal{U}_L\f$ for Taylor coefficients
        !> @param[in] L table with Taylor coefficients
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function UL(L)
          implicit none
          complex*16, intent(in) :: L(1:Lmax*(Lmax+1)/2+Lmax+1)
          complex*16, dimension(1:Lmax*(Lmax+1)/2+Lmax+1) :: UL

          integer :: ll,mm

          ! DEBUG
          !write(*,*) "L = ", L

          do ll = 0,Lmax
            do mm = 0,ll
              UL( tblinv(ll, mm) ) = tbl(L,ll,mm) / (2*ws+1)**(ll+1)
            end do
          end do

          ! DEBUG
          !write(*,*) "UL(L) = ", UL

        end function UL


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Formal summation of \f$M\f$ over NF, ie all (27) neighbouring boxes
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        complex*16 function MstarFunc(l,m)
          implicit none
          integer, intent(in) :: l, m

          integer :: i

          MstarFunc = 0

          do i = 1,num_neighbour_boxes
             MstarFunc = MstarFunc + Mvec(l, m, -lattice_vect(neighbour_boxes(:,i)) )
          end do

        end function MstarFunc


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Formal summation of \f$M\f$ over FF`, ie a lot of boxes
        !> with some overhead to avoid numerical elimination of small values
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        complex*16 function LstarFunc(l,m)
          use module_math_tools
          implicit none
          integer, intent(in) :: l, m
          real*8 :: rpart(num_neighbour_boxes*(num_neighbour_boxes-1)), ipart(num_neighbour_boxes*(num_neighbour_boxes-1)),rp,ip
          complex*16 :: tmp
          complex*16, parameter :: ic = (0,1)
          integer :: i, ii, k

          LstarFunc = 0
          k = 0

          do i = 1,num_neighbour_boxes-1 ! central box is being omitted in this loop
            do ii = 1,num_neighbour_boxes
              tmp = Lvec(l, m, -lattice_vect(neighbour_boxes(:,ii)) + (2*ws+1)*lattice_vect(neighbour_boxes(:,i)))
              k   = k+1
              ! we store the summands and order them before performing the sum
              ! to avoid numeric elimination
              rpart(k) =  real(tmp)
              ipart(k) = aimag(tmp)
            end do
          end do

          call sort_abs(rpart)
          call sort_abs(ipart)

          rp = 0.
          ip = 0.

          do i = k,1,-1
            rp = rp + rpart(i)
            ip = ip + ipart(i)
          end do

          LStarFunc = rp + ic*ip

        end function LstarFunc


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Converts cartesian coordinates to spherical system
        !> @param[in]  cartesian  cartesian vector [x, y, z]
        !> @param[out] spherical  spherical coordinates [ |r|, Cos(Theta), phi ]
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine cartesian_to_spherical(cartesian, spherical)
          implicit none
          real*8, intent(in)  :: cartesian(3)
          real*8, intent(out) :: spherical(3)

          spherical(1) = sqrt(dot_product(cartesian, cartesian))

          if (spherical(1) == 0) then
            spherical(2) = 1
          else
            spherical(2) = cartesian(3) / spherical(1)
          end if

          if ((cartesian(1) == 0) .and. (cartesian(2) == 0)) then
            spherical(3) = 0
          else
            spherical(3) = atan2(cartesian(2), cartesian(1))
          end if

        end subroutine cartesian_to_spherical


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Writes contents of table T to a output stream s in a structured way
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine PrintTable(s, T)
          implicit none
          complex*16, intent(in) :: T(1:Lmax*(Lmax+1)/2+Lmax+1)
          integer, intent(in) :: s

          integer :: ll,mm,idx

          idx = 0

          do ll = 0,Lmax
            do mm = 0,ll
              idx = idx + 1

              write(s,'(I6, I6, I6, D50.35, D50.35)') idx, ll, mm, tbl(T, ll, mm)
            end do
          end do


        end subroutine PrintTable


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Writes contents of table T to file s
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine WriteTableToFile(s, T)
          implicit none
          complex*16, intent(in) :: T(1:Lmax*(Lmax+1)/2+Lmax+1)
          character(len=*) :: s
          integer, parameter :: temp_file = 60

          open(temp_file, file=trim(s))

          write(temp_file,*) "# t_lattice_1 = ", t_lattice_1
          write(temp_file,*) "# t_lattice_2 = ", t_lattice_2
          write(temp_file,*) "# t_lattice_3 = ", t_lattice_3
          write(temp_file,*) "# periodicity_switches = ", periodicity_switches

          write(temp_file,*) "# Lmax    = ", Lmax
          write(temp_file,*) "# MaxIter = ", MaxIter
          write(temp_file,*) "##########################"
          write(temp_file,*) "  idx     l     m           real-part                                         imaginary-part"

          call PrintTable(temp_file, T)

          close(temp_file)

        end subroutine WriteTableToFile


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Dumps all parameters to the stream ifile
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine fmm_framework_param_dump(ifile)
          use treevars
          implicit none
          integer, intent(in) :: ifile

          if (me .ne. 0) return

          write(ifile,*) "LATTICE: ------------- Lattice fmm-framework switches ----------------"
          write(ifile,*) "LATTICE: Lmax = ", Lmax
          write(ifile,*) "LATTICE: LmaxL = ", LmaxL
          write(ifile,*) "LATTICE: MaxIter = ", MaxIter
          write(ifile,*) "LATTICE: ws = ", ws
          write(ifile,*) "LATTICE: t_lattice_1 = ", t_lattice_1
          write(ifile,*) "LATTICE: t_lattice_2 = ", t_lattice_2
          write(ifile,*) "LATTICE: t_lattice_3 = ", t_lattice_3
          write(ifile,*) "LATTICE: periodicity = ", periodicity
          write(ifile,*) "LATTICE: # neighbours = ", num_neighbour_boxes
          write (ifile,*)
        end subroutine fmm_framework_param_dump



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Shifts all particles that left the original box
        !> back into it
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    subroutine constrain_periodic(x, y, z, np_local)
	        use module_math_tools
	        use module_debug

	        implicit none

            real*8, intent(inout), dimension(1:np_local) :: x,y,z
            integer, intent(in) :: np_local

            integer :: p !< loop variable
	        real*8 :: lattice_coord(3), real_coord(3)

	        if (simplelattice) then

                do p = 1,np_local

                  if (periodicity(1)) x(p) = modulo(x(p), t_lattice_1(1))
                  if (periodicity(2)) y(p) = modulo(y(p), t_lattice_2(2))
                  if (periodicity(3)) z(p) = modulo(z(p), t_lattice_3(3))

                end do

 	        else
	          ! the lattice axes are not parallel to the outer cartesian axes
		        do p = 1,np_local

		            lattice_coord = matmul([x(p), y(p), z(p)], LatticeInv)

		            if (any(lattice_coord > 1.D+0) .or. any(lattice_coord < 0.D+0)) then

		              lattice_coord = modulo(lattice_coord, 1.D+0)

		              real_coord = matmul(lattice_coord, Lattice)

		              if (periodicity(1)) x(p) = real_coord(1)
		              if (periodicity(2)) y(p) = real_coord(2)
		              if (periodicity(3)) z(p) = real_coord(3)

		            end if
		        end do
		    end if

	    end subroutine


end module module_fmm_framework





