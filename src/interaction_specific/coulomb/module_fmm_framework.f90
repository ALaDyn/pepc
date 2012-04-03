! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2012 Juelich Supercomputing Centre, 
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!> Encapsulates calculation of the lattice contribution by means
!> of the FMM-approach to the lattice
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_fmm_framework
      use module_debug
      use module_mirror_boxes
      implicit none
      include 'mpif.h'
      private

      !! TODO: set dipole- and low-order stuff in MLattice to zero

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !> far- and near-field contribution to potential energy (has to be calculated in fields.p90)
      real*8, public :: potfarfield, potnearfield
      !> whether to do dipole correction or not, see [J.Chem.Phys. 107, 10131, eq. (19,20)]
      logical, public :: do_extrinsic_correction = .false.

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      public fmm_framework_init
      public fmm_framework_timestep
      public fmm_sum_lattice_force
      public lattice_vect
      public fmm_framework_param_dump

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! general stuff
      integer :: myrank
      integer :: MPI_COMM_fmm
      ! precision flags
      integer, parameter :: kfp                  = 8 ! numeric precision (kind value)
      integer, parameter :: MPI_REAL_fmm         = MPI_REAL8
      logical, parameter  :: chop_arrays         = .false.
      real(kfp), parameter :: prec = 1.E-16_kfp
      ! shortcut notations
      real(kfp), parameter :: zero = 0._kfp
      real(kfp), parameter :: one  = 1._kfp
      real(kfp), parameter :: two  = 2._kfp
      ! FMM-PARAMETERS
      integer, parameter :: Lmax_multipole = 50
      integer, parameter :: Lmax_taylor    = Lmax_multipole * 2
      integer, parameter :: MaxIter        = 16
      integer :: ws = 1
      logical, parameter :: use_pretabulated_lattice_coefficients = .false.
      ! FMM-VARIABLES
      integer, parameter :: fmm_array_length_multipole = Lmax_multipole*(Lmax_multipole+1)/2+Lmax_multipole+1
      integer, parameter :: fmm_array_length_taylor    = Lmax_taylor   *(Lmax_taylor   +1)/2+Lmax_taylor   +1
      ! internally calculated FMM variables
      complex(kfp) :: mu_cent(1:fmm_array_length_taylor)
      complex(kfp) :: omega_tilde(1:fmm_array_length_multipole)
      complex(kfp) :: MLattice(1:fmm_array_length_taylor)
      !> variables for extrinsic to intrinsic correction
      real(kfp) :: box_dipole(3) = zero
      real(kfp) :: quad_trace    = zero

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
        !>  - well-separation criterion is automatically set to module_mirror_boxes::mirror_box_layers
        !>  - module_mirror_boxes::calc_neighbour_boxes() must have been called before calling this routine
        !>
        !> @param[in] mpi_rank MPI rank of process for controlling debug output
        !> @param[in] mpi_comm MPI communicator to be used
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine fmm_framework_init(mpi_rank, mpi_comm)
          use module_debug
          use module_mirror_boxes, only : mirror_box_layers
          implicit none
          integer, intent(in) :: mpi_rank
          integer, intent(in) :: mpi_comm

          myrank       = mpi_rank
          MPI_COMM_fmm = mpi_comm
          ws           = mirror_box_layers

          LatticeCenter = 0.5*(t_lattice_1 + t_lattice_2 + t_lattice_3)

          do_periodic = any(periodicity(1:3))

           ! anything above has to be done in any case
          if (do_periodic) then
            MLattice = 0

            if (use_pretabulated_lattice_coefficients) then
              call load_lattice_coefficients(MLattice)
              DEBUG_WARNING('(a)', "Using pretabulated lattice coefficients from [J.Chem. Phys. 107, 10131]. These are only valid for unit-box simulations regions.")
            else
              call calc_lattice_coefficients(MLattice)
            endif

            if ((myrank == 0) .and. dbg(DBG_PERIODIC)) then
              call WriteTableToFile('MLattice.tab', MLattice, Lmax_taylor)
            end if
          end if

        end subroutine fmm_framework_init


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Refreshes Multipole information and Taylor coefficients,
        !> has to be called every timestep with particles that were used in tree buildup
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine fmm_framework_timestep(particles, nparticles)
          use module_pepc_types
          implicit none
          integer, intent(in) :: nparticles
          type(t_particle), dimension(nparticles), intent(in) :: particles

          if (do_periodic) then
            call calc_omega_tilde(particles, nparticles)
            call calc_mu_cent(omega_tilde, mu_cent)
            call calc_extrinsic_correction(particles, nparticles)
          endif
        end subroutine fmm_framework_timestep


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Calculates the lattice coefficients for computing mu_cent
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine calc_lattice_coefficients(ML)
          use module_debug
          implicit none

          complex(kfp), intent(out) :: ML(1:fmm_array_length_taylor)
          ! Variables for lattice coefficient calculation
          complex(kfp) :: Mstar(1:fmm_array_length_multipole)
          complex(kfp) :: Lstar(1:fmm_array_length_taylor)

          integer :: l, m, iter
          character(20) :: fn

          call pepc_status('LATTICE COEFFICIENTS: Starting calculation')

          Mstar = 0
          Lstar = 0

          ! pretabulation of necessary values
          do l = 0,Lmax_multipole
            do m = 0,l
              Mstar( tblinv(l, m, Lmax_multipole) ) = MstarFunc(l, m)
            end do
          end do

          do l = 0,Lmax_taylor
            do m = 0,l
              Lstar( tblinv(l, m, Lmax_taylor) )    = LstarFunc(l, m)
            end do
          end do

          call chop(Mstar)
          call chop(Lstar)

          if ((myrank == 0) .and. dbg(DBG_PERIODIC)) then
            call WriteTableToFile('Mstar.tab', Mstar, Lmax_multipole)
            call WriteTableToFile('Lstar.tab', Lstar, Lmax_taylor)
          end if

          ! zeroth step of iteration
          ML = Lstar

          do iter = 1,MaxIter

            ML = M2L( UL( ML ) , MStar ) + Lstar

            if ((myrank == 0) .and. dbg(DBG_PERIODIC)) then
              write(fn,'("MLattice.", I2.2, ".tab")') iter
              call WriteTableToFile(trim(fn), ML, Lmax_taylor)
            endif
          end do

          call chop(ML)

          ! ML(1:tblinv(3,3))=0
          call pepc_status('LATTICE COEFFICIENTS: finished calculation')

        end subroutine calc_lattice_coefficients


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Sets the lattice coefficients for computing mu_cent
        !> data taken from [Challacombe, White, Head-Gordon: J. Chem. Phys. 107, 10131]
        !> only works if Lmax = 20
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine load_lattice_coefficients(M)
          use module_debug
          implicit none

          complex(kfp), intent(inout) :: M(1:fmm_array_length_taylor)

          call pepc_status('LATTICE COEFFICIENTS: Loading')

          if (Lmax_taylor >= 4) then
            M(tblinv( 4,  0, Lmax_taylor)) =   2.8119304871888668206200481879919E+00_kfp
            M(tblinv( 4,  4, Lmax_taylor)) =   1.4059652435944334103100240939959E+01_kfp
          endif

          if (Lmax_taylor >= 6) then
            M(tblinv( 6,  0, Lmax_taylor)) =   5.4795908739321644069327702900830E-01_kfp
            M(tblinv( 6,  4, Lmax_taylor)) =  -3.8357136117525150848529392030580E+00_kfp
          endif

          if (Lmax_taylor >= 8) then
            M(tblinv( 8,  0, Lmax_taylor)) =   1.2156157302097918942115948482762E+02_kfp
            M(tblinv( 8,  4, Lmax_taylor)) =   1.2156157302097918942115948482762E+02_kfp
            M(tblinv( 8,  8, Lmax_taylor)) =   7.9015022463636473123753665137950E+03_kfp
          endif

          if (Lmax_taylor >= 10) then
            M(tblinv(10,  0, Lmax_taylor)) =   3.1179916736109123107822587274280E+02_kfp
            M(tblinv(10,  4, Lmax_taylor)) =  -6.8595816819440070837209692003420E+02_kfp
            M(tblinv(10,  8, Lmax_taylor)) =  -1.1661288859304812042325647640581E+04_kfp
          endif

          if (Lmax_taylor >= 12) then
            M(tblinv(12,  0, Lmax_taylor)) =   2.4245612747359092217199640516493E+05_kfp
            M(tblinv(12,  4, Lmax_taylor)) =   2.0375858264140266510162557633886E+05_kfp
            M(tblinv(12,  8, Lmax_taylor)) =   7.0682666545985000701644635107790E+05_kfp
            M(tblinv(12, 12, Lmax_taylor)) =   2.3702435984527078287639617913271E+08_kfp
          endif

          if (Lmax_taylor >= 14) then
            M(tblinv(14,  0, Lmax_taylor)) =   2.0954087119885542648713979286402E+06_kfp
            M(tblinv(14,  4, Lmax_taylor)) =  -2.6940969154138554834060830511089E+06_kfp
            M(tblinv(14,  8, Lmax_taylor)) =  -1.7062613797621084728238525990356E+07_kfp
            M(tblinv(14, 12, Lmax_taylor)) =  -6.5406686224214158124914349629700E+08_kfp
          endif

          if (Lmax_taylor >= 16) then
            M(tblinv(16,  0, Lmax_taylor)) =   5.4279858299650169624382885076980E+08_kfp
            M(tblinv(16,  4, Lmax_taylor)) =   2.2841041529105412872383486949334E+08_kfp
            M(tblinv(16,  8, Lmax_taylor)) =   1.2973301854895758582918144058332E+09_kfp
            M(tblinv(16, 12, Lmax_taylor)) =   2.5882484900055575638355343741650E+10_kfp
            M(tblinv(16, 16, Lmax_taylor)) =   6.9973653547984205656745320248030E+12_kfp
          endif

          if (Lmax_taylor >= 18) then
            M(tblinv(18,  0, Lmax_taylor)) =   1.4686049951258450810632984509870E+10_kfp
            M(tblinv(18,  4, Lmax_taylor)) =  -1.5376346487994712576061405817891E+10_kfp
            M(tblinv(18,  8, Lmax_taylor)) =  -2.4226921558569614141580552133902E+10_kfp
            M(tblinv(18, 12, Lmax_taylor)) =  -1.0692416604738659056152567751127E+12_kfp
            M(tblinv(18, 16, Lmax_taylor)) =  -3.9585194668444324327226917843820E+13_kfp
          endif

          if (Lmax_taylor >= 20) then
            M(tblinv(20,  0, Lmax_taylor)) =   2.9414124910043233182340700935067E+12_kfp
            M(tblinv(20,  4, Lmax_taylor)) =   5.0799363324581667612792095666680E+11_kfp
            M(tblinv(20,  8, Lmax_taylor)) =   4.3319375525806128280090124574150E+12_kfp
            M(tblinv(20, 12, Lmax_taylor)) =   4.7785845726839660008475961329560E+13_kfp
            M(tblinv(20, 16, Lmax_taylor)) =   1.6103883836731949966180674427718E+15_kfp
            M(tblinv(20, 20, Lmax_taylor)) =   5.0103139602723200237451484155740E+17_kfp
          endif

          if (Lmax_taylor >= 21) then
            DEBUG_WARNING(*, 'load_lattice_coefficients(): Lmax_taylor >= 21')
          endif

          call pepc_status('LATTICE COEFFICIENTS: finished')

        end subroutine load_lattice_coefficients


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Calculates the overall multipole expansion of the whole
        !> central box
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine calc_omega_tilde(particles, nparticles)
          use module_pepc_types
          use module_debug
          implicit none

          type(t_particle), dimension(nparticles), intent(in) :: particles
          integer, intent(in) :: nparticles

          integer :: ll, mm, p
          integer :: ierr
          real(kfp) :: S(3)
          real*8 :: R(3)

          omega_tilde = 0

          ! calculate multipole contributions of all local particles
          do p=1,nparticles
            R   = particles(p)%x - LatticeCenter
            S   = cartesian_to_spherical(R)

            do ll=0,Lmax_multipole
              do mm=0,ll
                omega_tilde( tblinv(ll, mm, Lmax_multipole) ) = omega_tilde( tblinv(ll, mm, Lmax_multipole) ) + omega(ll, mm, S, particles(p)%data%q)
              end do
            end do

          end do

          call chop(omega_tilde)

          ! sum multipole contributions from all processors - treat complex as two real numbers since e.g. complex*32 is not supported by mpi
          call MPI_ALLREDUCE(MPI_IN_PLACE, omega_tilde, 2*fmm_array_length_multipole, MPI_REAL_fmm, MPI_SUM, MPI_COMM_fmm, ierr)

          if ((myrank == 0) .and. dbg(DBG_PERIODIC)) then
            call WriteTableToFile('omega_tilde.tab', omega_tilde, Lmax_multipole)
          end if

          if (real(omega_tilde( tblinv(0, 0, Lmax_multipole))) > prec) then
            DEBUG_WARNING(*, 'The central box is not charge-neutral: Q_total=omega_tilde( tblinv(0, 0))=', omega_tilde( tblinv(0, 0, Lmax_multipole)), ' Ignoring, but this could lead to infinite energies etc.' )
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
        subroutine calc_extrinsic_correction(particles, nparticles)
          use module_debug
          use module_pepc_types
          implicit none

          type(t_particle), intent(in), dimension(nparticles) :: particles(:)
          integer, intent(in) :: nparticles

          real(kfp), parameter :: pi=acos(-one)
          real*8 :: r(3)

          integer :: p
          integer :: ierr

          box_dipole = 0.
          quad_trace = 0.

          if (do_extrinsic_correction) then

              ! calculate multipole contributions of all local particles
              do p=1,nparticles
                r = particles(p)%x - LatticeCenter

                box_dipole = box_dipole + particles(p)%data%q * r
                quad_trace = quad_trace + particles(p)%data%q * dot_product(r, r)
              end do

              ! sum contributions from all processors
              call MPI_ALLREDUCE(MPI_IN_PLACE, box_dipole, 3, MPI_REAL_fmm, MPI_SUM, MPI_COMM_fmm, ierr)
              call MPI_ALLREDUCE(MPI_IN_PLACE, quad_trace, 1, MPI_REAL_fmm, MPI_SUM, MPI_COMM_fmm, ierr)

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
        subroutine calc_mu_cent(omega, mu)
          implicit none
          complex(kfp), intent(in) :: omega(1:fmm_array_length_multipole)
          complex(kfp), intent(out) :: mu(1:fmm_array_length_taylor)

          ! contribution of all outer lattice cells, with regards to the centre of the original box
          mu = M2L(MLattice, omega)

          call chop(mu)

          if ((myrank == 0) .and. dbg(DBG_PERIODIC)) then
            call WriteTableToFile('mu_cent.tab', mu_cent, Lmax_taylor)
          end if

        end subroutine calc_mu_cent


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Calculates P^~ as given by [Challacombe et al, eq. (7)]
        !> This is consistent with [Kudin, eq(2)],
        !> but contradicts with [White, eq. (2a)]
        !>
        !> uses negative order relation as given by
        !>  http://mathworld.wolfram.com/AssociatedLegendrePolynomial.html
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        recursive real(kfp) function Ptilda(l, m, x) result(Pt)
          implicit none
          integer, intent(in) :: l, m
          real(kfp), intent(in) :: x

          if (m >= 0) then
            Pt = LegendreP(l, m, x) / factorial(l + m)
          else
            Pt = (-one)**m * Ptilda(l, abs(m), x)
          endif
        end function Ptilda

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Calculates P^~ as given by [Challocombe et al, eq. (8)]
        !> This is (almost) consistent with [Kudin, eq(3)],
        !> but contradicts with [White, eq. (2b)]
        !>
        !> uses negative order relation as given by
        !>  http://mathworld.wolfram.com/AssociatedLegendrePolynomial.html
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        recursive real(kfp) function P2tilda(l, m, x) result(P2t)
          implicit none
          integer, intent(in) :: l, m
          real(kfp), intent(in) :: x

          if (m >= 0) then
            P2t = LegendreP(l, m, x) * factorial(l - m)
          else
            P2t = (-one)**abs(m) * P2tilda(l, abs(m), x)
          endif

        end function P2tilda

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Calculates the chargeless moments of the multipole expansion
        !> for a certain particle (position given in spherical coordinates
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        complex(kfp) function OMultipole(l, m, s)
          implicit none
          integer, intent(in) :: l, m
          real(kfp), intent(in) :: s(3)

          complex(kfp), parameter :: i = (zero,one)

          if ((l<0) .or. (m<-l) .or. (m>l)) then
              OMultipole = 0 ! these are zero per definition
          else
              if (s(1) > 0 .or. l>0) then
                OMultipole = s(1)**real(l, kind=kfp) * &
                             Ptilda(l, m, s(2)) * exp(-i*real(m,kind=kfp)*s(3) )
              else
                ! avoid having to compute 0^0 = 1 since some runtimes do not like that
                OMultipole = Ptilda(l, m, s(2)) * exp(-i*real(m,kind=kfp)*s(3) )
              endif

          endif

        end function OMultipole

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Calculates the chargeless moments of the Taylor expansion
        !> for a certain particle (coordinates given in spherical system)
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        complex(kfp) function MTaylor(l, m, s)
          implicit none
          integer, intent(in) :: l, m
          real(kfp), intent(in) :: s(3)

          complex(kfp), parameter :: i = (zero,one)

          if ((l<0) .or. (m<-l) .or. (m>l)) then
              MTaylor = 0 ! these are zero per definition
          else
              MTaylor = one/(s(1)**real(l+1,kind=kfp)) * &
                        P2tilda(l, m, s(2)) * exp( i*real(m,kind=kfp)*s(3) )
          endif

        end function MTaylor

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Calculates the charged moments of the multipole expansion
        !> for a certain particle (position in spherical coordinates)
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        complex(kfp) function omega(l, m, s, q)
          implicit none
          integer, intent(in) :: l, m
          real(kfp), intent(in) :: s(3)
          real*8, intent(in) :: q

          omega = real(q, kind=kfp) * OMultipole(l, m, s)

        end function omega


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Calculates the charged moments of the Taylor expansion
        !> for a certain particle (coordinates in spherical system)
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        complex(kfp) function mu(l, m, s, q)
          implicit none
          integer, intent(in) :: l, m
          real(kfp), intent(in) :: s(3)
          real*8, intent(in) :: q

          mu = real(q, kind=kfp) * MTaylor(l, m, s)

        end function mu


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Calculates force at individual position that results
        !> from mirror boxes beyond the near field region,
        !> i.e. the lattice contribution
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine fmm_sum_lattice_force(pos, e_lattice, phi_lattice)
          implicit none

          real*8, intent(in) :: pos(3)
          real*8, intent(out) ::  e_lattice(3), phi_lattice

          complex(kfp) :: mu_shift(1:fmm_array_length_taylor), O_R(1:fmm_array_length_multipole)
          integer :: l, m
          real*8 :: R(3)
          real(kfp) :: S(3)

          if (.not. do_periodic) then
              e_lattice   = 0
              phi_lattice = 0
            else
              ! shift mu_cent to the position of our particle
              R        = pos - LatticeCenter
              S        = cartesian_to_spherical(R)

              do l = 0,Lmax_multipole
                do m = 0,l
                  O_R(tblinv(l, m, Lmax_multipole)) = OMultipole(l, m, S)
                end do
              end do

              mu_shift = L2L(O_R, mu_cent, 1)

              ! E = -grad(Phi)
              e_lattice  = -[  real(tbl(mu_shift,1,1, Lmax_taylor)), aimag(tbl(mu_shift,1,1, Lmax_taylor)), real(tbl(mu_shift,1,0, Lmax_taylor)) ]
              phi_lattice =    real(tbl(mu_shift,0,0, Lmax_taylor))

              if (do_extrinsic_correction) then    ! extrinsic correction
                e_lattice   = e_lattice   + box_dipole
                phi_lattice = phi_lattice - dot_product(R, box_dipole) + quad_trace
              endif

          end if

        end subroutine fmm_sum_lattice_force


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
        complex(kfp) function tbl(A, l, m, Lmax)
          implicit none
          integer, intent(in) :: l, m, Lmax
          complex(kfp), intent(in) :: A(*)

          if ((l<0)) then
            DEBUG_ERROR('("tbl(A,l,m) - invalid arguments. l=", I0, " m=", I0)', l, m)
          endif

          if ((l>Lmax) .or. (abs(m)>l)) then
            tbl = 0
          else
            tbl = A( tblinv(l, abs(m), Lmax) )

            if (m<0) tbl = (-one)**m * conjg(tbl)
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
        integer function tblinv(l,m,Lmax)
          use module_debug
          implicit none
          integer, intent(in) :: l, m,Lmax

          if ((l<0) .or. (m<0) .or. (m>l) .or. (l>Lmax)) then
            DEBUG_ERROR('("tblinv(l,m) - invalid arguments. l=", I0, " m=", I0)', l, m)
          endif

          tblinv = l*(l+1)/2 + 1 + m

        end function tblinv


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> M2M-Operator (denoted with \f$\triangleleft\f$ )
        !> eq. (6) in [Kudin]
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function M2M(O_a, O_b)
          implicit none
          complex(kfp), intent(in) :: O_a(1:fmm_array_length_multipole)
          complex(kfp), intent(in) :: O_b(1:fmm_array_length_multipole)
          complex(kfp), dimension(1:fmm_array_length_multipole) :: M2M

          complex(kfp) :: t
          integer :: l, m, j, k

          do l = 0,Lmax_multipole
            do m = 0,l

              t = zero

              do j = 0,l ! Attention, this sum only runs to l instead of infty|Lmax
                do k = -j,j
                  t = t + tbl(O_b, l-j, m-k, Lmax_multipole) * tbl(O_a, j, k, Lmax_multipole)
                end do
              end do

              M2M( tblinv(l, m, Lmax_multipole) ) = t
            end do
          end do

          call chop(M2M)

        end function M2M


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> M2L-Operator (denoted with \f$\otimes\f$ )
        !> eq. (7) in [Kudin]
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function M2L(M_b, O_a)
          implicit none
          complex(kfp), intent(in) :: M_b(1:fmm_array_length_taylor)
          complex(kfp), intent(in) :: O_a(1:fmm_array_length_multipole)
          complex(kfp), dimension(1:fmm_array_length_taylor) :: M2L

          complex(kfp) :: t
          integer :: l, m, j, k

          do l = 0,Lmax_taylor
            do m = 0,l

              t = 0_kfp

              do j = 0,Lmax_taylor ! should be infinity, but see last page of [Kudin]: there they state, that (p-l) is enough
                do k = -j,j
                  t = t + (-1)**j * tbl(M_b, j+l, k+m, Lmax_taylor) * tbl(O_a, j, k, Lmax_multipole) !TODO: evtl. (-1)**j hinzufuegen (vgl. S. 82 bei Ivo) ??
                end do
              end do

              M2L( tblinv(l, m, Lmax_taylor) ) = t
            end do
          end do

          call chop(M2L)

        end function M2L


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> L2L-Operator (denoted with \f$\triangleright\f$ )
        !> eq. (8) in [Kudin]
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function L2L(O_b, M_r, max_l)
          implicit none
          complex(kfp), intent(in) :: O_b(1:fmm_array_length_multipole)
          complex(kfp), intent(in) :: M_r(1:fmm_array_length_taylor)
          complex(kfp), dimension(1:fmm_array_length_taylor) :: L2L
          integer, intent(in), optional :: max_l

          complex(kfp) :: t
          integer :: l, m, j, k
          integer :: maxl

          if (present(max_l)) then
            maxl = max_l
          else
            maxl = Lmax_taylor
          endif

          do l = 0,maxl
            do m = 0,l

              t = 0_kfp

              do j = l,Lmax_multipole ! Attention, this sum starts at l since O_b(j-l<0,..) = 0 anyway
                do k = -j,j
                  t = t + tbl(O_b, j-l, k-m, Lmax_multipole) * tbl(M_r, j, k, Lmax_taylor)
                end do
              end do

              L2L( tblinv(l, m, Lmax_taylor) ) = t
            end do
          end do

          call chop(L2L)

        end function L2L


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Scaling Operator \f$\mathcal{U}_L\f$ for Taylor coefficients
        !> @param[in] L table with Taylor coefficients
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function UL(L)
          implicit none
          complex(kfp), intent(in) :: L(1:fmm_array_length_taylor)
          complex(kfp), dimension(1:fmm_array_length_taylor) :: UL
          real(kfp) :: scalefac

          integer :: ll, mm

          scalefac = real(2*ws+1, kind=kfp)

          do ll = 0,Lmax_taylor
            do mm = 0,ll
              UL( tblinv(ll, mm, Lmax_taylor) ) = L( tblinv(ll, mm, Lmax_taylor) ) / scalefac**real(ll+1, kind=kfp)
            end do
          end do

          call chop(UL)

        end function UL


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Formal summation of \f$M\f$ over NF, ie all (27) neighbouring boxes
        !> with some overhead to avoid numerical elimination of small values
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        complex(kfp) function MstarFunc(l,m)
          implicit none
          integer, intent(in) :: l, m
          real(kfp) :: rpart(num_neighbour_boxes), ipart(num_neighbour_boxes),rp,ip, s(3)
          complex(kfp) :: tmp
          complex(kfp), parameter :: ic = (zero,one)

          integer :: i

          do i = 1,num_neighbour_boxes
             s = cartesian_to_spherical(-lattice_vect(neighbour_boxes(:,i)))

             tmp      = OMultipole(l, m, s)
              ! we store the summands and order them before performing the sum
              ! to avoid numeric elimination
             rpart(i) =       real(tmp,  kind=kfp)
             ipart(i) = real(aimag(tmp), kind=kfp)
          end do

          call sort_abs(rpart)
          call sort_abs(ipart)

          rp = zero
          ip = zero

          do i = 1, num_neighbour_boxes,1 ! we sum up all contributions, starting from smallest
            rp = rp + rpart(i)
            ip = ip + ipart(i)
          end do

          MstarFunc = rp + ic*ip ! do not use cmplx()-function since it yields wrong results with complex*32 types

        end function MstarFunc


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Formal summation of \f$M\f$ over FF`, ie a lot of boxes
        !> with some overhead to avoid numerical elimination of small values
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        complex(kfp) function LstarFunc(l,m)
          implicit none
          integer, intent(in) :: l, m
          real(kfp) :: rpart(num_neighbour_boxes*(num_neighbour_boxes-1)), ipart(num_neighbour_boxes*(num_neighbour_boxes-1)),rp,ip, s(3)
          complex(kfp) :: tmp
          complex(kfp), parameter :: ic = (zero,one)
          integer :: i, ii, k

          k = 0

          ! sum over all boxes within FF' (cells in the far field of the central cell but in the near field of the central supercell 3x3x3 that embeds cell {0,0,0} in the center)
          do i = 1,num_neighbour_boxes-1 ! central box is being omitted in this loop
            do ii = 1,num_neighbour_boxes
              s   = cartesian_to_spherical((2*ws+1)*lattice_vect(neighbour_boxes(:,i)) + lattice_vect(neighbour_boxes(:,ii)))
              tmp = MTaylor(l, m, s)
              k   = k+1
              ! we store the summands and order them before performing the sum
              ! to avoid numeric elimination
              rpart(k) =       real(tmp,  kind=kfp)
              ipart(k) = real(aimag(tmp), kind=kfp)
            end do
          end do

          call sort_abs(rpart)
          call sort_abs(ipart)

          rp = zero
          ip = zero

          do i = 1,k ! we sum up all contributions, starting from smallest
            rp = rp + rpart(i)
            ip = ip + ipart(i)
          end do

          LStarFunc = rp + ic*ip ! do not use cmplx()-function since it yields wrong results with complex*32 types

        end function LstarFunc


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Converts cartesian coordinates to spherical system
        !> @param[in]  cartesian  cartesian vector [x, y, z]
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function cartesian_to_spherical(cartesian)
          implicit none
          real*8, intent(in)  :: cartesian(3)
          real(kfp), dimension(3) :: cartesian_to_spherical, c

          c = real(cartesian, kind=kfp)

          cartesian_to_spherical(1) = sqrt(dot_product(c, c))

          if (cartesian_to_spherical(1) == 0) then
            cartesian_to_spherical(2) = one
          else
            cartesian_to_spherical(2) = c(3) / cartesian_to_spherical(1)
          end if

          if ((c(1) == 0) .and. (c(2) == 0)) then
            cartesian_to_spherical(3) = zero
          else
            cartesian_to_spherical(3) = atan2(c(2), c(1))
          end if
        end function cartesian_to_spherical


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Writes contents of table T to a output stream s in a structured way
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine PrintTable(s, T, Lmax)
          implicit none
          complex(kfp), intent(in) :: T(:)
          integer, intent(in) :: Lmax
          integer, intent(in) :: s

          integer :: ll,mm,idx

          idx = 0

          do ll = 0,Lmax
            do mm = 0,ll
              idx = idx + 1

              write(s,'(I6, I6, I6, D50.35, D50.35)') idx, ll, mm, tbl(T, ll, mm, Lmax)
            end do
          end do

        end subroutine PrintTable


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Writes contents of table T to file s
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine WriteTableToFile(s, T, Lmax)
          implicit none
          complex(kfp), intent(in) :: T(:)
          integer, intent(in) :: Lmax
          character(len=*) :: s
          integer, parameter :: temp_file = 60

          open(temp_file, file=trim(s))

          write(temp_file,*) "# filename    = ", trim(s)
          write(temp_file,*) "# t_lattice_1 = ", t_lattice_1
          write(temp_file,*) "# t_lattice_2 = ", t_lattice_2
          write(temp_file,*) "# t_lattice_3 = ", t_lattice_3
          write(temp_file,*) "# periodicity_switches = ", periodicity_switches
          write(temp_file,*) "# Lmax_multipole = ", Lmax_multipole
          write(temp_file,*) "# Lmax_taylor    = ", Lmax_taylor

          write(temp_file,*) "# Lmax    = ", Lmax
          write(temp_file,*) "# MaxIter = ", MaxIter
          write(temp_file,*) "##########################"
          write(temp_file,*) "  idx     l     m           real-part                                         imaginary-part"

          call PrintTable(temp_file, T, Lmax)

          close(temp_file)

        end subroutine WriteTableToFile


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Dumps all parameters to the stream ifile
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine fmm_framework_param_dump(ifile)
          implicit none
          integer, intent(in) :: ifile

          if (myrank .ne. 0) return

          write(ifile,*) "LATTICE: ------------- Lattice fmm-framework switches ----------------"
          write(ifile,*) "LATTICE: Lmax_multipole = ", Lmax_multipole
          write(ifile,*) "LATTICE: Lmax_taylor    = ", Lmax_taylor
          write(ifile,*) "LATTICE: MaxIter        = ", MaxIter
          write(ifile,*) "LATTICE: ws             = ", ws
          write(ifile,*) "LATTICE: t_lattice_1    = ", t_lattice_1
          write(ifile,*) "LATTICE: t_lattice_2    = ", t_lattice_2
          write(ifile,*) "LATTICE: t_lattice_3    = ", t_lattice_3
          write(ifile,*) "LATTICE: periodicity    = ", periodicity
          write(ifile,*) "LATTICE: # neighbours    = ", num_neighbour_boxes
          write(ifile,*)
        end subroutine fmm_framework_param_dump


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Computes the associated Legendre polynomial \f$P_l_m (x)\f$.
        !> Here m and l are integers satisfying  \f$0 \leq m \leq l\f$,
        !> while x lies in the range \f$-1 \leq x \leq 1\f$.
        !>
        !> Code fragment for \f$P_l^m(x)\f$ taken from
        !>
        !> Numerical Recipes in Fortran 77: The Art of Scientific Computing
        !>              (ISBN 0-521-43064-X)
        !> pg. 246ff
        !>
        !> and modified to give \f$P_l_m (x)\f$:
        !> \f$ P_l_m (x) = (-1)^m P_l^m (x) \f$, see
        !>
        !> Abramowitz and Stegun: Handbook of Mathematical Functions
        !> Section 8. Legendre Functions (pg. 332)
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real(kfp) function LegendreP(l,m,x)
          use module_debug
          implicit none
          integer, intent(in) :: l, m
          real(kfp) ::x

          integer :: i,ll
          real(kfp) :: fact,pll,pmm,pmmp1,somx2

          pll = zero

          if ( (m < 0) .or. (m > l) .or. (abs(x) > 1) .or.  (l<0)) then
            DEBUG_ERROR(*,'Invalid arguments for LegendreP(',l,m,x,')')
          endif

          pmm = one     ! Compute P_m^m

          if (m > 0) then
            somx2 = sqrt((one-x)*(one+x))
            fact  = one

            do i = 1,m
               pmm  = -pmm * fact * somx2
               fact = fact+two
            enddo
          endif

          if (l == m) then
            LegendreP = pmm
          else
            pmmp1 = x*(two*m+one)*pmm  ! Compute P_m+1^m

            if (l == m+1) then
              LegendreP = pmmp1
            else                  ! Compute P_l^m , l > m + 1
              do ll = m+2,l
                pll   = (x*(two*ll-one)*pmmp1-(ll+m-one)*pmm)/real(ll-m,kind=kfp)
                pmm   = pmmp1
                pmmp1 = pll
              enddo

              LegendreP = pll
            endif
          endif

          LegendreP = (-one)**m * LegendreP

        end function LegendreP


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Calculates the factorial of the argument
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real(kfp) function factorial(n)
            use module_debug
            implicit none
            integer, intent(in) :: n

            if (n<0) then
              DEBUG_ERROR(*,"Tried to calculate factorial of negative argument.")
            end if

            if (n>170) then
              DEBUG_ERROR(*,"Tried to calculate factorial with n>170. This would lead to numeric overflow.")
            end if

            select case (n)
              case ( 0)
                factorial =                            one
              case ( 1)
                factorial =                            one
              case ( 2)
                factorial =                            two
              case ( 3)
                factorial =                            6._kfp
              case ( 4)
                factorial =                           24._kfp
              case ( 5)
                factorial =                          120._kfp
              case ( 6)
                factorial =                          720._kfp
              case ( 7)
                factorial =                         5040._kfp
              case ( 8)
                factorial =                        40320._kfp
              case ( 9)
                factorial =                       362880._kfp
              case (10)
                factorial =                      3628800._kfp
              case (11)
                factorial =                     39916800._kfp
              case (12)
                factorial =                    479001600._kfp
              case (13)
                factorial =                   6227020800._kfp
              case (14)
                factorial =                  87178291200._kfp
              case (15)
                factorial =                1307674368000._kfp
              case (16)
                factorial =               20922789888000._kfp
              case (17)
                factorial =              355687428096000._kfp
              case (18)
                factorial =             6402373705728000._kfp
              case (19)
                factorial =           121645100408832000._kfp
              case (20)
                factorial =          2432902008176640000._kfp
              case (21)
                factorial =         51090942171709440000._kfp
              case (22)
                factorial =       1124000727777607680000._kfp
              case (23)
                factorial =      25852016738884976640000._kfp
              case (24)
                factorial =     620448401733239439360000._kfp
              case default
                factorial = gamma(real(n+1, kfp))
            end select
        end function factorial


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Sorts the given values with a heap sort approach
        !> in order of ther absolute value
        !> compare (Numerical Recipes f90, p1171)
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          subroutine sort_abs(arr)
            implicit none
            real(kfp), intent(inout) :: arr(:)
            integer :: i,n

            n = size(arr)

            do i=n/2,1,-1                      ! Left range - hiring phase (heap creation)
               call sift_down(i,n)
            end do

            do i=n,2,-1                        ! Right range of sift-down is decr. from n-1 ->1
               ! during retirement/promotion (heap selection) phase.
               call swap( arr(1),arr(i) )      ! Clear space at end of array and retire top of heap into it
               call sift_down( 1,i-1)
            end do

          contains
            subroutine sift_down(l,r)        ! Carry out the sift-down on element arr(l) to maintain
              integer, intent(in) :: l,r     ! the heap structure
              integer :: j,jold    ! index
              real(kfp) :: a

              a    = arr(l)
              jold = l
              j    = l + l
              do                   ! do while j <= r
                 if (j > r) exit
                 if (j < r) then
                   if (abs(arr(j)) < abs(arr(j+1))) j = j+1
                 endif
                 if (abs(a) >= abs(arr(j))) exit       ! Found a`s level, so terminate sift-down

                 arr(jold) = arr(j)
                 jold      = j                    ! Demote a and continue
                 j         = j+j
              end do
              arr(jold) = a                  ! Put a into its slot

            end subroutine sift_down

            subroutine swap(p,q)
                real(kfp) :: p,q, dum
                dum = p
                p   = q
                q   = dum
            end subroutine swap

          end subroutine sort_abs


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !>  Sets all matrix entries that are smaller than 1.e-16 to 0.
        !> (separately for real and imaginary part)
        !> This is the same as Mathematicas Chop[]-function
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine chop(a)
          implicit none
          complex(kfp), intent(inout) :: a(:)
          integer :: i
          complex(kfp), parameter :: ic = (zero,one)
          real(kfp) :: re, im

          if (.not. chop_arrays) return

          DEBUG_WARNING(*, 'chopping some array')

          do i=1,size(a)
            re =       real(a(i),  kind=kfp)
            im = real(aimag(a(i)), kind=kfp)

            if (abs(re) < prec) re = zero
            if (abs(im) < prec) im = zero

            a(i) = re + ic*im ! do not use cmplx()-function here (see above)
          end do
        end subroutine chop


end module module_fmm_framework





