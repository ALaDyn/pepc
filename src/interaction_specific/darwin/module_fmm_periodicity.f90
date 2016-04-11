! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
!
! Copyright (C) 2002-2016 Juelich Supercomputing Centre,
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

#include "multipole.h"

!>
!> Encapsulates calculation of the lattice contribution by means
!> of the FMM-approach to the lattice
!>
module module_fmm_framework
      use module_pepc_kinds
      use module_interaction_specific_types, only : t_tree_node_interaction_data
      use module_shortcut
      implicit none
      include 'mpif.h'
      private

      ! public variable declarations

      !> far- and near-field contribution to potential energy (has to be calculated in fields.p90)
      real(kind_physics), public :: potfarfield, potnearfield
      !> whether to do dipole correction or not, see [Kudin 1998, eq. (2.8)]
      logical, public :: do_extrinsic_correction = .true.

      ! public subroutine declarations
      public fmm_framework_prepare
      public fmm_framework_timestep
      public fmm_sum_lattice_force
      public fmm_framework_param_dump

      ! private variable declarations

      ! general stuff
      integer :: myrank
      integer :: MPI_COMM_fmm
      ! precision flags
      ! for now, kfp must equal kind_physics, since the FMM operators (see module_multipole) are defined on kind_physics
      integer, parameter :: kfp                  = kind_physics ! numeric precision (kind value)
      integer, parameter :: MPI_REAL_fmm         = MPI_KIND_PHYSICS
      logical, parameter  :: chop_arrays         = .false.
!      real(kfp), parameter :: prec = 1.E-16_kfp
!      ! shortcut notations
!      real(kfp), parameter :: zero = 0._kfp
!      real(kfp), parameter :: one  = 1._kfp
!      real(kfp), parameter :: two  = 2._kfp
!      real(kfp), parameter :: half = one / two
      complex(kfp), parameter :: czero = (zero, zero)
      complex(kfp), parameter :: ic = (zero, one)
      ! FMM-PARAMETERS
      integer, parameter :: pMultipoleFMMP = 30
      integer, parameter :: qTaylorFMMP    = pMultipoleFMMP * 2
      integer, parameter :: MaxIter        = 32
      integer :: ws = 1
      ! internally calculated FMM variables
      complex(kfp) :: mu_cent(0:qTaylorFMMP)
      complex(kfp) :: omega_tilde(1:pMultipoleFMMP)
      complex(kfp) :: BLattice(0:qTaylorFMMP, 1:pMultipoleFMMP)
      !> fictitious charges
      type(t_tree_node_interaction_data) :: fictcharge(1:4)
      integer :: nfictcharge = 0

      ! subroutine-implementation

      interface WriteTableToFile
        module procedure WriteTableToFile1D, WriteTableToFile2D
      end interface WriteTableToFile

      interface chop
        module procedure chop1D, chop2D
      end interface chop

      interface matmul_careful
        module procedure matmul_carefulMM, matmul_carefulMV
      end interface matmul_careful

      contains

        !>
        !> Module Initialization, should be called on program startup
        !> after setting up all geometric parameters etc.
        !>  - well-separation criterion is automatically set to module_mirror_boxes::mirror_box_layers
        !>  - module_mirror_boxes::neighbour_boxes_prepare() must have been called before calling this routine
        !>
        !> @param[in] mpi_rank MPI rank of process for controlling debug output
        !> @param[in] mpi_comm MPI communicator to be used
        !>
        subroutine fmm_framework_prepare(mpi_rank, mpi_comm)
          use module_debug
          use module_mirror_boxes
          implicit none
          integer, intent(in) :: mpi_rank
          integer, intent(in) :: mpi_comm

          myrank       = mpi_rank
          MPI_COMM_fmm = mpi_comm
          ws           = mirror_box_layers

          if (periodicity(3)) then
            DEBUG_ERROR(*, 'Periodicity in the third coordinate with 2D backend.')
          end if
          if (.not. all(t_lattice_3 == [0,0,1])) then
            DEBUG_WARNING(*, 't_lattice_3 should be [0,0,1], resetting.')
            t_lattice_3 = [0,0,1]
            call neighbour_boxes_prepare()
          end if
          LatticeCenter(3) = 0

          BLattice = czero

          call calc_lattice_coefficients(BLattice)

          if ((myrank == 0) .and. dbg(DBG_PERIODIC)) then
            call WriteTableToFile('BLattice.tab', BLattice)
          end if
        end subroutine fmm_framework_prepare


        !>
        !> Refreshes Multipole information and Taylor coefficients,
        !> has to be called every timestep with particles that were used in tree buildup
        !>
        subroutine fmm_framework_timestep(particles)
          use module_pepc_types
          use module_mirror_boxes, only: check_lattice_boundaries
          use module_debug
          implicit none
          type(t_particle), intent(in) :: particles(:)

          if (.not. check_lattice_boundaries(particles)) then
            DEBUG_ERROR(*, 'Lattice contribution will be wrong. Aborting.')
          end if

          call calc_omega_tilde(particles)
          call calc_mu_cent(omega_tilde, mu_cent)
        end subroutine fmm_framework_timestep


        !>
        !> Calculates the lattice coefficients for computing mu_cent
        !>
        subroutine calc_lattice_coefficients(BL)
          use module_debug
          implicit none

          complex(kfp), intent(out) :: BL(0:qTaylorFMMP, 1:pMultipoleFMMP)
          ! Variables for lattice coefficient calculation
          complex(kfp) :: Astar(1:pMultipoleFMMP, 1:pMultipoleFMMP)
          complex(kfp) :: Bstar(0:qTaylorFMMP, 1:pMultipoleFMMP)

          integer :: k, l, iter
          character(20) :: fn

          call pepc_status('LATTICE COEFFICIENTS: Starting calculation')

          Astar = czero
          Bstar = czero

          ! pretabulation of necessary values
          do k = 1, pMultipoleFMMP
            do l = 1,k
              Astar(k, l) = AstarFunc(k, l)
            end do
          end do

          do k = 0, qTaylorFMMP
            do l = 1, pMultipoleFMMP
              Bstar(k, l)    = BstarFunc(k, l)
            end do
          end do

          call chop(Astar)
          call chop(Bstar)

          if ((myrank == 0) .and. dbg(DBG_PERIODIC)) then
            call WriteTableToFile('Astar.tab', Astar)
            call WriteTableToFile('Bstar.tab', Bstar)
          end if

          ! zeroth step of iteration
          BL = Bstar

          do iter = 1,MaxIter

            BL = matmul_careful( UB( BL ) , Astar ) + Bstar

            if ((myrank == 0) .and. dbg(DBG_PERIODIC)) then
              write(fn,'("BLattice.", I2.2, ".tab")') iter
              call WriteTableToFile(trim(fn), BL)
            endif
          end do

          call chop(BL)

          call pepc_status('LATTICE COEFFICIENTS: finished calculation')

        end subroutine calc_lattice_coefficients


        !>
        !> Calculates the overall multipole expansion of the whole
        !> central box
        !>
        subroutine calc_omega_tilde(particles)
          use module_pepc_types
          use module_mirror_boxes, only: LatticeOrigin, LatticeCenter, Lattice, periodicity
          use module_debug
          implicit none

          type(t_particle), intent(in) :: particles(:)

          integer :: k
          integer(kind_particle) :: p
          integer :: ierr
          real(kfp) :: qtotal

          qtotal = zero
          omega_tilde = czero

          ! calculate multipole contributions of all local particles
          do p=1,size(particles, kind=kind(p))
            call addparticle(particles(p)%x(1:2), particles(p)%data%q)
          end do

          call chop(omega_tilde)

          ! sum multipole contributions from all processors - treat complex as two real numbers since e.g. complex*32 is not supported by mpi
          call MPI_ALLREDUCE(MPI_IN_PLACE, qtotal, 1, MPI_REAL_fmm, MPI_SUM, MPI_COMM_fmm, ierr)
          call MPI_ALLREDUCE(MPI_IN_PLACE, omega_tilde, 2 * pMultipoleFMMP, MPI_REAL_fmm, MPI_SUM, MPI_COMM_fmm, ierr)

          call chop(omega_tilde)

          if (do_extrinsic_correction) then
            nfictcharge = 0

            if (periodicity(1)) then
              nfictcharge = nfictcharge + 1
              fictcharge(nfictcharge)%coc(1:2) = LatticeOrigin(1:2) + Lattice(1, 1:2)
              fictcharge(nfictcharge)%charge = real(omega_tilde(1), kind = kind_physics) / sqrt(dot_product(Lattice(1, 1:2), Lattice(1, 1:2)))

              nfictcharge = nfictcharge + 1
              fictcharge(nfictcharge)%coc(1:2) = LatticeOrigin(1:2)
              fictcharge(nfictcharge)%charge = -fictcharge(nfictcharge - 1)%charge
            end if

            if (periodicity(2)) then
              nfictcharge = nfictcharge + 1
              fictcharge(nfictcharge)%coc(1:2) = LatticeOrigin(1:2) + Lattice(2, 1:2)
              fictcharge(nfictcharge)%charge = real(aimag(omega_tilde(1)), kind = kind_physics) / sqrt(dot_product(Lattice(2, 1:2), Lattice(2, 1:2)))

              nfictcharge = nfictcharge + 1
              fictcharge(nfictcharge)%coc(1:2) = LatticeOrigin(1:2)
              fictcharge(nfictcharge)%charge = -fictcharge(nfictcharge - 1)%charge
            end if

            if (nfictcharge > 0) then
              do p=1,nfictcharge
                call addparticle(fictcharge(p)%coc(1:2), fictcharge(p)%charge)
              end do

              call chop(omega_tilde)
            end if
          end if

          if ((myrank == 0) .and. dbg(DBG_PERIODIC)) then
            call WriteTableToFile('omega_tilde.tab', omega_tilde)
          end if

          if (qtotal > 0.) then
            DEBUG_WARNING(*, 'The central box is not charge-neutral: Q_total=', qtotal, ' Setting to zero, resulting potentials might be wrong.' )
          end if

          contains
            subroutine addparticle(R, q)
              implicit none

              real(kind_physics), intent(in) :: R(2), q

              real(kfp) :: x0(2)
              complex(kfp) :: z0

              x0 = real(R(1:2) - LatticeCenter(1:2), kind = kfp)
              z0 = x0(1) + ic * x0(2)
              qtotal = qtotal + real(q, kind = kfp)

              do k=1,pMultipoleFMMP
                omega_tilde(k) = omega_tilde(k) + omega(k, z0, real(q, kind = kfp))
              end do
            end subroutine addparticle
        end subroutine calc_omega_tilde


        !>
        !> Calculates the lattice contribution with respect to the
        !> centre of the original box
        !>
        subroutine calc_mu_cent(omega, mu)
          use module_debug
          implicit none
          complex(kfp), intent(in) :: omega(1:pMultipoleFMMP)
          complex(kfp), intent(out) :: mu(0:qTaylorFMMP)

          ! contribution of all outer lattice cells, with regards to the centre of the original box
          mu = matmul_careful(BLattice, omega)

          call chop(mu)

          if ((myrank == 0) .and. dbg(DBG_PERIODIC)) then
            call WriteTableToFile('mu_cent.tab', mu_cent)
          end if

        end subroutine calc_mu_cent


        !>
        !> Calculates the charged moments of the multipole expansion
        !> for a certain particle
        !>
        complex(kfp) function omega(k, z0, q)
          use module_multipole
          implicit none
          integer, intent(in) :: k
          complex(kfp), intent(in) :: z0
          real(kfp), intent(in) :: q

          omega = - q * OMultipole(k, z0) / k

        end function omega


        !>
        !> Calculates force at individual position that results
        !> from mirror boxes beyond the near field region,
        !> i.e. the lattice contribution
        !>
        subroutine fmm_sum_lattice_force(pos, e_lattice, phi_lattice)
          use module_mirror_boxes, only: LatticeCenter, num_neighbour_boxes, lattice_vect, neighbour_boxes
          use module_multipole
          use module_debug
          implicit none

          real(kind_physics), intent(in) :: pos(3)
          real(kind_physics), intent(out) ::  e_lattice(2), phi_lattice

          integer :: k, ibox
          integer(kind_particle) :: p
          real(kfp) :: x0(2)
          complex(kfp) :: z0, cphi, ce
          real(kind_physics) :: delta(3), phitmp, etmp(2)

          x0 = real(pos(1:2) - LatticeCenter(1:2), kind = kfp)
          z0 = x0(1) + ic * x0(2)

          cphi = -mu_cent(0) ! OMultipole(0, z0) = 1
          ce   = czero       ! OmultipolePrime(0, z0) = 0

          do k = 1, qTaylorFMMP
            cphi = cphi - mu_cent(k) * OMultipole(k, z0)
            ce   = ce   - mu_cent(k) * OMultipolePrime(k, z0)
          end do

          ! E = -grad(Phi)
          e_lattice  = [ -real(ce, kind = kind_physics), real(aimag(ce), kind = kind_physics) ]
          phi_lattice = real(cphi, kind = kind_physics)

          if (do_extrinsic_correction) then    ! extrinsic correction
            do p=1,nfictcharge
              do ibox=1,num_neighbour_boxes
                delta = pos - lattice_vect(neighbour_boxes(:,ibox)) - fictcharge(p)%coc
                call log2d_kernel(fictcharge(p)%charge, delta(1:2), phitmp, etmp(1:2))
                e_lattice = e_lattice + etmp
                phi_lattice = phi_lattice + phitmp
              end do
            end do
          end if

          contains
            subroutine log2d_kernel(q, d, phi, e)
              implicit none

              real(kind_physics), intent(in) :: q, d(2)
              real(kind_physics), intent(out) :: phi, e(2)

              real(kfp) :: d2, rd2

              d2 = dot_product(d, d)
              rd2 = one / d2

              phi = -half * q * log(d2)
              e   = q * d * rd2
            end subroutine log2d_kernel


        end subroutine fmm_sum_lattice_force


        subroutine fmm_sum_lattice_force_darwin(pos,v, e_lattice, phi_lattice, B_lattice, A_lattice)
          use module_mirror_boxes, only: LatticeCenter, num_neighbour_boxes, lattice_vect, neighbour_boxes
          use module_multipole
          use module_debug
          implicit none

          real(kind_physics), intent(in)  ::  pos(3),v(3)
          real(kind_physics), intent(out) ::  e_lattice(2), phi_lattice, B_lattice(3), A_lattice(3)

          integer                         ::  k, ibox
          integer(kind_particle)          ::  p
          real(kfp)                       ::  x0(2)
          complex(kfp)                    ::  z0, cphi, ce, cA_log, cA_exp, cB_log
          real(kind_physics)              ::  delta(3), phitmp, etmp(2),Atmp(1:3),Btmp(1:3),Atmp1(1:2),Atmp2(1:2),Btmp1(1:2)

          x0 = real(pos(1:2) - LatticeCenter(1:2), kind = kfp)
          z0 = x0(1) + ic * x0(2)

          cphi     = -mu_cent(0)  ! OMultipole(0, z0) = 1
          ce       =  czero       ! OmultipolePrime(0, z0) = 0
          cA_log   =  czero       ! OmultipolePrime(0, z0) = 0
          cA_exp   =  czero       ! OmultipolePrime(0, z0) = 0
          cB_log   =  czero       ! OmultipolePrime(0, z0) = 0

          do k = 1, qTaylorFMMP
            cphi     = cphi     - mu_cent(k) * OMultipole(k, z0)
            ce       = ce       - mu_cent(k) * OMultipolePrime(k, z0)
            cA_log   = cA_log   - mu_cent(k) * OMultipolePrime(k, log(z0))
            cA_exp   = cA_exp   - mu_cent(k) * OMultipolePrime(k, exp(z0))
            cB_log   = cB_log   - mu_cent(k) * OMultipolePrime(k, z0)
          end do

          ! E = -grad(Phi)
          e_lattice         = [ -real(ce, kind = kind_physics), real(aimag(ce), kind = kind_physics) ]
          phi_lattice       = real(cphi, kind = kind_physics)
          Atmp1             = [ -real(cA_log, kind = kind_physics), real(aimag(cA_log), kind = kind_physics) ]
          Atmp2             = [ -real(cA_exp, kind = kind_physics), real(aimag(cA_exp), kind = kind_physics) ]
          Btmp1             = [ -real(cB_log, kind = kind_physics), real(aimag(cB_log), kind = kind_physics) ]

          A_lattice(1:3)    = phi_lattice*v
          A_lattice(1:2)    = A_lattice(1:2) + ( v(1)*Atmp1(1) + v(2)*Atmp1(2)  )*Atmp2(1:2)
          B_lattice(1:3)    = Btmp(2)*v(3)
          B_lattice(2)      =-Btmp(1)*v(3)
          B_lattice(3)      = Btmp(1)*v(2)   - Btmp(2)*v(1)

          if (do_extrinsic_correction) then    ! extrinsic correction
            do p=1,nfictcharge
              do ibox=1,num_neighbour_boxes
                delta = pos - lattice_vect(neighbour_boxes(:,ibox)) - fictcharge(p)%coc
                call log2d_kernel_darwin2D3V(fictcharge(p)%charge,fictcharge(p)%monoj(1:3), delta(1:2), phitmp, etmp(1:2),Atmp(1:3),Btmp(1:3))
                e_lattice       = e_lattice   + etmp
                phi_lattice     = phi_lattice + phitmp
                A_lattice       = A_lattice   + Atmp
                B_lattice       = B_lattice   + Btmp
              end do
            end do
          end if

          contains


            subroutine log2d_kernel_darwin2D3V(q, j, d, phi, e, A, B)
              use module_tool   , only: cross_product
              use module_shortcut, only: half,quarter,zero
              implicit none

              real(kind_physics), intent(in)     :: q, d(2),j(3)
              real(kind_physics), intent(out)    :: phi, e(2), A(3), B(3)

              real(kfp) :: d2, rd2,x(3)

              x(1:3)   = zero
              x(1:2)   = d(1:2)
              d2       = dot_product(d, d)
              rd2      = one / d2

              phi      = -half * q * log(d2)
              e        = q * d * rd2
              A(1:3)   = -quarter* q * log(d2)*v(1:3)
              A(1:2)   = A(1:2) + half*q*( j(1)*d(1) + j(2)*d(2) )*j(1:2)*rd2
              A(3)     = -half* q * log(d2)*j(3)
              B(1:3)   = -q*rd2 * cross_product(x,j)

            end subroutine log2d_kernel_darwin2D3V

        end subroutine fmm_sum_lattice_force_darwin




        !>
        !> Scaling Operator \f$\mathcal{U}_B\f$ for Taylor coefficients
        !> @param[in] B table with Taylor coefficients
        !>
        function UB(B)
          implicit none
          complex(kfp), intent(in) :: B(0:qTaylorFMMP, 1:pMultipoleFMMP)
          complex(kfp), dimension(0:qTaylorFMMP, 1:pMultipoleFMMP) :: UB
          real(kfp) :: scalefac

          integer :: k, l

          scalefac = real(2*ws+1, kind=kfp)

          do k = 0, qTaylorFMMP
            do l = 1, pMultipoleFMMP
              UB(k, l) = B(k, l) / scalefac**real(k + l, kind=kfp)
            end do
          end do

          call chop(UB)

        end function UB


        !>
        !> Formal summation of \f$A\f$ over NF, ie all (9) neighbouring boxes
        !> with some overhead to avoid numerical elimination of small values
        !>
        complex(kfp) function AstarFunc(k,l)
          use module_mirror_boxes, only: neighbour_boxes, num_neighbour_boxes,&
            lattice_vect
          use module_multipole
          implicit none
          integer, intent(in) :: k, l
          real(kfp) :: rpart(num_neighbour_boxes), ipart(num_neighbour_boxes), rp, ip, x0(3)
          complex(kfp) :: z0, tmp

          integer :: i

          do i = 1,num_neighbour_boxes
            x0 = -lattice_vect(neighbour_boxes(:, i))
            z0 = x0(1) + ic * x0(2)
            tmp      = ATranslate(k, l, z0)
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

          AstarFunc = rp + ic*ip ! do not use cmplx()-function since it yields wrong results with complex*32 types

        end function AstarFunc


        !>
        !> Formal summation of \f$B\f$ over FF`, ie a lot of boxes
        !> with some overhead to avoid numerical elimination of small values
        !>
        complex(kfp) function BstarFunc(k,l)
          use module_mirror_boxes, only: neighbour_boxes, num_neighbour_boxes,&
            lattice_vect
          use module_multipole
          implicit none
          integer, intent(in) :: k, l
          real(kfp) :: rpart(num_neighbour_boxes*(num_neighbour_boxes-1)), ipart(num_neighbour_boxes*(num_neighbour_boxes-1)), rp, ip, x0(3)
          complex(kfp) :: z0, tmp
          integer :: i, ii, j

          j = 0

          ! sum over all boxes within FF' (cells in the far field of the central cell but in the near field of the central supercell 3x3x3 that embeds cell {0,0,0} in the center)
          do i = 1,num_neighbour_boxes-1 ! central box is being omitted in this loop
            do ii = 1,num_neighbour_boxes
              x0  = (2*ws+1)*lattice_vect(neighbour_boxes(:,i)) + lattice_vect(neighbour_boxes(:,ii))
              z0  = x0(1) + ic * x0(2)
              tmp = BConvert(k, l, z0)
              j   = j+1
              ! we store the summands and order them before performing the sum
              ! to avoid numeric elimination
              rpart(j) =       real(tmp,  kind=kfp)
              ipart(j) = real(aimag(tmp), kind=kfp)
            end do
          end do

          call sort_abs(rpart)
          call sort_abs(ipart)

          rp = zero
          ip = zero

          do i = 1,j ! we sum up all contributions, starting from smallest
            rp = rp + rpart(i)
            ip = ip + ipart(i)
          end do

          BStarFunc = rp + ic*ip ! do not use cmplx()-function since it yields wrong results with complex*32 types

        end function BstarFunc


        !>
        !> Writes contents of table T to file s, 1 index
        !>
        subroutine WriteTableToFile1D(s, T)
          use module_mirror_boxes
          implicit none
          complex(kfp), intent(in) :: T(:)
          character(len=*) :: s
          integer, parameter :: temp_file = 60

          open(temp_file, file=trim(s))

          write(temp_file,*) "# filename    = ", trim(s)
          write(temp_file,*) "# t_lattice_1 = ", t_lattice_1
          write(temp_file,*) "# t_lattice_2 = ", t_lattice_2
          write(temp_file,*) "# t_lattice_3 = ", t_lattice_3
          write(temp_file,*) "# periodicity_switches = ", periodicity_switches
          write(temp_file,*) "# pMultipole  = ", pMultipoleFMMP
          write(temp_file,*) "# qTaylor     = ", qTaylorFMMP

          write(temp_file,*) "# MaxIter = ", MaxIter
          write(temp_file,*) "##########################"
          write(temp_file,*) "    k           real-part                                         imaginary-part"

          call PrintTable(temp_file, T)

          close(temp_file)

          contains
          subroutine PrintTable(s, T)
            implicit none
            complex(kfp), intent(in) :: T(:)
            integer, intent(in) :: s

            integer :: k

            do k = lbound(T, dim = 1), ubound(T, dim = 1)
              write(s,'(I6, D50.35, D50.35)') k, T(k)
            end do

          end subroutine PrintTable
        end subroutine WriteTableToFile1D


        !>
        !> Writes contents of table T to file s
        !>
        subroutine WriteTableToFile2D(s, T)
          use module_mirror_boxes
          implicit none
          complex(kfp), intent(in) :: T(:, :)
          character(len=*) :: s
          integer, parameter :: temp_file = 60

          open(temp_file, file=trim(s))

          write(temp_file,*) "# filename    = ", trim(s)
          write(temp_file,*) "# t_lattice_1 = ", t_lattice_1
          write(temp_file,*) "# t_lattice_2 = ", t_lattice_2
          write(temp_file,*) "# t_lattice_3 = ", t_lattice_3
          write(temp_file,*) "# periodicity_switches = ", periodicity_switches
          write(temp_file,*) "# pMultipole  = ", pMultipoleFMMP
          write(temp_file,*) "# qTaylor     = ", qTaylorFMMP

          write(temp_file,*) "# MaxIter = ", MaxIter
          write(temp_file,*) "##########################"
          write(temp_file,*) "    k     l           real-part                                         imaginary-part"

          call PrintTable(temp_file, T)

          close(temp_file)

          contains
          subroutine PrintTable(s, T)
            implicit none
            complex(kfp), intent(in) :: T(:, :)
            integer, intent(in) :: s

            integer :: k, l

            do k = lbound(T, dim = 1), ubound(T, dim = 1)
              do l = lbound(T, dim = 2),ubound(T, dim = 2)
                write(s,'(I6, I6, D50.35, D50.35)') k, l, T(k, l)
              end do
            end do

          end subroutine PrintTable
        end subroutine WriteTableToFile2D


        !>
        !> Dumps all parameters to the stream ifile
        !>
        subroutine fmm_framework_param_dump(ifile)
          use module_mirror_boxes
          implicit none
          integer, intent(in) :: ifile

          if (myrank .ne. 0) return

          write(ifile,*) "LATTICE: ------------- Lattice fmm-framework switches ----------------"
          write(ifile,*) "LATTICE: pMultipole   = ", pMultipoleFMMP
          write(ifile,*) "LATTICE: qTaylor      = ", qTaylorFMMP
          write(ifile,*) "LATTICE: MaxIter      = ", MaxIter
          write(ifile,*) "LATTICE: ws           = ", ws
          write(ifile,*) "LATTICE: t_lattice_1  = ", t_lattice_1
          write(ifile,*) "LATTICE: t_lattice_2  = ", t_lattice_2
          write(ifile,*) "LATTICE: t_lattice_3  = ", t_lattice_3
          write(ifile,*) "LATTICE: periodicity  = ", periodicity
          write(ifile,*) "LATTICE: # neighbours = ", num_neighbour_boxes
          write(ifile,*)
        end subroutine fmm_framework_param_dump


        function matmul_carefulMM(A, B) result(matmul_careful)
          implicit none

          complex(kfp), intent(in) :: A(:, :), B(:, :)
          complex(kfp) :: matmul_careful(size(A, dim = 1), size(B, dim = 2))

          complex(kfp) :: tmp
          real(kfp) :: rpart(size(A, dim = 2)), ipart(size(A, dim = 2)), rp, ip
          integer :: i, j, k

          do i = 1, size(A, dim = 1)
            do j = 1, size(B, dim = 2)

              rpart = zero
              ipart = zero

              do k = 1, size(A, dim = 2)
                tmp = A(i,k) * B(k, j)

                rpart(k) = real(tmp, kind = kfp)
                ipart(k) = real(aimag(tmp), kind = kfp)
              end do

              call sort_abs(rpart)
              call sort_abs(ipart)

              rp = zero
              ip = zero

              do k = 1, size(A, dim = 2)
                rp = rp + rpart(k)
                ip = ip + ipart(k)
              end do

              matmul_careful(i, j) = rp + ic * ip
            end do
          end do

        end function matmul_carefulMM


        function matmul_carefulMV(A, x) result(matmul_careful)
          implicit none

          complex(kfp), intent(in) :: A(:, :), x(:)
          complex(kfp) :: matmul_careful(size(A, dim = 1))

          complex(kfp) :: tmp
          real(kfp) :: rpart(size(A, dim = 2)), ipart(size(A, dim = 2)), rp, ip
          integer :: i, j

          do i = 1, size(A, dim = 1)
            rpart = zero
            ipart = zero

            do j = 1, size(A, dim = 2)
              tmp = A(i,j) * x(j)
              rpart(j) = real(tmp, kind = kfp)
              ipart(j) = real(aimag(tmp), kind = kfp)
            end do

            call sort_abs(rpart)
            call sort_abs(ipart)

            rp = zero
            ip = zero

            do j = 1, size(A, dim = 2)
              rp = rp + rpart(j)
              ip = ip + ipart(j)
            end do

            matmul_careful(i) = rp + ic * ip
          end do
        end function matmul_carefulMV


        !>
        !> Sorts the given values with a heap sort approach
        !> in order of ther absolute value
        !> compare (Numerical Recipes f90, p1171)
        !>
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


        !>
        !>  Sets all matrix entries that are smaller than 1.e-16 to 0.
        !> (separately for real and imaginary part)
        !> This is the same as Mathematicas Chop[]-function
        !>
        subroutine chop1D(a)
          use module_debug
          implicit none
          complex(kfp), intent(inout) :: a(:)
          integer :: i
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
        end subroutine chop1D


        subroutine chop2D(a)
          use module_debug
          implicit none
          complex(kfp), intent(inout) :: a(:,:)
          integer :: i, j
          real(kfp) :: re, im

          if (.not. chop_arrays) return

          DEBUG_WARNING(*, 'chopping some array')

          do i=lbound(a, dim = 2), ubound(a, dim = 2)
            do j=lbound(a, dim = 2), ubound(a, dim = 2)
              re =       real(a(i,j),  kind=kfp)
              im = real(aimag(a(i,j)), kind=kfp)

              if (abs(re) < prec) re = zero
              if (abs(im) < prec) im = zero

              a(i,j) = re + ic*im ! do not use cmplx()-function here (see above)
            end do
          end do
        end subroutine chop2D
end module module_fmm_framework
