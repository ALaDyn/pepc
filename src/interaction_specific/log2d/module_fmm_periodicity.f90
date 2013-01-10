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

#include "multipole.h"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!> Encapsulates calculation of the lattice contribution by means
!> of the FMM-approach to the lattice
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_fmm_periodicity
      use module_mirror_boxes, only: lattice_vect
      implicit none
      include 'mpif.h'
      private

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
      public fmm_periodicity_init
      public fmm_periodicity_timestep
      public fmm_periodicity_sum_lattice_force
      public lattice_vect
      public fmm_periodicity_param_dump

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
      integer, parameter :: pMultipoleFMMP = 20
      integer, parameter :: qTaylorFMMP    = pMultipoleFMMP * 2
      integer, parameter :: MaxIter        = 32
      integer :: ws = 1
      ! internally calculated FMM variables
      complex(kfp) :: mu_cent(0:qTaylorFMMP)
      complex(kfp) :: omega_tilde(1:pMultipoleFMMP)
      complex(kfp) :: BLattice(0:qTaylorFMMP, 1:pMultipoleFMMP)
      !> variables for extrinsic to intrinsic correction
      real(kfp) :: box_dipole(3) = zero
      real(kfp) :: quad_trace    = zero

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  subroutine-implementation  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      interface WriteTableToFile
        module procedure WriteTableToFile1D, WriteTableToFile2D
      end interface WriteTableToFile

      interface chop
        module procedure chop1D, chop2D
      end interface chop

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
        subroutine fmm_periodicity_init(mpi_rank, mpi_comm)
          use module_debug
          use module_mirror_boxes, only : periodicity, mirror_box_layers, &
            LatticeCenter, do_periodic
          implicit none
          integer, intent(in) :: mpi_rank
          integer, intent(in) :: mpi_comm

          myrank       = mpi_rank
          MPI_COMM_fmm = mpi_comm
          ws           = mirror_box_layers

          do_periodic = any(periodicity(1:3))

          ! anything above has to be done in any case
          if (do_periodic) then
            if (periodicity(3)) then
              DEBUG_ERROR(*, 'Periodicity along the z-axis with 2D backend.')
            end if
            LatticeCenter(3) = 0

            BLattice = 0

            call calc_lattice_coefficients(BLattice)

            if ((myrank == 0) .and. dbg(DBG_PERIODIC)) then
              call WriteTableToFile('BLattice.tab', BLattice)
            end if
          end if

        end subroutine fmm_periodicity_init


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Refreshes Multipole information and Taylor coefficients,
        !> has to be called every timestep with particles that were used in tree buildup
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine fmm_periodicity_timestep(particles, nparticles)
          use module_pepc_types
          use module_mirror_boxes, only: check_lattice_boundaries, do_periodic
          use module_debug
          implicit none
          integer, intent(in) :: nparticles
          type(t_particle), intent(in) :: particles(:)

          if (do_periodic) then
            if (.not. check_lattice_boundaries(particles, nparticles)) then
              DEBUG_ERROR(*, 'Lattice contribution will be wrong. Aborting.')
            endif
            
            call calc_omega_tilde(particles, nparticles)
            call calc_mu_cent(omega_tilde, mu_cent)
            call calc_extrinsic_correction(particles, nparticles)
          endif
        end subroutine fmm_periodicity_timestep


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Calculates the lattice coefficients for computing mu_cent
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

          Astar = 0
          Bstar = 0

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


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Calculates the overall multipole expansion of the whole
        !> central box
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine calc_omega_tilde(particles, nparticles)
          use module_pepc_types
          use module_mirror_boxes, only: LatticeCenter
          use module_debug
          implicit none

          type(t_particle), intent(in) :: particles(:)
          integer, intent(in) :: nparticles

          integer :: k, p
          integer :: ierr
          complex(kfp), parameter :: ic = (zero, one)
          complex(kfp) :: z0
          real*8 :: qtotal, x0(3)

          qtotal = 0
          omega_tilde = 0

          ! calculate multipole contributions of all local particles
          do p=1,nparticles
            x0 = particles(p)%x - LatticeCenter
            z0 = x0(1) + ic * x0(2)
            qtotal = qtotal + particles(p)%data%q

            do k=1,pMultipoleFMMP
              omega_tilde(k) = omega_tilde(k) + omega(k, z0, particles(p)%data%q)
            end do

          end do

          call chop(omega_tilde)

          ! sum multipole contributions from all processors - treat complex as two real numbers since e.g. complex*32 is not supported by mpi
          call MPI_ALLREDUCE(MPI_IN_PLACE, qtotal, 1, MPI_REAL8, MPI_SUM, MPI_COMM_fmm, ierr)
          call MPI_ALLREDUCE(MPI_IN_PLACE, omega_tilde, 2 * pMultipoleFMMP, MPI_REAL_fmm, MPI_SUM, MPI_COMM_fmm, ierr)

          call chop(omega_tilde)

          if ((myrank == 0) .and. dbg(DBG_PERIODIC)) then
            call WriteTableToFile('omega_tilde.tab', omega_tilde)
          end if

          if (qtotal > 0.) then
            DEBUG_WARNING(*, 'The central box is not charge-neutral: Q_total=', qtotal, ' Setting to zero, resulting potentials might be wrong.' )
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
        !>       ^ inside this publication, the volume factor is missing
        !>  [J. Chem. Phys. 101, 5024, eqn (5)] contains this volume
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine calc_extrinsic_correction(particles, nparticles)
          use module_debug
          use module_pepc_types
          use module_mirror_boxes, only : unit_box_volume
          implicit none

          type(t_particle), intent(in) :: particles(:)
          integer, intent(in) :: nparticles

          ! TODO: implement this

        end subroutine calc_extrinsic_correction


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Calculates the lattice contribution with respect to the
        !> centre of the original box
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine calc_mu_cent(omega, mu)
          use module_debug
          implicit none
          complex(kfp), intent(in) :: omega(1:pMultipoleFMMP)
          complex(kfp), intent(out) :: mu(0:qTaylorFMMP)

          ! contribution of all outer lattice cells, with regards to the centre of the original box
          mu = matmul(BLattice, omega)

          call chop(mu)

          if ((myrank == 0) .and. dbg(DBG_PERIODIC)) then
            call WriteTableToFile('mu_cent.tab', mu_cent)
          end if

        end subroutine calc_mu_cent


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Calculates the charged moments of the multipole expansion
        !> for a certain particle
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        complex(kfp) function omega(k, z0, q)
          use module_multipole
          implicit none
          integer, intent(in) :: k
          complex(kfp), intent(in) :: z0
          real*8, intent(in) :: q

          omega = - real(q, kind=kfp) * OMultipole(k, z0) / k

        end function omega


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Calculates force at individual position that results
        !> from mirror boxes beyond the near field region,
        !> i.e. the lattice contribution
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine fmm_periodicity_sum_lattice_force(pos, e_lattice, phi_lattice)
          use module_mirror_boxes, only: LatticeCenter, do_periodic
          use module_multipole
          use module_debug
          implicit none

          real*8, intent(in) :: pos(3)
          real*8, intent(out) ::  e_lattice(2), phi_lattice

          integer :: k
          real*8 :: x0(3)
          complex(kfp), parameter :: ic = (zero, one)
          complex(kfp) :: z0, cphi, ce

          if (.not. do_periodic) then
            e_lattice   = 0
            phi_lattice = 0
          else
            x0        = pos - LatticeCenter
            z0        = x0(1) + ic * x0(2)

            cphi = -mu_cent(0) ! OMultipole(0, z0) = 1
            ce   = 0           ! OmultipolePrime(0, z0) = 0

            do k = 1, qTaylorFMMP
              cphi = cphi - mu_cent(k) * OMultipole(k, z0)
              ce = ce + mu_cent(k) * OMultipolePrime(k, z0)
            end do

            ! E = -grad(Phi)
            e_lattice  = -[ real(ce, kind = 8), -real(aimag(ce), kind = 8) ]
            phi_lattice = real(cphi, kind = 8)

            if (do_extrinsic_correction) then    ! extrinsic correction
              ! TODO: implement
              DEBUG_WARNING(*, 'Extrinsic correction is not yet implemented.')
              !e_lattice   = e_lattice   + box_dipole
              !phi_lattice = phi_lattice - dot_product(R, box_dipole) + quad_trace
            end if

          end if

        end subroutine fmm_periodicity_sum_lattice_force


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Scaling Operator \f$\mathcal{U}_B\f$ for Taylor coefficients
        !> @param[in] B table with Taylor coefficients
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Formal summation of \f$A\f$ over NF, ie all (9) neighbouring boxes
        !> with some overhead to avoid numerical elimination of small values
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        complex(kfp) function AstarFunc(k,l)
          use module_mirror_boxes, only: neighbour_boxes, num_neighbour_boxes,&
            lattice_vect
          use module_multipole
          implicit none
          integer, intent(in) :: k, l
          real(kfp) :: rpart(num_neighbour_boxes), ipart(num_neighbour_boxes), rp, ip, x0(3)
          complex(kfp) :: z0, tmp
          complex(kfp), parameter :: ic = (zero,one)

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


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Formal summation of \f$B\f$ over FF`, ie a lot of boxes
        !> with some overhead to avoid numerical elimination of small values
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        complex(kfp) function BstarFunc(k,l)
          use module_mirror_boxes, only: neighbour_boxes, num_neighbour_boxes,&
            lattice_vect
          use module_multipole
          implicit none
          integer, intent(in) :: k, l
          real(kfp) :: rpart(num_neighbour_boxes*(num_neighbour_boxes-1)), ipart(num_neighbour_boxes*(num_neighbour_boxes-1)), rp, ip, x0(3)
          complex(kfp) :: z0, tmp
          complex(kfp), parameter :: ic = (zero,one)
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


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Writes contents of table T to file s, 1 index
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Writes contents of table T to file s
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Dumps all parameters to the stream ifile
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine fmm_periodicity_param_dump(ifile)
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
        end subroutine fmm_periodicity_param_dump


        function matmul_careful(A, B)
          implicit none

          complex(kfp), intent(in) :: A(:, :), B(:, :)
          complex(kfp) :: matmul_careful(size(A, dim = 1), size(B, dim = 2))

          complex(kfp), parameter :: ic = (zero, one)
          complex(kfp) :: tmp
          real(kfp) :: rpart(size(A, dim = 2)), ipart(size(A, dim = 2)), rp, ip
          integer :: i, j, k

          do i = 1, size(A, dim = 1)
            do j = 1, size(B, dim = 2)

              rpart = 0
              ipart = 0

              do k = 1, size(A, dim = 2)
                tmp = A(i,k) * B(k, j)

                rpart(k) = real(tmp, kind = kfp)
                ipart(k) = real(aimag(tmp), kind = kfp)
              end do

              call sort_abs(rpart)
              call sort_abs(ipart)

              rp = 0
              ip = 0

              do k = 1, size(A, dim = 2)
                rp = rp + rpart(k)
                ip = ip + ipart(k)
              end do

              matmul_careful(i, j) = rp + ic * ip
            end do
          end do

        end function matmul_careful


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
        subroutine chop1D(a)
          use module_debug
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
        end subroutine chop1D


        subroutine chop2D(a)
          use module_debug
          implicit none
          complex(kfp), intent(inout) :: a(:,:)
          integer :: i, j
          complex(kfp), parameter :: ic = (zero,one)
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


end module module_fmm_periodicity
