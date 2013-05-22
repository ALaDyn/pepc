! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2013 Juelich Supercomputing Centre, 
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>  Encapsulates anything necessary for particle setup,
!>  i.e. the old configure() and special_start() routines
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_particle_setup
      use physvars
      use module_mirror_boxes, only : LatticeOrigin
      implicit none
      save
      private

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      public update_particle_numbers
      public particle_setup

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real*8 :: BoxDimensions(3)

    contains

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> sets ne, ni, npart_total from total number of particles
        !> total number of particles must be divisible by Zion + 1
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine update_particle_numbers(particles, total_particles)
          use module_pepc_types
          implicit none
          type(t_particle), intent(inout), allocatable :: particles(:)
          integer(kind_particle), intent(in) :: total_particles
          integer(kind_particle) :: tmp

          ! we want the total particle number to be divisible by Zion+1
          ni          = floor(total_particles / (Zion + 1.))
          ne          = ni * Zion
          npart_total = ne + ni

          ! the initial number of particles on each individual mpi rank
          ! has to be divisible by two since electrons and ions are
          ! distributed pairwise
          tmp = (ni / n_cpu)

          if (((tmp * n_cpu) .ne. ni) .and. (mod(ni, int(n_cpu, kind=kind_particle)) > my_rank)) then
              tmp = tmp + 1
          endif

          if (allocated(particles)) deallocate(particles)
          allocate(particles(tmp * (Zion+1)))

        end subroutine


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> rescales coordinates to fit into x_plasma x y_plasma x z_plasma cuboid
        !> centered around the coordinate origin
        !> sets LatticeOrigin, BoxDimensions accordingly
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine rescale_coordinates_cuboid(particles)
          use module_pepc_types
          implicit none
          type(t_particle), intent(inout) :: particles(:)
          integer(kind_particle) :: i
          real*8 :: bb(3)

          bb = [x_plasma, y_plasma, z_plasma]

          do i = 1,size(particles, kind=kind(i))
            particles(i)%x = ((particles(i)%x - LatticeOrigin)/BoxDimensions - 0.5_8) * bb
          end do

          LatticeOrigin = -0.5*bb
          BoxDimensions =      bb
        end subroutine



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !>
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine particle_setup(particles, iconf)
          use module_velocity_setup
          use module_diagnostics
          use module_pepc_types
          implicit none
          include 'mpif.h'
          type(t_particle), intent(inout), allocatable :: particles(:)
          integer, intent(in) :: iconf
          integer(kind_particle) :: fences(-1:n_cpu-1)
          integer(kind_default) :: ierr

          call update_particle_numbers(particles, ne + ni)

          call MPI_SCAN(size(particles,kind=kind_particle), fences(my_rank), 1, MPI_KIND_PARTICLE, MPI_SUM, MPI_COMM_WORLD, ierr)
          call MPI_ALLGATHER(MPI_IN_PLACE, 1, MPI_KIND_PARTICLE, fences(0), 1, MPI_KIND_PARTICLE, MPI_COMM_WORLD, ierr)
          fences(-1) = 0

          select case (iconf)
            case (-1)
              if (my_rank == 0) write(*,*) "Using special start... case -1 (reading mpi-io checkpoint from timestep itime_in=", itime_in ,")"
              call read_particles(particles, itime_in)
              restart = .true.

            case(1)
              if (my_rank == 0) write(*,*) "Using special start... case 1 (homogeneous distribution)"
              call particle_setup_homogeneous(particles, fences)
              call rescale_coordinates_cuboid(particles)
              call init_generic(particles, fences)
              call init_fields_zero(particles)

         end select

        end subroutine




        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !>
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine particle_setup_homogeneous(particles,fences)
          use module_pepc_types
          implicit none
          type(t_particle), intent(inout) :: particles(:)
          integer(kind_particle), intent(in) :: fences(-1:n_cpu-1)
          integer(kind_pe) :: mpi_cnt
          integer(kind_particle) :: p
          real*8 :: xt, yt, zt

          do mpi_cnt = 0, n_cpu-1
            do p = 1, (fences(mpi_cnt) - fences(mpi_cnt-1))

              xt = par_rand()
              yt = par_rand()
              zt = par_rand()

              if ( my_rank == mpi_cnt ) then
                particles(p)%x(1) = xt
                particles(p)%x(2) = yt
                particles(p)%x(3) = zt
              end if
            end do
          end do

          LatticeOrigin = [ 0.0,  0.0,  0.0]
          BoxDimensions = [ 1.0,  1.0,  1.0]
        end subroutine


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> generic initialisation of charges, masses labels, work, labels
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine init_generic(particles,fences)
          use module_velocity_setup
          use module_units
          use module_pepc_types
          implicit none
          type(t_particle), intent(inout) :: particles(:)
          integer(kind_particle), intent(in) :: fences(-1:n_cpu-1)
          integer(kind_particle) :: nep, nip, i
          real*8, allocatable, dimension(:) :: tmp1, tmp2, tmp3
          real*8 :: vte_real, vti_real
          
          allocate(tmp1(size(particles, kind=kind_particle)), &
                   tmp2(size(particles, kind=kind_particle)), &
                   tmp3(size(particles, kind=kind_particle)))

          ! electrons and ions are distributed pairwise onto the mpi ranks
          nip = size(particles, kind=kind_particle) / (Zion + 1)
          nep = size(particles, kind=kind_particle) - nip

          particles(1:size(particles, kind=kind_particle)-1:2)%data%q  = qe        ! plasma electrons
          particles(2:size(particles, kind=kind_particle):2)%data%q    = qi        ! plasma ions
          particles(1:size(particles, kind=kind_particle)-1:2)%data%m  = mass_e    ! electron mass
          particles(2:size(particles, kind=kind_particle):2)%data%m    = mass_i    ! ion mass
          particles(1:size(particles, kind=kind_particle)-1:2)%label   = -fences(my_rank-1) - (/(i, i = 1, nep)/)      ! Electron labels
          particles(2:size(particles, kind=kind_particle):2)%label     =  fences(my_rank-1) + (/(i, i = 1, nep)/)      ! Ion labels

          vte_real = vte
          if (Te_initial_eV > 0.) vte_real = sqrt(3.*unit_kB*Te_initial_eV/unit_Ryd_in_eV/mass_e)
          vti_real = vti
          if (Ti_initial_eV > 0.) vti_real = sqrt(3.*unit_kB*Ti_initial_eV/unit_Ryd_in_eV/mass_i)

          if (my_rank==0) write(*,*) 'INIT_GENERIC: Initializing particle velocities to vte =',vte_real,' vti =',vti_real

          ! The following routines distribute their results to array slices
          if (vte_real > 0) then
             call maxwell3(tmp1,tmp2,tmp3,size(particles, kind=kind_particle),1_kind_particle,nep,vte_real)
          else
             call cold_start(tmp1(1:nep),tmp2(1:nep),tmp3(1:nep))
          endif

          if (vti_real > 0) then
             call maxwell3(tmp1,tmp2,tmp3,size(particles, kind=kind_particle),nep+1,nip,vti_real)
          else
             call cold_start(tmp1(nep+1:nip),tmp2(nep+1:nip),tmp3(nep+1:nip))
          endif
          ! we redistribute them to the pairwise ordering
          particles(1:size(particles,kind=kind_particle)-1:2)%data%v(1) = tmp1(    1:nep)
          particles(2:size(particles,kind=kind_particle)-0:2)%data%v(1) = tmp1(nep+1:nep+nip)
          particles(1:size(particles,kind=kind_particle)-1:2)%data%v(2) = tmp2(    1:nep)
          particles(2:size(particles,kind=kind_particle)-0:2)%data%v(2) = tmp2(nep+1:nep+nip)
          particles(1:size(particles,kind=kind_particle)-1:2)%data%v(3) = tmp3(    1:nep)
          particles(2:size(particles,kind=kind_particle)-0:2)%data%v(3) = tmp3(nep+1:nep+nip)

          particles(1:size(particles,kind=kind_particle))%work = 1.
          
          deallocate(tmp1, tmp2, tmp3)

        end subroutine




        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> initialisation of fields to zero
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine init_fields_zero(particles)
          use module_pepc_types
          implicit none
          type(t_particle), intent(inout) :: particles(:)
          integer(kind_particle) :: p

          do p=1,size(particles,kind=kind(p))
            particles(p)%data%b(1:3)    = 0.
            particles(p)%results%e(1:3) = 0.
            particles(p)%results%pot    = 0.
          end do

        end subroutine




        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> portable random number generator, see numerical recipes
        !> check for the random numbers:
        !> the first numbers should be 0.2853809, 0.2533582 and 0.0934685
        !> the parameter iseed is optional
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function par_rand(iseed)
          implicit none
          real :: par_rand
          integer, intent(in), optional :: iseed

          integer, parameter :: IM1  = 2147483563
          integer, parameter :: IM2  = 2147483399
          real,    parameter :: AM   = 1.0/IM1
          integer, parameter :: IMM1 = IM1-1
          integer, parameter :: IA1  = 40014
          integer, parameter :: IA2  = 40692
          integer, parameter :: IQ1  = 53668
          integer, parameter :: IQ2  = 52774
          integer, parameter :: IR1  = 12211
          integer, parameter :: IR2  = 3791
          integer, parameter :: NTAB = 32
          integer, parameter :: NDIV = 1+IMM1/NTAB
          real,    parameter :: eps_ = 1.2e-7 ! epsilon(eps_)
          real,    parameter :: RNMX = 1.0 - eps_

          integer :: j, k
          integer, volatile, save :: idum  = -1
          integer, volatile, save :: idum2 =  123456789
          integer, volatile, save :: iy    =  0
          integer, volatile, save :: iv(NTAB)


                  if (idum <=0 .or. present(iseed)) then
                    if (present(iseed)) then
                     idum = iseed
                    else
                      if (-idum < 1) then
                      idum = 1
                      else
                        idum = -idum
                      endif
                    endif

                    idum2 = idum

            do j = NTAB+7,0,-1
              k = idum/IQ1
              idum = IA1 * (idum-k*IQ1) - k*IR1
              if (idum < 0 ) idum = idum + IM1

              if (j<NTAB) iv(j+1) = idum

            end do
            iy = iv(1)
          end if

          k = idum/IQ1
          idum = IA1 * (idum-k*IQ1) - k*IR1
          if (idum < 0) idum = idum + IM1

          k = idum2/IQ2
          idum2 = IA2 * (idum2-k*IQ2) - k*IR2
          if (idum2 < 0) idum2 = idum2 + IM2

          j = iy/NDIV + 1
          iy = iv(j)-idum2
          iv(j) = idum

          if (iy < 1) iy = iy + IMM1
          par_rand = AM*iy
          if (par_rand > RNMX) par_rand = RNMX

        end function par_rand

end module module_particle_setup
