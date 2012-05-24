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
        !> sets ne, ni, np_local, npart_total from total number of particles
        !> total number of particles must be divisible by Zion + 1
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine update_particle_numbers(total_particles)
          implicit none
          integer, intent(in) :: total_particles
          integer :: tmp

          ! we want the total particle number to be divisible by Zion+1
          ni          = floor(total_particles / (Zion + 1.))
          ne          = ni * Zion
          npart_total = ne + ni

          ! the initial number of particles on each individual mpi rank
          ! has to be divisible by two since electrons and ions are
          ! distributed pairwise
          tmp = (ni / n_cpu)

          if (((tmp * n_cpu) .ne. ni) .and. (mod(ni, n_cpu) > my_rank)) then
              tmp = tmp + 1
          endif

          np_local = tmp * (Zion + 1)

        end subroutine


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Measures Actual Dimension and Origin of Simulation Box
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine  measure_dimensions(LatticeOrigin, BoxDimensions)
          implicit none
          include 'mpif.h'
          real*8, intent(out) :: LatticeOrigin(3), BoxDimensions(3)
          real*8 :: minc(3), maxc(3), h(3)
          integer :: ierr

          minc = [minval(particles(1:np_local)%x(1)), minval(particles(1:np_local)%x(2)), minval(particles(1:np_local)%x(3))]
          maxc = [maxval(particles(1:np_local)%x(1)), maxval(particles(1:np_local)%x(2)), maxval(particles(1:np_local)%x(3))]

          call MPI_ALLREDUCE(MPI_IN_PLACE, minc, 3, MPI_REAL8, MPI_MIN, MPI_COMM_PEPC, ierr)
          call MPI_ALLREDUCE(MPI_IN_PLACE, maxc, 3, MPI_REAL8, MPI_MAX, MPI_COMM_PEPC, ierr)

          h = maxc - minc

          where (h<=0.)
            h = 1.
          end where

          LatticeOrigin = minc
          BoxDimensions = h

        end subroutine


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> rescales coordinates to fit into x_plasma x y_plasma x z_plasma cuboid
        !> centered around the coordinate origin
        !> sets plasma_centre, LatticeOrigin, BoxDimensions accordingly
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine rescale_coordinates_cuboid()
          implicit none
          integer :: i
          real*8 :: bb(3)

          bb = [x_plasma, y_plasma, z_plasma]

          do i = 1,np_local
            particles(i)%x = ((particles(i)%x - LatticeOrigin)/BoxDimensions - 0.5_8) * bb
          end do

          plasma_centre =  (/ 0., 0., 0./) ! Centre of plasma
          LatticeOrigin = -0.5*bb
          BoxDimensions =      bb

          particle_shift_simunits = particle_shift * x_plasma
        end subroutine



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> rescales coordinates to fit into a sphere with radius r_plasma
        !> centered around the coordinate origin
        !> sets x_plasma = y_plasma = z_plasma := 2*r_sphere
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine rescale_coordinates_spherical()
          implicit none

          x_plasma = 2*r_sphere
          y_plasma = 2*r_sphere
          z_plasma = 2*r_sphere

          call rescale_coordinates_cuboid()
        end subroutine




        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !>
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine particle_setup(iconf)
          use module_icosahedron
          use module_velocity_setup
          use module_diagnostics
          implicit none
          include 'mpif.h'
          integer, intent(in) :: iconf
          integer :: fences(-1:n_cpu-1)
          integer :: ierr
          integer :: npp

          call update_particle_numbers(ne + ni)

          call MPI_SCAN(np_local, fences(my_rank), 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
          call MPI_ALLGATHER(MPI_IN_PLACE, 1, MPI_INTEGER, fences(0), 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
          fences(-1) = 0

          npp = fences(my_rank) - fences(my_rank-1)

          allocate ( particles(npp) )

          select case (iconf)
            case (-1)
              if (my_rank == 0) write(*,*) "Using special start... case -1 (reading mpi-io checkpoint from timestep itime_in=", itime_in ,")"
              call read_particles(itime_in)
              restart = .true.

            case(1)
              if (my_rank == 0) write(*,*) "Using special start... case 1 (homogeneous distribution)"
              call particle_setup_homogeneous(fences)
              call rescale_coordinates_cuboid()
              call init_generic(fences)
              call init_fields_zero()

            case(2)
              if (my_rank == 0) write(*,*) "Using special start... case 2 (one sphere benchmark)"
              call particle_setup_sphere(fences)
              call rescale_coordinates_spherical()
              call init_generic(fences)
              call init_fields_zero()

            case(3)
              if (my_rank == 0) write(*,*) "Using special start... case 3 (two sphere benchmark)"
              call particle_setup_doublesphere(fences)
              call rescale_coordinates_spherical()
              call init_generic(fences)
              call init_fields_zero()

            case(342)
              if (my_rank == 0) write(*,*) "Using special start... case 342 (42 spheres benchmark)"
              call particle_setup_manyspheres(fences)
              call rescale_coordinates_spherical()
              call init_generic(fences)
              call init_fields_zero()

            case(4)
             if (my_rank == 0) write(*,*) "Using special start... case 4: Plummer distribution (core cut)"
              call particle_setup_plummer(fences)
              call rescale_coordinates_spherical()
              call init_generic(fences)
              call init_fields_zero()

            case(5)
              if (my_rank == 0) write(*,*) "Using special start... case 5: 2D disc"
              call particle_setup_flatdisc(fences)
              z_plasma = 0
              call rescale_coordinates_spherical()
              call init_generic(fences)
              particles(1:np_local)%data%v(3) = 0.
              call init_fields_zero()

            case(6)
              if (my_rank == 0) write(*,*) "Using special start... case 6 (2D-homogeneous distribution)"
              call particle_setup_homogeneous(fences)
              z_plasma = 0
              call rescale_coordinates_cuboid()
              call init_generic(fences)
              particles(1:np_local)%data%v(3) = 0.
              call init_fields_zero()

            case(7)
              call update_particle_numbers((nint((npart_total/8)**(1./3.))**3)*8)

              if (my_rank == 0) then
                write(*,*) "Using special start... case 7 (3D Madelung Setup)"
                write(*,*) "Total particle number must be representable by 8*k^3. Set npart_total =", npart_total
              endif

              call particle_setup_madelung()
              !call rescale_coordinates_cuboid()
              call cold_start(particles(1:np_local)%data%v(1),particles(1:np_local)%data%v(2),particles(1:np_local)%data%v(3),np_local,1,np_local)
              call init_fields_zero()

            case(8)
              if (my_rank == 0) write(*,*) "Using special start... case 8 (fast homogeneous distribution)"
              call particle_setup_homogeneous_fast(fences)
              call rescale_coordinates_cuboid()
              call init_generic(fences)
              call init_fields_zero()

            case(12)
              npart_total   = get_nextlower_particles(npart_total/2)*2

              call update_particle_numbers(npart_total)
              if (my_rank == 0) then
                write(*,*) "Using special start... case 12 (Mackay Icosahedron)"
                write(*,*) "Total particle number must be two times a magic cluster number. Setting npart_total =", npart_total
              endif

              call particle_setup_icosahedron(fences)
              call rescale_coordinates_spherical()
              particles(1:np_local)%work = 1.
              call init_fields_zero()

            case(13)
              ni = nint((ne/Zion)**(1./3.))**3
              ne = ni * Zion
              call update_particle_numbers(ne + ni)

              if (my_rank == 0) then
                write(*,*) "Using special start... case 13 (homogeneous electron distribution, ionic lattice)"
                write(*,*) "Number of Ions per edge must be integer, setting ne =", ne, ", ni=", ni
              endif

              call particle_setup_ionlattice()
              call rescale_coordinates_cuboid()
              call cold_start(particles(1:np_local)%data%v(1),particles(1:np_local)%data%v(2),particles(1:np_local)%data%v(3),np_local,1,np_local)
              call init_fields_zero()

            case(14)
              if (my_rank == 0) write(*,*) "Using special start... case 14 (spherical cluster with linearly shifted electrons)"
              call particle_setup_sphere_shifted_electrons(fences, 1)
              call rescale_coordinates_spherical()
              call init_fields_zero()

            case(15)
              if (my_rank == 0) write(*,*) "Using special start... case 15 (spherical cluster with radially shifted electrons)"
              call particle_setup_sphere_shifted_electrons(fences, 2)
              call rescale_coordinates_spherical()
              call init_fields_zero()

         end select

         allocate(energy(1:3,np_local))

        end subroutine




        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !>
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine particle_setup_homogeneous(fences)
          implicit none
          integer, intent(in) :: fences(-1:n_cpu-1)
          integer :: mpi_cnt, p
          real*8 :: xt, yt, zt

          do mpi_cnt = 0, n_cpu-1
            do p = 1, (fences(mpi_cnt) - fences(mpi_cnt-1))

              xt = par_rand()
              yt = par_rand()
              zt = par_rand()

              if ( my_rank == mpi_cnt .and. p <= np_local ) then
                particles(p)%x(3) = zt
                particles(p)%x(2) = yt
                particles(p)%x(1) = xt
              end if
            end do
          end do

          LatticeOrigin = [-0.5, -0.5, -0.5]
          BoxDimensions = [ 1.0,  1.0,  1.0]
        end subroutine


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !>
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine particle_setup_sphere(fences)
          implicit none
          integer, intent(in) :: fences(-1:n_cpu-1)
          integer :: mpi_cnt, p
          real*8 :: xt, yt, zt

          do mpi_cnt = 0, n_cpu-1
            do p = 1, (fences(mpi_cnt) - fences(mpi_cnt-1))

              xt = 1.0_8
              yt = 1.0_8
              zt = 1.0_8

              do while ( (xt*xt + yt*yt + zt*zt) > 0.25_8)
                xt = par_rand() - 0.5_8
                yt = par_rand() - 0.5_8
                zt = par_rand() - 0.5_8
              end do

              if ( my_rank == mpi_cnt .and. p <= np_local ) then
                particles(p)%x(3) = zt
                particles(p)%x(2) = yt
                particles(p)%x(1) = xt
              end if
            end do
          end do

          LatticeOrigin = [-0.5, -0.5, -0.5]
          BoxDimensions = [ 1.0,  1.0,  1.0]
        end subroutine


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !>
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine particle_setup_sphere_shifted_electrons(fences, shiftmode)
          use module_units, only : pi
          use physvars, only : particle_shift
          implicit none
          integer, intent(in) :: fences(-1:n_cpu-1)
          integer, intent(in) :: shiftmode
          integer :: mpi_cnt, p
          real*8 :: xt, yt, zt, xe, ye, ze, scaleval

          ! "random" initialization of par_rand
          xt = par_rand(rngseed)

          do mpi_cnt = 0, n_cpu-1
            do p = 1, (fences(mpi_cnt) - fences(mpi_cnt-1)) / 2

              xt = 1.0_8
              yt = 1.0_8
              zt = 1.0_8

              do while ( (xt*xt + yt*yt + zt*zt) > 0.25_8)
                xt = par_rand() - 0.5_8
                yt = par_rand() - 0.5_8
                zt = par_rand() - 0.5_8
              end do

              if ( my_rank == mpi_cnt .and. p <= np_local/2 ) then

                ! distribute ion
                particles(p)%x(1:3)      = [xt, yt, zt]
                particles(p)%data%v(1:3) = [0., 0., 0.]
                particles(p)%data%q      =  qi
                particles(p)%data%m      =  mass_i
                particles(p)%label       =  fences(my_rank-1) + p

                ! distribute electron respectively
                select case (shiftmode)
                  case (1) ! linearly shifted electrons
                    xe = xt + particle_shift
                    ye = yt
                    ze = zt
                  case (2) ! radially shifted electrons
                    scaleval = sqrt(dot_product([xt, yt, zt], [xt, yt, zt]))
                    scaleval = (scaleval + particle_shift) / scaleval
                    xe = xt * scaleval
                    ye = yt * scaleval
                    ze = zt * scaleval
                  case default
                    xe = xt
                    ye = yt
                    ze = zt
                end select

                particles(np_local-p+1)%x(1:3) = [xe, ye, ze]
                particles(np_local-p+1)%data%q =  qe
                particles(np_local-p+1)%data%m =  mass_e
                particles(np_local-p+1)%label  = -fences(my_rank) + p

                ! chose random velocity
                xt = par_rand()*2*pi
                yt = par_rand()  *pi
                zt = par_rand()  *vte

                particles(p)%data%v(1:3) = [zt * cos(xt) * sin(yt), zt * cos(xt) * sin(yt), zt * cos(yt)]
              end if
            end do
          end do

          particles(1:np_local)%work = 1.

          LatticeOrigin = [-0.5, -0.5, -0.5]
          BoxDimensions = [ 1.0,  1.0,  1.0]
        end subroutine



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !>
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine particle_setup_doublesphere(fences)
          implicit none
          integer, intent(in) :: fences(-1:n_cpu-1)
          integer :: mpi_cnt, p
          real*8 :: xt, yt, zt

          do mpi_cnt = 0, n_cpu-1
            do p = 1, (fences(mpi_cnt) - fences(mpi_cnt-1))

              xt = 1.0_8
              yt = 1.0_8
              zt = 1.0_8

              do while ( (xt*xt + yt*yt + zt*zt) > 0.25_8)
                xt = par_rand() - 0.5_8
                yt = par_rand() - 0.5_8
                zt = par_rand() - 0.5_8
              end do

              if (par_rand() < 0.5) then
                xt = xt + 2.
              else
                xt = xt - 2.
              endif

              if ( my_rank == mpi_cnt .and. p <= np_local ) then
                particles(p)%x(3) = zt
                particles(p)%x(2) = yt
                particles(p)%x(1) = xt
              end if
            end do
          end do

          LatticeOrigin = [-2.5, -0.5, -0.5]
          BoxDimensions = [ 5.0,  1.0,  1.0]
        end subroutine


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !>
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine particle_setup_manyspheres(fences)
          implicit none
          integer, intent(in) :: fences(-1:n_cpu-1)
          integer :: mpi_cnt, p
          real*8 :: xt, yt, zt
          real*8 :: delta(3)

          do mpi_cnt = 0, n_cpu-1
            do p = 1, (fences(mpi_cnt) - fences(mpi_cnt-1))

              xt = 1.0_8
              yt = 1.0_8
              zt = 1.0_8

              do while ( (xt*xt + yt*yt + zt*zt) > 0.25_8)
                xt = par_rand() - 0.5_8
                yt = par_rand() - 0.5_8
                zt = par_rand() - 0.5_8
              end do

              delta = 20._8*GetSphereCenter(nint(42.*par_rand())) - 10._8

              if ( my_rank == mpi_cnt .and. p <= np_local ) then
                particles(p)%x(3) = zt + delta(1)
                particles(p)%x(2) = yt + delta(2)
                particles(p)%x(1) = xt + delta(3)
              end if
            end do
          end do

          LatticeOrigin = [-10.5, -10.5, -10.5]
          BoxDimensions = [ 21.0,  21.0,  21.0]
        end subroutine


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !>
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine particle_setup_plummer(fences)
          use module_units
          implicit none
          integer, intent(in) :: fences(-1:n_cpu-1)
          integer :: mpi_cnt, p
          real*8 :: xt, yt, zt
          real*8 :: r1, r2

          do mpi_cnt = 0, n_cpu-1
            do p = 1, (fences(mpi_cnt) - fences(mpi_cnt-1))

              r1 = 0.
              do while (.not.(r1 > 0.1 .and. r1 < 3))
                 r1 = (par_rand()**(-0.2D01/0.3D01)-0.1D01)**(-0.5D00)
              end do

              zt = (0.1D01-0.2D01*par_rand())*r1
              r2 = par_rand()
              xt = (r1**0.2D01-zt**0.2D01)**(0.5D00)*cos(0.2D01*pi*r2)
              yt = (r1**0.2D01-zt**0.2D01)**(0.5D00)*sin(0.2D01*pi*r2)

              if ( my_rank == mpi_cnt .and. p <= np_local ) then
                particles(p)%x(3) = zt
                particles(p)%x(2) = yt
                particles(p)%x(1) = xt
              end if
            end do
          end do

          call measure_dimensions(LatticeOrigin, BoxDimensions)

        end subroutine



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !>
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine particle_setup_flatdisc(fences)
          implicit none
          integer, intent(in) :: fences(-1:n_cpu-1)
          integer :: mpi_cnt, p
          real*8 :: xt, yt, zt

          do mpi_cnt = 0, n_cpu-1
            do p = 1, (fences(mpi_cnt) - fences(mpi_cnt-1))

              xt = 1.0_8
              yt = 1.0_8
              zt = 0.0_8

              do while ( (xt*xt + yt*yt) > 0.25_8)
                xt = par_rand() - 0.5_8
                yt = par_rand() - 0.5_8
              end do

              if ( my_rank == mpi_cnt .and. p <= np_local ) then
                particles(p)%x(3) = zt
                particles(p)%x(2) = yt
                particles(p)%x(1) = xt
              end if
            end do
          end do

          LatticeOrigin = [-0.5, -0.5, -0.5]
          BoxDimensions = [ 1.0,  1.0,  1.0]
        end subroutine


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !>
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine particle_setup_madelung()
          use module_mirror_boxes
          implicit none
          real*8 :: delta(3)
          integer :: i,j,k,n(3), myidx, globalidx

          n(1) = NINT((npart_total/8)**(1./3.))*2
          n(2) = NINT((npart_total/n(1)/4)**(1./2.))*2
          n(3) = npart_total/n(1)/n(2)

          if (my_rank == 0) write(*,*) "Using ", n, "particles per edge"

            delta(1) = sqrt(dot_product(t_lattice_1,t_lattice_1))/real(n(1))
            delta(2) = sqrt(dot_product(t_lattice_2,t_lattice_2))/real(n(2))
            delta(3) = sqrt(dot_product(t_lattice_3,t_lattice_3))/real(n(3))

            myidx     = 0
            globalidx = 0

            do i = 0, n(1)-1
              do j = 0, n(2)-1
                do k = 0, n(3)-1

                  globalidx = globalidx + 1

                  if ( mod((globalidx-1)/2,n_cpu) == my_rank) then ! distribute pairs of electron and ion, since np_local is constructed a bit weird
                    myidx = myidx + 1

                    particles(myidx)%x(1)  = (i + 0.5)*delta(1)
                    particles(myidx)%x(2)  = (j + 0.5)*delta(2)
                    particles(myidx)%x(3)  = (k + 0.5)*delta(3)
                    particles(myidx)%data%v(1) = 0
                    particles(myidx)%data%v(2) = 0
                    particles(myidx)%data%v(3) = 0

                    if (mod(i+j+k, 2) == 0) then
                      particles(myidx)%data%q       = qe
                      particles(myidx)%data%m       = mass_e
                      particles(myidx)%label = -globalidx
                    else
                      particles(myidx)%data%q       = qi
                      particles(myidx)%data%m       = mass_i
                      particles(myidx)%label = globalidx
                    end if

                  end if

                end do
              end do
            end do

            particles(1:np_local)%work = 1.

            if (myidx .ne. np_local) write(*,*) "ERROR in special_start(7): PE", my_rank, "set up", myidx, &
                "particles, but np_local=", np_local, "globalidx=", globalidx, "npart_total=",npart_total

            !write(*,*) "Particle positions: "
            !do i=1,np_local
            !  write(*,*) my_rank,i,x(i), y(i), z(i), q(i), pelabel(i)
            !end do

          LatticeOrigin = [ 0.0,  0.0,  0.0]
          BoxDimensions = [ 1.0,  1.0,  1.0]

        end subroutine



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !>
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine particle_setup_homogeneous_fast(fences)
          implicit none
          integer, intent(in) :: fences(-1:n_cpu-1)
          integer :: p

             ! initialize random number generator with some arbitrary seed
             particles(1)%x(1) = par_rand(my_rank + 13)

             do p = 1, (fences(my_rank) - fences(my_rank-1))
                particles(p)%x(1) = par_rand()
                particles(p)%x(2) = par_rand()
                particles(p)%x(3) = par_rand()
             end do

          LatticeOrigin = [-0.5, -0.5, -0.5]
          BoxDimensions = [ 1.0,  1.0,  1.0]
        end subroutine



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !>
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine particle_setup_icosahedron(fences)
          use module_units
          use module_icosahedron
          use physvars, only : particle_shift
          implicit none
          integer, intent(in) :: fences(-1:n_cpu-1)
          integer :: currlayer, particletype
          integer :: p
          real*8 :: r(3), xt, yt, zt

             ! "random" initialization of par_rand
             xt = par_rand(my_rank + rngseed)

             do p=1,np_local/2
               ! get particle position inside the cluster
               r = get_particle(p + fences(my_rank-1)/2-1, currlayer, particletype)

               ! put an ion there
               particles(np_local-p+1)%x(1)   = r(1)
               particles(np_local-p+1)%x(2)   = r(2)
               particles(np_local-p+1)%x(3)   = r(3)
               particles(np_local-p+1)%data%v(1)  = 0.
               particles(np_local-p+1)%data%v(2)  = 0.
               particles(np_local-p+1)%data%v(3)  = 0.
               particles(np_local-p+1)%data%q  = qi
               particles(np_local-p+1)%data%m  = mass_i
               !particles(np_local-p+1)%label  = 1 + p + fences(my_rank-1)/2-1
               particles(np_local-p+1)%label  = 1 + particletype

               ! and put an electron into near proximity
               xt = 2*pi*par_rand()
               yt =   pi*par_rand()
               zt = particle_shift*par_rand()

               particles(p)%x(1) = particles(np_local-p+1)%x(1) + zt * cos(xt) * sin(yt)
               particles(p)%x(2) = particles(np_local-p+1)%x(2) + zt * sin(xt) * sin(yt)
               particles(p)%x(3) = particles(np_local-p+1)%x(3) + zt * cos(yt)

               ! chose random velocity
               xt = par_rand()*2*pi
               yt = par_rand()  *pi
               zt = par_rand()  *vte

               particles(p)%data%v(1)  = zt * cos(xt) * sin(yt)
               particles(p)%data%v(2)  = zt * cos(xt) * sin(yt)
               particles(p)%data%v(3)  = zt * cos(yt)
               particles(p)%data%q  = qe
               particles(p)%data%m  = mass_e
               particles(p)%label  = -(1 + particletype)
             end do

             call measure_dimensions(LatticeOrigin, BoxDimensions)

        end subroutine



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !>
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine particle_setup_ionlattice()
          implicit none
          integer :: n, i, j, k
          real*8 :: delta(3)
          integer :: globalidx, myidx

             n = NINT(ni**(1./3.))
             delta(1) = 1./real(n)
             delta(2) = 1./real(n)
             delta(3) = 1./real(n)

             myidx     = 0
             globalidx = 0

             do i = 0, n-1
               do j = 0, n-1
                 do k = 0, n-1

                   globalidx = globalidx + 1

                   if ( mod((globalidx-1)/2,n_cpu) == my_rank) then ! distribute pairs of electron and ion, since np_local is constructed a bit weird
                     myidx = myidx + 1

                     particles(np_local-myidx+1)%x(1)  = (i + 0.5)*delta(1)
                     particles(np_local-myidx+1)%x(2)  = (j + 0.5)*delta(2)
                     particles(np_local-myidx+1)%x(3)  = (k + 0.5)*delta(3)
                     particles(np_local-myidx+1)%data%v(1) = 0
                     particles(np_local-myidx+1)%data%v(2) = 0
                     particles(np_local-myidx+1)%data%v(3) = 0

                     particles(np_local-myidx+1)%data%q = qi
                     particles(np_local-myidx+1)%data%m = mass_i
                   end if

                 end do
               end do
             end do

             do i = 1,ne
               globalidx = globalidx + 1

               if ( mod((globalidx-1)/2,n_cpu) == my_rank) then ! distribute pairs of electron and ion, since np_local is constructed a bit weird
                 myidx = myidx + 1

                 particles(np_local-myidx+1)%x(1)  = par_rand()
                 particles(np_local-myidx+1)%x(2)  = par_rand()
                 particles(np_local-myidx+1)%x(3)  = par_rand()

                 particles(np_local-myidx+1)%data%v(1) = 0
                 particles(np_local-myidx+1)%data%v(2) = 0
                 particles(np_local-myidx+1)%data%v(3) = 0

                 particles(np_local-myidx+1)%data%q = qe
                 particles(np_local-myidx+1)%data%m = mass_e
               end if
             end do

             if (myidx .ne. np_local) write(*,*) "ERROR in special_start(8): PE", my_rank, "set up", myidx, &
                   "particles, but np_local=", np_local, "globalidx=", globalidx, "npart_total=",npart_total
             if (globalidx .ne. npart_total) write(*,*) "ERROR in special_start(8): PE", my_rank, "set up globalidx=", globalidx, ", but npart_total=",npart_total

          LatticeOrigin = [ 0.0,  0.0,  0.0]
          BoxDimensions = [ 1.0,  1.0,  1.0]

        end subroutine




        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> generic initialisation of charges, masses labels, work, labels
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine init_generic(fences)
          use module_velocity_setup
          implicit none
          integer, intent(in) :: fences(-1:n_cpu-1)
          integer :: nep, nip
          integer :: i
          real*8, dimension(1:np_local) :: tmp1, tmp2, tmp3

          ! electrons and ions are distributed pairwise onto the mpi ranks
          nep = np_local / 2
          nip = nep

          particles(1:np_local-1:2)%data%q         = qe        ! plasma electrons
          particles(2:np_local:2)%data%q           = qi        ! plasma ions
          particles(1:np_local-1:2)%data%m         = mass_e    ! electron mass
          particles(2:np_local:2)%data%m           = mass_i    ! ion mass
          particles(1:np_local-1:2)%label   = -fences(my_rank-1) - (/(i, i = 1, nep)/)      ! Electron labels
          particles(2:np_local:2)%label     =  fences(my_rank-1) + (/(i, i = 1, nep)/)      ! Ion labels

          if (my_rank==0) write(*,*) 'INIT_GENERIC: Initializing particle velocities to vte =',vte,' vti =',vti

          ! The following routines distribute their results to array slices
          if (vte > 0) then
             call maxwell3(tmp1,tmp2,tmp3,np_local,1,nep,vte)
          else
             call cold_start(tmp1,tmp2,tmp3,np_local,1,nep)
          endif

          if (vti > 0) then
             call maxwell3(tmp1,tmp2,tmp3,np_local,nep+1,nip,vti)
          else
             call cold_start(tmp1,tmp2,tmp3,np_local,nep+1,nip)
          endif
          ! we redistribute them to the pairwise ordering
          particles(1:np_local-1:2)%data%v(1) = tmp1(    1:nep)
          particles(2:np_local-0:2)%data%v(1) = tmp1(nep+1:nep+nip)
          particles(1:np_local-1:2)%data%v(2) = tmp2(    1:nep)
          particles(2:np_local-0:2)%data%v(2) = tmp2(nep+1:nep+nip)
          particles(1:np_local-1:2)%data%v(3) = tmp3(    1:nep)
          particles(2:np_local-0:2)%data%v(3) = tmp3(nep+1:nep+nip)

          particles(1:np_local)%work = 1.

        end subroutine




        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> initialisation of fields to zero
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine init_fields_zero()
          implicit none

          particles(1:np_local)%results%e(1) = 0.
          particles(1:np_local)%results%e(2) = 0.
          particles(1:np_local)%results%e(3) = 0.
          particles(1:np_local)%results%pot = 0.

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


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> returns one (of 50 fixed) random 3D vectors
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function GetSphereCenter(idx)
          implicit none
          integer, intent(in) :: idx
          real*8, dimension(3) :: GetSphereCenter

          real*8, dimension(3, 50) :: random_vectors = reshape([  0.7011 ,  0.0942 ,  0.0012, &
            0.6663 ,  0.5985 ,  0.4624,&
            0.5391 ,  0.4709 ,  0.4243,&
            0.6981 ,  0.6959 ,  0.4609,&
            0.6665 ,  0.6999 ,  0.7702,&
            0.1781 ,  0.6385 ,  0.3225,&
            0.1280 ,  0.0336 ,  0.7847,&
            0.9991 ,  0.0688 ,  0.4714,&
            0.1711 ,  0.3196 ,  0.0358,&
            0.0326 ,  0.5309 ,  0.1759,&
            0.5612 ,  0.6544 ,  0.7218,&
            0.8819 ,  0.4076 ,  0.4735,&
            0.6692 ,  0.8200 ,  0.1527,&
            0.1904 ,  0.7184 ,  0.3411,&
            0.3689 ,  0.9686 ,  0.6074,&
            0.4607 ,  0.5313 ,  0.1917,&
            0.9816 ,  0.3251 ,  0.7384,&
            0.1564 ,  0.1056 ,  0.2428,&
            0.8555 ,  0.6110 ,  0.9174,&
            0.6448 ,  0.7788 ,  0.2691,&
            0.3763 ,  0.4235 ,  0.7655,&
            0.1909 ,  0.0908 ,  0.1887,&
            0.4283 ,  0.2665 ,  0.2875,&
            0.4820 ,  0.1537 ,  0.0911,&
            0.1206 ,  0.2810 ,  0.5762,&
            0.5895 ,  0.4401 ,  0.6834,&
            0.2262 ,  0.5271 ,  0.5466,&
            0.3846 ,  0.4574 ,  0.4257,&
            0.5830 ,  0.8754 ,  0.6444,&
            0.2518 ,  0.5181 ,  0.6476,&
            0.2904 ,  0.9436 ,  0.6790,&
            0.6171 ,  0.6377 ,  0.6358,&
            0.2653 ,  0.9577 ,  0.9452,&
            0.8244 ,  0.2407 ,  0.2089,&
            0.9827 ,  0.6761 ,  0.7093,&
            0.7302 ,  0.2891 ,  0.2362,&
            0.3439 ,  0.6718 ,  0.1194,&
            0.5841 ,  0.6951 ,  0.6073,&
            0.1078 ,  0.0680 ,  0.4501,&
            0.9063 ,  0.2548 ,  0.4587,&
            0.8797 ,  0.2240 ,  0.6619,&
            0.8178 ,  0.6678 ,  0.7703,&
            0.2607 ,  0.8444 ,  0.3502,&
            0.5944 ,  0.3445 ,  0.6620,&
            0.0225 ,  0.7805 ,  0.4162,&
            0.4253 ,  0.6753 ,  0.8419,&
            0.3127 ,  0.0067 ,  0.8329,&
            0.1615 ,  0.6022 ,  0.2564,&
            0.1788 ,  0.3868 ,  0.6135,&
            0.4229 ,  0.9160 ,  0.5822], shape(random_vectors))

          GetSphereCenter = 1._8*random_vectors(:, idx + 1)

        end function GetSphereCenter

end module module_particle_setup
