!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>  Encapsulates anything necessary for particle setup,
!>  i.e. the old configure() and special_start() routines
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_setup
      use physvars
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
        !> rescales coordinates to fit into x_plasma x y_plasma x z_plasma cuboid
        !> centered around the coordinate origin
        !> sets plasma_centre accordingly
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine rescale_coordinates_cuboid()
          implicit none
          real*8 :: minc, maxc

          minc = minval(x(1:np_local))
          maxc = maxval(x(1:np_local))
          x(1:np_local) = ((x(1:np_local) - minc)/(maxc-minc) - 0.5) * x_plasma
          minc = minval(y(1:np_local))
          maxc = maxval(y(1:np_local))
          y(1:np_local) = ((y(1:np_local) - minc)/(maxc-minc) - 0.5) * y_plasma
          minc = minval(z(1:np_local))
          maxc = maxval(z(1:np_local))
          z(1:np_local) = ((z(1:np_local) - minc)/(maxc-minc) - 0.5) * z_plasma

          plasma_centre =  (/ 0., 0., 0./) ! Centre of plasma
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
          integer :: ierr, nep, nip

		  call MPI_SCAN(np_local, fences(my_rank), 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
		  call MPI_ALLGATHER(MPI_IN_PLACE, 1, MPI_INTEGER, fences(0), 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
		  fences(-1) = 0

          ! electrons and ions are distributed pairwise onto the mpi ranks
		  nep = np_local/2
		  nip = nep

          select case (iconf)
            case (-1)
              if (my_rank == 0) write(*,*) "Using special start... case -1 (reading mpi-io checkpoint from timestep itime_in=", itime_in ,")"
              call read_particles(itime_in)
              work(1:np_local) = 1.

            case(1)
              if (my_rank == 0) write(*,*) "Using special start... case 1 (homogeneous distribution)"
              call particle_setup_homogeneous(fences)
              call rescale_coordinates_cuboid()
              call init_generic(fences, nep, nip)
              call init_fields_zero()

            case(2)
              if (my_rank == 0) write(*,*) "Using special start... case 2 (one sphere benchmark)"
              call particle_setup_sphere(fences)
              call rescale_coordinates_spherical()
              call init_generic(fences, nep, nip)
              call init_fields_zero()

            case(3)
              if (my_rank == 0) write(*,*) "Using special start... case 3 (two sphere benchmark)"
              call particle_setup_doublesphere(fences)
              call rescale_coordinates_spherical()
              call init_generic(fences, nep, nip)
              call init_fields_zero()

            case(342)
              if (my_rank == 0) write(*,*) "Using special start... case 342 (42 spheres benchmark)"
              call particle_setup_manyspheres(fences)
              call rescale_coordinates_spherical()
              call init_generic(fences, nep, nip)
              call init_fields_zero()

            case(4)
             if (my_rank == 0) write(*,*) "Using special start... case 4: Plummer distribution (core cut)"
              call particle_setup_plummer(fences)
              call rescale_coordinates_spherical()
              call init_generic(fences, nep, nip)
              call init_fields_zero()

            case(5)
              if (my_rank == 0) write(*,*) "Using special start... case 5: 2D disc"
              call particle_setup_flatdisc(fences)
              z_plasma = 0
              call rescale_coordinates_spherical()
              call init_generic(fences, nep, nip)
              uz(1:np_local) = 0.
              call init_fields_zero()

            case(6)
              if (my_rank == 0) write(*,*) "Using special start... case 6 (2D-homogeneous distribution)"
              call particle_setup_homogeneous(fences)
              z_plasma = 0
              call rescale_coordinates_cuboid()
              call init_generic(fences, nep, nip)
              uz(1:np_local) = 0.
              call init_fields_zero()

            case(7)
		      call update_particle_numbers((nint((npart_total/8)**(1./3.))**3)*8)

              if (my_rank == 0) then
                write(*,*) "Using special start... case 7 (3D Madelung Setup)"
                write(*,*) "Total particle number must be representable by 8*k^3. Set npart_total =", npart_total
              endif

              call particle_setup_madelung()
              call rescale_coordinates_cuboid()
              call cold_start(ux,uy,uz,nppm,1,np_local)
              call init_fields_zero()

            case(8)
              if (my_rank == 0) write(*,*) "Using special start... case 8 (fast homogeneous distribution)"
              call particle_setup_homogeneous_fast(fences)
              call rescale_coordinates_cuboid()
              call init_generic(fences, nep, nip)
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
              work(1:np_local) = 1.
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
              call cold_start(ux,uy,uz,nppm,1,np_local)
              call init_fields_zero()

         end select
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
                z(p) = zt
                y(p) = yt
                x(p) = xt
              end if
            end do
          end do

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

              do while ( (xt*xt + yt*yt + zt*zt) > 0.5_8)
                xt = par_rand() - 0.5_8
                yt = par_rand() - 0.5_8
                zt = par_rand() - 0.5_8
              end do

              if ( my_rank == mpi_cnt .and. p <= np_local ) then
                z(p) = zt
                y(p) = yt
                x(p) = xt
              end if
            end do
          end do

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

              do while ( (xt*xt + yt*yt + zt*zt) > 0.5_8)
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
                z(p) = zt
                y(p) = yt
                x(p) = xt
              end if
            end do
          end do

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

              do while ( (xt*xt + yt*yt + zt*zt) > 0.5_8)
                xt = par_rand() - 0.5_8
                yt = par_rand() - 0.5_8
                zt = par_rand() - 0.5_8
              end do

              delta = 20._8*GetSphereCenter(nint(42.*par_rand())) - 10._8

              if ( my_rank == mpi_cnt .and. p <= np_local ) then
                z(p) = zt + delta(1)
                y(p) = yt + delta(2)
                x(p) = xt + delta(3)
              end if
            end do
          end do

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
                z(p) = zt
                y(p) = yt
                x(p) = xt
              end if
            end do
          end do

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

              do while ( (xt*xt + yt*yt) > 0.5_8)
                xt = par_rand() - 0.5_8
                yt = par_rand() - 0.5_8
              end do

              if ( my_rank == mpi_cnt .and. p <= np_local ) then
                z(p) = zt
                y(p) = yt
                x(p) = xt
              end if
            end do
          end do

        end subroutine


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !>
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine particle_setup_madelung()
          use module_fmm_framework
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

		            x(myidx)  = (i + 0.5)*delta(1)
		            y(myidx)  = (j + 0.5)*delta(2)
		            z(myidx)  = (k + 0.5)*delta(3)
		            ux(myidx) = 0
		            uy(myidx) = 0
		            uz(myidx) = 0

		            if (mod(i+j+k, 2) == 0) then
		              q(myidx)       = qe
		              m(myidx)       = mass_e
		              pelabel(myidx) = -globalidx
		            else
		              q(myidx)       = qi
		              m(myidx)       = mass_i
                      pelabel(myidx) = globalidx
		            end if

		          end if

		        end do
		      end do
		    end do

		    work(1:np_local) = 1.

		    if (myidx .ne. np_local) write(*,*) "ERROR in special_start(7): PE", my_rank, "set up", myidx, &
		        "particles, but np_local=", np_local, "globalidx=", globalidx, "npart_total=",npart_total

		    !write(*,*) "Particle positions: "
		    !do i=1,np_local
		    !  write(*,*) my_rank,i,x(i), y(i), z(i), q(i), pelabel(i)
		    !end do

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
             x(1) = par_rand(my_rank + 13)

             do p = 1, (fences(my_rank) - fences(my_rank-1))
                x(p) = par_rand()
                y(p) = par_rand()
                z(p) = par_rand()
             end do

        end subroutine



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !>
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine particle_setup_icosahedron(fences)
          use module_units
          use module_icosahedron
          implicit none
          integer, intent(in) :: fences(-1:n_cpu-1)
          integer :: currlayer, particletype
          integer :: p
          real*8 :: r(3), xt, yt, zt

		     do p=1,np_local/2
		       ! get particle position inside the cluster
		       r = get_particle(p + fences(my_rank-1)/2-1, currlayer, particletype)

		       ! put an ion there
		       x(np_local-p+1)   = r(1)
		       y(np_local-p+1)   = r(2)
		       z(np_local-p+1)   = r(3)
		       ux(np_local-p+1)  = 0.
		       uy(np_local-p+1)  = 0.
		       uz(np_local-p+1)  = 0.
		        q(np_local-p+1)  = qi
		        m(np_local-p+1)  = mass_i
               !pelabel(np_local-p+1)  = 1 + p + fences(my_rank-1)/2-1
               pelabel(np_local-p+1)  = 1 + particletype

		       ! and put an electron into near proximity, "random initialization" of par_rand()
		       xt = 2*pi*par_rand(my_rank + 13)*2*pi
		       yt =   pi*par_rand()
		       zt = 0.01

		       x(p) = x(np_local-p+1) + zt * cos(xt) * sin(yt)
		       y(p) = y(np_local-p+1) + zt * sin(xt) * sin(yt)
		       z(p) = z(np_local-p+1) + zt * cos(yt)

		       ! chose random velocity
		       xt = par_rand()*2*pi
		       yt = par_rand()  *pi
		       zt = par_rand()  *vte

		       ux(p)  = zt * cos(xt) * sin(yt)
		       uy(p)  = zt * cos(xt) * sin(yt)
		       uz(p)  = zt * cos(yt)
		        q(p)  = qe
		        m(p)  = mass_e
		       pelabel(p)  = -(1 + particletype)
		     end do

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

		             x(np_local-myidx+1)  = (i + 0.5)*delta(1)
		             y(np_local-myidx+1)  = (j + 0.5)*delta(2)
		             z(np_local-myidx+1)  = (k + 0.5)*delta(3)
		             ux(np_local-myidx+1) = 0
		             uy(np_local-myidx+1) = 0
		             uz(np_local-myidx+1) = 0

		             q(np_local-myidx+1) = qi
		             m(np_local-myidx+1) = mass_i
		           end if

		         end do
		       end do
		     end do

		     do i = 1,ne
		       globalidx = globalidx + 1

		       if ( mod((globalidx-1)/2,n_cpu) == my_rank) then ! distribute pairs of electron and ion, since np_local is constructed a bit weird
		         myidx = myidx + 1

		         x(np_local-myidx+1)  = par_rand()
		         y(np_local-myidx+1)  = par_rand()
		         z(np_local-myidx+1)  = par_rand()

		         ux(np_local-myidx+1) = 0
		         uy(np_local-myidx+1) = 0
		         uz(np_local-myidx+1) = 0

		         q(np_local-myidx+1) = qe
		         m(np_local-myidx+1) = mass_e
		       end if
		     end do

		     if (myidx .ne. np_local) write(*,*) "ERROR in special_start(8): PE", my_rank, "set up", myidx, &
		           "particles, but np_local=", np_local, "globalidx=", globalidx, "npart_total=",npart_total
		     if (globalidx .ne. npart_total) write(*,*) "ERROR in special_start(8): PE", my_rank, "set up globalidx=", globalidx, ", but npart_total=",npart_total

        end subroutine




        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> generic initialisation of charges, masses labels, work, labels
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine init_generic(fences, nep, nip)
          use module_velocity_setup
          implicit none
          integer, intent(in) :: fences(-1:n_cpu-1)
          integer, intent(in) :: nep, nip
          integer :: i

		  q(1:nep)                  = qe        ! plasma electrons
		  q(nep + 1:np_local)       = qi        ! plasma ions (need Z* here)
		  m(1:nep)                  = mass_e    ! electron mass
		  m(nep + 1:np_local)       = mass_i    ! ion mass
		  pelabel(1:nep)            = fences(my_rank-1) + (/(i, i = 1, nep)/)      ! Electron labels
		  pelabel(nep + 1:np_local) = ne + fences(my_rank-1) + nep + (/(i, i = 1, (fences(my_rank) - fences(my_rank-1) - nep))/) ! Ion labels

		  if (my_rank==0) write(*,*) 'Initializing particle velocities to vte =',vte,' vti =',vti
		     if (vte > 0) then
		        call maxwell3(ux,uy,uz,nppm,1,nep,vte)
		     else
		        call cold_start(ux,uy,uz,nppm,1,nep)
		     endif

		     if (vti > 0) then
		        call maxwell3(ux,uy,uz,nppm,nep+1,nip,vti)
		     else
		        call cold_start(ux,uy,uz,nppm,nep+1,nip)
		     endif

		  work(1:np_local) = 1.

        end subroutine




        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> initialisation of fields to zero
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine init_fields_zero()
          implicit none

          ex(1:np_local) = 0.
          ey(1:np_local) = 0.
          ez(1:np_local) = 0.
          pot(1:np_local) = 0.

        end subroutine




		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!>
		!> portable random number generator, see numerical recipes
		!> check for the random numbers:
		!> the first numbers should be 0.2853809, 0.2533582 and 0.2533582
		!> the parameter iseed is optional
		!>
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		function par_rand(iseed)
		  implicit none

		  integer :: idum, idum2, iy, j, k
		  integer :: iv(32)

		  real :: par_rand
		  integer, intent(in), optional :: iseed

		  integer :: IM1, IM2, IMM1, IA1, IA2, IQ1, IQ2, IR1, IR2, NTAB, NDIV
		  real    :: AM, RNMX

		  save

		  data idum, idum2 /-1, 123456789/

		  IM1 = 2147483563
		  IM2 = 2147483399
		  AM  = 1.0/IM1
		  IMM1 = IM1-1
		  IA1 = 40014
		  IA2 = 40692
		  IQ1 = 53668
		  IQ2 = 52774
		  IR1 = 12211
		  IR2 = 3791
		  NTAB = 32
		  NDIV = 1+IMM1/NTAB
		  RNMX = 1.0 - 1.2e-7

		  if (idum < 0) then

		     if (present(iseed)) then
		       idum = iseed
		     else
		       idum = 1
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

end module module_setup