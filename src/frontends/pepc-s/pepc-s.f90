
! ==============================================================
!
!
!                  PEPC-S
!
!    Parallel Efficient Parallel Coulomb-solver: Single Call Version
!
!  ==============================================================
#include "fcs_fconfig.h"

module module_pepcs

  contains

    subroutine pepc(local_particles, local_max_particles, total_particles, &
                          positions, charges,  &
                          field, potentials,   &
                          virial_tensor,       &
                          lat_x, lat_y, lat_z, &
                          lat_period,          &
                          lat_corr,            &
                          eps, theta, db_level) bind(c,name='pepc')
        use treetypes
        use module_fmm_framework
        use module_mirror_boxes
        use module_pepcfields
        use module_setup
        use treevars
        implicit none
        include 'mpif.h'

        fcs_integer, intent(inout) :: local_particles, total_particles, local_max_particles
        fcs_real, intent(in) :: positions(3,local_max_particles), charges(local_max_particles)
        fcs_real, intent(out) :: field(3,local_max_particles), potentials(local_max_particles)
        fcs_real, intent(out), dimension(3,3) :: virial_tensor
        fcs_real, intent(in), dimension(3) :: lat_x, lat_y, lat_z
        fcs_integer, intent(in), dimension(3) :: lat_period
        fcs_integer, intent(in) :: lat_corr
        fcs_real, intent(in) :: theta, eps
        fcs_integer, intent(in) :: db_level

        integer, parameter :: itime = 1
        real, parameter :: np_mult_ = -45

        type(t_particle), allocatable :: particles(:)
        type(t_calc_force_params) :: cf_par
        integer :: i
        integer :: my_rank, n_cpu, ierr

        np_mult = np_mult_

        allocate(particles(local_particles))

        ! copy coordinates and charges to internal data structures
        do i=1,local_particles
            particles(i)%x(1:3) = positions(1:3,i)
            particles(i)%work   = 1._8
            particles(i)%label  = i ! this is just for debugging purposes
            particles(i)%data%q = charges(i)
        end do

        ! initialize calc force params
        cf_par%theta2      = theta**2
        cf_par%mac         = 0 ! Barnes-Hut MAC, everything else leads to an N^2 code
        cf_par%eps2        = eps**2
        cf_par%force_const = 1.
        cf_par%force_law   = 3 ! 3D-Coulomb

        ! =============================================================
        ! TODO: do this in some scafacos_init function instead
        ! of in every timestep
        ! Get the id number of the current task
        call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
        ! Get the number of MPI tasks
        call MPI_COMM_size(MPI_COMM_WORLD, n_cpu, ierr)
        ! pepc-initialization
        call libpepc_setup(my_rank,n_cpu,int(db_level))
        ! initialize framework for lattice contributions (is automatically ignored if periodicity = [false, false, false]
        t_lattice_1 = lat_x
        t_lattice_2 = lat_y
        t_lattice_3 = lat_z
        periodicity(1:3) = (lat_period(1:3) == 1)
        do_extrinsic_correction = (lat_corr == 1)
        call fmm_framework_init(my_rank, wellsep = 1)
        ! =============================================================

        call pepc_fields(local_particles, total_particles, particles, &
                cf_par, itime, num_neighbour_boxes, neighbour_boxes, .false., .false.)

        ! read fields and potentials from internal data structures
        do i=1,local_particles
            field(1:3,i)  = particles(i)%results%e
            potentials(i) = particles(i)%results%pot
        end do

        virial_tensor = 0. !TODO
        if (my_rank==0) write(*,*) "TODO: Virial unsupported in this version of pepc"

        ! =============================================================
        ! TODO: do this in some scafacos_finalize function instead
        ! of in every timestep
        ! finalize framework for lattice contributions
        call fmm_framework_finalize()
        ! cleanup of lpepc static data
        call libpepc_finalize()
        ! =============================================================


        deallocate(particles)

    end subroutine

end module
