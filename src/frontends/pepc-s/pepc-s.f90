
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
        use module_pepc_types
        use module_fmm_framework
        use module_mirror_boxes
        use module_pepc
        use treevars
        use module_calc_force, only : theta2, eps2, mac_select, force_law
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
        theta2      = theta**2
        mac_select  = 0 ! Barnes-Hut MAC, everything else leads to an N^2 code
        eps2        = eps**2
        force_law   = 3 ! 3D-Coulomb

        ! =============================================================
        ! TODO: do this in some scafacos_init function instead
        ! of in every timestep
        ! initialize framework for lattice contributions (is automatically ignored if periodicity = [false, false, false]
        t_lattice_1 = lat_x
        t_lattice_2 = lat_y
        t_lattice_3 = lat_z
        periodicity(1:3) = (lat_period(1:3) == 1)
        do_extrinsic_correction = (lat_corr == 1)
        call pepc_prepare()
        ! =============================================================

       call pepc_grow_and_traverse(local_particles, total_particles, particles, itime)

        ! read fields and potentials from internal data structures
        do i=1,local_particles
            field(1:3,i)  = particles(i)%results%e
            potentials(i) = particles(i)%results%pot
        end do

        virial_tensor = 0. !TODO
        if (my_rank==0) write(*,*) "TODO: Virial unsupported in this version of pepc"

        deallocate(particles)

    end subroutine

end module
