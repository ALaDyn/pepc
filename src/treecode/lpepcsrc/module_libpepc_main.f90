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

!>
!> main internal pepc routines
!>
module module_libpepc_main
    use module_debug, only : debug_level
    use treevars, only : np_mult, interaction_list_length_factor
    use module_spacefilling, only : curve_type
    use module_domains, only: weighted, force_cubic_domain
    use module_mirror_boxes, only: mirror_box_layers

    implicit none
    private

    public libpepc_restore_particles
    public libpepc_traverse_tree
    public libpepc_grow_tree
    public libpepc_read_parameters
    public libpepc_write_parameters

    namelist /libpepc/ debug_level, np_mult, curve_type, force_cubic_domain, weighted, interaction_list_length_factor, mirror_box_layers

    contains


    !>
    !> reads libpepc-specific parameters from file
    !>
    subroutine libpepc_read_parameters(filehandle)
        use module_debug, only: pepc_status
        implicit none
        integer, intent(in) :: filehandle

        call pepc_status("READ PARAMETERS, section libpepc")
        read(filehandle,NML=libpepc)

    end subroutine libpepc_read_parameters

    !>
    !> writes libpepc-specific parameters to file
    !>
    subroutine libpepc_write_parameters(filehandle)
        use module_debug, only: pepc_status
        implicit none
        integer, intent(in) :: filehandle

        write(filehandle, NML=libpepc)

    end subroutine

    !>
    !> Builds the tree from the given particles, redistributes particles
    !> to other MPI ranks if necessary (i.e. reallocates particles and changes np_local)
    !>
    subroutine libpepc_grow_tree(t, n, p, npl)
        use module_htable, only: htable_dump
        use module_pepc_types, only: t_particle, t_tree_node
        use module_timings
        use module_tree
        use module_domains
        use module_debug
        use module_comm_data, only: t_comm_data, comm_data_create
        use module_box, only: t_box, box_create
        use treevars, only: MPI_COMM_lpepc
        use module_spacefilling, only: compute_particle_keys
        implicit none
        include 'mpif.h'

        type(t_tree), intent(inout) :: t
        integer*8, intent(in) :: n !< total number of simulation particles (across all MPI ranks)
        type(t_particle), allocatable, intent(inout) :: p(:) !< input particle data, initializes %x, %data, %work appropriately (and optionally set %label) before calling this function
        integer, optional, intent(in) :: npl !< number of valid entries in p (local particles)

        type(t_particle) :: bp(2) 
        type(t_tree_node), pointer :: root_node
        type(t_comm_data) :: tree_comm_data
        integer*8, allocatable :: branch_keys(:)

        call pepc_status('GROW TREE')
        !call MPI_BARRIER( MPI_COMM_lpepc, ierr)  ! Wait for everyone to catch up
        call timer_start(t_all)
        call timer_start(t_fields_tree)

        call comm_data_create(tree_comm_data, MPI_COMM_lpepc)

        ! determine the bounding box of the particle configuration
        call box_create(t%bounding_box, p, tree_comm_data)
        call timer_start(t_domains_keys)
        ! assign SFC coordinate to each particle
        call compute_particle_keys(t%bounding_box, p)
        call timer_stop(t_domains_keys)

        ! Domain decomposition: allocate particle keys to PEs
        if (present(npl)) then
          call domain_decompose(t%decomposition, t%bounding_box, n, p, bp, tree_comm_data, npl = npl)
        else
          call domain_decompose(t%decomposition, t%bounding_box, n, p, bp, tree_comm_data)
        end if

        ! allocate the tree structure
        call tree_create(t, size(p), n, comm_data = tree_comm_data)

        ! build local part of tree
        call timer_start(t_local)
        call tree_build_from_particles(t, p, bp)
        ! build tree from local particle keys up to root
        call tree_build_upwards(t, p(:)%key_leaf)

        call tree_get_root(t, root_node, 'libpepc_grow_tree:root node')
        if (root_node%leaves .ne. size(p)) then
          call htable_dump(t%node_storage, p)
          DEBUG_ERROR(*, 'did not find all its particles inside the tree after local tree buildup: root_node%leaves =', root_node%leaves, ' but size(particles) =', size(p))
        end if
        call timer_stop(t_local)

        ! Should now have multipole information up to root list level(s) (only up to branch level, the information is correct)
        ! By definition, this is complete: each branch node is self-contained.
        ! This information has to be broadcast to the other PEs so that the top levels can be filled in.

        call timer_start(t_exchange_branches)
        call tree_exchange_branches(t, p, bp, branch_keys)
        call timer_stop(t_exchange_branches)

        ! build global part of tree
        call timer_start(t_global)
        call tree_build_upwards(t, branch_keys)
        call timer_stop(t_global)
        deallocate(branch_keys)

        call tree_get_root(t, root_node, 'libpepc_grow_tree:root node')
        if (root_node%leaves .ne. t%npart) then
          call htable_dump(t%node_storage, p)
          DEBUG_ERROR(*, 'did not find all particles inside the htable after global tree buildup: root_node%leaves =', root_node%leaves, ' but npart_total =', t%npart)
        endif

        call timer_stop(t_fields_tree)
        call pepc_status('TREE GROWN')

    end subroutine


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Traverses the complete tree for the given particles, i.e. computes
    !> the field values at their positions. Although missing information
    !> is automatically requested from remote MPI ranks, it is important
    !> that the particle coordinates fit to the local MPI ranks domain
    !> to avoid excessive communication
    !> If field values on some regular grid are needed, they can be
    !> generated using pepc_prepare_local_grid()
    !> Otherwise, it makes sense to provide the same particles as given/returned
    !> from to pepc_grow_tree()
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine libpepc_traverse_tree(t, p)
        use module_pepc_types
        use treevars
        use module_walk
        use module_mirror_boxes
        use module_walk_communicator
        use module_timings
        use module_interaction_specific
        use module_debug
        use module_tree, only: t_tree
        implicit none

        type(t_tree), target, intent(inout) :: t
        type(t_particle), target, intent(inout) :: p(:) !< input particle data, initializes %x, %data, %work appropriately (and optionally set %label) before calling this function

        integer :: ibox
        real*8 :: ttrav, ttrav_loc, tcomm(3) ! timing integrals

        call pepc_status('TRAVERSE TREE')

        call timer_reset(t_walk)
        call timer_reset(t_walk_local)
        call timer_reset(t_comm_total)
        call timer_reset(t_comm_recv)
        call timer_reset(t_comm_sendreqs)

        comm_loop_iterations = 0

        call timer_start(t_fields_passes)

        do ibox = 1,num_neighbour_boxes ! sum over all boxes within ws=1

            call debug_barrier() ! we have to synchronize the different walks to prevent problems with recognition of finished ranks by rank 0
                                 ! just for the case that some frontend calls traverse_tree() several times, all of them have to be
                                 ! synchronized individually - hence in any case there must be a barrier here

            ! tree walk finds interaction partners and calls interaction routine for particles on short list
            call tree_walk(t, p, ttrav, ttrav_loc, lattice_vect(neighbour_boxes(:,ibox)), tcomm)

            call timer_add(t_walk, ttrav)           ! traversal time (until all walks are finished)
            call timer_add(t_walk_local, ttrav_loc) ! traversal time (local)
            call timer_add(t_comm_total,    tcomm(TIMING_COMMLOOP))
            call timer_add(t_comm_recv,     tcomm(TIMING_RECEIVE))
            call timer_add(t_comm_sendreqs, tcomm(TIMING_SENDREQS))

        end do ! ibox = 1,num_neighbour_boxes

        ! add lattice contribution
        call timer_start(t_lattice)
        ! add lattice contribution and other per-particle-forces
        ! TODO: do not call calc_force_per_particle here!
        call calc_force_per_particle(p, size(p))
        call timer_stop(t_lattice)

        call timer_stop(t_fields_passes)

        call pepc_status('TRAVERSAL DONE')
    end subroutine libpepc_traverse_tree


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Restores the initial particle distribution (before calling pepc_grow_tree() ).
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine libpepc_restore_particles(t, np_local, particles)
        use module_pepc_types, only: t_particle
        use module_timings
        use module_debug, only : pepc_status
        use module_domains, only : domain_restore
        use module_tree, only: t_tree
        implicit none

        type(t_tree), intent(inout) :: t
        integer, intent(inout) :: np_local    !< number of particles on this CPU, i.e. number of particles in particles-array
        type(t_particle), allocatable, intent(inout) :: particles(:) !< input particle data on local MPI rank - is replaced by original particle data that was given before calling pepc_grow_tree()

        call pepc_status('RESTORE PARTICLES')

        ! restore initial particle order specified by calling routine to reassign computed forces
        call timer_start(t_restore)
        call domain_restore(t%decomposition, particles)
        call timer_stop(t_restore)
      
        np_local = size(particles) ! we have got the original particle order again

        call pepc_status('RESTORATION DONE')
    end subroutine libpepc_restore_particles


end module module_libpepc_main
