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

program pepc

    ! pepc modules
    use module_pepc
    use module_pepc_types
    use module_mirror_boxes
    use module_checkpoint
    use module_timings
    use module_debug
    use module_vtk

    ! frontend helper routines
    use helper
    use variables
    use particlehandling
    use integrator
    use output
    use diagnostics
    use module_cmdline
    use module_interaction_partners


    implicit none
    include 'mpif.h'

    ! timing variables
    real*8 :: timer(10)

    !!! initialize pepc library and MPI
    call pepc_initialize("pepc-f", my_rank, n_ranks, .true.)

    root = my_rank.eq.0


    !temporary aux variables
    diags=.false.
    interaction_partner_diags=.false.


    timer(1) = get_time()

    call read_args()

    if (do_resume)then
        call init_after_resume()
    else
        call set_default_parameters()
        call set_parameters()
        call init()

        call init_particles(plasma_particles)
        call init_wall_particles(wall_particles)          !wall particles and plasma particles are initialized seperately

        particles(1:npp)=plasma_particles(:)              !and than combined in one array for the pepc routines
        particles(npp+1:npp+nwp)=wall_particles(:)
    end if

    !probes for analysing interaction partner diags
    if (interaction_partner_diags) call init_probes(5)

    !get initial field configuration
    call pepc_particleresults_clear(particles, np)
    call pepc_grow_tree(np, tnp, particles)
    call pepc_traverse_tree(np, particles)
    if (diags) call pepc_tree_diagnostics()
    if (do_restore_particles) call pepc_restore_particles(np, particles)
    call pepc_timber_tree()


    !write initial configuration
    if(checkp_interval.ne.0) then
        call set_checkpoint()
    end if
    if(diag_interval.ne.0) then
        call write_particles(particles)
    end if

    timer(2) = get_time()
    if(root) write(*,'(a,es12.4)') " === init time [s]: ", timer(2) - timer(1)

    !MAIN LOOP ====================================================================================================
    DO step=startstep+1, nt+startstep
        timer(3) = get_time()

        ! Move particles according to electric field configuration
        call boris_nonrel(particles)
        timer(4) = get_time()

        ! Remove particles that left the simulation domain and reflux particles
        call hits_on_boundaries(particles)
        timer(5) = get_time()

        !pepc routines and timing
        call pepc_particleresults_clear(particles, np)
        call pepc_grow_tree(np, tnp, particles)
        call pepc_traverse_tree(np, particles)
        if (diags) call pepc_tree_diagnostics()
        if (do_restore_particles) call pepc_restore_particles(np, particles)
        if (interaction_partner_diags) call get_interaction_partners(5)
        call pepc_timber_tree()
        if (diags) call timings_GatherAndOutput(step)
        timer(6)=get_time()

        !diagnostics and checkpoints (positions and fields after current timestep)
        if(checkp_interval.ne.0) then
            if ((MOD(step,checkp_interval)==0).or.(step==nt+startstep)) then
                call set_checkpoint()
            end if
        end if
        if(diag_interval.ne.0) then
            if ((MOD(step,diag_interval)==0).or.(step==nt+startstep)) THEN
                call write_particles(particles)
            end if
        end if
        timer(7)=get_time()

        !output for root
        if(root) then
            write(*,*) " "
            write(*,'(a,i12)')    " ====== finished computing step  : ", step
            write(*,'(a,es12.4)') " == time in integrator [s]       : ", timer(4) - timer(3)
            write(*,'(a,es12.4)') " == time in particlehandling [s] : ", timer(5) - timer(4)
            write(*,'(a,es12.4)') " == time in pepc routines [s]    : ", timer(6) - timer(5)
            write(*,'(a,es12.4)') " == time in output routines [s]  : ", timer(7) - timer(6)
            write(*,'(a,es12.4)') " ====== simulation time          : ", step * dt
            write(*,'(a,es12.4)') " ====== run time                 : ", timer(3)-timer(1)
            write(*,*) " "
            write(*,*) " ================================================================================== "
            write(*,*) " ================================================================================== "
            write(*,*) " "
            write(*,*) " "
        end if


    end do
    !END OF MAIN LOOP ====================================================================================================

    deallocate(particles)
    timer(8) = get_time()

    if(root) then
        write(*,*)            " "
        write(*,'(a)')        " ===== finished pepc simulation"
        write(*,'(a,es12.4)') " ===== total run time [s]: ", timer(8) - timer(1)
    end if

    !!! cleanup pepc and MPI
    call pepc_finalize()

end program pepc

