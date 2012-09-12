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

    timer(1) = get_time()

    call read_args()

    if (do_resume)then
        call init_after_resume()
    else
        call set_parameter()
        call init()

        call init_particles(plasma_particles)
        call init_wall_particles(wall_particles)          !wall particles and plasma particles are initialized seperately

        particles(1:npp)=plasma_particles(:)              !and than combined in one array for the pepc routines
        particles(npp+1:npp+nwp)=wall_particles(:)
    end if

    !probes for analysign interaction partners
    call init_probes(6)

    timer(2) = get_time()
    if(root) write(*,'(a,es12.4)') " === init time [s]: ", timer(2) - timer(1)



    do step=startstep, nt-1+startstep
        if(root) then
            write(*,*) " "
            write(*,'(a,i12)')    " ====== computing step  :", step
            write(*,'(a,es12.4)') " ====== simulation time :", step * dt
        end if
    
        !pepc routines and timing
        timer(3) = get_time()
        call pepc_particleresults_clear(particles, np)
        call pepc_grow_tree(np, tnp, particles)
        call pepc_traverse_tree(np, particles)

        !tree traversal to get interaction partners of the probe particles
        call get_interaction_partners(6)

        if (dbg(DBG_STATS)) call pepc_statistics(step)
        call pepc_restore_particles(np, particles)
        call pepc_timber_tree()
        timer(4) = get_time()

        !diagnostics and checkpoints with timing
        !diagnostics are carried out, when fields are computed according to current positions
        !afterwards, they are moved and filtered (boundary-hits, ionization...)
        !checkpoint
        if(checkp_interval.ne.0) then
            if (MOD(step,checkp_interval)==0) then
                call set_checkpoint()
            end if
        end if
        !vtk output
        if(diag_interval.ne.0) then
            if ((MOD(step,diag_interval)==0).or.(step==nt-1+startstep)) THEN
                call write_particles(particles)
            end if
        end if
        !timing output for different processes
        !call timings_LocalOutput(step)
        call timings_GatherAndOutput(step)
        timer(5) = get_time()


        !integrator and filtering with timing
        !particles get their new velocities and positions, particles are filtered and new particles created
        call boris_nonrel(particles)

        if (open_sides) then
            !open sides (particles leave and are refluxed)
            call hits_on_boundaries_with_open_sides(particles)
        else
            !periodic sides (particles enter on other side after leaving domain)
            call hits_on_boundaries(particles)
        end if
        timer(6) = get_time()
    

        !output for root process
        if(root) write(*,'(a,es12.4)') " == time in pepc routines [s]                     : ", timer(4) - timer(3)
        if(root) write(*,'(a,es12.4)') " == time in output routines [s]                   : ", timer(5) - timer(4)
        if(root) write(*,'(a,es12.4)') " == time in integrator and particlehandling [s]   : ", timer(6) - timer(5)

    end do

    deallocate(particles)
    timer(7) = get_time()

    if(root) then
        write(*,*)            " "
        write(*,'(a)')        " ===== finished pepc simulation"
        write(*,'(a,es12.4)') " ===== total run time [s]: ", timer(7) - timer(1)
    end if

    !!! cleanup pepc and MPI
    call pepc_finalize()

end program pepc

