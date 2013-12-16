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
    use module_initialization
    use particlehandling
    use integrator
    use output
    use diagnostics
    use module_cmdline
    use module_geometry
    use module_species

    implicit none
    include 'mpif.h'

    ! timing variables
    real*8 :: timer(20)
    integer :: rc,ip,imp,irp
    !real*8 :: davor,danach


    !!! initialize pepc library and MPI
    call pepc_initialize("pepc-f", my_rank, n_ranks, .true.)

    root = my_rank.eq.0


    !temporary aux variables
    diags=.false.

    timer(1) = get_time()

    call read_args()

    call init_files()

    call set_default_parameters()
    call set_parameters()
    call init_species()
    call init_boundaries()

    call init_rng()
    call init_periodicity()
    call pepc_prepare(3_kind_dim)

    if (do_resume)then
        call init_after_resume()
    else
        call init()
    end if
    call get_number_of_particles(particles)

    call init_output_arrays()
    call write_parameters()

    !get initial field configuration
    call pepc_particleresults_clear(particles)
    call pepc_grow_tree(particles)
    call pepc_traverse_tree(particles)
    if (diags) call pepc_tree_diagnostics()
    call get_number_of_particles(particles)
    call pepc_timber_tree()


    !write initial configuration
    if(checkp_interval.ne.0) then
        call set_checkpoint()
    end if
    if(vtk_interval.ne.0) then
        call write_particles(particles)
    end if

    timer(2) = get_time()
    if(root) write(*,'(a,es12.4)') " === init time [s]: ", timer(2) - timer(1)
    if(root) write(*,*)""
    if(root) write(*,*)"spiegelladung=",spiegelladung


    !MAIN LOOP ====================================================================================================
    DO step=startstep+1, nt+startstep
        timer(3) = get_time()


        ! Move particles according to electric field configuration
        call boris_nonrel(particles)
        timer(4) = get_time()


        ! Remove particles that left the simulation domain and reflux particles
        call hits_on_boundaries(particles)
        call check_for_leakage(particles)
        timer(5) = get_time()


        IF (spiegelladung==1) THEN
            allocate(all_particles(size(particles)+sum(npps(1:2))),stat=rc)
            all_particles(1:size(particles))=particles
            imp=1
            DO ip=1,size(particles)
                IF (species(particles(ip)%data%species)%physical_particle .eqv. .false.) CYCLE
                all_particles(size(particles)+imp) = particles(ip)
                all_particles(size(particles)+imp)%data%q = -all_particles(size(particles)+imp)%data%q
                all_particles(size(particles)+imp)%x(1) = -all_particles(size(particles)+imp)%x(1) + 2*dx
                all_particles(size(particles)+imp)%data%species = -1
                imp = imp+1
            END DO
        END IF

        IF (spiegelladung==2) THEN
            allocate(all_particles(size(particles)+sum(npps(0:2))),stat=rc)
            all_particles(1:size(particles))=particles
            imp=1
            DO ip=1,size(particles)
                IF (particles(ip)%data%species > 2) CYCLE
                all_particles(size(particles)+imp) = particles(ip)
                all_particles(size(particles)+imp)%x(1) = -all_particles(size(particles)+imp)%x(1)
                all_particles(size(particles)+imp)%data%species = -1
                imp = imp+1
            END DO
        END IF

        !pepc routines
        call pepc_particleresults_clear(particles)
        IF (spiegelladung/=0) call pepc_particleresults_clear(all_particles)

        !grow tree
        if(root) write(*,'(a)') " == [main loop] grow tree"
        IF (spiegelladung==0)call pepc_grow_tree(particles)
        IF (spiegelladung==0)call get_number_of_particles(particles)
        IF (spiegelladung/=0)call pepc_grow_tree(all_particles)

        IF (spiegelladung/=0) THEN
            call get_number_of_particles(all_particles)
            deallocate(particles)
            allocate(particles(sum(npps)),stat=rc)
            irp=1
            DO ip=1,size(all_particles)
                IF (all_particles(ip)%data%species==-1) CYCLE
                particles(irp)=all_particles(ip)
                irp=irp+1
            END DO
        END IF

        timer(6)=get_time() !6-5: grow_tree
        !end grow tree

        !traverse tree
        if(root) write(*,'(a)') " == [main loop] traverse tree"
        call pepc_traverse_tree(particles)

        timer(7)=get_time() !7-6: traverse_tree
        !end traverse tree


        !diagnostics
        if (diags) call pepc_tree_diagnostics()
        IF (spiegelladung/=0)call get_number_of_particles(particles)

        timer(8)=get_time()
        !end diagnostics


        !timber tree
        if(root) write(*,'(a)') " == [main loop] timber tree"
        call pepc_timber_tree()


        timer(9)=get_time() !9-8: timber
        !end timber tree


        !further diagnostics
        if (diags) call timings_GatherAndOutput(step)


        timer(10)=get_time() !10-9 + 8-7: diag
        !end further diagnostics


        !vtk and checkpoints (positions and fields after current timestep)
        if(checkp_interval.ne.0) then
            if ((MOD(step,checkp_interval)==0).or.(step==nt+startstep)) then
                call set_checkpoint()
            end if
        end if
        if(vtk_interval.ne.0) then
            if ((MOD(step,vtk_interval)==0).or.(step==nt+startstep)) THEN
                IF (spiegelladung/=0) call write_particles(all_particles)
                IF (spiegelladung==0) call write_particles(particles)
            end if
        end if
        !end vtk and checkpoints


        !output
        call main_output(out)

        timer(11)=get_time()
        if(root) then
            write(*,*) " "
            write(*,'(a,i12)')    " ====== finished computing step  : ", step
            write(*,'(a,es12.4)') " ====== simulation time          : ", step * dt
            write(*,'(a,es12.4)') " ====== run time                 : ", timer(11)-timer(1)
            write(*,*) " "
            write(*,'(a)') " ================================================================================== "
            write(*,*) " "
            write(*,*) " "
            call timing_output(timer(4)-timer(3), timer(5)-timer(4), timer(6)-timer(5), timer(7)-timer(6),&
                               timer(10)-timer(9)+timer(8)-timer(7), timer(9)-timer(8), timer(11)-timer(10), out)
            call end_of_ts_output(step,out)
        end if
        if (MOD(step-startstep,10)==0) call flush_files()
        !end output

        IF (spiegelladung/=0) deallocate(all_particles)
    end do
    !END OF MAIN LOOP ====================================================================================================

    deallocate(particles)


    timer(12) = get_time()

    if(root) then
        write(*,*)            " "
        write(*,'(a)')        " ===== finished pepc simulation"
        write(*,'(a,es12.4)') " ===== total run time [s]: ", timer(12) - timer(1)
    end if

    call close_files()

    !!! cleanup pepc and MPI
    call pepc_finalize()

end program pepc

