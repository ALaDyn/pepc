!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>  Encapsulates ...
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_diagnostics
    implicit none
    private



    public write_total_momentum
    public write_particles
    public read_particles
    public cluster_diagnostics
    public verifydirect
    public compute_force_direct



contains

    subroutine compute_force_direct(nparticles, particles, nforceparticles)
        use physvars, only : my_rank, n_cpu, MPI_COMM_PEPC
        use module_pepc_types
        use module_timings
        use module_directsum
        use module_debug, only : pepc_status
        use module_interaction_specific_types, only: t_particle_results
        implicit none
        integer, intent(in) :: nparticles    !< number of particles on this CPU, i.e. number of particles in particles-array
        integer, intent(in) :: nforceparticles    !< number of particles to compute the force for, i.e. force is computed for particles(1:nforceparticles)
        type(t_particle), allocatable, intent(inout) :: particles(:) !< input particle data, initializes %x, %data appropriately (and optionally set %label) before calling this function

        integer :: i
        type(t_particle_results), allocatable :: directresults(:)

        call pepc_status('PEPC-MW: DIRECTSUM')

        call timer_start(t_all)

        call directforce(particles, nparticles, [(i,i=1,nforceparticles)], nforceparticles, directresults, my_rank, n_cpu, MPI_COMM_PEPC)
        particles(1:nforceparticles)%results = directresults(1:nforceparticles)

        deallocate(directresults)

        call timer_stop(t_all)

    end subroutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !>
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine cluster_diagnostics(itime, time_fs, momentum_acf)
        use physvars, only : MPI_COMM_PEPC, particles, energy, np_local, momentum_acf_from_timestep, my_rank, restart
        use module_acf
        implicit none
             include 'mpif.h'
        integer, intent(in) :: itime
        real*8, intent(in) :: time_fs
        type(acf) :: momentum_acf
        real*8 :: rsq, rclustersq
        integer :: nion, nboundelectrons, p, crit(2)
        logical, dimension(1:np_local) :: criterion
        real*8, dimension(1:np_local) :: distsq
        integer :: ierr
        real*8 :: mom(4)
        logical, save :: firstcall = .true.
        character(*), parameter :: filename = 'cluster.dat'

        ! calculate distance of particles from center [0., 0., 0.]
        do p = 1, np_local
            distsq(p) = dot_product(particles(p)%x, particles(p)%x)
        end do

        ! calculate rms radius of ion cluster
        rsq  = 0
        nion = 0
        do p = 1, np_local
            if (particles(p)%data%q > 0.) then
                rsq  = rsq  + distsq(p)
                nion = nion + 1
            endif
        end do

        call MPI_ALLREDUCE(MPI_IN_PLACE, nion, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_PEPC, ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE, rsq,  1, MPI_REAL8,   MPI_SUM, MPI_COMM_PEPC, ierr)

        rclustersq = 5./3.*rsq/(1.*nion) ! < square of rms cluster radius

        !> only select electrons that are inside the cluster or have negative total energy
        criterion = ((particles(1:np_local)%data%q < 0.) .and. ( (distsq(1:np_local) < rclustersq) .or. (energy(3,1:np_local) < 0.) ))

        crit(1) = count( (particles(1:np_local)%data%q < 0.) .and. (distsq(1:np_local) < rclustersq) )
        crit(2) = count( (particles(1:np_local)%data%q < 0.) .and. (energy(3,1:np_local) < 0.)  )

        call MPI_ALLREDUCE(MPI_IN_PLACE, crit, 2, MPI_INTEGER, MPI_SUM, MPI_COMM_PEPC, ierr)

        ! output total momentum of all negatively charged particles
        call write_total_momentum('momentum_electrons.dat', itime, time_fs, criterion, mom, nboundelectrons)

        if (itime > momentum_acf_from_timestep) then
            call momentum_acf%addval(mom(1:3)/nboundelectrons)
            call momentum_acf%to_file("momentum_electrons_Kt.dat")
        endif

        if (my_rank == 0) then
            if (firstcall .and. .not. restart) then
                firstcall = .false.
                open(88, FILE=trim(filename),STATUS='UNKNOWN', POSITION = 'REWIND')
                write(88,'("#",8(1x,a16))') "itime", "time_fs", "r_cluster^(rms)", "N_ion", "N_e^(free)", "N_e^(bound)", "N_e^(r<r_cl)", "N_e^(E<0)"
            else
                open(88, FILE=trim(filename),STATUS='UNKNOWN', POSITION = 'APPEND')
            endif

            write(88,'(" ",1(1x,i16),2(1x,g16.6),5(1x,i16))') itime, time_fs, sqrt(rclustersq), nion, nion-nboundelectrons, nboundelectrons, crit

            close(88)
        endif

    end subroutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !>
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine write_particles(allowcheckpoint)
        use physvars
        implicit none
        logical, intent(in) :: allowcheckpoint
        logical :: bin, asc, check, vtk
        integer, save :: lasttimestep = -1

        ! avoid calling this function several times per timestep
        if (lasttimestep .ne. itime) then

            bin = (idump_binary       > 0)
            if (bin) bin = ((mod(itime, idump_binary        ) == 0) .or. (itime == nt))
            asc = (idump              > 0)
            if (asc) asc = ((mod(itime, idump               ) == 0) .or. (itime == nt))
            check = (idump_checkpoint > 0) .and. (allowcheckpoint)
            if (check) check = ((mod(itime, idump_checkpoint) == 0) .or. (itime == nt))
            vtk = (idump_vtk          > 0)
            if (vtk) vtk = ((mod(itime, idump_vtk           ) == 0) .or. (itime == nt))

            call write_particles_type( bin, asc, check, vtk , itime==nt)
        endif

        lasttimestep = itime

    end subroutine write_particles


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !>
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine write_particles_type(binary, ascii, mpiio, vtk, final)
        use physvars
        use module_checkpoint
        use module_namelist
        implicit none
            include 'mpif.h'
        logical, intent(in) :: binary, ascii, mpiio, vtk, final
        integer*8 :: npart
        character(255) :: filename

        if (binary .or. ascii .or. mpiio) then

            npart = npart_total ! TODO: conversion real*4 --> real*8 will be unneccessary soon

            !!! write particle date as a binary file
            if (binary) then
                call write_particles_binary(my_rank, itime, np_local, particles, filename)
                call write_frontend_parameters_to_file(filename)
            endif

            !!! write particle date as a text file
            if (ascii) then
                call write_particles_ascii(my_rank, itime, np_local, particles, filename)
                call write_frontend_parameters_to_file(filename)
            endif

            !!! write particle checkpoint data using mpi-io
            if (mpiio) then
                call write_particles_mpiio(MPI_COMM_WORLD, my_rank, itime, np_local, npart, particles, filename)
                call write_frontend_parameters_to_file(filename)
            endif

        endif


        if (vtk) call write_particles_to_vtk(itime, final)

    end subroutine write_particles_type


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !>
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine read_particles(itime_in_)
        implicit none
        integer, intent(in) :: itime_in_

        call read_particles_type(itime_in_, .false., .false., .true.)

    end subroutine read_particles



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !>
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine read_particles_type(itime_in_, binary, ascii, mpiio)
        use physvars
        use module_checkpoint
        use module_namelist
        use module_prepare
        implicit none
            include 'mpif.h'
        logical, intent(in) :: binary, ascii, mpiio
        integer, intent(in) :: itime_in_
        integer*8 :: npart
        character(255) :: filename

        if (binary .or. ascii .or. mpiio) then

            !!! read particle date as a binary file
            if (binary) write(*,*) "read_particles(): binary mode unsupported" !call read_particles_binary(my_rank, itime, np_local, dp)

            !!! read particle date as a text file
            if (ascii)  write(*,*) "read_particles(): ascii mode unsupported" !call read_particles_ascii(my_rank, itime, np_local, dp)

            !!! read particle checkpoint data using mpi-io
            if (mpiio) then
                call read_particles_mpiio(itime_in_, MPI_COMM_WORLD, my_rank, n_cpu, itime, np_local, npart, particles, filename)
                call read_frontend_parameters_from_file(filename)
            endif

            npart_total = npart

            call pepcmw_prepare()

        endif

    end subroutine read_particles_type



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !>
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine write_particles_to_vtk(step, final)
        use physvars
        use module_vtk
        use module_units
        implicit none
        integer, intent(in) :: step
        logical, intent(in) :: final
        integer :: i
        type(vtkfile_unstructured_grid) :: vtk
        integer :: vtk_step

        if (step .eq. 0) then
            vtk_step = VTK_STEP_FIRST
        else if (step .eq. nt) then
            vtk_step = VTK_STEP_LAST
        else
            vtk_step = VTK_STEP_NORMAL
        endif

        call vtk%create_parallel("particles", step, my_rank, n_cpu, trun*unit_t0_in_fs, vtk_step)
        call vtk%write_headers(np_local, 0)
        call vtk%startpoints()
        call vtk%write_data_array("xyz", np_local, particles(:)%x(1), particles(:)%x(2), particles(:)%x(3))
        call vtk%finishpoints()
        call vtk%startpointdata()
        call vtk%write_data_array("velocity", np_local, particles(:)%data%v(1),    particles(:)%data%v(2),    particles(:)%data%v(3))
        call vtk%write_data_array("el_field", np_local, particles(:)%results%e(1), particles(:)%results%e(2), particles(:)%results%e(3))
        call vtk%write_data_array("el_pot", np_local, particles(:)%results%pot)
        call vtk%write_data_array("charge", np_local, particles(:)%data%q)
        call vtk%write_data_array("mass", np_local, particles(:)%data%m)
        call vtk%write_data_array("pelabel", np_local, particles(:)%label)
        call vtk%write_data_array("local index", np_local, [(i,i=1,np_local)])
        call vtk%write_data_array("processor", np_local, [(my_rank,i=1,np_local)])
        call vtk%write_data_array("work", np_local, particles(:)%work)
        call vtk%write_data_array("Epot", np_local, energy(1,1:np_local))
        call vtk%write_data_array("Ekin", np_local, energy(2,1:np_local))
        call vtk%write_data_array("Etot", np_local, energy(3,1:np_local))
        call vtk%finishpointdata()
        call vtk%dont_write_cells()
        call vtk%write_final()
        call vtk%close()

    end subroutine write_particles_to_vtk


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !>
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine write_total_momentum(filename, itime_, time_fs, selection, mom, ncontributions)
        use physvars
        use module_units
        implicit none
        include 'mpif.h'

        integer, intent(in) :: itime_
        character(*), intent(in) :: filename
        real*8, intent(in) :: time_fs
        logical, intent(in) :: selection(1:np_local)
        real*8, intent(out) :: mom(4)
        integer, intent(out) :: ncontributions
        real*8 :: tmp(4)
        real*8 :: r(4)
        integer :: p, ierr

        tmp = 0.
        ncontributions = 0

        do p = 1,np_local
            if (selection(p)) then
                r(1:3) = particles(p)%data%v
                r(4)   = sqrt(dot_product(r,r))
                tmp = tmp + r
                ncontributions = ncontributions + 1
            endif
        end do

        call MPI_ALLREDUCE(tmp, mom, 4, MPI_REAL8, MPI_SUM, MPI_COMM_PEPC, ierr )
        call MPI_ALLREDUCE(MPI_IN_PLACE, ncontributions, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_PEPC, ierr )

        if (my_rank == 0) then
            if (itime_ <= 1 .and. .not. restart) then
                open(87, FILE=trim(filename),STATUS='UNKNOWN', POSITION = 'REWIND')
            else
                open(87, FILE=trim(filename),STATUS='UNKNOWN', POSITION = 'APPEND')
            endif
            write(87,'(i10,5g25.12,i12)') itime_, time_fs, mom/ncontributions, ncontributions
            close(87)
        endif
    end subroutine write_total_momentum


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !>
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine verifydirect(particles, np_local, testidx, verbosity, my_rank, n_cpu, comm)
        use module_directsum
        use module_pepc_types
        implicit none
          include 'mpif.h'

        type(t_particle), intent(in) :: particles(1:np_local)
        integer, intent(in) :: verbosity !< verbosity level: 0 - only print max. relative deviations, 1 - additionally print all. relative deviations, 2 - additionally print all. calculated forces
        integer, intent(in) :: np_local !< number of local particles
        integer, dimension(:), intent(in) :: testidx !< field with particle indices that direct force has to be computed for
        integer :: ntest !< number of particles in testidx
        integer, intent(in) :: my_rank, n_cpu, comm
        real*8 :: deviation(4), deviation_max(4)
        real*8 :: field_abssum(4), field_average(4)

        integer :: i, ntest_total, ierr
        type(t_particle_results), target, dimension(:), allocatable :: res !< test results
        type(t_particle_results), pointer :: re
        integer :: p

        ntest = size(testidx)

        if (my_rank ==0) write(*,'("-- DIRECT VERIFICATION --")')

        call directforce(particles, np_local, testidx, ntest, res, my_rank, n_cpu, comm)

        deviation     = 0.
        deviation_max = 0.
        field_abssum  = 0.

        do i=1,ntest
            p = testidx(i)
            re=>res(i)
            deviation(1:3)       = abs( re%e - particles(p)%results%e )
            field_abssum(1:3)    = field_abssum(1:3)    + abs(re%e)
            deviation(4)         = abs( re%pot - particles(p)%results%pot )
            field_abssum(4)      = field_abssum(4)      + abs(re%pot)

            deviation_max  = max(deviation, deviation_max)

            if (verbosity > 1) then
                write(*,'("[",I6.6,":",I6.6,"]",3(x,F10.4), " | PEPC    ", 4(x,E20.13))') my_rank, p, particles(p)%x, particles(p)%results
                write(*,'("[",I6.6,":",I6.6,"]",33x,        " | DIRECT  ", 4(x,E20.13))') my_rank, p, re
            endif

            if (verbosity > 0) then
                write(*,'("[",I6.6,":",I6.6,"]",33x,        " | Abs.err ", 4(x,E20.13),"__")') my_rank, p, deviation
            endif
        end do

        call MPI_REDUCE(deviation_max,   deviation,             4, MPI_REAL8,   MPI_MAX, 0, comm, ierr)
        call MPI_REDUCE(field_abssum,    field_average,         4, MPI_REAL8,   MPI_SUM, 0, comm, ierr)
        call MPI_REDUCE(ntest,           ntest_total,           1, MPI_INTEGER, MPI_SUM, 0, comm, ierr)
        field_average         = field_average         / ntest_total

        if ((verbosity > -1) .and. (my_rank == 0)) then
            write(*,'("Maximum absolute deviation (ex, ey, ez, pot):             ", 4(x,E20.13))') deviation
            write(*,'("Average field values       (ex, ey, ez, pot):             ", 4(x,E20.13))') field_average
            write(*,'("Maximum relative deviation (ex, ey, ez, pot):             ", 4(x,F20.3))') deviation / field_average
        endif

        deallocate(res)

    end subroutine
end module module_diagnostics