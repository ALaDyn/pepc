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
    public periodic_system_diagnostics
    public verifydirect
    public compute_force_direct
    public fields_on_spherical_grid
    public add_grid_particles
    public dump_grid_particles
    public lattice_diagnostics



contains

    subroutine lattice_diagnostics()
      use module_fmm_framework
      use module_interaction_specific
      use module_mirror_boxes
      use physvars
      use module_pepc_types
      implicit none
      
      integer, parameter :: ntest = 15
      integer :: i
      real*8 :: mypos(3), e(3), p
      integer, parameter :: mynp = 8
      type(t_particle), allocatable, dimension(:) :: myp
      real*8, parameter, dimension(3,ntest) :: pos = reshape([&
                                0. ,  0. ,  0. , &
                               -0.5,  0. ,  0. , &
                                0. , -0.5,  0. , &
                                0. ,  0. , -0.5, &
                               +0.5,  0. ,  0. , &
                                0. , +0.5,  0. , &
                                0. ,  0. , +0.5, &
                               -0.5, -0.5, -0.5, &
                                0.5, -0.5, -0.5, &
                               -0.5,  0.5, -0.5, &
                                0.5,  0.5, -0.5, &
                               -0.5, -0.5,  0.5, &
                                0.5, -0.5,  0.5, &
                               -0.5,  0.5,  0.5, &
                                0.5,  0.5,  0.5 &
                              ], shape(pos))

      ! Defined with origin as lattice center (will shift correctly later)
      real*8, parameter, dimension(3,mynp) :: partpos = reshape([&
                               -0.25, -0.25, -0.25, &
                                0.25, -0.25, -0.25, &
                               -0.25,  0.25, -0.25, &
                                0.25,  0.25, -0.25, &
                               -0.25, -0.25,  0.25, &
                                0.25, -0.25,  0.25, &
                               -0.25,  0.25,  0.25, &
                                0.25,  0.25,  0.25 &
                              ], shape(partpos))
      allocate(myp(mynp))
      
      do i=1,mynp
        myp(i)%x   = partpos(1:3,i) + LatticeCenter
        if (i > mynp/2) then
          myp(i)%data%q =  1.
        else
          myp(i)%data%q = -1.
        endif
      end do

      ! we first output the lattice contribution from the original particles
      write(*,*) 'Lattice contribution from original particles onto test particles'      
      do i=1,ntest 
        mypos = LatticeCenter + pos(1:3,i) * [x_plasma, y_plasma, z_plasma]                                  
        call fmm_sum_lattice_force(mypos, e, p)
        write(*,*) mypos, "|", e
      end do
      
      ! and afterwards the lattice contribution that results from the artificial test case above
      write(*,*) 'Lattice contribution from artificial test case onto test particles'      
      call fmm_framework_timestep(myp(1:mynp))
                              
       do i=1,ntest 
        call fmm_sum_lattice_force(mypos, e, p)
        mypos = LatticeCenter + pos(1:3,i) * [x_plasma, y_plasma, z_plasma]                                  
        write(*,*) mypos, "|", e
      end do
      stop
    end subroutine


    subroutine compute_force_direct(nparticles, particles, nforceparticles)
        use physvars, only : my_rank, n_cpu, MPI_COMM_PEPC
        use module_pepc_types
        use module_timings
        use module_directsum
        use module_debug, only : pepc_status
        use module_interaction_specific_types, only: t_particle_results
        implicit none
        integer(kind_particle), intent(in) :: nparticles    !< number of particles on this CPU, i.e. number of particles in particles-array
        integer(kind_particle), intent(in) :: nforceparticles    !< number of particles to compute the force for, i.e. force is computed for particles(1:nforceparticles)
        type(t_particle), intent(inout) :: particles(:) !< input particle data, initializes %x, %data appropriately (and optionally set %label) before calling this function

        integer(kind_particle) :: i
        type(t_particle_results), allocatable :: directresults(:)

        call pepc_status('PEPC-MW: DIRECTSUM')

        call timer_start(t_all)

        call directforce(particles, [(i,i=1,nforceparticles)], nforceparticles, directresults, MPI_COMM_PEPC)
        particles(1:nforceparticles)%results = directresults(1:nforceparticles)

        deallocate(directresults)

        call timer_stop(t_all)

    end subroutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !>
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine cluster_diagnostics(itime, time_fs)
        use module_pepc_types
        use physvars, only : MPI_COMM_PEPC, particles, energy, np_local, my_rank, restart, spherical_grid_Nr, spherical_grid_Ntheta, spherical_grid_Nphi, rioncluster, relectroncluster
        implicit none
        include 'mpif.h'
        integer, intent(in) :: itime
        real*8, intent(in) :: time_fs
        real*8 :: rsq, rclustersq
        integer :: nion, nboundelectrons, crit(2)
        integer(kind_particle) :: p
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
        call write_spatially_resolved_data(itime, time_fs, criterion, spherical_grid_Nr, spherical_grid_Ntheta, spherical_grid_Nphi)
	
        rioncluster = sqrt(rclustersq)

        ! calculate rms radius of (bound-)electron cluster
        rsq  = 0
        do p = 1, np_local
            if (criterion(p)) then
                rsq  = rsq  + distsq(p)
            endif
        end do
	
        call MPI_ALLREDUCE(MPI_IN_PLACE, rsq,  1, MPI_REAL8,   MPI_SUM, MPI_COMM_PEPC, ierr)
        relectroncluster = sqrt(5./3.*rsq/(1.*nboundelectrons)) ! < rms electron cluster radius

        if (my_rank == 0) then
            if (firstcall .and. .not. restart) then
                firstcall = .false.
                open(88, FILE=trim(filename),STATUS='UNKNOWN', POSITION = 'REWIND')
                write(88,'("#",8(1x,a19))') "itime", "time_fs", "r_cluster^(rms,ion)", "r_cluster^(rms,el.)", "N_ion", "N_e^(free)", "N_e^(bound)", "N_e^(r<r_cl)", "N_e^(E<0)"
            else
                open(88, FILE=trim(filename),STATUS='UNKNOWN', POSITION = 'APPEND')
            endif

            write(88,'(" ",1(1x,i19),3(1x,g19.6),5(1x,i19))') itime, time_fs, rioncluster, relectroncluster, nion, nion-nboundelectrons, nboundelectrons, crit

            close(88)
        endif

    end subroutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !>
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine periodic_system_diagnostics(itime, time_fs)
        use physvars, only : particles, np_local
        implicit none
             include 'mpif.h'
        integer, intent(in) :: itime
        real*8, intent(in) :: time_fs
        logical, dimension(1:np_local) :: criterion
        integer :: nelectrons
        real*8 :: mom(4)
        character(*), parameter :: filename = 'cluster.dat'

        !> only select electrons
        criterion = (particles(1:np_local)%data%q < 0.)

        ! output total momentum of all negatively charged particles
        call write_total_momentum('momentum_electrons.dat', itime, time_fs, criterion, mom, nelectrons)

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
        use module_pepc_types
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
                call write_particles_binary(my_rank, itime, int(np_local, kind_default), particles, filename)
                call write_frontend_parameters_to_file(filename)
            endif

            !!! write particle date as a text file
            if (ascii) then
                call write_particles_ascii(my_rank, itime, int(np_local, kind_default), particles, filename)
                call write_frontend_parameters_to_file(filename)
            endif

            !!! write particle checkpoint data using mpi-io
            if (mpiio) then
                call write_particles_mpiio(MPI_COMM_WORLD, my_rank, itime, int(np_local, kind_default), npart, particles, filename)
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
        use module_pepc_types
        implicit none
        include 'mpif.h'
        logical, intent(in) :: binary, ascii, mpiio
        integer, intent(in) :: itime_in_
        character(255) :: filename
        integer(kind_default) :: nl

        if (binary .or. ascii .or. mpiio) then

            !!! read particle date as a binary file
            if (binary) write(*,*) "read_particles(): binary mode unsupported" !call read_particles_binary(my_rank, itime, np_local, dp)

            !!! read particle date as a text file
            if (ascii)  write(*,*) "read_particles(): ascii mode unsupported" !call read_particles_ascii(my_rank, itime, np_local, dp)

            !!! read particle checkpoint data using mpi-io
            if (mpiio) then
                call read_particles_mpiio(itime_in_, MPI_COMM_WORLD, my_rank, n_cpu, itime, nl, npart_total, particles, filename)
                np_local = nl
                call read_frontend_parameters_from_file(filename)
            endif

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
        integer(kind_particle) :: i
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
        call vtk%write_headers(np_local, 0_kind_particle)
        call vtk%startpoints()
        call vtk%write_data_array("xyz", particles(:)%x(1), particles(:)%x(2), particles(:)%x(3))
        call vtk%finishpoints()
        call vtk%startpointdata()
        call vtk%write_data_array("velocity", particles(:)%data%v(1),    particles(:)%data%v(2),    particles(:)%data%v(3))
        call vtk%write_data_array("el_field", particles(:)%results%e(1), particles(:)%results%e(2), particles(:)%results%e(3))
!        call vtk%write_data_array("el_field_near", particles(:)%results%e_near(1), particles(:)%results%e_near(2), particles(:)%results%e_near(3))
!        call vtk%write_data_array("el_field_far", particles(:)%results%e_far(1), particles(:)%results%e_far(2), particles(:)%results%e_far(3))
        call vtk%write_data_array("el_pot", particles(:)%results%pot)
!        call vtk%write_data_array("el_pot_near", particles(:)%results%pot_near)
!        call vtk%write_data_array("el_pot_far", particles(:)%results%pot_far)
        call vtk%write_data_array("charge", particles(:)%data%q)
        call vtk%write_data_array("mass", particles(:)%data%m)
        call vtk%write_data_array("pelabel", particles(:)%label)
        call vtk%write_data_array("local index", [(i,i=1,np_local)])
        call vtk%write_data_array("processor", [(my_rank,i=1,np_local)])
        call vtk%write_data_array("work", particles(:)%work)
        call vtk%write_data_array("Epot", energy(1,1:np_local))
        call vtk%write_data_array("Ekin", energy(2,1:np_local))
        call vtk%write_data_array("Etot", energy(3,1:np_local))
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
        integer(kind_particle) :: p
        integer(kind_default) :: ierr

        tmp = 0.
        ncontributions = 0
        r   = 0.

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


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !>
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine write_spatially_resolved_data(itime_, time_fs, selection, NR, NTheta, NPhi)
        use physvars
        use module_units
        implicit none
        include 'mpif.h'

        integer, intent(in) :: itime_
        real*8, intent(in) :: time_fs
        logical, intent(in) :: selection(1:np_local)

        integer, intent(in) :: NR, NTheta, NPhi

        real*8 :: maxR

        integer, parameter :: RAWDATA_PX = 1
        integer, parameter :: RAWDATA_PY = 2
        integer, parameter :: RAWDATA_PZ = 3
        integer, parameter :: RAWDATA_P  = 4
        integer, parameter :: RAWDATA_N  = 5

        real*8,  dimension(RAWDATA_PX:RAWDATA_N, 1:NR, 1:NTheta, 1:NPhi) :: rawdata
        character*50, parameter :: filename = 'spatially_resolved.dat'
       
        integer(kind_particle) :: p
        integer :: ierr, iR, iTheta, iPhi, idata
        real*8 :: rspherical(3), deltaR, deltaPhi, deltaTheta

        rawdata = 0.
        maxR    = r_sphere

        deltaR     = maxr / NR
        deltaPhi   = 2*pi / NPhi
        deltaTheta =   pi / NTheta

        do p = 1,np_local
            if (selection(p)) then
                ! determine box indices
                rspherical    = cartesian_to_spherical(particles(p)%x)
                ! cartesian_to_spherical returns cos(theta) as second argument
                rspherical(2) = acos(rspherical(2))

                iR     =        floor( rspherical(1) / deltaR     )          + 1
                iPhi   = modulo(floor( rspherical(2) / deltaTheta ), NTheta) + 1
                iTheta = modulo(floor( rspherical(3) / deltaPhi   ), NPhi  ) + 1

                if (iR <= NR) then ! other particles lie outside the grid
                  rawdata(RAWDATA_PX:RAWDATA_PZ, iR, iTheta, iPhi) =             particles(p)%data%v
                  rawdata(RAWDATA_P,             iR, iTheta, iPhi) = dot_product(particles(p)%data%v, particles(p)%data%v)
                  rawdata(RAWDATA_N,             iR, iTheta, iPhi) = rawdata(RAWDATA_N, iR, iTheta, iPhi) + 1
                endif

            endif
        end do

        ! collect data globally
        if (my_rank == 0) then
          call MPI_REDUCE(MPI_IN_PLACE, rawdata, 5*NR*NTheta*NPhi, MPI_REAL8, MPI_SUM, 0, MPI_COMM_PEPC, ierr )
        else
          call MPI_REDUCE(rawdata, MPI_IN_PLACE, 5*NR*NTheta*NPhi, MPI_REAL8, MPI_SUM, 0, MPI_COMM_PEPC, ierr )
        endif

        ! output data to file
        if (my_rank == 0) then
            if (itime_ <= 1 .and. .not. restart) then
                open(87, FILE=trim(filename),STATUS='UNKNOWN', POSITION='REWIND', FORM='unformatted')
                ! write header
                write(87) nt, maxR, rioncluster, relectroncluster, NR, NTheta, NPhi
            else
                open(87, FILE=trim(filename),STATUS='UNKNOWN', POSITION='APPEND', FORM='unformatted')
            endif

            write(87) itime_, time_fs
            do iR = 1,NR
              do iTheta = 1,NTheta
                do iPhi = 1,NPhi
                  do idata = 1, 5
                    write (87) rawdata(idata, iR, iTheta, iPhi)
                  end do
                end do
              end do
            end do

            close(87)
        endif

    end subroutine write_spatially_resolved_data


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !> Converts cartesian coordinates to spherical system
        !> @param[in]  cartesian  cartesian vector [x, y, z]
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function cartesian_to_spherical(cartesian)
          implicit none
          integer, parameter :: kfp = 8
          real(kfp), parameter :: zero = 0._kfp
          real(kfp), parameter :: one  = 1._kfp
          real(kfp), parameter :: two  = 2._kfp
          real*8, intent(in)  :: cartesian(3)
          real(kfp), dimension(3) :: cartesian_to_spherical, c

          c = real(cartesian, kind=kfp)

          cartesian_to_spherical(1) = sqrt(dot_product(c, c))

          if (cartesian_to_spherical(1) == 0) then
            cartesian_to_spherical(2) = one
          else
            cartesian_to_spherical(2) = c(3) / cartesian_to_spherical(1)
          end if

          if ((c(1) == 0) .and. (c(2) == 0)) then
            cartesian_to_spherical(3) = zero
          else
            cartesian_to_spherical(3) = atan2(c(2), c(1))
          end if
        end function cartesian_to_spherical

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

        type(t_particle), intent(inout) :: particles(1:np_local)
        integer, intent(in) :: verbosity !< verbosity level: 0 - only print max. relative deviations, 1 - additionally print all. relative deviations, 2 - additionally print all. calculated forces
        integer(kind_particle), intent(in) :: np_local !< number of local particles
        integer(kind_particle), dimension(:), intent(in) :: testidx !< field with particle indices that direct force has to be computed for
        integer(kind_particle) :: ntest !< number of particles in testidx
        integer, intent(in) :: my_rank, n_cpu, comm
        real*8 :: deviation(4), deviation_max(4)
        real*8 :: field_abssum(4), field_average(4)

        integer :: ntest_total, ierr
        type(t_particle_results), target, dimension(:), allocatable :: res !< test results
        type(t_particle_results), pointer :: re
        integer(kind_particle) :: p, i

        ntest = size(testidx)

        if (my_rank ==0) write(*,'("-- DIRECT VERIFICATION --")')

        call directforce(particles, testidx, ntest, res, comm)

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

    
    subroutine add_grid_particles(particles, np_local, npart_total, rmax, my_rank, num_pe)
      use module_pepc_types
      use physvars, only : spherical_grid_Nr, spherical_grid_Ntheta, spherical_grid_Nphi, ngrid_local, ngrid_global, grid_rmax
      implicit none
      integer(kind_particle), intent(inout) :: np_local, npart_total
      type(t_particle), allocatable, intent(inout), dimension(:) :: particles
      real*8, intent(in) :: rmax
      integer, intent(in) :: my_rank, num_pe

      type(t_particle), dimension(:), allocatable :: grid, tmp
      integer :: i
      
      grid_rmax = rmax

      call create_spherical_grid(grid, ngrid_local, ngrid_global, grid_rmax, spherical_grid_Nr, spherical_grid_Ntheta, spherical_grid_Nphi, my_rank, num_pe)      
      
      npart_total = npart_total + ngrid_global
      
      if (my_rank .eq. 0) then
      
        write(*,*) 'add_grid_particles() - adding ", ngrid_local, " grid particles on rank 0'
        
        do i=1,ngrid_local
          grid(i)%label     = i
          grid(i)%data%q    = 0
          grid(i)%data%m    = huge(1._8)
          grid(i)%data%v(:) = 0.
        end do
      
        allocate(tmp(np_local))
        tmp(1:np_local) = particles(1:np_local)
        deallocate(particles)
        allocate(particles(np_local+ngrid_local))
        
        particles(            1:ngrid_local         ) = grid(:)
        particles(ngrid_local+1:ngrid_local+np_local) = tmp(:)
        
        deallocate(grid, tmp)

        np_local = np_local + ngrid_local
        
      endif
      
    end subroutine
    
    subroutine dump_grid_particles(my_rank, filename, particles, itime, time_fs)
      use module_pepc_types
      use physvars, only : spherical_grid_Nr, spherical_grid_Ntheta, spherical_grid_Nphi, ngrid_local, ngrid_global, grid_rmax
      implicit none
      type(t_particle), intent(in), dimension(:) :: particles
      integer, intent(in) :: itime, my_rank
      real*8, intent(in) :: time_fs
      character(*), intent(in) :: filename
      
      call dump_spherical_grid(my_rank, filename, itime, time_fs, particles, ngrid_local, grid_rmax, spherical_grid_Nr, spherical_grid_Ntheta, spherical_grid_Nphi)
    end subroutine

    
    
    subroutine create_spherical_grid(grid, ngrid_local, ngrid_global, rmax, Nr, Ntheta, Nphi, my_rank, num_pe)
      use module_pepc_types
      use module_units
      implicit none
      type(t_particle), dimension(:), allocatable :: grid
      integer, intent(out) :: ngrid_local, ngrid_global
      real*8, intent(in) :: rmax
      integer, intent(in) :: Nr, Ntheta, Nphi, my_rank, num_pe
      
      integer :: ir, itheta, iphi, idx
      real*8 :: r, theta, phi
      
      ngrid_global = (Nr+1)*(Ntheta+1)*(Nphi+1)
      
      if (my_rank .ne. 0) then
        ngrid_local = 0
        allocate(grid(ngrid_local))
      else
        ngrid_local = ngrid_global
        idx = 0
        allocate(grid(ngrid_local))
      
        do ir=0,Nr
          do itheta=0,Ntheta
            do iphi=0,Nphi
              idx   = idx + 1
              theta =      pi / Ntheta * itheta
              phi   = 2._8*pi / Nphi   * iphi
              r     =    rmax / Nr     * ir
	      
              grid(idx)%x       = r * [ cos(phi)*sin(theta) , sin(phi)*sin(theta), cos(theta) ]
              grid(idx)%work    = 1._8
              grid(idx)%data%q  = 1._8
	    
            end do
          end do
        end do
      end if
      
    end subroutine create_spherical_grid
    
    subroutine fields_on_spherical_grid(itime, time_fs, filename_dump, rmax, my_rank, num_pe)
      use module_pepc
      use module_pepc_types
      use physvars, only : spherical_grid_Nr, spherical_grid_Ntheta, spherical_grid_Nphi
      implicit none
      character(*), intent(in) :: filename_dump
      real*8, intent(in) :: rmax
      integer, intent(in) :: itime
      real*8, intent(in) :: time_fs
      integer, intent(in) :: my_rank, num_pe
      type(t_particle), dimension(:), allocatable :: grid
      integer :: ngrid_local, ngrid_global
      
      call create_spherical_grid(grid, ngrid_local, ngrid_global, rmax, spherical_grid_Nr, spherical_grid_Ntheta, spherical_grid_Nphi, my_rank, num_pe) 
      call pepc_particleresults_clear(grid)
      call pepc_traverse_tree(grid)
      call dump_spherical_grid(my_rank, filename_dump, itime, time_fs, grid, ngrid_local, rmax, spherical_grid_Nr, spherical_grid_Ntheta, spherical_grid_Nphi)
      
      deallocate(grid)   
      
    end subroutine fields_on_spherical_grid
    


    subroutine dump_spherical_grid(my_rank, filename, itime, time_fs, grid, ngrid, rmax, Nr, Ntheta, Nphi)
        use module_pepc_types
        use physvars, only : nt, restart, rioncluster, relectroncluster
        implicit none

        integer, intent(in) :: itime, my_rank
        real*8, intent(in) :: time_fs
        integer, intent(in) :: Nr, Ntheta, Nphi
        real*8, intent(in) :: rmax
        type(t_particle), dimension(:), intent(in) :: grid
        integer, intent(in) :: ngrid
        character(*), intent(in) :: filename
        logical, save :: dumpedgrid = .false.

        integer :: p, iR, iTheta, iPhi, idata, idx
        real*8 :: rspherical(3), deltaR, deltaPhi, deltaTheta

        if (my_rank .ne. 0) return

        ! output data to file
        if (itime <= 1 .and. .not. restart) then
            open(87, FILE=trim(filename),STATUS='UNKNOWN', POSITION = 'REWIND', FORM='unformatted')
            ! write header
            write(87) nt, rmax, rioncluster, relectroncluster, Nr, Ntheta, Nphi
        else
            open(87, FILE=trim(filename),STATUS='UNKNOWN', POSITION = 'APPEND', FORM='unformatted')
        endif
	
        idx = 0

        write(87) itime, time_fs
        do iR = 0,NR
          do iTheta = 0,NTheta
            do iPhi = 0,NPhi
	      idx = idx + 1
              write (87) grid(idx)%results%e(1:3), grid(idx)%results%pot
            end do
          end do
        end do

        close(87)
	
        if (.not. dumpedgrid) then
          dumpedgrid = .true.
          open(87, FILE=trim(filename)//"_grid.dat",STATUS='UNKNOWN', POSITION = 'REWIND', FORM='unformatted')
          ! write header
          write(87) nt, rmax, Nr, Ntheta, Nphi
          idx = 0
	  
          do iR = 0,NR
            do iTheta = 0,NTheta
              do iPhi = 0,NPhi
	        idx = idx + 1
                write (87) grid(idx)%x
              end do
            end do
          end do
	  
        endif

    end subroutine dump_spherical_grid

      
end module module_diagnostics
