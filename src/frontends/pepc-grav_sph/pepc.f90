! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
!
! Copyright (C) 2002-2023 Juelich Supercomputing Centre,
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

! ==============================================================
!
!
!                  PEPC-SPH
!
!    Parallel Efficient Parallel Coulomb-solver: Electrostatics
!
!   $ Revision $
!
!   Driver code for Coulomb-solver library lpepc
!
!  Major changes:
!   June 2005: Development begun
!
!   See README.compile for summary of program units
!  ==============================================================

program pepce

   ! TODO: use the openmp module only, when compiling with openmp: !$ use omp_lib
   ! module_deps has to be changed to remove "!$" when using openmp
   ! TODO: use omp_lib, only: ...
   use omp_lib

   use treevars, only: &
      num_threads

   use module_interaction_specific, only: &
      mac_select, force_law

   use module_nn, only: &
      nn_prepare_particleresults

   use module_interaction_specific_types, only: &
      PARTICLE_TYPE_FIXED

   use physvars
! , only: &
!        trun, &
!        particles, &
!        nt, &
!        npart_total, &
!        np_local, &
!        n_cpu, &
!        my_rank, &
!        itime, &
!        ispecial, &
!        dt, &
!        idim

   use module_neighbour_test, only: &
      validate_n_nearest_neighbour_list, &
      draw_neighbours

   use module_sph, only: &
      sph, &
      sph_kernel_tests

   use module_timings, only: &
      timer_start, &
      timer_stop, &
      timings_LocalOutput, &
      timings_GatherAndOutput, &
      t_tot

   use module_mirror_boxes, only: &
      neighbour_boxes, &
      num_neighbour_boxes, &
      do_periodic, &
      constrain_periodic, &
      periodicity

   use module_pepc, only: &
      pepc_initialize, &
      pepc_read_parameters_from_first_argument, &
      pepc_finalize, &
      pepc_grow_tree, &
      pepc_traverse_tree, &
      pepc_prepare, &
      pepc_timber_tree, &
      global_tree, &
      pepc_particleresults_clear

   use files, only: &
      openfiles, &
      closefiles, &
      write_particles, &
      process_commandline_arguments, &
      run_type

   use particle_pusher, only: &
      velocities, &
      push

   use module_particle_setup, only: &
      particle_setup

   use mpi
   implicit none

   integer :: omp_thread_num
   integer :: ierr, ifile, i

   logical            :: param_file_avail
   character(len=255) :: param_file_name
   logical            :: data_file_avail
   character(len=255) :: data_file_name

   ! Allocate array space for tree
   call pepc_initialize("pepc-sph", my_rank, n_cpu, .true.)

   call process_commandline_arguments()

   ! Set up O/P files
   call openfiles

   ! Time stamp
   if (my_rank .eq. 0) call stamp(6, 1)
   if (my_rank .eq. 0) call stamp(15, 1)

   ! Each CPU gets copy of initial data
   ! Read in parameter file and allocate particles array
   call pepc_setup()

   ! Set the number of openmp threads.
   ! Set this only, when compiling with openmp (with !$)
   ! Set number of openmp threads to the same number as pthreads used in the walk
!$ call omp_set_num_threads(num_threads)

   ! Inform the user that openmp is used, and with how many threads
!$OMP PARALLEL PRIVATE(omp_thread_num)
!$ omp_thread_num = OMP_GET_THREAD_NUM()
!$ if ((my_rank .eq. 0) .and. (omp_thread_num .eq. 0)) write (*, *) 'Using OpenMP with', OMP_GET_NUM_THREADS(), 'threads.'
!$OMP END PARALLEL

   if ((run_type .eq. 3) .or. (run_type .eq. 2)) ispecial = -1
   ! Set up particles
   call particle_setup(ispecial)

!  call special_start(ispecial)

   ! initialize calc force params
   mac_select = 3 ! NN MAC
   force_law = 5 ! NN "Interaction"

   particles(:)%work = 1._8

   if (run_type .eq. 3) then ! convert
      call write_particles(itime, 0.)
   else
      call write_particles(itime - 1, -0.1)
   end if

   call pepc_prepare(int(idim, 1))

   ! TODO: where should do_periodic be set?
   do_periodic = periodicity(1) .or. periodicity(2) .or. periodicity(3)

   if (my_rank .eq. 0) then
      call sph_kernel_tests(idim)
   end if

   write (*, *) "pepc run_type:", run_type
!  call flush()

   if (run_type .eq. 1) then ! normal simulation from default or with initial particle distribution
      write (*, *) "run_type = 1"

      !
      if (.not. initialized_v_minus_half) then

         call pepc_grow_tree(particles)
         np_local = size(particles, kind=kind(np_local))

         call pepc_particleresults_clear(particles) ! initialize neighbour lists etc
         call nn_prepare_particleresults(global_tree, particles) ! improve neighbour list to speedup neighbour search

         if (do_gravity) then
            ! summing gravitational forces
            mac_select = 1 ! nn-mac
            force_law = 3  ! neighbour list force law
            call pepc_traverse_tree(particles)
         end if

         if (do_sph) then
            ! search neighbours for sph
            mac_select = 3 ! nn-mac
            force_law = 5  ! neighbour list force law

            call pepc_traverse_tree(particles)

            !     call validate_n_nearest_neighbour_list(np_local, particles, itime, num_neighbour_boxes, neighbour_boxes)

            ! call draw_neighbours(np_local, particles, itime)

            call sph(np_local, particles, itime, num_neighbour_boxes, neighbour_boxes, idim)

         end if

         ! quick n dirty
         do i = 1, np_local
            particles(i)%results%sph_force = particles(i)%results%sph_force - 39.4784176 * particles(i)%results%e
         end do

         ! Integrate velocity for t-dt/2
         do i = 1, np_local
            if (btest(particles(i)%data%type, PARTICLE_TYPE_FIXED)) then
               ! don`t move this particle
            else
               particles(i)%data%v_minus_half = particles(i)%data%v - particles(i)%results%sph_force * dt / 2.
               ! this is not perfect, but better than nothing
            end if
         end do

         ! the temperature is currently not integrated here, knowing that that produces an energy error
         ! do i = 1, np_local
         !    particles(i)%data%temperature = particles(i)%data%temperature + particles(i)%results%temperature_change * dt
         !    if(particles(i)%data%temperature < 0.000000001 ) then
         !       particles(i)%data%temperature = 0.000000001
         !    end if
         ! end do

         ! TODO: calculate work for sph
         particles(:)%work = 1._8

         call pepc_timber_tree()
      end if

   else if (run_type .eq. 2) then
      write (*, *) "run_type = 2"

   else if (run_type .eq. 3) then
      write (*, *) "run_type = 3"
      itime = nt + 1

   else
      write (*, *) "Error: unknown run_type in pepc.f90."
      call MPI_ABORT(MPI_COMM_WORLD, 13, ierr)

   end if

   ! Loop over all timesteps. <= nt to write out the last timestep
   do while (itime .le. nt)

      if (my_rank .eq. 0) then
         ifile = 6
         write (ifile, '(//a,i8,(3x,a,f12.3))') &
            ' Timestep ', itime &
            , ' total run time = ', trun
      end if

      call timer_start(t_tot)

      call pepc_grow_tree(particles)
      np_local = size(particles, kind=kind(np_local))

      ! TODO: remove this debug output
      ! write(*,*) 'num_neighbour_boxes:', num_neighbour_boxes
      ! write(*,*) 'neigbour_boxes:', neighbour_boxes

      call pepc_particleresults_clear(particles) ! initialize neighbour lists etc
      call nn_prepare_particleresults(global_tree, particles) ! improve neighbour list to speedup neighbour search

      if (do_gravity) then
         ! summing gravitational forces
         mac_select = 1 ! nn-mac
         force_law = 3  ! neighbour list force law
         call pepc_traverse_tree(particles)
      end if

      if (do_sph) then
         ! search neighbours for sph
         mac_select = 3 ! nn-mac
         force_law = 5  ! neighbour list force law

         call pepc_traverse_tree(particles)

         !     call validate_n_nearest_neighbour_list(np_local, particles, itime, num_neighbour_boxes, neighbour_boxes)

         ! call draw_neighbours(np_local, particles, itime)

         call sph(np_local, particles, itime, num_neighbour_boxes, neighbour_boxes, idim)

      end if

      ! quick n dirty
      do i = 1, np_local
         particles(i)%results%sph_force = particles(i)%results%sph_force - 39.4784176 * particles(i)%results%e
      end do

      ! Because the density and force is computed according to the current particle positions and velocities, everything is written out here.
      call write_particles(itime, trun)

      ! Now the time-step is incremented and the particles are moved.
      ! The result is not written out again, so the particles_...000000...dat/vtk contains the initial particle configuration with the densities and forces.
      ! To test the integration based on the computed forces, at least 2 time-steps have to be computed.
      itime = itime + 1
      trun = trun + dt

      ! Integrator
      call velocities(1, np_local, dt)
      call push(1, np_local, dt)  ! update positions

      do i = 1, np_local
         particles(i)%data%temperature = particles(i)%data%temperature + particles(i)%results%temperature_change * dt
         if (particles(i)%data%temperature .lt. 0.000000001) then
            particles(i)%data%temperature = 0.000000001
         end if
      end do

      write (*, *) "do_periodic:", do_periodic

      ! periodic systems demand periodic boundary conditions
      if (do_periodic) call constrain_periodic(particles)

      ! timings dump
      call timer_stop(t_tot) ! total loop time without diags

      call pepc_timber_tree()
      call timings_LocalOutput(itime)
      call timings_GatherAndOutput(itime)

      ! TODO: calculate work for sph
      particles(:)%work = 1._8

   end do

   ! The final particle positions could be written out here, but the densities and forces would belong to the previous timestep.
   ! This could lead to misinterpretations in the postprocessing and is explicitly NOT done.
   ! So we knowingly throw away the result of the last time-step.
   ! To get nt timesteps written out, the time-loop is while itime <= nt

   ! deallocate array space for particles
   call cleanup(my_rank, n_cpu)

   ! Time stamp
   if (my_rank .eq. 0) call stamp(6, 2)
   if (my_rank .eq. 0) call stamp(15, 2)

   ! Tidy up O/P files
   call closefiles

   ! cleanup of lpepc static data
   call pepc_finalize()

end program pepce
