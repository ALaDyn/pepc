
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

  use module_pepc
  use module_walk, only: num_walk_threads

  use module_calc_force, only: &
       mac_select, force_law

  use physvars, only: &
       trun, &
       particles, &
       weighted, &
       nt, &
       npart_total, &
       np_mult, &
       np_local, &
       n_cpu, &
       my_rank, &
       itime, &
       ispecial, &
       dt, &
       db_level, &
       curve_type

  use module_neighbour_test, only: &
       validate_n_nearest_neighbour_list, &
       draw_neighbours

  use module_sph, only: &
       sph
  
  use module_timings, only: &
       timer_start, &
       timer_stop, &
       timings_LocalOutput, &
       timings_GatherAndOutput, &
       t_tot
       
  use module_mirror_boxes, only: &
       neighbour_boxes, &
       num_neighbour_boxes

  use files, only: &
       openfiles, &
       closefiles


  implicit none

  integer :: omp_thread_num
  integer :: ierr, ifile

  ! Allocate array space for tree
  call pepc_initialize("pepc-sph", my_rank, n_cpu, .true.)

  ! Set up O/P files
  call openfiles

  ! Time stamp
  if (my_rank==0) call stamp(6,1)
  if (my_rank==0) call stamp(15,1)

  ! Each CPU gets copy of initial data
  call pepc_setup()


  ! Set the number of openmp threads.
  ! Set this only, when compiling with openmp (with !$)
  ! Set number of openmp threads to the same number as pthreads used in the walk
  !$ call omp_set_num_threads(num_walk_threads)

  ! Inform the user that openmp is used, and with how many threads
  !$OMP PARALLEL PRIVATE(omp_thread_num)
  !$ omp_thread_num = OMP_GET_THREAD_NUM()
  !$ if( (my_rank .eq. 0) .and. (omp_thread_num .eq. 0) ) write(*,*) 'Using OpenMP with', OMP_GET_NUM_THREADS(), 'threads.'
  !$OMP END PARALLEL

  ! Set up particles
  call special_start(ispecial)

  ! initialize calc force params
  mac_select  = 1 ! NN MAC
  force_law   = 5 ! NN "Interaction"

  particles(:)%work = 1._8


  ! Loop over all timesteps
  do while (itime < nt)
     itime = itime + 1
     trun  = trun  + dt

     if (my_rank==0 ) then
        ifile=6
           write(ifile,'(//a,i8,(3x,a,f12.3))') &
                ' Timestep ',itime &
                ,' total run time = ',trun 
     endif
     
     call timer_start(t_tot)
     
     call pepc_grow_and_traverse(np_local, npart_total, particles, itime, .true., .true., .true.)

     
     ! timings dump
     call timer_stop(t_tot) ! total loop time without diags

     call timings_LocalOutput(itime)
     call timings_GatherAndOutput(itime)

!     call draw_neighbours(np_local, particles, itime)


     ! do i=1, np_local
     !   write(37+my_rank,*) i, "|", particle_results(i)%neighbour_nodes(:)
     !   flush(6)
     ! end do

     call validate_n_nearest_neighbour_list(np_local, particles, &
          itime, num_neighbour_boxes, neighbour_boxes)

     call sph(np_local, particles, itime, num_neighbour_boxes, neighbour_boxes)


  end do

  ! deallocate array space for particles
  call cleanup(my_rank,n_cpu)
  
  ! Time stamp
  if (my_rank==0) call stamp(6,2)
  if (my_rank==0) call stamp(15,2)

  ! Tidy up O/P files
  call closefiles

  ! cleanup of lpepc static data
  call pepc_finalize()

end program pepce
