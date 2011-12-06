
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

  use module_setup
  use module_tree_walk, only: num_walk_threads

  use treetypes, only: &
       t_calc_force_params

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
       sph_density, &
       update_particle_props
  
  use timings, only: &
       timer_start, &
       timer_stop, &
       timings_LocalOutput, &
       timings_GatherAndOutput, &
       t_tot
       
  use module_mirror_boxes, only: &
       neighbour_boxes, &
       num_neighbour_boxes

  use module_pepcfields, only: &
       pepc_fields

  use files, only: &
       openfiles, &
       closefiles


  implicit none
  include 'mpif.h'

  integer :: omp_thread_num
  integer :: ierr, ifile
  type(t_calc_force_params) ::cf_par

  ! Allocate array space for tree
  call libpepc_setup("pepc-sph", my_rank, n_cpu)

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
  cf_par%mac         = 1 ! NN MAC
  cf_par%force_law   = 5 ! NN "Interaction"

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
     
     call MPI_BARRIER( MPI_COMM_WORLD, ierr)  ! Wait for everyone to catch up
     call timer_start(t_tot)
     
     call pepc_fields(np_local, npart_total, particles, &
          cf_par, itime, num_neighbour_boxes, neighbour_boxes, .true., .true.)
     
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

!     call sph_density(np_local, particles, itime, num_neighbour_boxes, neighbour_boxes)

     call update_particle_props(np_local, particles)


  end do

  ! deallocate array space for particles
  call cleanup(my_rank,n_cpu)
  
  ! Time stamp
  if (my_rank==0) call stamp(6,2)
  if (my_rank==0) call stamp(15,2)

  ! Tidy up O/P files
  call closefiles

  ! cleanup of lpepc static data
  call libpepc_finalize()

end program pepce
