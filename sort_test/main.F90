
program main

  use physvars
  implicit none
  include 'mpif.h'

  ! timing stuff

  integer :: ierr, lvisit_active

  ! Initialize the MPI system
  call MPI_INIT(ierr)

  ! Get the id number of the current task
  call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)

  ! Get the number of MPI tasks
  call MPI_COMM_size(MPI_COMM_WORLD, n_cpu, ierr)

  call openfiles       ! Set up O/P files

  call setup           ! Each CPU gets copy of initial data

  call configure       ! Set up particles

!  ------------------ VISIT ------------------
#ifdef VISIT_NBODY
  if (my_rank ==0 .and. vis_on) then
     write(*,*) "Initialzing VISIT..."

     call flvisit_nbody2_init ! Start up VISIT interface to xnbody
     call flvisit_nbody2_check_connection(lvisit_active)

     call vis_parts_nbody(1)

     call flvisit_nbody2_close ! Tidy up VISIT interface to xnbody
#else
#endif
!  ------------------ VISIT ------------------

  call sorting

  call closefiles      ! Tidy up O/P files

  ! End the MPI run
  call MPI_FINALIZE(ierr)


end program main




