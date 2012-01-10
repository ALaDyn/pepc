program pepc

  ! pepc modules
  use module_pepc
  use module_pepc_types
  
  ! frontend helper routines
  use helper
  implicit none
    
  ! timing variables
  real*8 :: timer(5)
    
  ! loop variables
  integer :: it
  
  ! return code variable
  integer :: rc
  
  ! particle data (position, velocity, mass, charge)
  type(t_particle), allocatable :: particles(:)
 
  timer(1) = get_time()
 
  !!! initialize pepc library and MPI
  call pepc_initialize("pepc-mini", my_rank, n_ranks, .true.)

  call set_parameter()

  ! set initially number of local particles
  np = tnp / n_ranks
  if(my_rank.eq.0) np = np + MOD(tnp, n_ranks)

  allocate(particles(np), stat=rc)
  if(rc.ne.0) write(*,*) " === particle allocation error!"
  
  call init_particles(particles)

  timer(2) = get_time()

  if(my_rank.eq.0) write(*,*) " === init time [s]: ", timer(2) - timer(1)
 
  do it=0, nt
    write(*,*) " == computing step ", it
    
    timer(3) = get_time()
    call pepc_grow_tree(np, tnp, particles)
    call pepc_traverse_tree(np, particles)
    call pepc_timber_tree()
    
    ! get new particle number, no resorting
    ! reduce size by two, due to the boundary partices
    np = size(particles) - 2
    
    call push_particles(particles)
    
    if(vtk_output) call write_particles(particles, it, (it.eq.nt))

    timer(4) = get_time()
    if(my_rank.eq.0) write(*,*) " == time in step [s]: ", timer(4) - timer(3)
    
  end do 
 
  deallocate(particles)

  timer(5) = get_time()

  if(my_rank.eq.0) write(*,*) " ===== finished pepc simulation"
  if(my_rank.eq.0) write(*,*) " ===== total run time [s]: ", timer(5) - timer(1)

  !!! cleanup pepc and MPI
  call pepc_finalize()

end program pepc

