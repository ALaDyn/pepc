program pepc

  ! pepc modules
  use module_pepc
  use module_pepc_types
  
  ! frontend helper routines
  use helper
  implicit none
    
  ! timing variables
  real*8 :: timer(5)
      
  timer(1) = get_time()
 
  !!! initialize pepc library and MPI
  call pepc_initialize("pepc-mini", my_rank, n_ranks, .true.)

  call set_parameter()
  
  call init_particles(particles)

  call print_main("after init")

  timer(2) = get_time()

  if(my_rank.eq.0) write(*,*) " === init time [s]: ", timer(2) - timer(1)
 
  do step=0, nt
    if(my_rank.eq.0) write(*,*) " "
    if(my_rank.eq.0) write(*,*) " ====== computing step ", step
    
    timer(3) = get_time()
    
    call print_main("before clean")
    call pepc_particleresults_clear(particles, np)
    call print_main("after clean")    
    call pepc_grow_tree(np, tnp, particles)
    call print_main("after grow")
    call pepc_traverse_tree(np, particles)
    
    call print_main("after traverse")
    
    if(domain_output) call write_domain(particles)
    
    call pepc_timber_tree()
    !call pepc_restore_particles(np, particles)
    
    if(particle_direct .gt. 0.0) call test_particles()  
    
    if(particle_output) call write_particles(particles)

    if(particle_filter) call filter_particles(particles)
        
    call push_particles(particles)    
    
    call print_main("after push")
    
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

