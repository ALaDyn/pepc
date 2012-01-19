program pepc

  ! pepc modules
  use module_pepc
  use module_pepc_types
  
  ! frontend helper routines
  use helper
  implicit none
    
  ! timing variables
  real*8 :: timer(5)
      
  !!! initialize pepc library and MPI
  call pepc_initialize("pepc-mini", my_rank, n_ranks, .true.)

  root = my_rank.eq.0

  timer(1) = get_time()

  call set_parameter()
  
  call init_particles(particles)

  timer(2) = get_time()

  if(root) write(*,'(a,es12.4)') " === init time [s]: ", timer(2) - timer(1)
 
  do step=0, nt
    if(root) then
      write(*,*) " "
      write(*,'(a,i12)')    " ====== computing step  :", step
      write(*,'(a,es12.4)') " ====== simulation time :", step * dt
    end if
    
    timer(3) = get_time()
    
    call pepc_particleresults_clear(particles, np)
    call pepc_grow_tree(np, tnp, particles)
    call pepc_traverse_tree(np, particles)
    
    if(domain_output) call write_domain(particles)
    
    if(particle_probe) call compute_field()
    
    call pepc_timber_tree()
    !call pepc_restore_particles(np, particles)
    
    if(particle_test) call test_particles()  
    
    if(particle_output) call write_particles(particles)

    if(particle_filter) call filter_particles(particles)
        
    call push_particles(particles)    
    
    timer(4) = get_time()
    if(root) write(*,'(a,es12.4)') " == time in step [s]                              : ", timer(4) - timer(3)
    
  end do 
 
  deallocate(particles)

  timer(5) = get_time()

  if(root) then
    write(*,*)            " "
    write(*,'(a)')        " ===== finished pepc simulation"
    write(*,'(a,es12.4)') " ===== total run time [s]: ", timer(5) - timer(1)
  end if

  !!! cleanup pepc and MPI
  call pepc_finalize()

end program pepc

