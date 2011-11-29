!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_neighbour_test
  implicit none
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!  public type declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  integer, parameter :: num_neighbour_particles = 50
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!  subroutine-implementation  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> \brief validate list of n next neighbours as produced by pthreads based tree-walk with modified mac
  !> 
  !> compute distances of particles to all particles from all processes, blockwise and 
  !> if neccesary for all periodic copies of the sim-box.
  !> sort these distance to get n next neighbours (excluding self) and compare this list to result of
  !> the tree-based neighbour search.
  !> writes total number of not matching neighbours in lists to nn_validatation_result

  !> \author Andreas Breslau
  !> \date 2011.11.29

  !   param[in,out]  Name      Description
  !> \param[in]      npshort   number of particles in this chunk to validate neighbour lists for 
  !> \param[in]      pshort    array containing the indexes of the local particle arrays x, y, z, m ... for the particles in the actual chunk
  !> \param[in]      pass      number of actual chunk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine validate_n_nearest_neighbour_list(np_local, npart_total, particles, particle_results, &
       itime, num_neighbour_boxes, neighbour_boxes)
    
    use treetypes
    
    use physvars, only: &
         n_cpu, &
         my_rank
    
    use module_htable, only: & 
         htable, &
         key2addr
    
    implicit none
    include 'mpif.h'


    
    integer, intent(in) :: np_local    !< # particles on this CPU
    integer, intent(in) :: npart_total !< total # simulation particles
    type(t_particle), allocatable, intent(in) :: particles(:)
    type(t_particle_results), intent(in), allocatable :: particle_results(:)
    integer, intent(in) :: itime  ! timestep
    integer, intent(in) :: num_neighbour_boxes !< number of shift vectors in neighbours list (must be at least 1 since [0, 0, 0] has to be inside the list)
    integer, intent(in) :: neighbour_boxes(3, num_neighbour_boxes) ! list with shift vectors to neighbour boxes that shall be included in interaction calculation, at least [0, 0, 0] should be inside this list
    
    integer :: ierr
    integer :: actual_pe
    integer*8, allocatable :: key_buffer(:)
    real*8, allocatable :: x_buffer(:)
    real*8, allocatable :: y_buffer(:)
    real*8, allocatable :: z_buffer(:)

    integer :: local_particle_index
    integer :: remote_particle_index
    real*8, allocatable :: distances2(:,:)
    integer*8, allocatable :: keys(:,:)
    integer*8 :: tmp_key
    real*8 :: tmp_dist2
    logical :: found
    integer :: not_found
    integer, dimension(1) :: tmp_loc

    integer :: maxdist
    integer :: index_in_test_neighbour_list
    integer :: index_in_result_neighbour_list
    real*8 :: dist2
    integer :: actual_node
    integer :: actual_address

    character(100) :: filename


    integer, dimension(n_cpu) :: all_np_local
    integer :: max_np_local
    

    ! num_neighbour_particles+1 because the particle itself is stored within the lists first and removed later
    allocate( distances2( num_neighbour_particles +1, np_local), keys( num_neighbour_particles +1, np_local), STAT= ierr )

    distances2(1:num_neighbour_particles+1, 1:np_local) = huge(0._8)
    keys(1:num_neighbour_particles+1, 1:np_local) = -13
    maxdist = 1



    if( ierr .ne. 0 ) then
       write (*,*) 'allocate of fields for keys and distances to next neighbours in validate_n_nearest_neighbour_list failed in module_neighbour_test.f90'
    end if

    call MPI_ALLGATHER( np_local, 1, MPI_INTEGER, all_np_local, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr )

    max_np_local = maxval(all_np_local)
    
    allocate( key_buffer( max_np_local ), x_buffer( max_np_local ), y_buffer( max_np_local ), z_buffer( max_np_local ), STAT=ierr )
    
    if( ierr .ne. 0 ) then
       write (*,*) 'allocate of buffers in validate_n_nearest_neighbour_list failed in module_neighbour_test.f90'
    end if

    
    do actual_pe =0 , n_cpu-1                                  ! each process sends data

       if( actual_pe .eq. my_rank ) then
          
          key_buffer( 1:np_local ) = particles( 1:np_local )%key 
          x_buffer( 1:np_local )   = particles( 1:np_local )%x(1)
          y_buffer( 1:np_local )   = particles( 1:np_local )%x(2)
          z_buffer( 1:np_local )   = particles( 1:np_local )%x(3)

          call MPI_BCAST( key_buffer, all_np_local( actual_pe+1 ), MPI_INTEGER8, actual_pe, MPI_COMM_WORLD, ierr )
          call MPI_BCAST( x_buffer,   all_np_local( actual_pe+1 ), MPI_REAL8,    actual_pe, MPI_COMM_WORLD, ierr )
          call MPI_BCAST( y_buffer,   all_np_local( actual_pe+1 ), MPI_REAL8,    actual_pe, MPI_COMM_WORLD, ierr )
          call MPI_BCAST( z_buffer,   all_np_local( actual_pe+1 ), MPI_REAL8,    actual_pe, MPI_COMM_WORLD, ierr )

          write(*,*), all_np_local( actual_pe +1)
          write(*,'(30(O30,x))'), particles( 1:np_local )%key
          write(*,'(30(O30,x))'), key_buffer(1:np_local )

          write(*,*) 'testing key2addr'
          call flush
          write(*,*) 'key2addr: ', key2addr( particles(1)%key, 'key2addr test' )

          
          do local_particle_index = 1, np_local
             
             do remote_particle_index = 1, all_np_local( actual_pe+1 )
                
                dist2 = ( x_buffer( remote_particle_index ) - particles( local_particle_index )%x(1) ) **2 &
                     + ( y_buffer( remote_particle_index )  - particles( local_particle_index )%x(2) ) **2 &
                     + ( z_buffer( remote_particle_index )  - particles( local_particle_index )%x(3) ) **2
                
                if( dist2 < distances2(maxdist, local_particle_index ) ) then
                   distances2(maxdist, local_particle_index ) = dist2
                   keys(maxdist, local_particle_index ) = key_buffer( remote_particle_index )
                   tmp_loc = maxloc(distances2(1:num_neighbour_particles+1, local_particle_index ) )
                   maxdist = tmp_loc(1)        !< this is needed because maxloc returns an array
                end if
                
             end do
             
          end do
          
       end if

    end do

    ! now distances2 and keys contain the num_neighbour_particles closest neighbours and the particle itself

    ! move particle self to end of list and ignore in the following
    do local_particle_index = 1, np_local
       ! find closest particle, should be self
       tmp_loc = minloc(distances2(1:num_neighbour_particles+1, local_particle_index ) )
       ! move this particle to the end of the list
       tmp_dist2 = distances2( tmp_loc(1), local_particle_index )
       distances2( tmp_loc(1), local_particle_index ) = distances2( num_neighbour_particles+1, local_particle_index )
       distances2( num_neighbour_particles+1, local_particle_index ) = tmp_dist2

       tmp_key = keys( tmp_loc(1), local_particle_index )
       keys( tmp_loc(1), local_particle_index ) = keys( num_neighbour_particles+1, local_particle_index )
       keys( num_neighbour_particles+1, local_particle_index ) = tmp_key
    end do

    ! now distances2 and keys contain the num_neighbour_particles closest neighbours, ignore last element
    
    do local_particle_index = 1, np_local
       if(distances2( num_neighbour_particles+1, local_particle_index ) .ne. 0. ) then
          write(*,*) 'particle self was not removed from neighbour list in validate_n_nearest_neighbour_list failed in module_neighbour_test.f90', local_particle_index
       end if
    end do

   
!     IF( tree_nn_debug ) THEN

!    write( NN_filename, '(a,i6.6,a,i6.6,a)' ) "nn_validate_", itime-1, "_", my_rank, ".list"

!        ! \bug ab: with ACCESS='APPEND' compilation failes on jugene
! !       OPEN(99, FILE=NN_filename, ACCESS='APPEND')
!    open(99, FILE= NN_filename, POSITION='APPEND')
       
    DO local_particle_index = 1, np_local
       
       write( 99, '(i6.6,3(E12.5),i6,O30)' ) local_particle_index, particles(local_particle_index)%x(1), particles(local_particle_index)%x(2), particles(local_particle_index)%x(3), particles(local_particle_index)%label, particles(local_particle_index)%key
       write(99,*) keys(1:num_neighbour_particles, local_particle_index)
    
!           DO j = 1, num_stored( i )
!              WRITE( 99, '(O30,E12.5)' ) stored_keys( i , j ), sqrt( stored_distances( i , j ) )
!           END DO
          
    !END DO
       
!       CLOSE(99)
       
    END DO

    close(99)

    ! now compare list with result of neighbour search

    not_found = 0
    
    do local_particle_index = 1, np_local
       do index_in_test_neighbour_list = 1, num_neighbour_particles   ! ignore particle self (num_neighbour_particles+1)

          write(*,'(a,O30)') 'key: ', keys(index_in_test_neighbour_list, local_particle_index)
          call flush
          actual_address = key2addr( keys(index_in_test_neighbour_list, local_particle_index), 'neighbour test' )
          actual_node = htable(actual_address)%node

          found = .false.

          do index_in_result_neighbour_list = 1, num_neighbour_particles
             if( actual_node .eq. particle_results(local_particle_index)%neighbour_nodes(index_in_result_neighbour_list) ) then
                found = .true.
                exit
             end if
          end do

          if( .not. found ) then
             
             write( filename, '(a,i6.6,a,i6.6,a)' ) "validation_", itime-1, "_", my_rank, ".errors"
             
             ! \bug ab: with ACCESS='APPEND' compilation failes on jugene
             !             OPEN( 76, FILE=filename, ACCESS='APPEND')
             open(76, FILE= filename, POSITION='Append')
             ! use format O30 for keys for octal
             write( 76 , '(a,i6,a,i19.19,a,i6)' ) "for particle ", particles(local_particle_index)%label, " the neighbour with key ", &
                  keys(index_in_test_neighbour_list, local_particle_index), " was not found in the result neighbour list"
!             WRITE( 76, *) '  ', nn_keys(1:n_nn)
!             WRITE( 76, *) '  ', xcoc( nodelist(1:n_nn,i) )
!             WRITE( 76, *) '  ', sqrt(dist2_list(1:n_nn,i))
!             WRITE( 76, *) '  ', nterm(i), r_nn( pshort( i ) ), tree_x( pshort( i))
!             WRITE( 76, '(O30,E12.5E2)') stored_keys(i, j), x_from_key
!             WRITE( 76, *) '  ', stored_keys( i , 1:n_nn)
             close( 76 )

             not_found = not_found +1

          end if

       end do

    end do

  end subroutine validate_n_nearest_neighbour_list


end module module_neighbour_test
