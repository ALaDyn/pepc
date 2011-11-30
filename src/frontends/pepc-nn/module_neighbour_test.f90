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

    use treevars, only: &
         xmin, &
         ymin, &
         boxsize, &
         tree_nodes
    
    use physvars, only: &
         n_cpu, &
         my_rank
    
    use module_htable, only: & 
         htable, &
         key2addr
    
    use module_spacefilling, only: &
         coord_to_key_lastlevel, &
         key_to_coord

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
    integer*8 :: actual_node
    integer*8 :: node_key
    integer :: actual_address

    character(100) :: filename


    integer, dimension(n_cpu) :: all_np_local
    integer :: max_np_local
    

    ! variables for gle output
    character(40) :: cfile
    character(50) :: outfile    
    integer :: actual_neighbour
    integer*8 :: actual_key
    real*8 :: neighbour_x
    real*8 :: neighbour_y
    real*8 :: neighbour_z
    real*8 :: smoothing_length
    character(12), parameter :: colors(0:15) = (/"orange      ", "cyan        ", "magenta     ", "blue        ", "green       ", &
         "red         ","yellow      ","grey20      ","brown       ", "yellowgreen ", "violet      ", "royalblue   ", &
         "plum        ", "goldenrod   ",  "powderblue  ", "lime        "/) 

    ! end variable declaration


    ! num_neighbour_particles+1 because the particle itself is stored within the lists first and removed later
    allocate( distances2( num_neighbour_particles +1, np_local), keys( num_neighbour_particles +1, np_local), STAT= ierr )

    if( ierr .ne. 0 ) then
       write (*,*) 'allocate of fields for keys and distances to next neighbours in validate_n_nearest_neighbour_list failed in module_neighbour_test.f90'
    end if

    distances2(1:num_neighbour_particles+1, 1:np_local) = huge(0._8)
    keys(1:num_neighbour_particles+1, 1:np_local) = -13_8
    maxdist = 1

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

       end if
       
       call MPI_BCAST( key_buffer, all_np_local( actual_pe+1 ), MPI_INTEGER8, actual_pe, MPI_COMM_WORLD, ierr )
       call MPI_BCAST( x_buffer,   all_np_local( actual_pe+1 ), MPI_REAL8,    actual_pe, MPI_COMM_WORLD, ierr )
       call MPI_BCAST( y_buffer,   all_np_local( actual_pe+1 ), MPI_REAL8,    actual_pe, MPI_COMM_WORLD, ierr )
       call MPI_BCAST( z_buffer,   all_np_local( actual_pe+1 ), MPI_REAL8,    actual_pe, MPI_COMM_WORLD, ierr )
       
!       write(*,*), all_np_local( actual_pe +1)
!       write(*,'(30(O30,x))'), particles( 1:np_local )%key
!       write(*,'(30(O30,x))'), key_buffer(1:np_local )

       do local_particle_index = 1, np_local

          tmp_loc = maxloc(distances2(1:num_neighbour_particles+1, local_particle_index ) )
          maxdist = tmp_loc(1)        !< this is needed because maxloc returns an array
          
          ! actual_pe+1 because actual_pe is from 0 to n_cpu-1 and fortran arrays are normally from 1 to ...
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!  draw neighbours  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    write(cfile, '(a,i6.6)') 'test_particles_neighbours_', itime
    
    !  Header file written out by root PE: does box and includes particle O/P from all PEs
    if ( my_rank .eq. 0 ) then
       
       outfile = TRIM(cfile) // "_header.gle"
       open(60,file=outfile)
       
       !  initialise graphics filter
       write (60,'(a,4(/a),2(/a,2f13.4))'), &
            'size 18 18', &
            'set font rm', &
            'set lwidth 0.001 lstyle 1', &
            'psize=0.005', &
            'begin translate 0.5 0.5', &
            'begin scale ', 17./boxsize, 17./boxsize, &
            'begin translate ', -xmin, -ymin
       
       !     write (60,'(a,2f13.4)') 'amove', xmin, ymin
       !     write (60,'(a,2f13.4)') 'box ',boxsize,boxsize
       
       close(60)
    endif

    
    !  Now do particles belonging to each processor domain
    write (outfile,'(2a,i3.3,a)') TRIM(cfile), "_dom", my_rank, ".gle"
    open (60,file=outfile) 
    
    do local_particle_index=1,np_local
       write (60,'(a,a)') 'set color ',colors( mod(my_rank,8) )
       write (60,'(a,2f13.4)') 'amove ', particles(local_particle_index)%x(1), particles(local_particle_index)%x(2)
       write (60,'(2a)') 'circle psize fill ',colors( mod(my_rank,8))
    end do
    
    close(60)


    ! Now write one gle file for each particle, which includes above written files with all particles and overplots the actual particle
    ! in black, the neighbours of this particle with a black circle and a big grey circle for the maximum neighbour radius
    do local_particle_index =1, np_local
  
       write (outfile, '(a,a,i6.6,a)') TRIM(cfile), '_', particles(local_particle_index)%label, '.gle'
       open (60,file=outfile)

       ! include the header file
       write (60,'(3a)') 'include ', TRIM(cfile), "_header.gle"

       ! print smoothing-length circle
       smoothing_length = sqrt(maxval(distances2(1:num_neighbour_particles, local_particle_index)))

       write (60, '(a)') 'set color black'
       write (60, '(a,2f13.4)') 'amove ', particles(local_particle_index)%x(1), particles(local_particle_index)%x(2)
       ! TODO change g30.10 back to f13.4
       write (60, '(a,g30.10,a)') 'circle ', smoothing_length, ' fill lightgray'

       ! include files with local particles from all domains
       do actual_pe =0, n_cpu-1
          write (60,'(3a,i3.3,a)') 'include ', TRIM(cfile), "_dom", actual_pe, ".gle"
       end do

       ! overplot neighbours with a black circle
       do actual_neighbour = 1, num_neighbour_particles
          actual_key = keys(actual_neighbour, local_particle_index)
          call key_to_coord(actual_key, neighbour_x, neighbour_y, neighbour_z)
          
          write (60, '(a)') 'set color black'
          write (60, '(a,2f13.4)') 'amove ', neighbour_x, neighbour_y
          !        write (60, '(a,2f13.4)') 'amove ', xcoc(next_neighbours(j,p)), ycoc(next_neighbours(j,i))
          write (60, '(a)') 'circle psize'
       end do
       
       ! highlight current particle
       write (60, '(a)') 'set color black'
       write (60, '(a,2f13.4)') 'amove ', particles(local_particle_index)%x(1), particles(local_particle_index)%x(2)
       write (60, '(a)') 'circle psize fill black'

       ! footer
       write (60,'(a/a/a)') 'end translate','end scale','end translate'
       close(60)
       
    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!  end draw neighbours  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    ! now compare list with result of neighbour search

    not_found = 0
    
    do local_particle_index = 1, np_local
       do index_in_result_neighbour_list = 1, num_neighbour_particles

          actual_node = particle_results(local_particle_index)%neighbour_nodes(index_in_result_neighbour_list)
          node_key = coord_to_key_lastlevel(tree_nodes(actual_node)%coc(1), tree_nodes(actual_node)%coc(2), tree_nodes(actual_node)%coc(3))

          found = .false.
          
          do index_in_test_neighbour_list = 1, num_neighbour_particles   ! ignore particle self (num_neighbour_particles+1)
             
             if( node_key .eq. keys(index_in_test_neighbour_list, local_particle_index) ) then
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
                  node_key, " was not found in the result neighbour list"
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


    write(*,*) my_rank, 'not found: ', not_found

  end subroutine validate_n_nearest_neighbour_list





! ======================
!
!   DRAW all particles colored by domain and neighbours
!   for postprocessing by GLE 

!   call from fields.f90 after 
!
! ======================


  subroutine draw_neighbours(np_local, npart_total, particles, particle_results, itime)
    
    use treetypes
    
    use treevars, only: &
         xmin, &
         ymin, &
         boxsize, &
         tree_nodes
    
    use physvars, only: &
         n_cpu, &
         my_rank
    
    implicit none
    include 'mpif.h'
    
    
    integer, intent(in) :: np_local    !< # particles on this CPU
    integer, intent(in) :: npart_total !< total # simulation particles
    type(t_particle), allocatable, intent(in) :: particles(:)
    type(t_particle_results), intent(in), allocatable :: particle_results(:)
    integer, intent(in) :: itime  ! timestep
    
    integer :: ierr
    character(30) :: cfile
    character(40) :: outfile    
    integer :: actual_pe
    integer :: local_particle_index
    integer :: actual_neighbour
    integer*8 :: actual_node
    
    character(12), parameter :: colors(0:15) = (/"orange      ", "cyan        ", "magenta     ", "blue        ", "green       ", &
         "red         ","yellow      ","grey20      ","brown       ", "yellowgreen ", "violet      ", "royalblue   ", &
         "plum        ", "goldenrod   ",  "powderblue  ", "lime        "/) 
    
    write(cfile, '(a,i6.6)') 'particles_neighbours_', itime
    
    !  Header file written out by root PE: does box and includes particle O/P from all PEs
    if ( my_rank .eq. 0 ) then
       
       outfile = TRIM(cfile) // "_header.gle"
       open(60,file=outfile)
       
       !  initialise graphics filter
       write (60,'(a,4(/a),2(/a,2f13.4))'), &
            'size 18 18', &
            'set font rm', &
            'set lwidth 0.001 lstyle 1', &
            'psize=0.005', &
            'begin translate 0.5 0.5', &
            'begin scale ', 17./boxsize, 17./boxsize, &
            'begin translate ', -xmin, -ymin
       
       !     write (60,'(a,2f13.4)') 'amove', xmin, ymin
       !     write (60,'(a,2f13.4)') 'box ',boxsize,boxsize
       
       close(60)
    endif

    
    !  Now do particles belonging to each processor domain
    write (outfile,'(2a,i3.3,a)') TRIM(cfile), "_dom", my_rank, ".gle"
    open (60,file=outfile) 
    
    do local_particle_index=1,np_local
       write (60,'(a,a)') 'set color ',colors( mod(my_rank,8) )
       write (60,'(a,2f13.4)') 'amove ', particles(local_particle_index)%x(1), particles(local_particle_index)%x(2)
       write (60,'(2a)') 'circle psize fill ',colors( mod(my_rank,8))
    end do
    
    close(60)


    ! Now write one gle file for each particle, which includes above written files with all particles and overplots the actual particle
    ! in black, the neighbours of this particle with a black circle and a big grey circle for the maximum neighbour radius
    do local_particle_index =1, np_local
  
       write (outfile, '(a,a,i6.6,a)') TRIM(cfile), '_', particles(local_particle_index)%label, '.gle'
       open (60,file=outfile)

       ! include the header file
       write (60,'(3a)') 'include ', TRIM(cfile), "_header.gle"

       ! print smoothing-length circle
       write (60, '(a)') 'set color black'
       write (60, '(a,2f13.4)') 'amove ', particles(local_particle_index)%x(1), particles(local_particle_index)%x(2)
       write (60, '(a,f13.4,a)') 'circle ', sqrt(particle_results(local_particle_index)%maxdist2), ' fill lightgray'

       ! include files with local particles from all domains
       do actual_pe =0, n_cpu-1
          write (60,'(3a,i3.3,a)') 'include ', TRIM(cfile), "_dom", actual_pe, ".gle"
       end do

       ! overplot neighbours with a black circle
       do actual_neighbour = 1, num_neighbour_particles
          actual_node = particle_results(local_particle_index)%neighbour_nodes(actual_neighbour)
          
          write (60, '(a)') 'set color black'
          write (60, '(a,2f13.4)') 'amove ', tree_nodes(actual_node)%coc(1), tree_nodes(actual_node)%coc(2)
          !        write (60, '(a,2f13.4)') 'amove ', xcoc(next_neighbours(j,p)), ycoc(next_neighbours(j,i))
          write (60, '(a)') 'circle psize'
       end do
       
       ! highlight current particle
       write (60, '(a)') 'set color black'
       write (60, '(a,2f13.4)') 'amove ', particles(local_particle_index)%x(1), particles(local_particle_index)%x(2)
       write (60, '(a)') 'circle psize fill black'

       ! footer
       write (60,'(a/a/a)') 'end translate','end scale','end translate'
       close(60)
       
    end do
    
  end subroutine draw_neighbours
  

end module module_neighbour_test
