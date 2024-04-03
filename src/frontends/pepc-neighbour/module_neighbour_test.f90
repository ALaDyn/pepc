! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
!
! Copyright (C) 2002-2021 Juelich Supercomputing Centre,
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_neighbour_test

   use module_interaction_specific_types, only: &
      num_neighbour_particles

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

   subroutine validate_n_nearest_neighbour_list(np_local, particles, &
                                                itime, num_neighbour_boxes, neighbour_boxes)

      ! TODO: use the openmp module only, when compiling with openmp: !$ use omp_lib
      ! module_deps has to be changed to remove "!$" when using openmp
      use omp_lib

      use module_pepc_types, only: &
         t_particle, &
         t_tree_node

      use module_pepc, only: &
         t => global_tree

      use physvars, only: &
         n_cpu, &
         my_rank

      use module_spacefilling, only: &
         coord_to_key, &
         key_to_coord

      use module_mirror_boxes, only: &
         lattice_vect

      use mpi
      implicit none

      integer, intent(in) :: np_local    !< # particles on this CPU
      type(t_particle), intent(in) :: particles(:)
      integer, intent(in) :: itime  ! timestep
      integer, intent(in) :: num_neighbour_boxes !< number of shift vectors in neighbours list (must be at least 1 since [0, 0, 0] has to be inside the list)
      integer, intent(in) :: neighbour_boxes(3, num_neighbour_boxes) ! list with shift vectors to neighbour boxes that shall be included in interaction calculation, at least [0, 0, 0] should be inside this list

      integer :: ierr
      integer :: actual_pe
      real*8, allocatable :: coord_buffer(:, :)
!    real*8, allocatable :: y_buffer(:)
!    real*8, allocatable :: z_buffer(:)

      integer :: local_particle_index
      integer :: remote_particle_index
      real*8, allocatable :: distances2(:, :)
      real*8, allocatable :: positions(:, :, :)
      real*8 :: tmp_real8
      logical :: found
      integer :: not_found
      integer :: all_not_found
      integer, dimension(1) :: tmp_loc
      logical :: tree_nn_debug
      logical :: draw_neighbour_test

      integer :: i
      integer :: maxdist
      integer :: index_in_test_neighbour_list
      integer :: index_in_result_neighbour_list
      real*8 :: dist2
      type(t_tree_node), pointer :: actual_node
      integer*8 :: node_key
      integer*8 :: neighbour_key

      character(100) :: filename

      integer, dimension(n_cpu) :: all_np_local
      integer :: max_np_local
      integer :: ibox

      real*8, dimension(3) :: vbox

      ! variables for gle output
      character(40) :: cfile
      character(50) :: outfile
      integer :: actual_neighbour
      real*8 :: neighbour_x
      real*8 :: neighbour_y
      real*8 :: neighbour_z
      real*8 :: smoothing_length
      character(12), parameter :: colors(0:15) = (/"orange      ", "cyan        ", "magenta     ", "blue        ", "green       ", &
                                                   "red         ", "yellow      ", "grey20      ", "brown       ", "yellowgreen ", "violet      ", "royalblue   ", &
                                                   "plum        ", "goldenrod   ", "powderblue  ", "lime        "/)

      ! end variable declaration

      tree_nn_debug = .false.
      draw_neighbour_test = .false.

      ! num_neighbour_particles+1 because the particle itself is stored within the lists first and removed later
      allocate (distances2(num_neighbour_particles + 1, np_local), positions(3, num_neighbour_particles + 1, np_local), STAT=ierr)

      if (ierr .ne. 0) then
         write (*, *) 'allocate of fields for keys and distances to next neighbours in validate_n_nearest_neighbour_list failed in module_neighbour_test.f90'
      end if

      distances2(1:num_neighbour_particles + 1, 1:np_local) = huge(0._8)
      positions(1:3, 1:num_neighbour_particles + 1, 1:np_local) = -13._8

      call MPI_ALLGATHER(np_local, 1, MPI_INTEGER, all_np_local, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

      max_np_local = maxval(all_np_local)

      allocate (coord_buffer(3, max_np_local), STAT=ierr)

      if (ierr .ne. 0) then
         write (*, *) 'allocate of buffers in validate_n_nearest_neighbour_list failed in module_neighbour_test.f90'
      end if

      do actual_pe = 0, n_cpu - 1                                  ! each process sends data

         if (actual_pe .eq. my_rank) then

            coord_buffer(1, 1:np_local) = particles(1:np_local)%x(1)
            coord_buffer(2, 1:np_local) = particles(1:np_local)%x(2)
            coord_buffer(3, 1:np_local) = particles(1:np_local)%x(3)

         end if

         call MPI_BCAST(coord_buffer, all_np_local(actual_pe + 1) * 3, MPI_REAL8, actual_pe, MPI_COMM_WORLD, ierr)
!       call MPI_BCAST( y_buffer,   all_np_local( actual_pe+1 ), MPI_REAL8,    actual_pe, MPI_COMM_WORLD, ierr )
!       call MPI_BCAST( z_buffer,   all_np_local( actual_pe+1 ), MPI_REAL8,    actual_pe, MPI_COMM_WORLD, ierr )

!       write(*,*), all_np_local( actual_pe +1)
!       write(*,'(30(O30,x))'), particles( 1:np_local )%key

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(local_particle_index, tmp_loc, maxdist, remote_particle_index, dist2, ibox, vbox)
         do local_particle_index = 1, np_local

            tmp_loc = maxloc(distances2(1:num_neighbour_particles + 1, local_particle_index))
            maxdist = tmp_loc(1)        !< this is needed because maxloc returns an array

            ! loop over all periodic neighbour boxes
            do ibox = 1, num_neighbour_boxes

               vbox = lattice_vect(neighbour_boxes(:, ibox))

               ! actual_pe+1 because actual_pe is from 0 to n_cpu-1 and fortran arrays are normally from 1 to ...
               do remote_particle_index = 1, all_np_local(actual_pe + 1)

                  dist2 = sum((coord_buffer(1:3, remote_particle_index) - particles(local_particle_index)%x + vbox)**2)

!( x_buffer( remote_particle_index ) - particles( local_particle_index )%x(1) ) **2 &
!                     + ( y_buffer( remote_particle_index )  - particles( local_particle_index )%x(2) ) **2 &
!                     + ( z_buffer( remote_particle_index )  - particles( local_particle_index )%x(3) ) **2

                  if (dist2 .lt. distances2(maxdist, local_particle_index)) then
                     distances2(maxdist, local_particle_index) = dist2
                     positions(1:3, maxdist, local_particle_index) = coord_buffer(1:3, remote_particle_index)

!                   positions(1, maxdist, local_particle_index ) = x_buffer( remote_particle_index )
!                   positions(2, maxdist, local_particle_index ) = y_buffer( remote_particle_index )
!                   positions(3, maxdist, local_particle_index ) = z_buffer( remote_particle_index )
                     tmp_loc = maxloc(distances2(1:num_neighbour_particles + 1, local_particle_index))
                     maxdist = tmp_loc(1)        !< this is needed because maxloc returns an array
                  end if

               end do

            end do

         end do
!$OMP END PARALLEL DO

      end do

      ! now distances2 and keys contain the num_neighbour_particles closest neighbours and the particle itself

      ! move particle self to end of list and ignore in the following
      do local_particle_index = 1, np_local
         ! find closest particle, should be self
         tmp_loc = minloc(distances2(1:num_neighbour_particles + 1, local_particle_index))
         ! move this particle to the end of the list
         tmp_real8 = distances2(tmp_loc(1), local_particle_index)
         distances2(tmp_loc(1), local_particle_index) = distances2(num_neighbour_particles + 1, local_particle_index)
         distances2(num_neighbour_particles + 1, local_particle_index) = tmp_real8

         do i = 1, 3
            tmp_real8 = positions(i, tmp_loc(1), local_particle_index)
            positions(i, tmp_loc(1), local_particle_index) = positions(i, num_neighbour_particles + 1, local_particle_index)
            positions(i, num_neighbour_particles + 1, local_particle_index) = tmp_real8
         end do

      end do

      ! now distances2 and keys contain the num_neighbour_particles closest neighbours, ignore last element

      do local_particle_index = 1, np_local
         if (distances2(num_neighbour_particles + 1, local_particle_index) .ne. 0.) then
            write (*, *) 'particle self was not removed from neighbour list in validate_n_nearest_neighbour_list failed in module_neighbour_test.f90', local_particle_index
         end if
      end do

      if (tree_nn_debug) then

         write (filename, '(a,i6.6,a,i6.6,a)') "nn_validate_", itime, "_", my_rank, ".list"

         !        ! \bug ab: with ACCESS='APPEND' compilation failes on jugene
         ! !       OPEN(99, FILE=NN_filename, ACCESS='APPEND')
         open (99, file=filename, position='APPEND')

         do local_particle_index = 1, np_local

            write (99, '(i6.6,3(E12.5),i6,O30)') local_particle_index, particles(local_particle_index)%x(1), particles(local_particle_index)%x(2), particles(local_particle_index)%x(3), particles(local_particle_index)%label, particles(local_particle_index)%key
            do actual_neighbour = 1, num_neighbour_particles
               write (99, *) coord_to_key(t%bounding_box, positions(1:3, actual_neighbour, local_particle_index))
            end do

         end do

         close (99)

      end if ! tree_nn_debug

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !!!!!!!!!!!!!!  draw neighbours  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! This block is basically a copy of the draw_neighbours subroutine.
      ! Because the smoothing-length and the coodinates of the neighbours cannot be accessed with particle_results,
      ! these two parts below have to be adjusted

      if (draw_neighbour_test) then

         write (cfile, '(a,i6.6)') 'test_particles_neighbours_', itime

         !  Header file written out by root PE: does box and includes particle O/P from all PEs
         if (my_rank .eq. 0) then

            outfile = TRIM(cfile)//"_header.gle"
            open (60, file=outfile)

            !  initialise graphics filter
            write (60, '(a,4(/a),2(/a,2f13.4))') &
               'size 18 18', &
               'set font rm', &
               'set lwidth 0.001 lstyle 1', &
               'psize=0.005', &
               'begin translate 0.5 0.5', &
               'begin scale ', 17./t%bounding_box%boxsize(1), 17./t%bounding_box%boxsize(2), &
               'begin translate ', -t%bounding_box%boxmin(1), -t%bounding_box%boxmin(2)

            !     write (60,'(a,2f13.4)') 'amove', boxmin(1), boxmin(2)
            !     write (60,'(a,2f13.4)') 'box ',boxsize(1),boxsize(2)

            close (60)
         end if

         !  Now do particles belonging to each processor domain
         write (outfile, '(2a,i3.3,a)') TRIM(cfile), "_dom", my_rank, ".gle"
         open (60, file=outfile)

         do local_particle_index = 1, np_local
            write (60, '(a,a)') 'set color ', colors(mod(my_rank, 8))
            write (60, '(a,2f13.4)') 'amove ', particles(local_particle_index)%x(1), particles(local_particle_index)%x(2)
            write (60, '(2a)') 'circle psize fill ', colors(mod(my_rank, 8))
         end do

         close (60)

         ! Now write one gle file for each particle, which includes above written files with all particles and overplots the actual particle
         ! in black, the neighbours of this particle with a black circle and a big grey circle for the maximum neighbour radius
         do local_particle_index = 1, np_local

            write (outfile, '(a,a,i6.6,a)') TRIM(cfile), '_', particles(local_particle_index)%label, '.gle'
            open (60, file=outfile)

            ! include the header file
            write (60, '(3a)') 'include ', TRIM(cfile), "_header.gle"

            ! print smoothing-length circle
            smoothing_length = sqrt(maxval(distances2(1:num_neighbour_particles, local_particle_index)))

            write (60, '(a)') 'set color black'
            write (60, '(a,2f13.4)') 'amove ', particles(local_particle_index)%x(1), particles(local_particle_index)%x(2)
            ! TODO change g30.10 back to f13.4
            write (60, '(a,g30.10,a)') 'circle ', smoothing_length, ' fill lightgray'

            ! include files with local particles from all domains
            do actual_pe = 0, n_cpu - 1
               write (60, '(3a,i3.3,a)') 'include ', TRIM(cfile), "_dom", actual_pe, ".gle"
            end do

            ! overplot neighbours with a black circle
            do actual_neighbour = 1, num_neighbour_particles
               neighbour_x = positions(1, actual_neighbour, local_particle_index)
               neighbour_y = positions(2, actual_neighbour, local_particle_index)
               neighbour_z = positions(3, actual_neighbour, local_particle_index)

               write (60, '(a)') 'set color black'
               write (60, '(a,2f13.4)') 'amove ', neighbour_x, neighbour_y
               write (60, '(a)') 'circle psize'
            end do

            ! highlight current particle
            write (60, '(a)') 'set color black'
            write (60, '(a,2f13.4)') 'amove ', particles(local_particle_index)%x(1), particles(local_particle_index)%x(2)
            write (60, '(a)') 'circle psize fill black'

            ! footer
            write (60, '(a/a/a)') 'end translate', 'end scale', 'end translate'
            close (60)

         end do

      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!  end draw neighbours  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! now compare list with result of neighbour search

      not_found = 0

      do local_particle_index = 1, np_local
         do index_in_result_neighbour_list = 1, num_neighbour_particles

            actual_node => t%nodes(particles(local_particle_index)%results%neighbour_nodes(index_in_result_neighbour_list))
            node_key = coord_to_key(t%bounding_box, actual_node%interaction_data%coc)

            found = .false.

            do index_in_test_neighbour_list = 1, num_neighbour_particles   ! ignore particle self (num_neighbour_particles+1)

               neighbour_key = coord_to_key(t%bounding_box, positions(1:3, index_in_test_neighbour_list, local_particle_index))

               if (node_key .eq. neighbour_key) then
                  found = .true.
                  exit
               end if
            end do

            if (.not. found) then

               write (filename, '(a,i6.6,a,i6.6,a)') "validation_", itime, "_", my_rank, ".errors"

               ! \bug ab: with ACCESS='APPEND' compilation failes on jugene
               !             OPEN( 76, FILE=filename, ACCESS='APPEND')
               open (76, FILE=filename, POSITION='Append')
               ! use format O30 for keys for octal
               write (76, '(a,i6,a,i19.19,a,i6)') "for particle ", particles(local_particle_index)%label, " the neighbour with key ", &
                  node_key, " was not found in the result neighbour list"
!             WRITE( 76, *) '  ', nn_keys(1:n_nn)
!             WRITE( 76, *) '  ', xcoc( nodelist(1:n_nn,i) )
!             WRITE( 76, *) '  ', sqrt(dist2_list(1:n_nn,i))
!             WRITE( 76, *) '  ', nterm(i), r_nn( pshort( i ) ), tree_x( pshort( i))
!             WRITE( 76, '(O30,E12.5E2)') stored_keys(i, j), x_from_key
!             WRITE( 76, *) '  ', stored_keys( i , 1:n_nn)
               close (76)

               not_found = not_found + 1

            end if

         end do

      end do

      ! TODO: remove this debugging output
      ! write(*,*) my_rank, 'not found: ', not_found

      deallocate (distances2, positions)
      deallocate (coord_buffer)

      call MPI_REDUCE(not_found, all_not_found, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

      if (my_rank .eq. 0) then

         open (55, file='nn_validation_results.dat', status='UNKNOWN', position='APPEND')

         write (55, *) itime, "nn_validate found ", all_not_found, "differences in nn-lists (", (t%npart * num_neighbour_particles), "neighbours total, ~ ", &
            int(all_not_found * 100./(t%npart * num_neighbour_particles)), "% )"

         close (55)

      end if

   end subroutine validate_n_nearest_neighbour_list

! ======================
!
!   DRAW all particles colored by domain and neighbours
!   for postprocessing by GLE
!
! ======================

   subroutine draw_neighbours(np_local, particles, itime)

      use module_pepc, only: &
         t => global_tree

      use module_pepc_types, only: &
         t_particle, &
         t_tree_node

      use physvars, only: &
         n_cpu, &
         my_rank

      use module_mirror_boxes, only: &
         num_neighbour_boxes, &
         neighbour_boxes, &
         lattice_vect

      use mpi
      implicit none

      integer, intent(in) :: np_local    !< # particles on this CPU
      type(t_particle), intent(in) :: particles(:)
      integer, intent(in) :: itime  ! timestep

      character(30) :: cfile
      character(40) :: outfile
      integer :: actual_pe
      integer :: local_particle_index
      integer :: actual_neighbour
      type(t_tree_node), pointer :: actual_node
      integer :: ibox
      real*8, dimension(3) :: vbox

      character(12), parameter :: colors(0:15) = ["orange      ", "cyan        ", "magenta     ", "blue        ", "green       ", &
                                                  "red         ", "yellow      ", "grey20      ", "brown       ", "yellowgreen ", "violet      ", "royalblue   ", &
                                                  "plum        ", "goldenrod   ", "powderblue  ", "lime        "]

      write (cfile, '(a,i6.6)') 'particles_neighbours_', itime

      !  Header file written out by root PE: does box and includes particle O/P from all PEs
      if (my_rank .eq. 0) then

         outfile = TRIM(cfile)//"_header.gle"
         open (60, file=outfile)

         !  initialise graphics filter
         write (60, '(a,4(/a),2(/a,2f13.4))') &
            'size 18 18', &
            'set font rm', &
            'set lwidth 0.001 lstyle 1', &
            'psize=0.005', &
            'begin translate 0.5 0.5', &
            'begin scale ', 17./t%bounding_box%boxsize(1), 17./t%bounding_box%boxsize(2), &
            'begin translate ', -t%bounding_box%boxmin(1), -t%bounding_box%boxmin(2)

         !     write (60,'(a,2f13.4)') 'amove', boxmin(1), boxmin(2)
         !     write (60,'(a,2f13.4)') 'box ',boxsize(1),boxsize(2)

         close (60)
      end if

      !  Now do particles belonging to each processor domain
      write (outfile, '(2a,i3.3,a)') TRIM(cfile), "_dom", my_rank, ".gle"
      open (60, file=outfile)

      do local_particle_index = 1, np_local
         write (60, '(a,a)') 'set color ', colors(mod(my_rank, 8))
         write (60, '(a,2f13.4)') 'amove ', particles(local_particle_index)%x(1), particles(local_particle_index)%x(2)
         write (60, '(2a)') 'circle psize fill ', colors(mod(my_rank, 8))
      end do

      close (60)

      ! Now write one gle file for each particle, which includes above written files with all particles and overplots the actual particle
      ! in black, the neighbours of this particle with a black circle and a big grey circle for the maximum neighbour radius
      do local_particle_index = 1, np_local

         write (outfile, '(a,a,i6.6,a)') TRIM(cfile), '_', particles(local_particle_index)%label, '.gle'
         open (60, file=outfile)

         ! include the header file
         write (60, '(3a)') 'include ', TRIM(cfile), "_header.gle"

         ! print smoothing-length circle
         ! for all periodic neighbour boxes
         do ibox = 1, num_neighbour_boxes

            vbox = lattice_vect(neighbour_boxes(:, ibox))

            write (60, '(a)') 'set color black'
            write (60, '(a,2f13.4)') 'amove ', particles(local_particle_index)%x(1) + vbox(1), particles(local_particle_index)%x(2) + vbox(2)
            write (60, '(a,f13.4,a)') 'circle ', sqrt(particles(local_particle_index)%results%maxdist2), ' fill lightgray'

         end do

         ! include files with local particles from all domains
         do actual_pe = 0, n_cpu - 1
            write (60, '(3a,i3.3,a)') 'include ', TRIM(cfile), "_dom", actual_pe, ".gle"
         end do

         ! overplot neighbours with a black circle
         do actual_neighbour = 1, num_neighbour_particles
            actual_node => t%nodes(particles(local_particle_index)%results%neighbour_nodes(actual_neighbour))

            write (60, '(a)') 'set color black'
            write (60, '(a,2f13.4)') 'amove ', actual_node%interaction_data%coc(1), actual_node%interaction_data%coc(2)
            !        write (60, '(a,2f13.4)') 'amove ', xcoc(next_neighbours(j,p)), ycoc(next_neighbours(j,i))
            write (60, '(a)') 'circle psize'
         end do

         ! highlight current particle
         write (60, '(a)') 'set color black'
         write (60, '(a,2f13.4)') 'amove ', particles(local_particle_index)%x(1), particles(local_particle_index)%x(2)
         write (60, '(a)') 'circle psize fill black'

         ! footer
         write (60, '(a/a/a)') 'end translate', 'end scale', 'end translate'
         close (60)

      end do

   end subroutine draw_neighbours

end module module_neighbour_test

