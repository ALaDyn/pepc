!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!> Encapsulates domain decomposition and restoration of original particle order
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_domains
    implicit none
    save
    private

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    integer, public :: weighted = 1 !< set to 0 to disable load balancing, 1 to enable load balancing



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    public tree_domains
    public restore


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !>  Domain decomposition:
    !>  Share particle keys amoung PEs
    !>  - weighting according to load incurred on previous timestep
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine tree_domains(particles, nppm, indxl,irnkl,islen,irlen,fposts,gposts,npnew,npold,neighbour_pe_particles)

        use treevars
        use module_interaction_specific
        use module_utils
        use module_timings
        use module_spacefilling
        use module_branching
        use module_debug
        implicit none
        include 'mpif.h'

        type(t_particle), allocatable, intent(inout) :: particles(:)

        integer, intent(in) :: nppm !< maximum allowed number of particles on this PE
        integer, intent(out) :: indxl(nppm),irnkl(nppm)
        integer, intent(out) :: islen(num_pe),irlen(num_pe)
        integer, intent(out) :: fposts(num_pe+1),gposts(num_pe+1)
        integer :: npnew,npold
        integer, intent(out) :: neighbour_pe_particles !< number of particles that have been fetched from neighbouring PEs - they are stored in particles(npp+1:npp+neighbour_pe_particles)

        integer :: i, prev, next, handle(4)

        real*8 :: s
        real*8 :: xmin_local, xmax_local, ymin_local, ymax_local, zmin_local, zmax_local

        integer :: state(MPI_STATUS_SIZE), ierr
        logical :: fixedduplicate

        ! arrays for parallel sort

        type (t_particle) :: ship_parts(nppm), get_parts(nppm)

        integer*8 :: w1(nppm)

        real*8 :: xboxsize, yboxsize, zboxsize

        real*8 imba

        type (t_particle) :: ship_props, get_props


        interface
            subroutine slsort_keys(nin,nmax,keys,workload,balance_weight,max_imbalance,nout,indxl,irnkl,scounts,rcounts,sdispls,rdispls,keys2,irnkl2,size,rank)
                integer,intent(in) :: nin,nmax,balance_weight,size,rank
                real*8,intent(in) :: max_imbalance
                integer,intent(out) :: nout,indxl(*),irnkl(*),scounts(*),rcounts(*),sdispls(*),rdispls(*),irnkl2(*)
                integer*8,intent(out) :: keys2(*)
                integer*8,intent(inout) :: keys(*)
                real*8,intent(inout) :: workload(*)
            end subroutine slsort_keys
        end interface

        real*8 work2(nppm)
        integer irnkl2(nppm)
        integer*8 :: local_keys(nppm)

        call timer_start(t_domains)
        call timer_start(t_domains_keys)

        call pepc_status('DOMAIN DECOMPOSITION')

        ! Find limits of local simulation region
        xmin_local = minval(particles(1:npp)%x(1))
        xmax_local = maxval(particles(1:npp)%x(1))
        ymin_local = minval(particles(1:npp)%x(2))
        ymax_local = maxval(particles(1:npp)%x(2))
        zmin_local = minval(particles(1:npp)%x(3))
        zmax_local = maxval(particles(1:npp)%x(3))

        ! Find global limits
        call MPI_ALLREDUCE(xmin_local, xmin, 1, MPI_REAL8, MPI_MIN,  MPI_COMM_WORLD, ierr )
        call MPI_ALLREDUCE(xmax_local, xmax, 1, MPI_REAL8, MPI_MAX,  MPI_COMM_WORLD, ierr )
        call MPI_ALLREDUCE(ymin_local, ymin, 1, MPI_REAL8, MPI_MIN,  MPI_COMM_WORLD, ierr )
        call MPI_ALLREDUCE(ymax_local, ymax, 1, MPI_REAL8, MPI_MAX,  MPI_COMM_WORLD, ierr )
        call MPI_ALLREDUCE(zmin_local, zmin, 1, MPI_REAL8, MPI_MIN,  MPI_COMM_WORLD, ierr )
        call MPI_ALLREDUCE(zmax_local, zmax, 1, MPI_REAL8, MPI_MAX,  MPI_COMM_WORLD, ierr )

        xboxsize = xmax-xmin
        yboxsize = ymax-ymin
        zboxsize = zmax-zmin

        ! Safety margin - put buffer region around particles
        xmax = xmax + xboxsize/10000.
        xmin = xmin - xboxsize/10000.
        ymax = ymax + yboxsize/10000.
        ymin = ymin - yboxsize/10000.
        zmax = zmax + zboxsize/10000.
        zmin = zmin - zboxsize/10000.

        boxsize = max(xmax-xmin, ymax-ymin, zmax-zmin)

        if (dbg(DBG_DOMAIN)) then
          DEBUG_WARNING('(4(a15,f12.4/))',
            'xmin = ',xmin,'xmax = ',xmax,
            'ymin = ',ymin,'ymax = ',ymax,
            'zmin = ',zmin,'zmax = ',zmax,
            'boxsize = ',boxsize )
        endif

        s=boxsize/2**nlev       ! refinement length

        call compute_particle_keys(particles)

        if (dbg(DBG_DOMAIN)) call print_particle_list(particles, npp, &
                                     'Particle list before key sort (see t_particle in module_pepc_types.f90 for meaning of the columns)')

        call timer_stop(t_domains_keys)
        call timer_start(t_domains_sort)

        ! Define wraps for ring network  0 -> 1 -> 2 -> ... ... -> num_pe-1 -> 0 ...
        if (me == 0) then
            prev = num_pe - 1
        else
            prev = me-1
        endif

        if (me == num_pe-1 ) then
            next = 0
        else
            next = me+1
        endif

        imba = 0.01

        npold = npp
        npnew = npp

        call timer_start(t_domains_add_sort)

        ! start permutation of local key list
        work2(1:npp) = particles(1:npp)%work

        call timer_start(t_domains_sort_pure)

        local_keys(1:npold) = particles(1:npold)%key

        ! perform index sort on keys !TODO: remove the "-2", compare other cases with "+2" and "npp+1" etc.
        call slsort_keys(npold,nppm-2,local_keys,work2,weighted,imba,npnew,indxl,irnkl,islen,irlen,fposts,gposts,w1,irnkl2,num_pe,me)

        ! FIXME: every processor has to have at least one particle
        if (npnew < 2) then
            DEBUG_ERROR('("rank less than two particles after sorting (had ", I8, " before) - currently this can lead to errors --> aborting")', npold)
        endif

        call timer_stop(t_domains_sort_pure)

        call timer_stop(t_domains_sort)
        call timer_start(t_domains_ship)

        ! Now permute particle properties
        ! Set up particle structure
        call timer_start(t_domains_add_pack)

        do i=1,npold
            ship_parts(i) = particles( indxl(i) )
        enddo

        call timer_stop(t_domains_add_pack)

        deallocate(particles) ! has size npold until here, i.e. npp == npold

        call timer_start(t_domains_add_alltoallv)

        ! perform permute
        call MPI_alltoallv(  ship_parts, islen, fposts, mpi_type_particle, &
        get_parts, irlen, gposts, mpi_type_particle, &
        MPI_COMM_WORLD,ierr)

        call timer_stop(t_domains_add_alltoallv)

        allocate(particles(npnew+2)) ! TODO: the limit particles from neighbouring PEs are put into the final two places - this is for branching and correct insertion into the tree and should be done there with local variables instead
        npp = npnew

        call timer_start(t_domains_add_unpack)

        do i=1,npp
            particles( irnkl(i) ) = get_parts(i)
        enddo

        call timer_stop(t_domains_add_unpack)

        if (npp > nppm) then
            DEBUG_ERROR('("More than nppm particles after sorting: nppm = ", I0, " < npp = ",I0,". All local particle fields are too shirt. Aborting.")', nppm, npp)
        endif

        call timer_stop(t_domains_ship)
        call timer_stop(t_domains_add_sort)
        call timer_start(t_domains_bound)

        if (dbg(DBG_DOMAIN)) call print_particle_list(particles, npp, &
                                     'Particle list after key sort (see t_particle in module_pepc_types.f90 for meaning of the columns)')

        particles(1:npp)%pid = me  ! new owner

        ! Each PE now has sorted segment of particles of length npp
        ! Note that now npp /= npart/num_pe, only approx depending on key distribution, or target shape.

        ! check for duplicate keys
        ! we expect the sorting library to avoid duplicate keys across processor boundaries
        ! hence, we should not modify the first and the last particles keys to avoid damaging this restriction
        do i=1,npp-2
            if (particles(i)%key == particles(i+1)%key) then
                DEBUG_INFO('("Identical keys found for i = ", I0, " and its successor, key = ", O0, ", labels = ", I0,x, I0)', i, particles(i)%key, particles(i)%label, particles(i+1)%label)

                fixedduplicate = .false.

                ! first, we try to shift down the lower key
                if (i >= 2) then ! do not merge both if-statements, since in case of compiler optimization rearranging both conditions, an acces to particles(i-1=0) might occur
                  if (particles(i)%key - 1 .ne. particles(i-1)%key) then
                    particles(i)%key = particles(i)%key - 1
                    ! adjust position to fit key
                    call key_to_coord_dim(particles(i)%key, particles(i)%x, idim, particles(i)%x)
                    DEBUG_INFO('("shifting (i)-th particles key down to ", O0)', particles(i)%key)
                    fixedduplicate = .true.
                  endif
                endif

                if (.not. fixedduplicate) then
                  ! we have to shift up the upper key - if the keys are dense, this might propagate further
                  ! upwards until a gap, i.e. particles(i+1)%key - particles(i)%key > 1, exists
                  particles(i+1)%key = particles(i+1)%key + 1
                  ! adjust position
                  call key_to_coord_dim(particles(i+1)%key, particles(i+1)%x, idim, particles(i+1)%x)
                  DEBUG_INFO('("shifting (i+1)-th particles key up to ", O0)', particles(i+1)%key)
                endif
            endif
        end do
        ! we did not check the pair (npp-1) (npp) for duplicate keys since we did not want to touch processor boundaries --> work downwards again
        do i=npp,2,-1
          if (particles(i)%key == particles(i-1)%key) then
            DEBUG_INFO('("Identical keys found in second pass for i = ", I0, " and its predecessor, key = ", O0, ", labels = ", I0,x, I0)', i, particles(i)%key, particles(i)%label, particles(i+1)%label)

            if (i == 2) then
              DEBUG_ERROR(*, "Obviously, keys are dense on this PE, i.e. there is no gap between any keys to resolve the identical-key-conflicts. Aborting.")
            endif

            ! we have to shift down the lower key
            particles(i-1)%key = particles(i-1)%key - 1
            call key_to_coord_dim(particles(i-1)%key, particles(i-1)%x, idim, particles(i-1)%x)
            DEBUG_INFO('("shifting (i-1)-th particles key down to ", O0)', particles(i-1)%key)
          else
            exit ! from this loop - no more duplicate keys left
          endif
        end do


        ! Copy boundary particles to adjacent PEs to ensure proper tree construction
        !  - if we do not do this, can get two particles on separate PEs 'sharing' a leaf
        ship_props = particles( 1 )

        ! Ship 1st particle data to end of list of LH neighbour PE
        neighbour_pe_particles = 0

        if (me /= 0 ) then
            call MPI_ISEND( ship_props, 1, mpi_type_particle, prev, 1, MPI_COMM_WORLD, handle(1), ierr )
            call MPI_REQUEST_FREE(handle(1),ierr)
        endif

        ! Place incoming data at end of array
        if ( me /= num_pe-1) then
            call MPI_RECV( get_props, 1, mpi_type_particle, next, 1,  MPI_COMM_WORLD, state, ierr )
            neighbour_pe_particles = neighbour_pe_particles + 1
            particles(npp+neighbour_pe_particles) = get_props
        endif

        ! Ship  end particle data to start of list of RH neighbour PE
        ship_props = particles( npp )

        if (me /= num_pe-1 ) then
            call MPI_ISEND( ship_props, 1, mpi_type_particle, next, 2, MPI_COMM_WORLD, handle(3), ierr )
            call MPI_REQUEST_FREE(handle(3),ierr)
        endif

        ! Place incoming data at end of array

        if ( me /= 0) then
            call MPI_RECV( get_props, 1, mpi_type_particle, prev, 2,  MPI_COMM_WORLD, state, ierr )
            neighbour_pe_particles = neighbour_pe_particles + 1
            particles(npp + neighbour_pe_particles) = get_props
        endif

        ! Initialize VLD-stuff
        call branches_initialize_VLD(particles)

        if (dbg(DBG_DOMAIN)) call print_particle_list(particles, npp + neighbour_pe_particles, &
                                     'Particle list after boundary swap (see t_particle in module_pepc_types.f90 for meaning of the columns)')

        call timer_stop(t_domains_bound)
        call timer_stop(t_domains)

    end subroutine tree_domains


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !>  Restore initial particle order
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine restore(npnew,npold,nppm_ori,indxl,irnkl,islen,irlen,fposts,gposts,&
                             particles)
        use module_interaction_specific
        use treevars, only : num_pe, npp
        use module_pepc_types
        use module_debug, only : pepc_status
        implicit none
        include 'mpif.h'

        integer, intent(in) :: npnew,npold,nppm_ori
        integer, intent(in) :: indxl(nppm_ori),irnkl(nppm_ori)
        integer, intent(in) :: islen(num_pe),irlen(num_pe)
        integer, intent(in) :: fposts(num_pe+1),gposts(num_pe+1)
        type(t_particle),         intent(inout), allocatable :: particles(:)

        integer :: i, ierr

        type (t_particle)         :: get_parts(npold), ship_parts(npnew)

        call pepc_status('RESTORE DOMAINS')

        do i=1,npnew
          ship_parts(i) = particles(indxl(i))
        enddo

        deallocate(particles) ! had size npnew

        ! perform permute
        call MPI_alltoallv(  ship_parts, islen, fposts, MPI_TYPE_particle, &
              get_parts, irlen, gposts, MPI_TYPE_particle, &
              MPI_COMM_WORLD,ierr )

        allocate(particles(npold))
        npp = npold

        do i=1,npold
            particles(irnkl(i)) = get_parts(i)
        enddo

    end subroutine restore


    subroutine print_particle_list(particles, npart, callinfo)
      use module_pepc_types
      use module_debug
      implicit none
      type(t_particle), intent(in) :: particles(:)
      integer, intent(in) :: npart
      character(*), intent(in) :: callinfo

      integer :: j

      call debug_ipefile_open()
      write (debug_ipefile,'(/a/)') callinfo
      do j=1,npart
        write(debug_ipefile,'(i10)',advance='no') j
        write(debug_ipefile,*)                     particles(j)
        write(debug_ipefile,'(/)')
      end do
      write(debug_ipefile,'(/)')
      call debug_ipefile_close()

   end subroutine

end module module_domains
