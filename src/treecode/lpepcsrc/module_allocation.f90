!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!> Contains all lpepcsrc-global allocation and deallocation routines
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_allocation
    implicit none

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
    !>
    !> Allocates all tree- and multipole-specific arrays,
    !> Initializes several estimations on maximum used memory
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine allocate_tree()
        use treevars
        use module_timings
        use module_debug
        use module_htable
        use module_branching
        use module_interaction_specific, only : get_number_of_interactions_per_particle
        implicit none

        integer :: nbaddr ! number of bits in the hashing function

        call timer_start(t_allocate)

        if (allocated(htable)) call deallocate_tree()

        call pepc_status('ALLOCATE TREE')

        call get_number_of_interactions_per_particle(npart, nintmax)
        nintmax = interaction_list_length_factor * nintmax

        !  Space for # table and tree arrays
        !  TODO: need good estimate for max # branches
        size_tree = max(30*nintmax+4*npp,10000)

        if (np_mult>0) then
            nbaddr = int(max(log(1.*size_tree)/log(2.),15.))
            maxaddress = size_tree
        else
            maxaddress = int(abs(np_mult)*10000)
            nbaddr = int(max(log(1.*maxaddress)/log(2.) ,15.))
        endif

        if (num_pe==1) then
            maxleaf = npart
            maxtwig = maxaddress/2
        else if (num_pe.lt.1024) then
            maxleaf = maxaddress/3
            maxtwig = 2*maxleaf
        else
            !  required # twigs increases with P because of branches
            maxleaf = int(maxaddress/(1.+log(1.*num_pe)/3.))
            maxtwig = maxaddress-maxleaf
        endif

        if (maxleaf < npp) then
          DEBUG_WARNING('("maxleaf = ", I0, " < npp+2 = ", I0, ".",/,"Setting maxleaf = npp+2 for now, but expect that to fail during walk. You should increase np_mult.")',maxleaf, npp+2)
          maxleaf = npp + 2
        endif

        hashconst = 2**nbaddr-1
        allocate ( htable(0:maxaddress), free_addr(maxaddress), point_free(0:maxaddress), &
        branch_owner(branch_max_global), pebranch(branch_max_global) )

        ! Allocate memory for tree node properties
        allocate ( tree_nodes(-maxtwig:maxleaf) )
        ! allocate memory for storing number of branches per PE
        allocate ( nbranches(num_pe+2) )

        call timer_stop(t_allocate)

    end subroutine allocate_tree



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !> Deallocates all tree- and multipole-specific arrays
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine deallocate_tree()
        use treevars
        use module_debug
        use module_htable
        implicit none

        call pepc_status('DEALLOCATE TREE')

        deallocate ( htable, free_addr, point_free, branch_key, branch_owner, pebranch )

        ! multipole moments
        deallocate ( tree_nodes )
        ! number of branches per PE
        deallocate ( nbranches )

    end subroutine deallocate_tree



end module module_allocation
