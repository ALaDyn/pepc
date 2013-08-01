module pmg_helper

    use pfasst

    type pmg_nml_t
        integer, dimension(:), allocatable :: m           ! number of inner nodes (one dim only)
        integer, dimension(:), allocatable :: order       ! FD order
        integer :: echo_pmg    ! Verbose level for PMG
        integer :: pnodes      ! number of nodes for prolongation
        real(pfdp) :: lambda
    end type pmg_nml_t

    real(pfdp), dimension(19), parameter :: lap_4c = &
    1.0_pfdp/6.0_pfdp* (/ -24.0_pfdp, 2.0_pfdp, 2.0_pfdp, 2.0_pfdp, 2.0_pfdp, 2.0_pfdp, 2.0_pfdp,&
    1.0_pfdp, 1.0_pfdp, 1.0_pfdp, 1.0_pfdp, 1.0_pfdp, 1.0_pfdp,&
    1.0_pfdp, 1.0_pfdp, 1.0_pfdp, 1.0_pfdp, 1.0_pfdp, 1.0_pfdp /)

    real(pfdp), dimension(7),  parameter :: lap_2 = &
    1.0_pfdp/1.0_pfdp* (/ -6.0_pfdp, 1.0_pfdp, 1.0_pfdp, 1.0_pfdp, 1.0_pfdp, 1.0_pfdp, 1.0_pfdp /)

    real(pfdp), dimension(19), parameter :: mass_4c = &
    1.0_pfdp/12.0_pfdp*(/ 6.0_pfdp, 1.0_pfdp, 1.0_pfdp, 1.0_pfdp, 1.0_pfdp, 1.0_pfdp, 1.0_pfdp, &
    0.0_pfdp, 0.0_pfdp, 0.0_pfdp, 0.0_pfdp, 0.0_pfdp, 0.0_pfdp, &
    0.0_pfdp, 0.0_pfdp, 0.0_pfdp, 0.0_pfdp, 0.0_pfdp, 0.0_pfdp /)
    real(pfdp), dimension(7),  parameter :: mass_2 = &
    1.0_pfdp/1.0_pfdp* (/ 1.0_pfdp, 0.0_pfdp, 0.0_pfdp, 0.0_pfdp, 0.0_pfdp, 0.0_pfdp, 0.0_pfdp /)



contains
  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!
    !! Initialize PMG for a given spatial MPI communicator -> returns cartesian MPI structure
    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine init_pmg(pmg_comm, pmg_nml, nlevels, MPI_COMM_SPACE)
        use pf_mod_mpi
        use encap
        implicit none

        type(pmg_comm_t), intent(out) :: pmg_comm
        type(pmg_nml_t), intent(out) :: pmg_nml
        integer, intent(in) :: nlevels, MPI_COMM_SPACE

        logical, dimension(1:3) :: mpi_periods
        integer :: mpi_err

        call MPI_COMM_RANK(MPI_COMM_SPACE, pmg_comm%mpi_rank, mpi_err)
        call MPI_COMM_SIZE(MPI_COMM_SPACE, pmg_comm%mpi_size, mpi_err)

        ! create cartesian communicator
        pmg_comm%mpi_dims = 0
        mpi_periods = .false.
        call MPI_DIMS_CREATE(pmg_comm%mpi_size, 3, pmg_comm%mpi_dims, mpi_err)
        call MPI_CART_CREATE(MPI_COMM_SPACE, 3, pmg_comm%mpi_dims, mpi_periods, .true., &
        pmg_comm%mpi_comm, mpi_err)

        ! get rank and coordinates in cartesian communicator
        call MPI_COMM_RANK(pmg_comm%mpi_comm, pmg_comm%mpi_rank, mpi_err)
        call MPI_COMM_SIZE(pmg_comm%mpi_comm, pmg_comm%mpi_size, mpi_err)
        call MPI_CART_COORDS(pmg_comm%mpi_comm, pmg_comm%mpi_rank, 3, pmg_comm%mpi_coords, mpi_err)

        call read_in_pmg_params(pmg_nml, pmg_comm%mpi_rank, pmg_comm%mpi_comm, nlevels)
        
    end subroutine init_pmg


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!
    !! Sets up PMG parameter on each level and defines initial RHS (as "old" u value)
    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine setup_solver(wk, pmg_comm, pmg_nml, nlevels, pf_tol)
        use pf_mod_mpi
        use encap
        implicit none
        
        type(app_data_t), pointer, intent(inout) :: wk(:)
        type(pmg_comm_t), intent(in) :: pmg_comm
        type(pmg_nml_t), intent(inout) :: pmg_nml
        integer, intent(in) :: nlevels
        real(pfdp), intent(in) :: pf_tol
        
        integer :: i, nc, nn, i1, i2, i3, den, num, mpi_err
        real(pfdp), allocatable :: pmat_tmp(:,:)
        type(pmg_pars_t), pointer :: F

        allocate(wk(nlevels))

        ! Compute weights for prolongation
        nc = pmg_nml%pnodes-1
        nn = pmg_nml%pnodes
        allocate(pmat_tmp(nc,nn))
        pmat_tmp = 0

        do i1 = 1,nc
            do i2 = 1,nn
                den = 1
                num = 1
                do i3 = 1,nn
                    if (i2.ne.i3) then
                        den = den*(2*(i2-1) - 2*(i3-1))
                        num = num*(2*(i1-1)+1 - 2*(i3-1))
                    end if
                end do
                pmat_tmp(i1,i2) = 1.0*num/den
            end do
        end do

        do i = 1, nlevels

            allocate(wk(i)%pmg_par)

            F => wk(i)%pmg_par

            F%pmg_comm = pmg_comm

            F%pnodes = pmg_nml%pnodes

            allocate(F%pmat(nc,nn))
            F%pmat = pmat_tmp

            ! solver parameters
            F%nu1 = 2
            F%nu2 = 2
            F%omega = 0.6D0
            F%maxiter = 50
            F%periodic_integer = 0
            F%echo_pmg = pmg_nml%echo_pmg
            F%lambda = pmg_nml%lambda

            ! global number of nodes per direction, scaled by level (TODO: non-cubic domains?)
            F%m = pmg_nml%m(i)
            F%n = pmg_nml%m(i)
            F%o = pmg_nml%m(i)

            ! calculate local grid coordinates
            F%m_start = (F%m/pmg_comm%mpi_dims(1))*pmg_comm%mpi_coords(1) + min(pmg_comm%mpi_coords(1),mod(F%m,pmg_comm%mpi_dims(1)))
            F%m_end = (F%m/pmg_comm%mpi_dims(1))*(pmg_comm%mpi_coords(1)+1) + min(pmg_comm%mpi_coords(1)+1,mod(F%m,pmg_comm%mpi_dims(1))) - 1
            F%n_start = (F%n/pmg_comm%mpi_dims(2))*pmg_comm%mpi_coords(2) + min(pmg_comm%mpi_coords(2),mod(F%n,pmg_comm%mpi_dims(2)))
            F%n_end = (F%n/pmg_comm%mpi_dims(2))*(pmg_comm%mpi_coords(2)+1) + min(pmg_comm%mpi_coords(2)+1,mod(F%n,pmg_comm%mpi_dims(2))) - 1
            F%o_start = (F%o/pmg_comm%mpi_dims(3))*pmg_comm%mpi_coords(3) + min(pmg_comm%mpi_coords(3),mod(F%o,pmg_comm%mpi_dims(3)))
            F%o_end = (F%o/pmg_comm%mpi_dims(3))*(pmg_comm%mpi_coords(3)+1) + min(pmg_comm%mpi_coords(3)+1,mod(F%o,pmg_comm%mpi_dims(3))) - 1

            if (F%m_end-F%m_start+1 .le. F%pnodes/2+1) then
                write(*,*) 'ERROR: grid to small on level',i,' and direction m:',F%m_start,F%m_end,'.... aborting'
                call MPI_ABORT(MPI_COMM_WORLD,0,mpi_err)
            end if
            if (F%n_end-F%n_start+1 .le. F%pnodes/2+1) then
                write(*,*) 'ERROR: grid to small on level',i,' and direction n:',F%n_start,F%n_end,'.... aborting'
                call MPI_ABORT(MPI_COMM_WORLD,0,mpi_err)
            end if
            if (F%o_end-F%o_start+1 .le. F%pnodes/2+1) then
                write(*,*) 'ERROR: grid to small on level',i,' and direction o:',F%o_start,F%o_end,'.... aborting'
                call MPI_ABORT(MPI_COMM_WORLD,0,mpi_err)
            end if

            ! local number of nodes
            F%nvar = (F%m_end-F%m_start+1) * (F%n_end-F%n_start+1) * (F%o_end-F%o_start+1)

            ! grid spacing
            F%h = 1.0D0/dble(min(F%m+1,F%n+1,F%o+1))

            ! Set FD order on the levels
            F%order = pmg_nml%order(i)

            ! Set stencil depending of the prescribed FD order of the current level
            select case(F%order)

                case(41) ! compact 4th order stencil
                    F%stencil_size = 19
                    allocate(F%xoff(F%stencil_size),F%yoff(F%stencil_size),F%zoff(F%stencil_size))
                    F%xoff(1:F%stencil_size) = (/ 0, 1, -1, 0,  0, 0,  0,  1, -1, 1, -1, 1,  1, -1, -1, 0,  0,  0,  0 /)
                    F%yoff(1:F%stencil_size) = (/ 0, 0,  0, 1, -1, 0,  0, -1, -1, 1,  1, 0,  0,  0,  0, 1,  1, -1, -1 /)
                    F%zoff(1:F%stencil_size) = (/ 0, 0,  0, 0,  0, 1, -1,  0,  0, 0,  0, 1, -1,  1, -1, 1, -1,  1, -1 /)
                    allocate(F%stencil_lap(F%stencil_size))
                    F%stencil_lap = lap_4c
                    allocate(F%stencil_mass(F%stencil_size))
                    F%stencil_mass = mass_4c
                    F%use_mass = .true.
                case default ! 2nd order stencil
                    F%stencil_size = 7
                    allocate(F%xoff(F%stencil_size),F%yoff(F%stencil_size),F%zoff(F%stencil_size))
                    F%xoff(1:F%stencil_size) = (/ 0, 1, -1, 0,  0, 0,  0 /)
                    F%yoff(1:F%stencil_size) = (/ 0, 0,  0, 1, -1, 0,  0 /)
                    F%zoff(1:F%stencil_size) = (/ 0, 0,  0, 0,  0, 1, -1 /)
                    allocate(F%stencil_lap(F%stencil_size))
                    F%stencil_lap = lap_2
                    allocate(F%stencil_mass(F%stencil_size))
                    F%stencil_mass = mass_2
                    F%use_mass = .false.
            end select

            ! fix tolerance for PMG
            F%tol = pf_tol/10

            ! define halo extension
            F%xmin = min(minval(F%xoff),-F%pnodes/2)
            F%xmax = max(maxval(F%xoff),F%pnodes/2)
            F%ymin = min(minval(F%yoff),-F%pnodes/2)
            F%ymax = max(maxval(F%yoff),F%pnodes/2)
            F%zmin = min(minval(F%zoff),-F%pnodes/2)
            F%zmax = max(maxval(F%zoff),F%pnodes/2)

            call f_cuboid_alloc(F%m_end - F%m_start + 1, F%n_end - F%n_start + 1, F%o_end - F%o_start + 1, F%mem_in_ptr)
            call f_cuboid_alloc(F%m_end - F%m_start + 1, F%n_end - F%n_start + 1, F%o_end - F%o_start + 1, F%mem_out_ptr)

            allocate(wk(i)%array(F%m_start+F%xmin:F%m_end+F%xmax, F%n_start+F%ymin:F%n_end+F%ymax, F%o_start+F%zmin:F%o_end+F%zmax))

        end do

        deallocate(pmg_nml%m, pmg_nml%order, pmat_tmp)

    end subroutine setup_solver


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!
    !! Destroy PMG data structure
    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine pmg_finalize(wk, nlevels)
        use encap
        implicit none

        type(app_data_t), intent(inout), pointer :: wk(:)
        integer, intent(in) :: nlevels

        integer :: i

        do i = 1, nlevels

            call f_cuboid_free(wk(i)%pmg_par%m_end - wk(i)%pmg_par%m_start + 1, &
                               wk(i)%pmg_par%n_end - wk(i)%pmg_par%n_start + 1, &
                               wk(i)%pmg_par%o_end - wk(i)%pmg_par%o_start + 1, &
                               wk(i)%pmg_par%mem_in_ptr)
            call f_cuboid_free(wk(i)%pmg_par%m_end - wk(i)%pmg_par%m_start + 1, &
                               wk(i)%pmg_par%n_end - wk(i)%pmg_par%n_start + 1, &
                               wk(i)%pmg_par%o_end - wk(i)%pmg_par%o_start + 1, &
                               wk(i)%pmg_par%mem_out_ptr)

            deallocate(wk(i)%array)

        end do

        deallocate(wk)

    end subroutine pmg_finalize


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!
    !! Exchange halo for a given array (must be of u-shape, though, i.e. PFASST data needs conversion first)
    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine exchange_halo(sol)
        use pf_mod_mpi
        use encap
        implicit none

        type(app_data_t), intent(inout) :: sol

        integer :: left_rank, right_rank, down_rank, up_rank, front_rank, back_rank
        integer :: mpi_err, cnt

        type(pmg_pars_t) :: F

        F = sol%pmg_par

        call MPI_CART_SHIFT(F%pmg_comm%mpi_comm, 0, 1, left_rank, right_rank, mpi_err)
        call MPI_CART_SHIFT(F%pmg_comm%mpi_comm, 1, 1, down_rank, up_rank, mpi_err)
        call MPI_CART_SHIFT(F%pmg_comm%mpi_comm, 2, 1, back_rank, front_rank, mpi_err)


        !TODO: This code only works if the compiler creates temp arrays for the MPI calls (as gfortran does)

        ! 1. step: send and receive halo in x direction, first left, then right
        cnt = F%xmax * (F%o_end - F%o_start + 1) * (F%n_end - F%n_start + 1)
        call MPI_SENDRECV(sol%array(F%m_start:F%m_start+F%xmax-1, F%n_start:F%n_end, F%o_start:F%o_end), cnt, MPI_DOUBLE_PRECISION, left_rank, 1, &
        sol%array(F%m_end+1:F%m_end+F%xmax, F%n_start:F%n_end, F%o_start:F%o_end), cnt, MPI_DOUBLE_PRECISION, right_rank, 1, F%pmg_comm%mpi_comm, MPI_STATUS_IGNORE, mpi_err)
        cnt = abs(F%xmin) * (F%o_end - F%o_start + 1) * (F%n_end - F%n_start + 1)
        call MPI_SENDRECV(sol%array(F%m_end+F%xmin+1:F%m_end, F%n_start:F%n_end, F%o_start:F%o_end), cnt, MPI_DOUBLE_PRECISION, right_rank, 1, &
        sol%array(F%m_start+F%xmin:F%m_start-1, F%n_start:F%n_end, F%o_start:F%o_end), cnt, MPI_DOUBLE_PRECISION, left_rank, 1, F%pmg_comm%mpi_comm, MPI_STATUS_IGNORE, mpi_err)

        ! 2. step: send and receive halo in y direction (now incl. x-halo), first down, then up
        cnt = (F%m_end+F%xmax - (F%m_start+F%xmin) + 1) * F%ymax * (F%o_end-F%o_start+1)
        call MPI_SENDRECV(sol%array(F%m_start+F%xmin:F%m_end+F%xmax, F%n_start:F%n_start+F%ymax-1, F%o_start:F%o_end), cnt, MPI_DOUBLE_PRECISION, down_rank, 1, &
        sol%array(F%m_start+F%xmin:F%m_end+F%xmax, F%n_end+1:F%n_end+F%ymax, F%o_start:F%o_end), cnt, MPI_DOUBLE_PRECISION, up_rank, 1, F%pmg_comm%mpi_comm, MPI_STATUS_IGNORE, mpi_err)
        cnt = (F%m_end+F%xmax - (F%m_start+F%xmin) + 1) * abs(F%ymin) * (F%o_end-F%o_start+1)
        call MPI_SENDRECV(sol%array(F%m_start+F%xmin:F%m_end+F%xmax, F%n_end+F%ymin+1:F%n_end, F%o_start:F%o_end), cnt, MPI_DOUBLE_PRECISION, up_rank, 1, &
        sol%array(F%m_start+F%xmin:F%m_end+F%xmax, F%n_start+F%ymin:F%n_start-1, F%o_start:F%o_end), cnt, MPI_DOUBLE_PRECISION, down_rank, 1, F%pmg_comm%mpi_comm, MPI_STATUS_IGNORE, mpi_err)

        ! 3. step: send and receive halo in z direction (now incl. x- and y-halo), first back, then front
        cnt = (F%m_end+F%xmax - (F%m_start+F%xmin) + 1) * (F%n_end+F%ymax - (F%n_start+F%ymin) + 1) * F%zmax
        call MPI_SENDRECV(sol%array(F%m_start+F%xmin:F%m_end+F%xmax, F%n_start+F%ymin:F%n_end+F%ymax, F%o_start:F%o_start+F%zmax-1), cnt, MPI_DOUBLE_PRECISION, back_rank, 1, &
        sol%array(F%m_start+F%xmin:F%m_end+F%xmax, F%n_start+F%ymin:F%n_end+F%ymax, F%o_end+1:F%o_end+F%zmax), cnt, MPI_DOUBLE_PRECISION, front_rank, 1, F%pmg_comm%mpi_comm, MPI_STATUS_IGNORE, mpi_err)
        cnt = (F%m_end+F%xmax - (F%m_start+F%xmin) + 1) * (F%n_end+F%ymax - (F%n_start+F%ymin) + 1) * abs(F%zmin)
        call MPI_SENDRECV(sol%array(F%m_start+F%xmin:F%m_end+F%xmax, F%n_start+F%ymin:F%n_end+F%ymax, F%o_end+F%zmin+1:F%o_end), cnt, MPI_DOUBLE_PRECISION, front_rank, 1, &
        sol%array(F%m_start+F%xmin:F%m_end+F%xmax, F%n_start+F%ymin:F%n_end+F%ymax, F%o_start+F%zmin:F%o_start-1), cnt, MPI_DOUBLE_PRECISION, back_rank, 1, F%pmg_comm%mpi_comm, MPI_STATUS_IGNORE, mpi_err)

    end subroutine exchange_halo

    subroutine dump_grid(sol,iter)
        use iso_c_binding
        use encap
        implicit none
        type(app_data_t), intent(in) :: sol
        integer, intent(in) :: iter

        type(pmg_pars_t) :: F
        integer :: i, j, k
        character(50) :: resfile

        F = sol%pmg_par
        write(resfile,'(a,i2.2,a,i6.6,a)') "pmg_data/results_", iter, "_", F%pmg_comm%mpi_rank,".dat"
        open(11,file=resfile)
        do k = F%o_start, F%o_end
            do j = F%n_start, F%n_end
                do i = F%m_start, F%m_end
                    write(11,'(3i4,e42.30)') i,j,k,sol%array(i,j,k)
                end do
            end do
        end do
        close(11)

    end subroutine dump_grid


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!
    !! Simple read-in routine, using the first command line argument
    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine read_in_pmg_params(pmg_namelist, rank, comm, nlevels)
        use pf_mod_mpi
        implicit none
        
        integer, intent(in) :: rank, comm, nlevels
        type(pmg_nml_t), intent(out) :: pmg_namelist

        logical :: available
        character(len=255) :: file_name
        integer :: ierr
        integer, parameter :: para_file_id = 10

        ! variables for pmg namelist
        integer, dimension(:), allocatable :: m, order
        integer :: echo_pmg = 0
        integer :: pnodes = 0
        real(pfdp) :: lambda = 0.0
        namelist /pmg_nml/ m, order, echo_pmg, pnodes, lambda

        ! rank 0 reads in first command line argument
        available = .false.
        if (rank .eq. 0) then
            if( COMMAND_ARGUMENT_COUNT() .ne. 0 ) then
                call GET_COMMAND_ARGUMENT(1, file_name)
                available = .true.
                if(rank .eq. 0) write(*,*) "found parameter file: ", file_name
            end if
        end if

        ! broadcast file name, read actual inputs from namelist file (name ist fixed, sorry!)
        call MPI_BCAST( available, 1, MPI_LOGICAL, 0, comm, ierr )
        if (available) then

            allocate(order(nlevels), m(nlevels), pmg_namelist%m(nlevels), pmg_namelist%order(nlevels))

            call MPI_BCAST( file_name, 255, MPI_CHARACTER, 0, comm, ierr )
            open(para_file_id,file=trim(file_name),action='read')
            rewind(para_file_id)
            read(para_file_id, NML=pmg_nml)
            close(para_file_id)

            pmg_namelist%m = m
            pmg_namelist%order = order
            pmg_namelist%pnodes = pnodes
            pmg_namelist%echo_pmg = echo_pmg
            pmg_namelist%lambda = lambda

            deallocate(order, m)

        end if

    end subroutine read_in_pmg_params


end module pmg_helper
