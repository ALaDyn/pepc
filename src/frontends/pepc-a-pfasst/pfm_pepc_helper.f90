module pfm_pepc_helper
    use pfasst
    implicit none
    
contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!
    !! Sets up parameters on each level and defines initial RHS (as "old" u value)
    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine setup_solver(wk, nlevels)
        use pf_mod_mpi
        use pfm_encap
        implicit none
        
        type(app_data_t), pointer, intent(inout) :: wk(:)
        integer, intent(in) :: nlevels
        
        integer :: i

        allocate(wk(nlevels))

        do i = 1, nlevels
          associate (F=>wk(i))
            ! TODO
            ! local number of nodes
            F%nvar = 2500
            allocate(F%particles(F%nvar))
            F%theta = 0.3 + 0.4*(i-1)/(nlevels-1)
            F%mpi_rank = 0
          end associate
        end do

    end subroutine setup_solver


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!
    !! Destroy PMG data structure
    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine pfh_pepc_finalize(wk, nlevels)
        use pfm_encap
        implicit none

        type(app_data_t), intent(inout), pointer :: wk(:)
        integer, intent(in) :: nlevels

        integer :: i

        do i = 1, nlevels
            deallocate(wk(i)%particles)
        end do

        deallocate(wk)

    end subroutine pfh_pepc_finalize

end module pfm_pepc_helper
