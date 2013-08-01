module pfm_feval
    use pfm_encap
    use pf_mod_mpi
    use module_debug
    implicit none

contains

    !> Initialize feval, i.e. transfer initial u to PFASST data object
    subroutine feval_init(y0, yend, t, ctx)
        use iso_c_binding
        implicit none

        type(app_data_t), intent(out) :: y0, yend
        real(pfdp), intent(in) :: t
        type(c_ptr), intent(in) :: ctx

        type(app_data_t), pointer :: wk

        call pepc_status('|------> feval_init()')
        call c_f_pointer(ctx,wk)

       !TODO
       ! allocate(y0%array(F%m_start+F%xmin:F%m_end+F%xmax,F%n_start+F%ymin:F%n_end+F%ymax,F%o_start+F%zmin:F%o_end+F%zmax))
       ! allocate(yend%array(F%m_start+F%xmin:F%m_end+F%xmax,F%n_start+F%ymin:F%n_end+F%ymax,F%o_start+F%zmin:F%o_end+F%zmax))
        !y0%array = 0.0_pfdp
        !yend%array = 0.0_pfdp

        !call exact(y0,t)

    end subroutine feval_init


    !> Eval. exact solution of PDE
    subroutine exact(y,t)
        implicit none

        type(app_data_t), intent(inout) :: y
        real(pfdp), intent(in) :: t

        call pepc_status('|------> exact()')
        !TODO
        ! y%array(i,j,k) = sin(pi*(i+1)*F%h)*sin(pi*(j+1)*F%h)*sin(pi*(k+1)*F%h)*cos(t)

    end subroutine exact


    !> End feval
    subroutine feval_finalize(y0, yend)
        implicit none

        type(app_data_t), intent(inout) :: y0, yend

        call pepc_status('|------> feval_finalize()')
        ! TODO
        !deallocate(y0%array)
        !deallocate(yend%array)

    end subroutine feval_finalize



    subroutine mult_by_stencil(datin, datout, stencil)
        implicit none
        type(app_data_t), intent(inout) :: datin
        type(app_data_t), intent(inout) :: datout
        real(pfdp), dimension(1:), intent(in) :: stencil

        call pepc_status('|------> mult_by_stencil()')
        ! TODO
        !call exchange_halo(datin)
        !
        !datout%array = 0
        !do k = F%o_start, F%o_end
        !    do j = F%n_start, F%n_end
        !        do i = F%m_start, F%m_end
        !            do l = 1, size(stencil)
        !                datout%array(i,j,k) = datout%array(i,j,k) + stencil(l) * datin%array(i+F%xoff(l),j+F%yoff(l),k+F%zoff(l))
        !            end do
        !        end do
        !    end do
        !end do
        !
        !call MPI_BARRIER(F%pmg_comm%mpi_comm, ierr)

    end subroutine mult_by_stencil


      ! Evaluate the explicit function at y, t.
    subroutine eval_f1(yptr, t, level, ctx, f1ptr)
        implicit none
        type(c_ptr), intent(in), value :: yptr, f1ptr, ctx
        real(pfdp),  intent(in)        :: t
        integer,     intent(in)        :: level

        type(app_data_t), pointer :: y, f1

        call pepc_status('|------> eval_f1()')
        call c_f_pointer(yptr,y)
        call c_f_pointer(f1ptr,f1)

        !TODO
        ! f1%array(i,j,k) = -sin(pi*(i+1)*F%h)*sin(pi*(j+1)*F%h)*sin(pi*(k+1)*F%h)*(sin(t)-3*F%lambda*pi**2*cos(t))
        !call MPI_BARRIER(F%pmg_comm%mpi_comm, ierr)

    end subroutine eval_f1


    !> Evaluates implicit part of f, i.e. calculates Au
    subroutine eval_f2(yptr, t, level, ctx, f2ptr)
        type(c_ptr), intent(in), value :: yptr, f2ptr, ctx
        real(pfdp),  intent(in)        :: t
        integer,     intent(in)        :: level

        integer :: ierr
        type(app_data_t), pointer :: y, f2

        call pepc_status('|------> eval_f2()')
        call c_f_pointer(yptr,y)
        call c_f_pointer(f2ptr,f2)

        ! TODO
        !call mult_by_stencil(y,f2,F%stencil_lap)
        !f2%array = F%lambda/(F%h**2) * f2%array
        !call MPI_BARRIER(F%pmg_comm%mpi_comm, ierr)

    end subroutine eval_f2


    !> Computes y_new by solving the SDC equation, returns also f2(y_new)
    subroutine comp_f2(yptr, t, dt, rhsptr, level, ctx, f2ptr)
        type(c_ptr), intent(in), value :: yptr, rhsptr, f2ptr, ctx
        real(pfdp),  intent(in)        :: t, dt
        integer,     intent(in)        :: level

        integer :: ierr
        type(app_data_t), pointer :: y, f2, rhs

        call pepc_status('|------> comp_f2()')
        call c_f_pointer(yptr,y)
        call c_f_pointer(f2ptr,f2)
        call c_f_pointer(rhsptr,rhs)

        ! Solve system with rhs for u

        !TODO
        !call mg_fortran(y%array, rhs%array, y%pmg_par%mem_out_ptr, rhs%pmg_par%mem_in_ptr, &
        !                F%maxiter, F%tol, F%m, F%n, F%o, F%m_start, F%m_end, F%n_start, F%n_end, F%o_start, F%o_end,&
        !                F%periodic_integer, F%nu1, F%nu2, F%omega, F%stencil_size, F%stencil_mass - dt*F%lambda/(F%h**2)*F%stencil_lap, &
        !                F%xoff, F%yoff, F%zoff, F%xmin, F%ymin, F%zmin, F%xmax, F%ymax, F%zmax, F%echo_pmg, F%pmg_comm%mpi_comm)
        !
        !
        !if (F%use_mass) then
        !
        !    call mult_by_stencil(y,f2,F%stencil_mass)
        !    ! Use rhs for faster final calculation of f2
        !    f2%array = (f2%array - rhs%array) / dt
        !
        !else
        !
        !    ! Use rhs for faster final calculation of f2
        !    f2%array = (y%array - rhs%array) / dt
        !
        !end if
        !
        !call MPI_BARRIER(F%pmg_comm%mpi_comm, ierr)

    end subroutine comp_f2

end module pfm_feval
