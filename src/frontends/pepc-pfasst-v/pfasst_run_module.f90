module pfasst_run_module


contains
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine run_parallel(y_start, dt, tend)

    use pfasst_helper_module
    use pfasst_calc_module
    use pfasst_transfer_module

    implicit none

    ! arguments
    real(kind=8), intent(in) :: y_start(NvarF), dt, tend

    ! local stuff
    integer :: nblock, i, j, k, step
    real(kind=8) :: fF(NvarF), t0


    !!!! set initial conditions
    y0F = y_start
    call start_timer(TTOTAL)


    !!!! time "block" loop
    nblock = int(ceiling(tend/dt/num_space_instances))
    do i = 1, nblock

       step = my_rank_pfasst + (i-1)*num_space_instances

       t0   = step * Dt

       y0newF = y0F

       call restrict(y0F, y0G, NvarF, NvarG, 1)
       y0newG = y0G

       !  store initial condition as provisional solution for all F nodes
       ySDC_F  = spread(y0F, DIM=2, NCOPIES=nnodesF)
       call eval_f1(y0F, t0, NvarF, 1, fF)
       fSDC_F(:,:,1) = spread(fF, DIM=2, NCOPIES=nnodesF)
       call eval_f2(y0F, t0, NvarF, 1, fF)
       fSDC_F(:,:,2) = spread(fF, DIM=2, NCOPIES=nnodesF)

       !!!! predictor loop

       ! restrict ySDC_F to get first ySDC_G and compute FAS correction
       call restrict_time_space_fas(t0, dt, &
            ySDC_F, ySDC_G, fSDC_F, fSDC_G, NvarF, NvarG, 1, tau_G)

       call start_timer(TPREDICTOR)
       do k = 1, my_rank_pfasst + 1

          ! get new initial value (skip on first iteration)
          if (k > 1) &
               call receive(y0newG, NvarG, my_rank_pfasst-1, k-1)

          ! sweep with new initial value
          do j = 1, num_coarse_sweeps
             if (use_fas > 0) then
                call sdcsweep(y0newG, t0, Dt, ySDC_G, fSDC_G, NvarG, &
                     NnodesG, 2, tau_G)
             else
                call sdcsweep(y0newG, t0, Dt, ySDC_G, fSDC_G, NvarG, &
                     NnodesG, 2)
             end if
          end do
          yendG = ySDC_G(:,NnodesG)

          call send(yendG, NvarG, my_rank_pfasst+1, k)
          !call echo_error(y0newG, t0, dt, yendG, NvarG, step, -k)
       end do
       call end_timer(TPREDICTOR,echo_timings=1)


       !!!! pfasst iterations

       do k = 1, Niter

          call start_timer(TITERATION)

          ! interpolate G to F
          call start_timer(TINTERPOLATE)
          call interpolate_time_space(t0, Dt, ySDC_F, ySDC_G, fSDC_F, fSDC_G, 1)
          call end_timer(TINTERPOLATE, k, 1)

          !interpolation_order = save_interpolation_order

          if (k == 1) then
             y0newF = ySDC_F(:,1)
             yendF = ySDC_F(:,NnodesF)
             yendG = ySDC_G(:,NnodesG)

             !call dump(DFINE,   step, 0, yendF)
             !call dump(DCOARSE, step, 0, yendG)
             !call echo_error(y0F, t0-dt, dt, y0newF, NvarF, step, 0)
          end if

          ! fine sweep
          call start_timer(TFINE)
          do j = 1, num_fine_sweeps
             call sdcsweep(y0newF, t0, dt, ySDC_F, fSDC_F, NvarF, NnodesF, 1)
          end do
          call end_timer(TFINE, k, 1)

          yendF = ySDC_F(:,NnodesF)

          !call start_timer(TIO)
          !call dump(DFINE,   step, k, yendF)
          !call end_timer(TIO,1)

          ! call compute_residual(y0newF, Dt, fSDC_F, ySDC_F(:,NnodesF), res)
          ! if (res < pfasst_residual_tol .and. dy0 < pfasst_u0_tol) then
          !    print '("pfasst no-forward conditions met: step: ",i3,"; iter: ",i3)', step+1, k
          ! else
          ! end if

          !call start_timer(TIO)
          !call dump_residual(k, step, Dt, y0newF, fSDC_F, yendF)
          !call end_timer(TIO,1)

          ! restrict ySDC_F to get ySDC_G guess and compute FAS correction
          call start_timer(TRESTRICT)
          call restrict_time_space_fas(t0, dt, &
               ySDC_F, ySDC_G, fSDC_F, fSDC_G, NvarF, NvarG, 1, tau_G)
          call end_timer(TRESTRICT, k,1)

          !call start_timer(TIO)
          !call dump(DFAS, step, k-1, tau_G(:,1))
          !call end_timer(TIO,1)

          ! get new initial value
          call start_timer(TRECEIVE)
          if (my_rank_pfasst > 0) then
             call receive(y0newF, NvarF, my_rank_pfasst-1, 100+k)
             call restrict(y0newF, y0newG, NvarF, NvarG, 1)
          end if
          call end_timer(TRECEIVE, k, 1)

          !if (my_rank_pfasst > 0) &
          !     dy0 = maxval(abs(y0newF-ySDC_F(:,1)))

          ! coarse sweep
          call start_timer(TCOARSE)
          do j = 1, num_coarse_sweeps
             if (use_fas > 0) then
                call sdcsweep(y0newG, t0, dt, ySDC_G, fSDC_G, NvarG, NnodesG, 2, tau_G)
             else
                call sdcsweep(y0newG, t0, dt, ySDC_G, fSDC_G, NvarG, NnodesG, 2)
             end if
          end do
          call end_timer(TCOARSE, k, 1)

          yendG = ySDC_G(:,NnodesG)

          !call start_timer(TIO)
          !call dump(DCOARSE, step, k, yendG)
          !call end_timer(TIO)

          ! interpolate the G correction to F
          call interpolate(yendF, yendG, NvarF, NvarG, 1)

          ! send new value forward
          call start_timer(TSEND)
          call send(yendF, NvarF, my_rank_pfasst+1, 100+k)
          call end_timer(TSEND, k, 1)

          !call echo_error(y0F, t0, dt, yendF, NvarF, step, k)
          call end_timer(TITERATION, k, 1)
       end do

       ! broadcast yendF (non-pipelined time loop)
       if (nblock > 1) &
            call broadcast(yendF, NvarF, n_cpu_pfasst-1)

       y0F = yendF
    end do

    call end_timer(TTOTAL, echo_timings=1)

    !call dump_timings(runtimes, my_rank_pfasst)

  end subroutine run_parallel

  end module pfasst_run_module
