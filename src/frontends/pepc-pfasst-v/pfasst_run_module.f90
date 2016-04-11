! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2016 Juelich Supercomputing Centre, 
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

module pfasst_run_module


contains

  subroutine run_rk4(y_start, delta_t, nsteps)

    use pfasst_helper_module
    use pfasst_calc_module
    use pfasst_transfer_module
    use physvars
    use manipulate_particles
    use files, only : dump

    implicit none

    real(kind=8), intent(in) :: y_start(NvarF), delta_t
    integer, intent(in) :: nsteps

    integer :: i, k, step
    real(kind=8) :: t0, res, err
    real(kind=8) :: fF(NvarF), y_tmp(NvarF)

    !!!! set initial conditions
    y0F = y_start
    call start_timer(TTOTAL)

    y_tmp = 0.

    !!!! time step loop
    do i = 1,nsteps

       call start_timer(TIO)

       step = (i-1)
       t0   = step * delta_t

       call eval_f1(y0F, t0, NvarF, 1, fF)
       y0newF = y0F + (1.0/2.0)*delta_t*fF
       y_tmp = y0F + 1.0/6.0*delta_t*fF

       call eval_f1(y0newF, t0, NvarF, 1, fF)
       y0newF = y0F + (1.0/2.0)*delta_t*fF
       y_tmp = y_tmp + 1.0/3.0*delta_t*fF

       call eval_f1(y0newF, t0, NvarF, 1, fF)
       y0newF = y0F + (1.0/1.0)*delta_t*fF
       y_tmp = y_tmp + 1.0/3.0*delta_t*fF

       call eval_f1(y0newF, t0, NvarF, 1, fF)
       y_tmp = y_tmp + 1.0/6.0*delta_t*fF

       y0F = y_tmp
       y0newF = y_tmp

       call end_timer(TIO, step, echo_timings=echo_timings)

       call pfasst_to_pepc(vortex_particles(1:np), np, y0F(1:NvarF))

       call dump(step,real(t0))

    end do

    call end_timer(TTOTAL, echo_timings=echo_timings)

  end subroutine run_rk4


  subroutine run_rk3(y_start, delta_t, nsteps)

    use pfasst_helper_module
    use pfasst_calc_module
    use pfasst_transfer_module
    use physvars
    use manipulate_particles
    use files, only : dump

    implicit none

    real(kind=8), intent(in) :: y_start(NvarF), delta_t
    integer, intent(in) :: nsteps

    integer :: i, k, step
    real(kind=8) :: t0, res, err
    real(kind=8) :: fF(NvarF)

    !!!! set initial conditions
    y0F = y_start
    call start_timer(TTOTAL)

    !!!! time step loop
    do i = 1,nsteps

       call start_timer(TIO)

       step = (i-1)
       t0   = step * delta_t

       call eval_f1(y0F, t0, NvarF, 1, fF)
       y0newF = y0F + (1.0/3.0)*delta_t*fF

       call eval_f1(y0newF, t0, NvarF, 1, fF)
       y0newF = y0F + (1.0/2.0)*delta_t*fF

       call eval_f1(y0newF, t0, NvarF, 1, fF)
       y0F = y0F + delta_t*fF

       y0newF = y0F

       call end_timer(TIO, step, echo_timings=echo_timings)

       call pfasst_to_pepc(vortex_particles(1:np), np, y0F(1:NvarF))

       call dump(step,real(t0))

    end do

    call end_timer(TTOTAL, echo_timings=echo_timings)

  end subroutine run_rk3

  subroutine run_serial(y_start, delta_t, nsteps)

    use pfasst_helper_module
    use pfasst_calc_module
    use pfasst_transfer_module
    use physvars
    use manipulate_particles
    use files, only : dump

    implicit none

    real(kind=8), intent(in) :: y_start(NvarF), delta_t
    integer, intent(in) :: nsteps

    integer :: i, k, step
    real(kind=8) :: t0, res, err
    real(kind=8) :: fF(NvarF)

    !!!! set initial conditions
    y0F = y_start
    call start_timer(TTOTAL)

    !write(*,*) tend, delta_t, int(tend/delta_t), int(ceiling(tend/delta_t))

    !write(*,*) nsteps
    !!!! time step loop
    do i = 1, nsteps

       !step = my_rank_pfasst + (i-1)*Nproc
       call start_timer(TIO)

       step = i-1
       t0   = step * delta_t

       ! store initial condition as provisional solution
       call start_timer(TINITIAL)
       ySDC_F  = spread(y0F, DIM=2, NCOPIES=nnodesF)
       call eval_f1(y0F, t0, NvarF, 1, fF)
       fSDC_F(:,:,1) = spread(fF, DIM=2, NCOPIES=nnodesF)
       call eval_f2(y0F, t0, NvarF, 1, fF)
       fSDC_F(:,:,2) = spread(fF, DIM=2, NCOPIES=nnodesF)
       y0newF = y0F
       call end_timer(TINITIAL,echo_timings=echo_timings)

       !call start_timer(TIO)
       !call dump(DFINE, step, 0, y0newF)
       !call end_timer(TIO,echo_timings=echo_timings)

       !call echo_error(y0F, t0-delta_t, delta_t, y0newF, NvarF, step, 0)

       ! cycle through iterations
       do k = 1, Niter
          call start_timer(TITERATION)

          ! compute the fine propagator with old y0
          call start_timer(TFINE)
          call sdcsweep(y0newF, t0, delta_t, ySDC_F, fSDC_F, NvarF, NnodesF, 1)
          call end_timer(TFINE, k, 0)

          yendF = ySDC_F(:,NnodesF)

          !call start_timer(TIO)
          !call dump(DFINE, step, k, yendF)
          !call end_timer(TIO,echo_timings=0)

          !call echo_error(y0F, t0, delta_t, yendF, NvarF, step, k)

          !call start_timer(TIO)
          call dump_residual(k, step, delta_t, y0newF, fSDC_F, yendF, my_rank_space)
          !call end_timer(TIO,echo_timings=0)

          !if (parallel > -1) then
          !   call compute_residual(y0newF, delta_t, fSDC_F, yendF, res)
          !   if (res < serial_residual_tol) then
          !      print '("serial residual condition met after: ",i3," iterations",i3)', k
          !      exit
          !   end if
          !end if

          call end_timer(TITERATION, k, echo_timings=echo_timings)
       end do

       y0F = yendF

       call end_timer(TIO, step, echo_timings=echo_timings)

       call pfasst_to_pepc(vortex_particles(1:np), np, y0F(1:NvarF))

       if ((rem_freq .gt. 0) .and. (mod(i,rem_freq)==0)) then

          call remeshing()
          ! TODO: reallocate variables for pfasst
          call pepc_to_pfasst_part(vortex_particles(1:np), np, y0F(1:NvarF))

       end if

       call dump(step,real(t0))

    end do

    call end_timer(TTOTAL, echo_timings=echo_timings)

  end subroutine run_serial


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine run_parallel(y_start, delta_t, tend)

    use pfasst_helper_module
    use pfasst_calc_module
    use pfasst_transfer_module
    use physvars
    use files, only : dump

    implicit none

    ! arguments
    real(kind=8), intent(in) :: y_start(NvarF), delta_t, tend

    ! local stuff
    integer :: nblock, i, j, k, step
    real(kind=8) :: fF(NvarF), t0


    !!!! set initial conditions
    y0F = y_start
    call start_timer(TTOTAL)


    !!!! time "block" loop
    nblock = int(ceiling(tend/delta_t/num_space_instances))
    do i = 1, nblock

       call start_timer(TIO)
       step = my_rank_pfasst + (i-1)*num_space_instances

       t0   = step * delta_t

       y0newF = y0F

       call restrict(y0F, y0G, NvarF, NvarG, 1)
       y0newG = y0G

       !  store initial condition as provisional solution for all F nodes
       call start_timer(TINITIAL)
       ySDC_F  = spread(y0F, DIM=2, NCOPIES=nnodesF)
       call eval_f1(y0F, t0, NvarF, 1, fF)
       fSDC_F(:,:,1) = spread(fF, DIM=2, NCOPIES=nnodesF)
       call eval_f2(y0F, t0, NvarF, 1, fF)
       fSDC_F(:,:,2) = spread(fF, DIM=2, NCOPIES=nnodesF)
       call end_timer(TINITIAL,echo_timings=echo_timings)
       !!!! predictor loop

       ! restrict ySDC_F to get first ySDC_G and compute FAS correction
       call restrict_time_space_fas(t0, delta_t, &
            ySDC_F, ySDC_G, fSDC_F, fSDC_G, NvarF, NvarG, 1, tau_G)

       call start_timer(TPREDICTOR)
       do k = 1, my_rank_pfasst + 1

          ! get new initial value (skip on first iteration)
          if (k > 1) then
               call receive(y0newG, NvarG, my_rank_pfasst-1, k-1)
          end if

          ! sweep with new initial value
          do j = 1, num_coarse_sweeps
             if (use_fas > 0) then
                call sdcsweep(y0newG, t0, delta_t, ySDC_G, fSDC_G, NvarG, &
                     NnodesG, 2, tau_G)
             else
                call sdcsweep(y0newG, t0, delta_t, ySDC_G, fSDC_G, NvarG, &
                     NnodesG, 2)
             end if
          end do
          yendG = ySDC_G(:,NnodesG)

          call send(yendG, NvarG, my_rank_pfasst+1, k)

          !call echo_error(y0newG, t0, delta_t, yendG, NvarG, step, -k)
       end do
       call end_timer(TPREDICTOR,echo_timings=echo_timings)


       !!!! pfasst iterations

       do k = 1, Niter

          call start_timer(TITERATION)

          ! interpolate G to F
          call start_timer(TINTERPOLATE)
          call interpolate_time_space(t0, delta_t, ySDC_F, ySDC_G, fSDC_F, fSDC_G, 1)
          call end_timer(TINTERPOLATE, k, 0)

          !interpolation_order = save_interpolation_order

          if (k == 1) then
             y0newF = ySDC_F(:,1)
             yendF = ySDC_F(:,NnodesF)
             yendG = ySDC_G(:,NnodesG)

             !call dump(DFINE,   step, 0, yendF)
             !call dump(DCOARSE, step, 0, yendG)
             !call echo_error(y0F, t0-delta_t, delta_t, y0newF, NvarF, step, 0)
          end if

          ! fine sweep
          call start_timer(TFINE)
          do j = 1, num_fine_sweeps
             call sdcsweep(y0newF, t0, delta_t, ySDC_F, fSDC_F, NvarF, NnodesF, 1)
          end do
          call end_timer(TFINE, k, 0)

          yendF = ySDC_F(:,NnodesF)

          !call start_timer(TIO)
          !call dump(DFINE,   step, k, yendF)
          !call end_timer(TIO,1)

          ! call compute_residual(y0newF, delta_t, fSDC_F, ySDC_F(:,NnodesF), res)
          ! if (res < pfasst_residual_tol .and. dy0 < pfasst_u0_tol) then
          !    print '("pfasst no-forward conditions met: step: ",i3,"; iter: ",i3)', step+1, k
          ! else
          ! end if

          !call start_timer(TIO)
          call dump_residual(k, step, delta_t, y0newF, fSDC_F, yendF, my_rank_space)
          !call end_timer(TIO,1)

          ! restrict ySDC_F to get ySDC_G guess and compute FAS correction
          call start_timer(TRESTRICT)
          call restrict_time_space_fas(t0, delta_t, &
               ySDC_F, ySDC_G, fSDC_F, fSDC_G, NvarF, NvarG, 1, tau_G)
          call end_timer(TRESTRICT, k,0)

          !call start_timer(TIO)
          !call dump(DFAS, step, k-1, tau_G(:,1))
          !call end_timer(TIO,1)

          ! get new initial value
          call start_timer(TRECEIVE)
          if (my_rank_pfasst > 0) then
             call receive(y0newF, NvarF, my_rank_pfasst-1, 100+k)
             call restrict(y0newF, y0newG, NvarF, NvarG, 1)
          end if
          call end_timer(TRECEIVE, k, 0)

          !if (my_rank_pfasst > 0) &
          !     dy0 = maxval(abs(y0newF-ySDC_F(:,1)))

          ! coarse sweep
          call start_timer(TCOARSE)
          do j = 1, num_coarse_sweeps
             if (use_fas > 0) then
                call sdcsweep(y0newG, t0, delta_t, ySDC_G, fSDC_G, NvarG, NnodesG, 2, tau_G)
             else
                call sdcsweep(y0newG, t0, delta_t, ySDC_G, fSDC_G, NvarG, NnodesG, 2)
             end if
          end do
          call end_timer(TCOARSE, k, 0)

          yendG = ySDC_G(:,NnodesG)

          !call start_timer(TIO)
          !call dump(DCOARSE, step, k, yendG)
          !call end_timer(TIO)

          ! interpolate the G correction to F
          call interpolate(yendF, yendG, NvarF, NvarG, 1)

          ! send new value forward
          call start_timer(TSEND)
          call send(yendF, NvarF, my_rank_pfasst+1, 100+k)
          call end_timer(TSEND, k, 0)

          !call echo_error(y0F, t0, delta_t, yendF, NvarF, step, k)
          call end_timer(TITERATION, k, echo_timings=echo_timings)
       end do

       ! broadcast yendF (non-pipelined time loop)
       if (nblock > 1) &
            call broadcast(yendF, NvarF, n_cpu_pfasst-1)

       y0F = yendF

       call end_timer(TIO, step, echo_timings=echo_timings)

       call pfasst_to_pepc(vortex_particles(1:np), np, y0F(1:NvarF))

       call dump(step,real(t0))

    end do

    call end_timer(TTOTAL, echo_timings=echo_timings)

    !call dump_timings(runtimes, my_rank_pfasst)

  end subroutine run_parallel

  end module pfasst_run_module
