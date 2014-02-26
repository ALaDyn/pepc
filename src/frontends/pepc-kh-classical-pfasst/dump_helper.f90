module dump_helper
  implicit none

  contains

  subroutine perform_all_dumps(step, pepc_pars, physics_pars, time_pars, field_grid, particles)
    use encap
    use physics_helper
    use checkpoint_helper
    use field_helper
    use pepc_helper
    implicit none
    type(physics_pars_t), intent(in) :: physics_pars
    type(pepc_pars_t), intent(in) :: pepc_pars
    type(time_pars_t), intent(in) :: time_pars
    type(t_particle), intent(in) :: particles(:)
    type(field_grid_t), intent(inout) :: field_grid
    integer(kind = 4) :: step

    logical :: do_pdump, do_fdump, do_cdump
    real(kind = 8) :: timer_pio, timer_fcomp, timer_fio, timer_chkpt

    do_pdump = .false.
    do_fdump = .false.
    do_cdump = .false.

    if (pepc_pars%pdump .ne. 0) do_pdump = mod(step, pepc_pars%pdump) .eq. 0
    if (pepc_pars%fdump .ne. 0) do_fdump = mod(step, pepc_pars%fdump) .eq. 0
    if (pepc_pars%cdump .ne. 0) do_cdump = mod(step, pepc_pars%cdump) .eq. 0

    if (do_cdump) then
      call t_start(timer_chkpt)
      call write_checkpoint(pepc_pars, time_pars, step, physics_pars, field_grid, particles)
      call t_stop(timer_chkpt)
    end if

    if(do_fdump) then
      call t_start(timer_fcomp)
!TODO      call compute_field(pepc_pars, field_grid, particles)
      call t_stop(timer_fcomp)

      call t_start(timer_fio)
      call write_field_on_grid(pepc_pars%pepc_comm, time_pars, step, &
          physics_pars, field_grid)
      call t_stop(timer_fio)
    end if

    if (do_pdump) then
       call t_start(timer_pio)
       call write_domain(time_pars, step, particles)
       call write_particles(pepc_pars, time_pars, step, particles)
       call t_stop(timer_pio)
    end if

    call physics_dump(pepc_pars, physics_pars, time_pars, step, particles)

    if (pepc_pars%pepc_comm%root_stdio) then
       if (do_cdump) &
          write(*,'(a,es12.4)') " == time in checkpoint I/O (s)                    : ", timer_chkpt
       if (do_pdump) &
          write(*,'(a,es12.4)') " == time in particle I/O (s)                      : ", timer_pio
       if (do_fdump) then
          write(*,'(a,es12.4)') " == time in field computation (s)                 : ", timer_fcomp
          write(*,'(a,es12.4)') " == time in field I/O (s)                         : ", timer_fio
       end if
    end if

  end subroutine

  subroutine t_start(t)
    use pepc_helper
    implicit none
    real(kind = 8), intent(out) :: t
    t = -get_time()
  end subroutine t_start

  subroutine t_stop(t)
    use pepc_helper
    implicit none
    real(kind = 8), intent(inout) :: t
    t = t + get_time()
  end subroutine t_stop

  subroutine t_resume(t)
    use pepc_helper
    implicit none
    real(kind = 8), intent(inout) :: t
    t = t - get_time()
  end subroutine t_resume

end module dump_helper
