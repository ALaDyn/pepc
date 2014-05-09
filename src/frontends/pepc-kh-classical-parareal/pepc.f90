! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
!
! Copyright (C) 2002-2014 Juelich Supercomputing Centre,
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

program pepc
  use module_pepc
  use module_pepc_types
  use module_timings

  use constants
  use encap
  use pepc_helper
  use time_helper
  use field_helper
  use physics_helper
  use checkpoint_helper
  use module_rng

  implicit none

  type(pepc_pars_t) :: pepc_pars
  type(time_pars_t) :: time_pars
  type(physics_pars_t) :: physics_pars
  type(field_grid_t) :: field_grid
  type(t_particle), dimension(:), allocatable :: p

  integer(kind = 4) :: step
  real(kind = 8) :: timer_total, timer_init, timer_step, timer_pcomp, &
    timer_pio, timer_fcomp, timer_fio, timer_dynamics, timer_chkpt

  logical :: root, do_pdump, do_fdump, do_cdump

  character(20) :: mode
  character(16384) :: params_file_name

  call get_command_argument(1, mode)
  call get_command_argument(2, params_file_name)

  ! initialize pepc library and MPI
  call pepc_setup(params_file_name, pepc_pars)

  call t_start(timer_total)
  call t_start(timer_init)

  root = pepc_pars%pepc_comm%mpi_rank.eq.0

  call rng_init(pepc_pars%pepc_comm%mpi_rank + 1)

  if (root) then
    print *, "== [", FRONTEND_NAME, "]"
    print *, "   running on", pepc_pars%pepc_comm%mpi_size, " MPI ranks."
    print *, ""
  end if

  call t_stop(timer_init)

  if (trim(mode) == "init" .or. trim(mode) == "step") then
    call t_resume(timer_init)

    call setup_time(params_file_name, pepc_pars%pepc_comm, time_pars)
    call setup_physics(params_file_name, physics_pars)
    call setup_field_grid(params_file_name, pepc_pars%pepc_comm, field_grid)

    call pepc_prepare(2_kind_dim)
    call t_stop(timer_init)

    if (trim(mode) == "init") then
      call t_resume(timer_init)
      call setup_particles(pepc_pars, physics_pars, p)
      call t_stop(timer_init)

      if (root) write(*,'(a,es12.4)') " == time in setup (s)                            : ", timer_init

      if(root) then
        write(*,*) " "
        write(*,'(a,i12)')    " ====== computing step  :", 0
        write(*,'(a,es12.4)') " ====== simulation time :", 0.0_8
        write(*,*) ""
      end if

      call t_start(timer_pcomp)

      call pepc_particleresults_clear(p)

      call pepc_grow_tree(p)
      call pepc_traverse_tree(p)
      p(:)%results%e(1) = p(:)%results%e(1) * force_const
      p(:)%results%e(2) = p(:)%results%e(2) * force_const
      p(:)%results%pot  = p(:)%results%pot  * force_const

      call pepc_timber_tree()
      call t_stop(timer_pcomp)

      call t_start(timer_chkpt)
      call write_checkpoint(pepc_pars, time_pars, 0, physics_pars, field_grid, p)
      call t_stop(timer_chkpt)

      call physics_dump(pepc_pars, physics_pars, time_pars, step, p)

      call timings_GatherAndOutput(step, 0, step == 0)

      if (root) then
        write(*,'(a,es12.4)') " == time in force computation (s)                 : ", timer_pcomp
        write(*,'(a,es12.4)') " == time in checkpoint I/O (s)                    : ", timer_chkpt
      end if

    else if (trim(mode) == "step") then
      call t_resume(timer_init)
      if (root) then
        print *, ""
        print *, "   pdump   = ", pepc_pars%pdump
        print *, "   fdump   = ", pepc_pars%fdump
        print *, "   cdump   = ", pepc_pars%cdump
        print *, ""
      end if

      call load_particles(pepc_pars, time_pars%nresume, p)
      call t_stop(timer_init)

      if (root) write(*,'(a,es12.4)') " == time in setup (s)                            : ", timer_init

      do_pdump = .false.
      do_fdump = .false.
      do_cdump = .false.

      do step = time_pars%nresume + 1, time_pars%nresume + time_pars%nsteps
        if(root) then
          write(*,*) " "
          write(*,'(a,i12)')    " ====== computing step  :", step
          write(*,'(a,es12.4)') " ====== simulation time :", time_of_step(step, time_pars)
          write(*,*) ""
        end if

        call t_start(timer_step)

        call t_start(timer_dynamics)
        call push_particles(time_pars, physics_pars, p)
        call constrain_particles(physics_pars, p)
        call t_stop(timer_dynamics)

        call t_start(timer_pcomp)

        call pepc_particleresults_clear(p)

        call pepc_grow_tree(p)
        call pepc_traverse_tree(p)
        p(:)%results%e(1) = p(:)%results%e(1) * force_const
        p(:)%results%e(2) = p(:)%results%e(2) * force_const
        p(:)%results%pot  = p(:)%results%pot  * force_const

        call t_stop(timer_pcomp)

        do_pdump = (pepc_pars%pdump .ne. 0) .and. (mod(step, pepc_pars%pdump) .eq. 0)
        do_fdump = (pepc_pars%fdump .ne. 0) .and. (mod(step, pepc_pars%fdump) .eq. 0)
        do_cdump = (pepc_pars%cdump .ne. 0) .and. (mod(step, pepc_pars%cdump) .eq. 0)

        if (do_cdump) then
          call t_start(timer_chkpt)
          call write_checkpoint(pepc_pars, time_pars, step, physics_pars, field_grid, p)
          call t_stop(timer_chkpt)
        end if

        if(do_fdump) then
          call t_start(timer_fcomp)
          call compute_field(pepc_pars, field_grid, p)
          call t_stop(timer_fcomp)

          call t_start(timer_fio)
          call write_field_on_grid(pepc_pars%pepc_comm, time_pars, step, &
            physics_pars, field_grid)
          call t_stop(timer_fio)
        end if

        if (do_pdump) then
          call t_start(timer_pio)
          call write_domain(time_pars, step, p)
          call write_particles(pepc_pars, time_pars, step, p)
          call t_stop(timer_pio)
        end if

        call t_resume(timer_pcomp)
        call pepc_timber_tree()
        call t_stop(timer_pcomp)

        call physics_dump(pepc_pars, physics_pars, time_pars, step, p)

        call timings_GatherAndOutput(step, 0, step == 0)

        call t_stop(timer_step)

        if (root) then
          write(*,'(a,es12.4)') " == time in step (s)                              : ", timer_step
          if (do_cdump) &
              write(*,'(a,es12.4)') " == time in checkpoint I/O (s)                    : ", timer_chkpt
          write(*,'(a,es12.4)') " == time in force computation (s)                 : ", timer_pcomp
          if (do_pdump) &
              write(*,'(a,es12.4)') " == time in particle I/O (s)                      : ", timer_pio
          if (do_fdump) then
              write(*,'(a,es12.4)') " == time in field computation (s)                 : ", timer_fcomp
              write(*,'(a,es12.4)') " == time in field I/O (s)                         : ", timer_fio
          end if
          write(*,'(a,es12.4)') " == time in pusher (s)                            : ", timer_dynamics
        end if

      end do
    end if ! mode == "step"

    deallocate(p)

    call t_stop(timer_total)

    if(root) then
      write(*,*)            " "
      write(*,'(a)')        " ===== finished pepc simulation"
      write(*,'(a,es12.4)') " ===== total run time without setup (s): ", timer_total - timer_init
      write(*,'(a,es12.4)') " ===== total run time with setup (s):    ", timer_total
    end if

    !!! cleanup pepc and MPI
    call pepc_finalize()
  else if (trim(mode) == "correct") then
    print *, "Parareal correction not implemented!"

  else if (trim(mode) == "check") then
    print *, "Parareal convergence check not implemented!"

  else
    print *, "Unknown mode of operation: ", trim(mode)

  end if

contains

  subroutine t_start(t)
    implicit none

    real(kind = 8), intent(out) :: t

    t = -get_time()
  end subroutine t_start

  subroutine t_stop(t)
    implicit none

    real(kind = 8), intent(inout) :: t

    t = t + get_time()
  end subroutine t_stop

  subroutine t_resume(t)
    implicit none

    real(kind = 8), intent(inout) :: t

    t = t - get_time()
  end subroutine t_resume

end program pepc

