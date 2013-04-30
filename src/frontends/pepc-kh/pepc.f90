! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2013 Juelich Supercomputing Centre, 
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
  ! pepc modules
   use module_pepc
   use module_pepc_types
   use module_timings
  
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

   call t_start(timer_total)
   call t_start(timer_init)

   ! initialize pepc library and MPI
   call pepc_setup(p, pepc_pars)

   root = pepc_pars%pepc_comm%mpi_rank.eq.0

   call rng_init(pepc_pars%pepc_comm%mpi_rank + 1)

   if (root) then
      print *, "== [pepc-kh]"
      print *, "   running on", pepc_pars%pepc_comm%mpi_size, " MPI ranks."
      print *, "   pdump   = ", pepc_pars%pdump
      print *, "   fdump   = ", pepc_pars%fdump
      print *, "   cdump   = ", pepc_pars%cdump
      print *, ""
   end if

   call setup_time(time_pars, pepc_pars%pepc_comm)
   call setup_physics(physics_pars, time_pars, p, pepc_pars)
   call setup_field_grid(field_grid, pepc_pars%pepc_comm)

   call pepc_prepare(2_kind_dim)

   call t_stop(timer_init)

   if (root) write(*,'(a,es12.4)') " == time in setup (s)                            : ", timer_init

   do_pdump = .false.
   do_fdump = .false.
   do_cdump = .false.

   do step = time_pars%nresume, time_pars%nsteps
    if(root) then
      write(*,*) " "
      write(*,'(a,i12)')    " ====== computing step  :", step
      write(*,'(a,es12.4)') " ====== simulation time :", step * time_pars%dt
      write(*,*) ""
    end if
    
    call t_start(timer_step)

    if (pepc_pars%pdump .ne. 0) do_pdump = mod(step, pepc_pars%pdump) .eq. 0
    if (pepc_pars%fdump .ne. 0) do_fdump = mod(step, pepc_pars%fdump) .eq. 0
    if (pepc_pars%cdump .ne. 0) do_cdump = mod(step, pepc_pars%cdump) .eq. 0

    if (do_cdump) then
      call t_start(timer_chkpt)
      call write_checkpoint(pepc_pars, time_pars, step, physics_pars, field_grid, p)
      call t_stop(timer_chkpt)
    end if

    call t_start(timer_pcomp)

    call pepc_particleresults_clear(p)

    call pepc_grow_tree(p)
    pepc_pars%npp = size(p, kind=kind(pepc_pars%npp))
    call pepc_traverse_tree(p)
    p(:)%results%e(1) = p(:)%results%e(1) * force_const
    p(:)%results%e(2) = p(:)%results%e(2) * force_const
    p(:)%results%pot  = p(:)%results%pot  * force_const

    call apply_external_field()

    call t_stop(timer_pcomp)

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

    call pepc_timber_tree()
    !call pepc_restore_particles(pepc_pars%npp, p)

    call physics_dump(pepc_pars, physics_pars, time_pars, step, p)
        
    call t_start(timer_dynamics)
    call push_particles(pepc_pars, time_pars, physics_pars, p)
    call constrain_particles(pepc_pars, time_pars, physics_pars, p)
    call t_stop(timer_dynamics)

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

