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
  use module_debug
  ! frontend helper routines
  use pepcboris_helper
  use pepcboris_integrator

  use pf_mod_verlet, only: pf_verlet_create, pf_verlet_destroy
  use pfm_helper
  use pfm_encap
  use pfm_feval
  use pfm_hooks
  use pf_mod_comm_mpi, only: pf_mpi_create, pf_mpi_destroy, pf_mpi_setup
  use pf_mod_pfasst, only: pf_pfasst_create, pf_pfasst_setup, pf_pfasst_destroy
  use pf_mod_parallel, only: pf_pfasst_run
  use pf_mod_hooks

  implicit none

  ! variables for pfasst
  integer(kind_default):: MPI_COMM_SPACE, MPI_COMM_TIME, mpi_err
  type(pf_pfasst_t) :: pf
  type(pf_nml_t) :: pf_nml
  type(pf_comm_t) :: tcomm
  type(pf_encap_t), target :: encap
  type(pf_sweeper_t), target :: sweeper
  type(level_params_t), pointer :: level_params(:)
  type(app_data_t), pointer :: y0, yend

  ! particle data (position, velocity, mass, charge)
  type(t_particle), allocatable, target :: particles(:)
  integer :: step, i

  ! Take care of communication stuff
  call pfm_init_pfasst(pf_nml, MPI_COMM_SPACE, MPI_COMM_TIME)
  ! initialize pepc library and MPI
  pepcboris_nml%comm=MPI_COMM_SPACE
  call pepc_initialize('pepc-boris-pfasst', pepcboris_nml%rank, pepcboris_nml%nrank, .false., db_level_in=DBG_STATUS, comm=pepcboris_nml%comm)
  ! this is not an MPI-parallel application
  if (pepcboris_nml%nrank > 1) then
    DEBUG_ERROR(*, "This is not an MPI-parallel application")
  endif
  ! frontend parameter initialization, particle configuration etc.
  call pepcboris_init(pepcboris_nml, particles, dt=pf_nml%tend/pf_nml%nsteps/(pf_nml%nnodes(pf_nml%nlevels)-1), nt=pf_nml%nsteps*(pf_nml%nnodes(pf_nml%nlevels)-1))  ! we use the finest (i.e. highest) level here
  ! commit all internal pepc variables
  call pepc_prepare(dim)
  ! prepare table with level-dependent parameters
  call pfm_setup_solver_level_params(particles, level_params, pf_nml, dim, pepcboris_nml%rank, MPI_COMM_SPACE)
  ! initial potential will be needed for energy computation
  call eval_force(particles, level_params(pf_nml%nlevels), pepcboris_nml, step=0, comm=MPI_COMM_SPACE, clearresults=.true.) ! again, use parameters of finest level

  select case (pepcboris_nml%workingmode)
    case (WM_PFASST)
      ! Set up PFASST object
      call pf_mpi_create(tcomm, MPI_COMM_TIME)
      call pf_pfasst_create(pf, tcomm, pf_nml%nlevels)

      call pfm_encap_init(encap, particles)
      call pf_verlet_create(sweeper, calc_Efield, build_rhs, impl_solver)
      call pfm_fill_pfasst_object(pf, encap, sweeper, pf_nml, level_params)

      ! Initial conditions for pfasst
      call feval_init(y0, yend, pf_nml%nlevels, pf%levels(pf_nml%nlevels)%levelctx, encap%encapctx, particles)

      call pf_mpi_setup(tcomm, pf)
      call pf_pfasst_setup(pf)

      ! Add user-defined calls, e.g. diagnostics, here
      call pf_add_hook(pf, pf_nml%nlevels, PF_POST_STEP, dump_particles_hook)
      call pf_add_hook(pf, pf_nml%nlevels, PF_PRE_STEP, dump_particles_hook)

      ! Here we go       pfasst-object, initial value, dt, t_end, number of steps, final solution
      call pf_pfasst_run(pf, c_loc(y0), pf_nml%tend/pf_nml%nsteps, pf_nml%tend, nsteps=pf_nml%nsteps, qend=c_loc(yend))

      ! Remove everything (hopefully)
      call pf_mpi_destroy(tcomm)
      call feval_finalize(y0,yend)
      call pf_pfasst_destroy(pf)
      call pf_verlet_destroy(sweeper)
      call pfm_finalize_solver_level_params(level_params, pf_nml%nlevels)
    case (WM_VERLET)
      associate (dt => pepcboris_nml%dt, &
                 nt => pepcboris_nml%nt)

        do step=0,nt
          if(pepcboris_nml%rank==0) then
            write(*,*) " "
            write(*,'(a,i12,"/",i0)')   ' ====== computing step  :', step, nt
            write(*,'(a,f12.4)')        ' ====== simulation time :', step*dt
          end if

          do i=1,size(particles,kind=kind(i))
            write(47,*) step*dt, i, particles(i)%x, particles(i)%data%v
          end do

          call push_particles_velocity_verlet_boris(particles, dt)
          call eval_force(particles, level_params(pf_nml%nlevels), pepcboris_nml, step, MPI_COMM_SPACE, clearresults=.true.)
          call update_velocities_velocity_verlet_boris(particles, dt)

        end do
      end associate
    case (WM_ANALYTIC)
      associate (dt => pepcboris_nml%dt, &
                 nt => pepcboris_nml%nt, &
                 params => pepcboris_nml%setup_params)

        block
          complex*16 :: u, Rp, Rm
          real*8 :: Omegap, Omegam, Rscrm, Rscrp, Iscrm, Iscrp, Omegasq
          complex*16, parameter :: ic = (0._8,1._8)
          real*8, parameter :: sqrttwo = sqrt(2._8)

          Omegasq = sqrt((params(PARAMS_OMEGAB)**2)/4._8 - params(PARAMS_OMEGAE)**2)
          Omegap  = params(PARAMS_OMEGAB)/2._8 + Omegasq
          Omegam  = params(PARAMS_OMEGAB)/2._8 - Omegasq

          Rscrm = (params(PARAMS_VX0) + Omegap*params(PARAMS_X0)) / (Omegap - Omegam)
          Rscrp =  params(PARAMS_X0)  - Rscrm
          Iscrm = (params(PARAMS_VY0) + Omegap*params(PARAMS_Y0)) / (Omegap - Omegam)
          Iscrp =  params(PARAMS_Y0)  - Iscrm

          Rp = Rscrp + ic*Iscrp
          Rm = Rscrm + ic*Iscrm

          do step=0,nt
            if(pepcboris_nml%rank==0) then
              write(*,*) " "
              write(*,'(a,i12,"/",i0)')   ' ====== computing step  :', step, nt
              write(*,'(a,f12.4)')        ' ====== simulation time :', step*dt
            end if

            u = Rp * exp(-ic*Omegap*step*dt) + Rm * exp(-ic*Omegam*step*dt)

            particles(1)%x(1) = real(u)
            particles(1)%x(2) = aimag(u)
            particles(1)%x(3) = params(PARAMS_Z0) * cos(sqrttwo*params(PARAMS_OMEGAE)*step*dt) + &
              params(PARAMS_VZ0)/(sqrttwo*params(PARAMS_OMEGAE)) * sin(sqrttwo*params(PARAMS_OMEGAE)*step*dt)

            do i=1,size(particles,kind=kind(i))
              write(49,*) step*dt, i, particles(i)%x, particles(i)%data%v
            end do
          end do
        end block
      end associate
    case default
        DEBUG_ERROR(*,'Invalid working mode:', pepcboris_nml%workingmode)
  end select

  call pepc_finalize()
  call MPI_COMM_FREE(MPI_COMM_TIME, mpi_err)
  call MPI_COMM_FREE(MPI_COMM_SPACE, mpi_err)
  call MPI_FINALIZE( mpi_err )

end program pepc

