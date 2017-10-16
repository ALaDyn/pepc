! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
!
! Copyright (C) 2002-2017 Juelich Supercomputing Centre,
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

!>
!> helper module
!>
module pepcboris_diagnostics
  use module_pepc_kinds
  implicit none
  private
  save

  public dump_particles
  public dump_iterations
  public dump_energy
  public backup_velocities

  real(kind_physics), allocatable, public :: vold(:,:)

  contains
  
  logical function dontdump(step)
    use pepcboris_helper
    implicit none
    integer, intent(in) :: step
    dontdump = .not. ((step == 1) .or. (step == pepcboris_nml%nt) .or. (mod(step, pepcboris_nml%dumpstep) == 0) )
  end function

  subroutine backup_velocities(particles)
    use pf_mod_dtype
    use module_pepc_kinds
    use module_pepc_types
    implicit none
    type(t_particle), intent(in) :: particles(:)
    integer(kind_particle) :: p

    if (.not. allocated(vold)) allocate(vold(3,size(particles)))

    do p = 1, size(particles, kind=kind_particle)
      vold(1:3,p) = particles(p)%data%v(1:3)
    end do
  end subroutine

  subroutine dump_iterations(step, dt, hook, niter, residual)
    use pf_mod_dtype
    use pepcboris_helper
    use pepcboris_paralleldump
    implicit none
    integer, intent(in) :: step, niter, hook
    real(pfdp), intent(in) :: residual
    real(pfdp), intent(in) :: dt

    character(len=PARALLELDUMP_MAXLEN) :: line

    if (dontdump(step)) return

    write(line,*) step*dt, step, hook, niter, residual
    call paralleldump_dump(pepcboris_nml%workingmode + IFILE_SUMMAND_NITER, line)
  end subroutine

  subroutine dump_particles(vtk_step, step, dt, particles, comm, do_average)
    use pf_mod_dtype
    use module_pepc_kinds
    use module_pepc_types
    use module_pepc
    use module_debug
    use module_vtk_helpers
    use module_vtk
    use pepcboris_paralleldump
    use pepcboris_helper, only: IFILE_SUMMAND_PARTICLES, IFILE_SUMMAND_PARTICLES_AVG, pepcboris_nml
    implicit none
    integer, intent(in) :: vtk_step
    integer, intent(in) :: step
    real(pfdp), intent(in) :: dt
    type(t_particle), intent(in) :: particles(:)
    integer, intent(in) :: comm
    logical, intent(in) :: do_average
    integer(kind_particle) :: p
    integer :: istream, dumptype
    real(kind_physics) :: avg(3,4), delt(3)

    character(len=PARALLELDUMP_MAXLEN) :: line

    if (dontdump(step)) return

    istream  = pepcboris_nml%workingmode + IFILE_SUMMAND_PARTICLES
    dumptype = pepcboris_nml%dumptype

    select case (dumptype)
      case (1)
        ! linear output of all particles to fort.istream file
        if (do_average) then
          DEBUG_ASSERT(allocated(vold))
          ! we have to average over old and new velocities
          do p=1,size(particles,kind=kind(p))
            write(line,*) step*dt, p, particles(p)%x, (particles(p)%data%v + vold(:,p))/2._kind_physics
            call paralleldump_dump(istream, line)
          end do
        else
          do p=1,size(particles,kind=kind(p))
            write(line,*) step*dt, p, particles(p)%x, particles(p)%data%v
            call paralleldump_dump(istream, line)
          end do
        endif
    end select

    select case (dumptype)
      case (2,4)
        ! vtk output
        ! FIXME: I am rather sure that this will not work with time-parallel execution
        call vtk_write_particles("particles", comm, step, real(step*dt, kind(1._8)), vtk_step, particles, particle_output_data)
    end select

    select case (dumptype)
      case (3,4)
        avg = 0._kind_physics
        ! average position and velocity
        do p=1,size(particles,kind=kind(p))
          avg(:,1) = avg(:,1) +  particles(p)%x
          avg(:,2) = avg(:,2) +  particles(p)%data%v
        end do
        ! compute average absolute distance vector from center of mass
        do p=1,size(particles,kind=kind(p))
          delt(:)  = particles(p)%x - avg(:,1)
          avg(:,3) = avg(:,3) +  abs( delt )

          avg(1,4) =     avg(1,4) + sqrt(dot_product(delt, delt))
          avg(2,4) = max(avg(2,4),  sqrt(dot_product(delt, delt)))
        end do

        if (do_average) then
          DEBUG_ASSERT(allocated(vold))
          ! we have to average over old and new velocities
          avg(:,2) = ( avg(:,2) + sum(vold,dim=2) ) / 2._kind_physics
        endif

        avg = avg / size(particles,kind=kind(p))

        write(line,*) step*dt, 1, avg(:,1), avg(:,2), avg(:,3), avg(1,4), avg(2,4), particles(1)%x
        call paralleldump_dump(pepcboris_nml%workingmode + IFILE_SUMMAND_PARTICLES_AVG, line)
    end select

    contains

    subroutine particle_output_data(d, r, vtkf)
      use module_vtk
      use module_interaction_specific_types
      implicit none

      type(t_particle_data), intent(in) :: d(:)
      type(t_particle_results), intent(in) :: r(:)
      type(vtkfile_unstructured_grid), intent(inout) :: vtkf

      call vtk_write_particle_data_results(d, r, vtkf)
      if (do_average) then
        ! we have to average over old and new velocities
        DEBUG_ASSERT(allocated(vold))
        call vtkf%write_data_array("v_avg", ( d(:)%v(1)+vold(1,:) )/2._kind_physics, &
                                            ( d(:)%v(2)+vold(2,:) )/2._kind_physics, &
                                            ( d(:)%v(3)+vold(3,:) )/2._kind_physics)
      else
        call vtkf%write_data_array("v_avg", d(:)%v(1), d(:)%v(2), d(:)%v(3))
      endif
    end subroutine
  end subroutine

  subroutine dump_energy(step, t, particles, level_params, comm, do_average)
    use pf_mod_dtype
    use module_pepc_kinds
    use module_pepc_types
    use pepcboris_paralleldump
    use pfm_encap, only : level_params_t
    use pfm_feval, only : eval_force
    use pepcboris_helper, only : IFILE_SUMMAND_ENERGY, pepcboris_nml
    use module_debug
    implicit none
    integer, intent(in) :: step
    real(pfdp), intent(in) :: t
    type(t_particle), intent(in) :: particles(:)
    type(level_params_t), intent(in) :: level_params
    integer(kind_default), intent(in) :: comm
   logical, intent(in) :: do_average

    integer(kind_particle) :: p
    type(t_particle), allocatable :: particles_tmp(:)
    real(kind_physics) :: epot, ekin, etot, vtmp(3)
    character(len=PARALLELDUMP_MAXLEN) :: line

    if (dontdump(step)) return

    epot = 0._kind_physics
    ekin = 0._kind_physics
    ! we need a temporary copy here to prevent destruction the results stored in the particles-array
    allocate(particles_tmp(size(particles)))
    particles_tmp(:) = particles(:)
    call eval_force(particles_tmp, level_params, pepcboris_nml, comm, clearresults=.true.)
    if (do_average) then
      DEBUG_ASSERT(allocated(vold))
      ! we have to average over old and new velocities
      do p=1,size(particles_tmp,kind=kind(p))
        epot = epot + particles_tmp(p)%data%q * particles_tmp(p)%results%pot
        vtmp = (particles_tmp(p)%data%v + vold(:,p) ) / 2._kind_physics
        ekin = ekin + particles_tmp(p)%data%m/2._kind_physics * dot_product(vtmp, vtmp)
      end do
    else
      do p=1,size(particles_tmp,kind=kind(p))
        epot = epot + particles_tmp(p)%data%q * particles_tmp(p)%results%pot
        ekin = ekin + particles_tmp(p)%data%m/2._kind_physics * dot_product(particles_tmp(p)%data%v,particles_tmp(p)%data%v)
      end do
    endif


    deallocate(particles_tmp)

    etot = epot + ekin
    write(line,*) t, epot, ekin, etot
    call paralleldump_dump(pepcboris_nml%workingmode + IFILE_SUMMAND_ENERGY, line)

  end subroutine


end module
