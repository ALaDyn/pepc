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

!>
!> helper module
!>
module pepcboris_diagnostics
  implicit none
  private
  save

  public dump_particles
  public dump_energy
  public backup_velocities

  real*8, allocatable, public :: vold(:,:)

  contains

   subroutine backup_velocities(particles)
    use module_pepc_types
    implicit none
    type(t_particle), intent(in) :: particles(:)
    integer(kind_particle) :: p

    if (.not. allocated(vold)) allocate(vold(3,size(particles)))

    do p = 1, size(particles, kind=kind_particle)
      vold(1:3,p) = particles(p)%data%v(1:3)
    end do
  end subroutine

  subroutine dump_particles(vtk_step, step, dt, particles, comm, do_average)
    use module_pepc_types
    use module_pepc
    use module_debug
    use module_vtk_helpers
    use module_vtk
    use pepcboris_helper, only: IFILE_SUMMAND_PARTICLES, IFILE_SUMMAND_PARTICLES_AVG, pepcboris_nml
    implicit none
    integer, intent(in) :: vtk_step
    integer, intent(in) :: step
    real*8, intent(in) :: dt
    type(t_particle), intent(in) :: particles(:)
    integer, intent(in) :: comm
    logical, intent(in) :: do_average
    integer(kind_particle) :: p
    integer :: istream, dumptype
    type(t_particle) :: pavg

    istream  = pepcboris_nml%workingmode + IFILE_SUMMAND_PARTICLES
    dumptype = pepcboris_nml%dumptype

    select case (dumptype)
      case (0)
        ! linear output to fort.istream file
        if (do_average) then
          DEBUG_ASSERT(allocated(vold))
          ! we have to average over old and new velocities
          do p=1,size(particles,kind=kind(p))
            write(istream,*) step*dt, p, particles(p)%x, (particles(p)%data%v + vold(:,p))/2._8
          end do
        else
          do p=1,size(particles,kind=kind(p))
            write(istream,*) step*dt, p, particles(p)%x, particles(p)%data%v
          end do
        endif
      case (1,3)
        ! vtk output
        call vtk_write_particles("particles", comm, step, step*dt, vtk_step, particles, particle_output_data)

        if (dumptype==3) then
          pavg%x      = 0.
          pavg%data%v = 0.

          if (do_average) then
            DEBUG_ASSERT(allocated(vold))
            ! we have to average over old and new velocities
            do p=1,size(particles,kind=kind(p))
              pavg%x      = pavg%x      +  particles(p)%x
              pavg%data%v = pavg%data%v + (particles(p)%data%v + vold(:,p))/2._8
            end do
          else
            do p=1,size(particles,kind=kind(p))
              pavg%x      = pavg%x      +  particles(p)%x
              pavg%data%v = pavg%data%v +  particles(p)%data%v
            end do
          endif

          pavg%x      = pavg%x      / size(particles,kind=kind(p))
          pavg%data%v = pavg%data%v / size(particles,kind=kind(p))

          write(pepcboris_nml%workingmode + IFILE_SUMMAND_PARTICLES_AVG,*) step*dt, 1, pavg%x, pavg%data%v
        endif


      case default
        DEBUG_ERROR(*, 'dump_particles() - invalid value for dumptype:', dumptype)
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
        call vtkf%write_data_array("v_avg", ( d(:)%v(1)+vold(1,:) )/2._8, &
                                            ( d(:)%v(2)+vold(2,:) )/2._8, &
                                            ( d(:)%v(3)+vold(3,:) )/2._8)
      else
        call vtkf%write_data_array("v_avg", d(:)%v(1), d(:)%v(2), d(:)%v(3))
      endif
    end subroutine
  end subroutine

  subroutine dump_energy(t, particles, level_params, comm, do_average)
    use module_pepc_types
    use pfm_encap, only : level_params_t
    use pfm_feval, only : eval_force
    use pepcboris_helper, only : IFILE_SUMMAND_ENERGY, pepcboris_nml
    use module_debug
    implicit none
    real*8, intent(in) :: t
    type(t_particle), intent(in) :: particles(:)
    type(level_params_t), intent(in) :: level_params
    integer(kind_default), intent(in) :: comm
   logical, intent(in) :: do_average

    integer(kind_particle) :: p
    type(t_particle), allocatable :: particles_tmp(:)
    real*8 :: epot, ekin, etot, vtmp(3)
    integer :: istream

    istream  = pepcboris_nml%workingmode + IFILE_SUMMAND_ENERGY

    epot = 0.
    ekin = 0.

    ! we need a temporary copy here to prevent destruction the results stored in the particles-array
    allocate(particles_tmp(size(particles)))
    particles_tmp(:) = particles(:)
    call eval_force(particles_tmp, level_params, pepcboris_nml, -1, comm, clearresults=.true.)
    if (do_average) then
      DEBUG_ASSERT(allocated(vold))
      ! we have to average over old and new velocities
      do p=1,size(particles_tmp,kind=kind(p))
        epot = epot + particles_tmp(p)%data%q * particles_tmp(p)%results%pot
        vtmp = (particles_tmp(p)%data%v + vold(:,p) ) / 2._8
        ekin = ekin + particles_tmp(p)%data%m/2._8 * dot_product(vtmp, vtmp)
      end do
    else
      do p=1,size(particles_tmp,kind=kind(p))
        epot = epot + particles_tmp(p)%data%q * particles_tmp(p)%results%pot
        ekin = ekin + particles_tmp(p)%data%m/2._8 * dot_product(particles_tmp(p)%data%v,particles_tmp(p)%data%v)
      end do
    endif


    deallocate(particles_tmp)

    etot = epot + ekin
    write(istream,*) t, epot, ekin, etot

  end subroutine


end module
