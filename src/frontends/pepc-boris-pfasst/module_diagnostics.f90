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

  subroutine dump_particles(t, particles, istream, do_average)
    use module_pepc_types
    use module_debug
    implicit none
    real*8, intent(in) :: t
    type(t_particle), intent(in) :: particles(:)
    integer, intent(in) :: istream
    logical, intent(in) :: do_average
    integer(kind_particle) :: p

    if (do_average) then
      DEBUG_ASSERT(allocated(vold))
      ! we have to average over old and new velocities
      do p=1,size(particles,kind=kind(p))
        write(istream,*) t, p, particles(p)%x, (particles(p)%data%v + vold(:,p))/2._8
      end do
    else
      do p=1,size(particles,kind=kind(p))
        write(istream,*) t, p, particles(p)%x, particles(p)%data%v
      end do
    endif
  end subroutine

  subroutine dump_energy(t, particles, istream, level_params, nml, comm, do_average)
    use module_pepc_types
    use pfm_encap, only : level_params_t
    use pfm_feval, only : eval_force
    use pepcboris_helper, only : pepcboris_nml_t
    use module_debug
    implicit none
    real*8, intent(in) :: t
    type(t_particle), intent(in) :: particles(:)
    integer, intent(in) :: istream
    type(level_params_t), intent(in) :: level_params
    type(pepcboris_nml_t), intent(in) :: nml
    integer(kind_default), intent(in) :: comm
   logical, intent(in) :: do_average

    integer(kind_particle) :: p
    type(t_particle), allocatable :: particles_tmp(:)
    real*8 :: epot, ekin, etot, vtmp(3)

    epot = 0.
    ekin = 0.

    ! we need a temporary copy here to prevent destruction the results stored in the particles-array
    allocate(particles_tmp(size(particles)))
    particles_tmp(:) = particles(:)
    call eval_force(particles_tmp, level_params, nml, -1, comm, clearresults=.true.)
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
