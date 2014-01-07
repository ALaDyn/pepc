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
!> diagnostics and output for pepc-andreev
!>
module pepca_diagnostics
  implicit none
  private
  save

    public write_particles_vtk
    public gather_and_write_densities
    public write_domain
    public diagnose_energy
    public dumpnow
    public compare_particles_to_checkpoint
    public compare_particles_to_particles
    public get_checkpoint_id

    integer, public, parameter :: E_KIN_E = 1
    integer, public, parameter :: E_POT_E = 2
    integer, public, parameter :: E_TOT_E = 3
    integer, public, parameter :: E_KIN_I = 4
    integer, public, parameter :: E_POT_I = 5
    integer, public, parameter :: E_TOT_I = 6
    integer, public, parameter :: E_KIN   = 7
    integer, public, parameter :: E_POT   = 8
    integer, public, parameter :: E_TOT   = 9
    integer, public, parameter :: E_MAXIDX = E_TOT

    integer, public, parameter :: NORM_MAXIMUM = 1
    integer, public, parameter :: NORM_L2      = 2
    integer, public, parameter :: norm         = NORM_L2

  contains

  logical function dumpnow(stepmod, step, nt)
    implicit none
    integer, intent(in) :: stepmod, step, nt

    dumpnow =        (stepmod>0) &
              .and. ((mod(step, stepmod)==0) .or. (step==nt))
  end function

  integer function get_checkpoint_id(physical_time)
    use pepca_units
    implicit none
    real*8, intent(in) :: physical_time

    get_checkpoint_id = aint(physical_time*unit_time_as_per_simunit)
  end function


  integer function vtk_step_of_step(step, nt) result(vtk_step)
    use module_vtk
    use pepca_helper
    implicit none

    integer, intent(in) :: step, nt

    if (step == 0) then
      vtk_step = VTK_STEP_FIRST
    else if (step == nt ) then
      vtk_step = VTK_STEP_LAST
    else
      vtk_step = VTK_STEP_NORMAL
    endif
  end function vtk_step_of_step


  subroutine write_particles_vtk(p, step, nt, realtime, comm)
    use module_vtk_helpers
    use module_pepc_types
    use pepca_units
    implicit none

    include 'mpif.h'

    type(t_particle), intent(in) :: p(:)
    real*8, intent(in) :: realtime
    integer, intent(in) :: step, nt
    integer(kind_default), intent(in) :: comm

    integer :: vtk_step

    vtk_step = vtk_step_of_step(step, nt)
    call vtk_write_particles("particles", comm, step, realtime, vtk_step, p, vtk_results, unit_length_micron_per_simunit)

    contains

    subroutine vtk_results(d, r, vtkf)
      use module_vtk
      use module_interaction_specific_types
      implicit none

      type(t_particle_data), intent(in) :: d(:)
      type(t_particle_results), intent(in) :: r(:)
      type(vtkfile_unstructured_grid), intent(inout) :: vtkf

      call vtk_write_particle_data_results(d, r, vtkf)
    end subroutine
  end subroutine write_particles_vtk


  subroutine gather_and_write_densities(p, Ngrid, step, nt, realtime, rank)
    use module_vtk_helpers
    use module_pepc_types
    use pepca_units
    use module_pepc, only: global_tree
    use pepca_helper, only: dim
    implicit none

    include 'mpif.h'

    type(t_particle), intent(in) :: p(:)
    real*8, intent(in) :: realtime
    integer, intent(in) :: step, nt
    integer, intent(in) :: Ngrid(1:3)
    integer(kind_pe), intent(in) :: rank

    type t_coordarray
      real*8, pointer :: coords(:)
    end type

    integer :: vtk_step, i, g
    integer :: cell(3)
    integer(kind_particle) :: ip
    integer, dimension(2,3) :: dims
    type(t_coordarray) :: grid(3)
    real*8, allocatable :: dens_el(:,:,:), dens_ion(:,:,:)
    integer(kind_default) :: ierr
    real*8 :: cell_volume

    vtk_step = vtk_step_of_step(step, nt)

    allocate (dens_el(Ngrid(1), Ngrid(2), Ngrid(3)), dens_ion(Ngrid(1), Ngrid(2), Ngrid(3)))
    dens_el    = 0.
    dens_ion   = 0.
    dims(:, :) = 0
    dims(2, 1:dim) = Ngrid(1:dim)

    do g=1,3
      allocate(grid(g)%coords(dims(1,g):dims(2,g)))
      do i=dims(1,g),dims(2,g)
        grid(g)%coords(i) = global_tree%bounding_box%boxmin(g) + global_tree%bounding_box%boxsize(g)/Ngrid(g)*i
      end do
    end do

    cell(dim+1:) = 1
    do ip = 1, size(p,kind=kind_particle)
      cell(1:dim) = nint((Ngrid(1:dim)-1)*(p(ip)%x(1:dim) - global_tree%bounding_box%boxmin(1:dim)) / global_tree%bounding_box%boxsize(1:dim)) + 1

      if (p(ip)%data%q < 0) then
        dens_el( cell(1), cell(2), cell(3)) = dens_el( cell(1), cell(2), cell(3)) + 1.
      else if (p(ip)%data%q > 0) then
        dens_ion(cell(1), cell(2), cell(3)) = dens_ion(cell(1), cell(2), cell(3)) + 1.
      endif
    end do


    if (rank==0) then
      ! we perform this output on a single rank to prevent having to think about the reduction above
      call MPI_REDUCE(MPI_IN_PLACE, dens_el,  size(dens_el,  kind=kind_default), MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      call MPI_REDUCE(MPI_IN_PLACE, dens_ion, size(dens_ion, kind=kind_default), MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

      ! compute real physical densities, i.e. particles per micron**2(3)
      cell_volume = product((global_tree%bounding_box%boxsize(1:dim) * unit_length_micron_per_simunit) / Ngrid(1:dim))
      dens_el(:,:,:)  = dens_el /cell_volume
      dens_ion(:,:,:) = dens_ion/cell_volume

      call vtk_write_densities_on_grid("densities", step, realtime, vtk_step, dims, dims, &
        grid(1)%coords, grid(2)%coords, grid(3)%coords, &
        dens_el,  'n_el' , &
        dens_ion, 'n_ion', &
        MPI_COMM_SELF, coord_scale=unit_length_micron_per_simunit)
    else
      call MPI_REDUCE(dens_el,  dens_el,  size(dens_el,  kind=kind_default), MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      call MPI_REDUCE(dens_ion, dens_ion, size(dens_ion, kind=kind_default), MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    endif

    do g=1,3
      deallocate(grid(g)%coords)
    end do
    deallocate (dens_el, dens_ion)

  end subroutine gather_and_write_densities


  subroutine write_domain(p, step, nt, realtime)
    use module_vtk
    use module_vtk_helpers
    use module_pepc_types
    use pepca_units
    use module_pepc, only: global_tree
    implicit none

    type(t_particle), allocatable, intent(in) :: p(:)
    real*8, intent(in) :: realtime
    integer, intent(in) :: step, nt

    integer :: vtk_step

    ! output of tree diagnostics
    vtk_step = vtk_step_of_step(step, nt)
    call vtk_write_branches(  step, realtime, vtk_step, global_tree, coord_scale=unit_length_micron_per_simunit)
    call vtk_write_leaves(    step, realtime, vtk_step, global_tree, coord_scale=unit_length_micron_per_simunit)
    call vtk_write_spacecurve(step, realtime, vtk_step, p,           coord_scale=unit_length_micron_per_simunit)
  end subroutine write_domain


  subroutine diagnose_energy(p, energies, step, realtime, comm, do_dump)
    use module_pepc_types
    use pepca_units
    use module_debug
    implicit none
    include 'mpif.h'

    type(t_particle), intent(in) :: p(:)
    real*8, intent(in) :: realtime
    integer, intent(in) :: step
    integer(kind_default), intent(in) :: comm
    logical, intent(in) :: do_dump
    integer(kind_default) :: ierr
    real*8 :: energies(E_MAXIDX)

    integer, parameter :: file_energies = 42

    integer(kind_particle) :: ip
    real*8 :: gam, ekin, epot

    call pepc_status('------------- diagnose_energy')

    energies = 0.

    do ip=1, size(p, kind=kind_particle)
      ! v represents momentum p/mc = gamma*v/c
      ! gamma = sqrt(1 + (p/mc)^2)
      gam  = sqrt( 1.0 + dot_product(p(ip)%data%v, p(ip)%data%v) )
      ekin = (gam-1._8) * p(ip)%data%m
      epot = p(ip)%data%q * p(ip)%results%pot / 2._8

      if (p(ip)%data%q > 0) then
        energies(E_KIN_I) = energies(E_KIN_I) + ekin
        energies(E_POT_I) = energies(E_POT_I) + epot
      else if (p(ip)%data%q < 0) then
          energies(E_KIN_E) = energies(E_KIN_E) + ekin
          energies(E_POT_E) = energies(E_POT_E) + epot
      else
        write(*,*) 'unexpected species in energy computation'
      endif
    end do

    energies(E_TOT_E) = energies(E_KIN_E) + energies(E_POT_E)
    energies(E_TOT_I) = energies(E_KIN_I) + energies(E_POT_I)
    energies(E_KIN  ) = energies(E_KIN_E) + energies(E_KIN_I)
    energies(E_POT  ) = energies(E_POT_E) + energies(E_POT_I)
    energies(E_TOT  ) = energies(E_TOT_E) + energies(E_TOT_I)

    ! we perform this output on a single rank to prevent having to think about the reduction above
    call MPI_ALLREDUCE(MPI_IN_PLACE, energies,  size(energies,  kind=kind_default), MPI_REAL8, MPI_SUM, comm, ierr)

    if (do_dump) then
      if (step == 0) then
        open(unit=file_energies, file='energy.dat', status='unknown', position='rewind',action='write')
        write(file_energies, '("#",a8,x,a10,  9(x,a16))') 'step', 'time (fs)', 'E_kin(e)', 'E_pot(e)', 'E_tot(e)', 'E_kin(i)', 'E_pot(i)', 'E_tot(i)', 'E_kin', 'E_pot', 'E_tot'
      else
        open(unit=file_energies, file='energy.dat', status='old', position='append',action='write')
      endif

      write(file_energies,   '( x ,I8,x,f10.4,9(x,g16.10))') step, realtime, energies

      close(file_energies)
    endif
  end subroutine diagnose_energy


  subroutine compare_particles_to_checkpoint(particles, itime_in, comm, xerr, verr)
    use module_pepc_types
    use module_checkpoint
    use module_debug
    implicit none
    type(t_particle), intent(in) :: particles(:)
    integer(kind_default), intent(in) :: itime_in
    integer(kind_default), intent(in) :: comm
    real*8, intent(out) :: xerr, verr

    type(t_particle), allocatable :: particles_ref(:)
    integer(kind_default) :: itime
    integer(kind_particle) :: n_total
    character(100) :: filename
    logical :: success

    call pepc_status('------------- compare_particles_to_checkpoint')

    call read_particles_mpiio(itime_in, comm, itime, n_total, particles_ref, filename, size(particles), file_exists=success, noparams=.true.)

    if (.not. success) then
      DEBUG_WARNING('("trying to read checkpoint with itime=", I0, " --> filename =", a," but the file does not exist. Please run the executable with use_pfasst=.false. in params to generate reference data first.")', itime_in, trim(filename))
      xerr = -1.
      verr = -1.
    else
      call compare_particles_to_particles(particles, particles_ref, comm, xerr, verr)
    endif

  end subroutine


  subroutine compare_particles_to_particles(particles, particles_ref, comm, xerr, verr)
    use module_debug
    use module_pepc_types
    implicit none
    include 'mpif.h'
    type(t_particle), intent(in) :: particles(:), particles_ref(:)
    integer(kind_default), intent(in) :: comm
    real*8, intent(out) :: xerr, verr

    integer(kind_particle) :: i
    integer(kind_default) :: ierr
    real*8 :: xsum2, vsum2, mydiff(3)

    call pepc_status('------------- compare_particles_to_particles')

    DEBUG_ASSERT(size(particles) == size(particles_ref))

    xerr  = 0
    verr  = 0
    xsum2 = 0
    vsum2 = 0

    select case (norm)

      case (NORM_MAXIMUM)
        do i=1,size(particles, kind=kind_particle)
          xerr = max(xerr, maxval(abs((particles(i)%x     -particles_ref(i)%x     )/particles_ref(i)%x     )))
          verr = max(verr, maxval(abs((particles(i)%data%v-particles_ref(i)%data%v)/particles_ref(i)%data%v)))
        end do

        call MPI_ALLREDUCE(MPI_IN_PLACE, xerr, 1, MPI_REAL8, MPI_MAX, comm, ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE, verr, 1, MPI_REAL8, MPI_MAX, comm, ierr)

      case (NORM_L2)
        do i=1,size(particles, kind=kind_particle)
          xsum2 = xsum2 + dot_product(particles_ref(i)%x,     particles_ref(i)%x)
          vsum2 = vsum2 + dot_product(particles_ref(i)%data%v,particles_ref(i)%data%v)

          mydiff = particles(i)%x-particles_ref(i)%x
          xerr = xerr + dot_product(mydiff, mydiff)

          mydiff = particles(i)%data%v-particles_ref(i)%data%v
          verr = verr + dot_product(mydiff, mydiff)
        end do

        call MPI_ALLREDUCE(MPI_IN_PLACE, xerr,  1, MPI_REAL8, MPI_SUM, comm, ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE, verr,  1, MPI_REAL8, MPI_SUM, comm, ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE, xsum2, 1, MPI_REAL8, MPI_SUM, comm, ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE, vsum2, 1, MPI_REAL8, MPI_SUM, comm, ierr)

        xerr = sqrt(xerr) / sqrt(xsum2)
        verr = sqrt(verr) / sqrt(vsum2)

      case default
        DEBUG_ERROR(*, 'Invalid norm selected')
    end select

  end subroutine

end module














