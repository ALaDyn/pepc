! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2012 Juelich Supercomputing Centre, 
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
module helper
  use module_pepc_types
  use module_timings
  use module_units
  implicit none

  ! timing variables
  integer, parameter :: t_user_total       = t_userdefined_first
  integer, parameter :: t_user_init        = t_userdefined_first + 1
  integer, parameter :: t_user_step        = t_userdefined_first + 2
  integer, parameter :: t_user_directsum   = t_userdefined_first + 3
  integer, parameter :: t_user_particleio  = t_userdefined_first + 4
  
  ! MPI variables
  integer(kind_pe) :: my_rank, n_ranks
  logical :: root

  ! time variables
  real*8 :: dt
  integer :: step

  ! control variables
  integer :: nt                   ! number of timesteps
  logical :: particle_output      ! turn vtk output on/off
  logical :: domain_output        ! turn vtk output on/off
  integer :: diag_interval        ! number of timesteps between all diagnostics and IO
  
  integer, parameter :: particle_direct = -1 ! number of particle for direct summation

  contains

  subroutine set_parameter()
      
    use module_pepc
    use module_interaction_specific, only : theta2, eps2
    implicit none
      
    integer, parameter :: fid = 12
    character(255)     :: para_file
    logical            :: read_para_file

    namelist /pepcandreev/ nt, dt, particle_output, domain_output, diag_interval
    
    ! set default parameter values
    nt                = 500
    dt                = 0.5e-3
    particle_output   = .true.
    domain_output     = .true.
    diag_interval     = 10    
 
    ! read in namelist file
    call pepc_read_parameters_from_first_argument(read_para_file, para_file)

    if (read_para_file) then
      if(root) write(*,'(a)') " == reading parameter file, section pepcandreev: ", para_file
      open(fid,file=para_file)
      read(fid,NML=pepcandreev)
      close(fid)
    else
      if(root) write(*,*) " == no param file, using default parameter "
    end if    

    if(root) then
      write(*,'(a,i12)')       " == number of time steps      : ", nt
      write(*,'(a,es12.4)')    " == time step                 : ", dt
      write(*,'(a,i12)')       " == diag & IO interval        : ", diag_interval
      write(*,'(a,l12)')       " == particle output           : ", particle_output
      write(*,'(a,l12)')       " == domain output             : ", domain_output
    end if

    call pepc_prepare(2_kind_dim)
  end subroutine set_parameter


  subroutine read_particles(p, file_el, nel, file_ion, nion)
    implicit none
    
    type(t_particle), allocatable, intent(inout) :: p(:)
    integer(kind_particle), intent(in) :: nel, nion
    character(*), intent(in) :: file_el, file_ion
    
    integer, parameter :: filehandle = 47
    integer(kind_particle) :: ip

    if(root) then
      write(*,'(a, 2(x,a))') " == [load] reading particles from files", file_el, file_ion
    
      ! we read all particles the root rank
      allocate(p(nel + nion))

      ! Script for file preprocessing: sed -i 's/\([0-9]\)-/\1 -/g' filename
      !call system("sed -i 's/\([0-9]\)-/\1 -/g' " // trim(file_el))
      open(filehandle, file=trim(file_el), STATUS='OLD', POSITION = 'REWIND', ACTION='READ')
      do ip=1,nel
        read(filehandle, *)   p(ip)%x(1:2), p(ip)%data%v(1:2)
        p(ip)%x(3) = 0.
        p(ip)%v(3) = 0.
        ! other stuff
        p(ip)%label       = -ip
        p(ip)%data%q      =  unit_qe
        p(ip)%data%m      =  unit_me
        p(ip)%work        =  1.0
        ! rescale units
        p(ip)%x = p(ip)%x * 1.e6 / unit_abohr_in_m ! micron --> aB
        p(ip)%v = p(ip)%v / p(ip)%data%m / unit_c ! m*c -> simunits
      end do  
      close(filehandle)

      !call system("sed -i 's/\([0-9]\)-/\1 -/g' " // trim(file_ion))
      open(filehandle, file=trim(file_ion), STATUS='OLD', POSITION = 'REWIND', ACTION='READ')
      do ip=nel+1,nel+nion
        read(filehandle, *)   p(ip)%x(1:2), p(ip)%data%v(1:2)
        p(ip)%x(3) = 0.
        p(ip)%v(3) = 0.
        ! other stuff
        p(ip)%label       =  ip-nel 
        p(ip)%data%q      =  unit_qp
        p(ip)%data%m      =  unit_mp
        p(ip)%work        =  1.0
        ! rescale units
        p(ip)%data%x = p(ip)%data%x * 1.e6 / unit_abohr_in_m ! micron --> aB
        p(ip)%data%v = p(ip)%data%v / p(ip)%data%m / unit_c ! m*c -> simunits
      end do  
      close(filehandle)
    else
      allocate(p(0))
    endif

  end subroutine read_particles

    
  subroutine push_particles(p)
    use module_mirror_boxes
    implicit none
    
    type(t_particle), allocatable, intent(inout) :: p(:)
    integer(kind_particle) :: ip
    real*8  :: acc(1:3), gam

    if(root) write(*,'(a)') " == [pusher] push particles "

    do ip=1, size(p, kind=kind_particle)
      acc(1:2) = p(ip)%data%q * p(ip)%results%e(1:2) / p(ip)%data%m
      acc(3)   = 0.
      
      p(ip)%data%v = p(ip)%data%v + dt * acc
      gam          = sqrt( 1.0 + ( dot_product(p(ip)%data%v, p(ip)%data%v) ) / unit_c2 )
      p(ip)%x      = p(ip)%x + p(ip)%data%v / gam * dt
    end do
  end subroutine push_particles

  
  integer function vtk_step_of_step(step) result(vtk_step)
    use module_vtk
    implicit none

    integer, intent(in) :: step

    if (step .eq. 0) then
      vtk_step = VTK_STEP_FIRST
    else if (nt - 1 - step < diag_interval ) then
      vtk_step = VTK_STEP_LAST
    else
      vtk_step = VTK_STEP_NORMAL
    endif
  end function vtk_step_of_step


  subroutine write_particles(p)
    use module_vtk_helpers
    implicit none

    include 'mpif.h'
    
    type(t_particle), intent(in) :: p(:)

    integer :: vtk_step
    
    call timer_start(t_user_particleio)
    vtk_step = vtk_step_of_step(step)
    call vtk_write_particles("particles", MPI_COMM_WORLD, step, dt * step, vtk_step, p, vtk_coulomb_results)
    call timer_stop(t_user_particleio)
    if(root) write(*,'(a,es12.4)') " == [write particles] time in vtk output [s]      : ", timer_read(t_user_particleio)

    contains

    subroutine vtk_coulomb_results(d, r, vtkf)
      use module_vtk
      use module_interaction_specific_types
      implicit none

      type(t_particle_data), intent(in) :: d(:)
      type(t_particle_results), intent(in) :: r(:)
      type(vtkfile_unstructured_grid), intent(inout) :: vtkf
      
      call vtk_write_particle_data_results(d, r, vtkf)
    end subroutine
  end subroutine write_particles

  
  subroutine write_domain(p)
    use module_vtk
    use module_vtk_helpers
    implicit none
  
    type(t_particle), allocatable, intent(in) :: p(:)

    integer :: vtk_step
  
    ! output of tree diagnostics
    vtk_step = vtk_step_of_step(step)
    call vtk_write_branches(step,  dt * step, vtk_step)
    call vtk_write_leaves(step, dt * step, vtk_step)
    call vtk_write_spacecurve(step, dt * step, vtk_step, p)
  end subroutine write_domain
  
end module
