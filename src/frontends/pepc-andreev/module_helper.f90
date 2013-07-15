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
module pepca_helper
  use module_pepc_types
  use module_timings
  use pepca_units
  use pepca_globals
  implicit none

  ! timing variables
  integer, parameter :: t_user_total       = t_userdefined_first
  integer, parameter :: t_user_init        = t_userdefined_first + 1
  integer, parameter :: t_user_step        = t_userdefined_first + 2
  
  ! MPI variables
  integer(kind_pe) :: my_rank, n_ranks

  ! time variables
  real*8 :: dt
  ! interaction cutoff parameter
  real*8 :: eps

  ! control variables
  integer :: nt                   ! number of timesteps
  integer :: particle_output_interval      ! turn vtk output on/off
  integer :: domain_output_interval        ! turn vtk output on/off
  

  contains

  subroutine set_parameter()
      
    use module_pepc
    use module_interaction_specific, only : theta2, eps2
    use treevars, only : num_threads, np_mult
    use pepca_globals
    implicit none
      
    integer, parameter :: fid = 12
    character(255)     :: para_file
    logical            :: read_para_file

    namelist /pepcandreev/ nt, dt, particle_output_interval, domain_output_interval, eps
    
    ! set default parameter values
    nt                        = 5000
    dt                        = 0.2 ! here, dt is still in fs
    particle_output_interval  = 25
    domain_output_interval    =  0
    eps                       = 1.e-5 ! in microns
    theta2                    = 0.36
 
    num_threads = 8
    np_mult = -500

    ! read in namelist file
    call pepc_read_parameters_from_first_argument(read_para_file, para_file)
    
    dt   =  dt / unit_time_fs_per_simunit ! now, dt is in sim units
    eps2 = (eps/unit_length_micron_per_simunit)**2

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
      write(*,'(a,es12.4)')    " == time step (fs)            : ", dt*unit_time_fs_per_simunit
      write(*,'(a,es12.4)')    " == final time (ns)           : ", dt*nt*unit_time_fs_per_simunit
      write(*,'(a,l12)')       " == particle output interval  : ", particle_output_interval
      write(*,'(a,l12)')       " == domain output interval    : ", domain_output_interval
    end if

    call pepc_prepare(dim)
  end subroutine set_parameter


  subroutine read_particles(p, file_el, nel, file_ion, nion)
    use pepca_globals
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
        read(filehandle, *)   p(ip)%x(1:dim), p(ip)%data%v(1:dim)
        p(ip)%x(dim+1:)      = 0.
        p(ip)%data%v(dim+1:) = 0.
        ! other stuff
        p(ip)%label       = -ip
        p(ip)%data%q      =  unit_qe
        p(ip)%data%m      =  unit_me
        p(ip)%work        =  1.0
        
        p(ip)%x = p(ip)%x / unit_length_micron_per_simunit
      end do  
      close(filehandle)

      !call system("sed -i 's/\([0-9]\)-/\1 -/g' " // trim(file_ion))
      open(filehandle, file=trim(file_ion), STATUS='OLD', POSITION = 'REWIND', ACTION='READ')
      do ip=nel+1,nel+nion
        read(filehandle, *)   p(ip)%x(1:dim), p(ip)%data%v(1:dim)
        p(ip)%x(dim+1:)        = 0.
        p(ip)%data%v(dim+1:)   = 0.
        ! other stuff
        p(ip)%label       =  ip-nel 
        p(ip)%data%q      =  unit_qp
        p(ip)%data%m      =  unit_mp
        p(ip)%work        =  1.0
        
        p(ip)%x = p(ip)%x / unit_length_micron_per_simunit
      end do  
      close(filehandle)
    else
      allocate(p(0))
    endif

  end subroutine read_particles
  
  
end module
