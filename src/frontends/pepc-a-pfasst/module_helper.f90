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
  integer, parameter :: t_user_diag        = t_userdefined_first + 3

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

    namelist /pepcandreev/ nt, dt, particle_output_interval, domain_output_interval, eps, Ngrid
    
    ! set default parameter values
    nt                        = 20000
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
      if(root_space) write(*,'(a)') ' == reading parameter file, section pepcandreev: ', para_file
      open(fid,file=para_file)
      read(fid,NML=pepcandreev)
      close(fid)
    else
      if(root_space) write(*,*) ' == no param file, using default parameters '
    end if    

    if(root_space) then
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
    integer, intent(in) :: nel, nion
    character(*), intent(in) :: file_el, file_ion
    
    integer, parameter :: filehandle = 47
    integer(kind_particle) :: ip

    if(root_space) then
      write(*,'(a, 2(x,a))') " == [load] reading particles from files", file_el, file_ion
    
      ! we read all particles the root rank

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
    endif

  end subroutine read_particles

  
  subroutine generate_particles(p, nel, nion)
    use pepca_globals
    implicit none
    
    type(t_particle), allocatable, intent(inout) :: p(:)
    integer, intent(in) :: nel, nion
    
    real*8 :: pos(1:3), vel(1:3)
    integer :: i,j,k,l

    ! stupid parallel random number generation
    if(root_space) write(*,'(a, 2(x,a))') " == [generate] generating particles"
      
    pos = 0
    vel = 0
    l   = 0
    do i=0,nrank_space-1
      do j=1,nel+nion
        l = l+1
        do k=1,dim
          pos(k) = par_rand()
          vel(k) = par_rand()
        end do

        if (i==rank_space) then
          p(j)%x(1:dim)      = pos(1:dim)
          p(j)%data%v(1:dim) = vel(1:dim)
          if (j<=nel) then
            p(j)%label       = -l
            p(j)%data%q      =  unit_qe
            p(j)%data%m      =  unit_me
          else
            p(j)%label       = -l
            p(j)%data%q      =  unit_qe
            p(j)%data%m      =  unit_me
          endif

          p(j)%work          =  1.0
        endif
      end do      
    end do

  end subroutine generate_particles
  
  !>
  !> portable random number generator, see numerical recipes
  !> check for the random numbers:
  !> the first numbers should be 0.2853809, 0.2533582 and 0.0934685
  !> the parameter iseed is optional
  !>
  function par_rand(iseed)
    implicit none
    real :: par_rand
    integer, intent(in), optional :: iseed
    
    integer, parameter :: IM1  = 2147483563
    integer, parameter :: IM2  = 2147483399
    real,    parameter :: AM   = 1.0/IM1
    integer, parameter :: IMM1 = IM1-1
    integer, parameter :: IA1  = 40014
    integer, parameter :: IA2  = 40692
    integer, parameter :: IQ1  = 53668
    integer, parameter :: IQ2  = 52774
    integer, parameter :: IR1  = 12211
    integer, parameter :: IR2  = 3791
    integer, parameter :: NTAB = 32
    integer, parameter :: NDIV = 1+IMM1/NTAB
    real,    parameter :: eps_ = 1.2e-7 ! epsilon(eps_)
    real,    parameter :: RNMX = 1.0 - eps_
    
    integer :: j, k
    integer, volatile, save :: idum  = -1
    integer, volatile, save :: idum2 =  123456789
    integer, volatile, save :: iy    =  0
    integer, volatile, save :: iv(NTAB)
    
    
    if (idum <=0 .or. present(iseed)) then
       if (present(iseed)) then
          idum = iseed
       else
          if (-idum < 1) then
             idum = 1
          else
             idum = -idum
          endif
       endif
       
       idum2 = idum
       
       do j = NTAB+7,0,-1
          k = idum/IQ1
          idum = IA1 * (idum-k*IQ1) - k*IR1
          if (idum < 0 ) idum = idum + IM1
          
          if (j<NTAB) iv(j+1) = idum
          
       end do
       iy = iv(1)
    end if
    
    k = idum/IQ1
    idum = IA1 * (idum-k*IQ1) - k*IR1
    if (idum < 0) idum = idum + IM1
    
    k = idum2/IQ2
    idum2 = IA2 * (idum2-k*IQ2) - k*IR2
    if (idum2 < 0) idum2 = idum2 + IM2
    
    j = iy/NDIV + 1
    iy = iv(j)-idum2
    iv(j) = idum
    
    if (iy < 1) iy = iy + IMM1
    par_rand = AM*iy
    if (par_rand > RNMX) par_rand = RNMX
  end function par_rand
  
end module
