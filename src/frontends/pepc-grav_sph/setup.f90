! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2016 Juelich Supercomputing Centre, 
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

!  ================================
!
!         SETUP
!
!   $Revision$
!
!     Initialise constants and 
!      simulation variables
!
!  ================================


subroutine pepc_setup()

  use physvars, only: &
       particles, &
       theta, &
       nt, &
       np_local, &
       npart_total, &
       trun, &
       n_cpu, &
       my_rank, &
       mac, &
       ispecial, &
       idim, &
       eps, &
       dt, &
       dump_time, &
       cp_time, &
       do_sph, &
       do_gravity

  use module_mirror_boxes, only: &
       t_lattice_1, &
       t_lattice_2, &
       t_lattice_3, &
       periodicity

  use module_pepc, only: &
       pepc_get_para_file

  use module_interaction_specific_types, only: &
       num_neighbour_particles

  use files, only: &
       parameter_file_name, &
       parameter_file_available

  implicit none
  include 'mpif.h'

  integer :: ierr, npart_tmp

  namelist /pepcsph/ npart_total, &
       mac, theta, eps, ispecial, &
       idim, nt, dt, &
       t_lattice_1, t_lattice_2, t_lattice_3, periodicity, &
       dump_time, cp_time, num_neighbour_particles, &
       do_gravity, do_sph


  !  Default input set
  ispecial        =   6

  ! particle number
  npart_total  = 10000 ! default particle number

  ! physics stuff
  mac         = 0
  theta       = 0.6

  eps           = 0.01

  ! control
  nt           = 1
  dt           = 0.01
  trun         = 0.

  
  if (parameter_file_available) then
     if(my_rank .eq. 0) write(*,*) "reading parameter file, section pepcsph: ", parameter_file_name
     open(10,file=parameter_file_name)
     read(10,NML=pepcsph)
     close(10)
  else
     if(my_rank .eq. 0) write(*,*) "##### using default parameter #####"
  end if


  write(*,*) 'periodicity:', periodicity
  write(*,*) 't_lattice_1:', t_lattice_1


  ! Derived parameters

  np_local = npart_total/n_cpu
  if (np_local*n_cpu .ne. npart_total .and. mod(npart_total,n_cpu) > my_rank)  np_local = npart_total/n_cpu+1


  allocate ( particles(np_local) )


  if (my_rank == 0) then
     write(*,*) "Starting PEPC-SPH with",n_cpu," Processors, simulating",np_local, &
                         " Particles on each Processor in",nt,"timesteps..."
  end if

end subroutine pepc_setup




