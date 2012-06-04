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


program pepc

  ! pepc modules
  use module_pepc
  use module_pepc_types
  use module_mirror_boxes
  
  ! frontend helper routines
  use helper
  use variables
  use particlehandling
  use integrator

   
  implicit none
    
  ! timing variables
  real*8 :: timer(10)
  integer :: j

      
  !!! initialize pepc library and MPI
  call pepc_initialize("pepc-f", my_rank, n_ranks, .true.)

  root = my_rank.eq.0

  timer(1) = get_time()

  call set_parameter()
  call init()

  open(unit=100, file="/lustre/jwork/fsnafzj/fsafzj08/test0.dat", status='unknown', position='rewind')
  open(unit=101, file="/lustre/jwork/fsnafzj/fsafzj08/test1.dat", status='unknown', position='rewind')
  open(unit=102, file="/lustre/jwork/fsnafzj/fsafzj08/test2.dat", status='unknown', position='rewind')
  open(unit=103, file="/lustre/jwork/fsnafzj/fsafzj08/test3.dat", status='unknown', position='rewind')
  open(unit=104, file="/lustre/jwork/fsnafzj/fsafzj08/test4.dat", status='unknown', position='rewind')
  open(unit=105, file="/lustre/jwork/fsnafzj/fsafzj08/test5.dat", status='unknown', position='rewind')
  open(unit=106, file="/lustre/jwork/fsnafzj/fsafzj08/test6.dat", status='unknown', position='rewind')
  open(unit=107, file="/lustre/jwork/fsnafzj/fsafzj08/test7.dat", status='unknown', position='rewind')

  call init_particles(plasma_particles)
  call init_wall_particles(wall_particles)

  particles(1:npp)=plasma_particles(:)
  particles(npp+1:npp+nwp)=wall_particles(:)

  timer(2) = get_time()

      !do j=1,np
      !    write(0,*)particles(j)%label,particles(j)%x,particles(j)%data%v
      !end do

  if(root) write(*,'(a,es12.4)') " === init time [s]: ", timer(2) - timer(1)

 
  DO step=0, nt-1
    if(root) then
      write(*,*) " "
      write(*,'(a,i12)')    " ====== computing step  :", step
      write(*,'(a,es12.4)') " ====== simulation time :", step * dt
    end if
    
    timer(3) = get_time()


    call pepc_particleresults_clear(particles, np)
    call pepc_grow_tree(np, tnp, particles)
    call pepc_traverse_tree(np, particles)

    call pepc_timber_tree()
    call pepc_restore_particles(np, particles)

    timer(4) = get_time()

    call boris_nonrel(particles)    

    call hits_on_boundaries(particles)
!do j=1,np
!if (particles(j)%label==1092) THEN
! write(*,*) step, " in pepc",particles(j)%data%B,particles(j)%x,particles(j)%data%q,particles(j)%data%v
!endif   
!end do  
do j=1,np
if (IsNaN(particles(j)%data%B(1))) THEN
 write(*,*) step, " NAN in pepc",particles(j)%label,particles(j)%data%B,particles(j)%x,particles(j)%data%q,particles(j)%data%v
endif   
end do  


      !do j=1,np
      !    write(step+1,*)particles(j)%label,particles(j)%x,particles(j)%data%v
      !end do

    timer(5) = get_time()
    
    if(root) write(*,'(a,es12.4)') " == time in pepc routines [s]                              : ", timer(4) - timer(3)
    if(root) write(*,'(a,es12.4)') " == time in integrator and particlehandling [s]            : ", timer(5) - timer(4)
    
  END DO 
 
  

  deallocate(particles)

  timer(6) = get_time()

  if(root) then
    write(*,*)            " "
    write(*,'(a)')        " ===== finished pepc simulation"
    write(*,'(a,es12.4)') " ===== total run time [s]: ", timer(6) - timer(1)
  end if

  !!! cleanup pepc and MPI
  call pepc_finalize()

end program pepc

