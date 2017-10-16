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


program pepc

    ! pepc modules
    use module_pepc
    use module_pepc_types
    use module_checkpoint
  
    ! frontend helper routines
    use helper
    use variables
    use fnpy

    implicit none
    include 'mpif.h'

    real(KIND=8) :: v2sum(3),vsum(3)


    !!! initialize pepc library and MPI
    call pepc_initialize("pepc-numpy", my_rank, n_ranks, .true.)

    root = my_rank.eq.0
  
    timer(1) = get_time()
   
    call init_after_resume()

    timer(2) = get_time()
    if (vtk) then
        call write_particles_vtk(particles)
        timer(3) = get_time()
        deallocate(particles)
        if (root) write(*,'(a)')        " ===== finished file conversion to VTK format"

    else
        if (root) allocate(all_particles(npart))
        if (root) allocate(particles_npy(0:19*npart-1))
        if (root) allocate(particles_npy_reshaped(0:18,0:npart-1))

        call MPI_GATHER(particles,np,MPI_TYPE_PARTICLE,all_particles,npart,MPI_TYPE_PARTICLE,0,MPI_COMM_WORLD,ierr)



        v2sum=0.0
        vsum=0.0
        DO j=1,npart
            IF(particles(j)%data%species==2) THEN
                v2sum=v2sum+particles(j)%data%v*particles(j)%data%v
                vsum=vsum+particles(j)%data%v
            END IF
        END DO
        v2sum=v2sum/1000000
        vsum=vsum/1000000
        write(*,*)v2sum*0.5*mp/e
        write(*,*)vsum


        if (root) then
        do j=0,npart-1
            particles_npy(j)=all_particles(j+1)%x(1)
            particles_npy(j+1*npart)=all_particles(j+1)%x(2)
            particles_npy(j+2*npart)=all_particles(j+1)%x(3)
            particles_npy(j+3*npart)=all_particles(j+1)%data%v(1)
            particles_npy(j+4*npart)=all_particles(j+1)%data%v(2)
            particles_npy(j+5*npart)=all_particles(j+1)%data%v(3)
            particles_npy(j+6*npart)=sqrt(particles_npy(j+5*npart)**2+particles_npy(j+3*npart)**2+particles_npy(j+4*npart)**2)
            particles_npy(j+7*npart)=all_particles(j+1)%results%e(1)
            particles_npy(j+8*npart)=all_particles(j+1)%results%e(2)
            particles_npy(j+9*npart)=all_particles(j+1)%results%e(3)
            particles_npy(j+10*npart)=sqrt(particles_npy(j+9*npart)**2+particles_npy(j+8*npart)**2+particles_npy(j+7*npart)**2)
            particles_npy(j+11*npart)=all_particles(j+1)%results%pot
            particles_npy(j+12*npart)=all_particles(j+1)%data%m
            particles_npy(j+13*npart)=all_particles(j+1)%data%q
            particles_npy(j+14*npart)=j+1
            particles_npy(j+15*npart)=all_particles(j+1)%label
            particles_npy(j+16*npart)=all_particles(j+1)%work
            !particles_npy(j+17*npart)=all_particles(j+1)%pid
            particles_npy(j+17*npart)=-1
            particles_npy(j+18*npart)=all_particles(j+1)%data%species
        end do
        endif
        particles_npy_reshaped=reshape(particles_npy,[19_8,npart],order=[2,1])
        timer(3) = get_time()

        call save_double(file_out, shape(particles_npy_reshaped), particles_npy_reshaped)

        deallocate(particles,all_particles,particles_npy_reshaped,particles_npy)
        if (root) write(*,'(a)')        " ===== finished file conversion to npy format"
    end if


    timer(4) = get_time()

    if(root) then
        write(*,*)            " "
        write(*,'(a,es12.4)') " ===== time for initialization [s]: ", timer(2) - timer(1)
        if (vtk) then
            write(*,'(a,es12.4)') " ===== time for writing data [s]:   ", timer(3) - timer(2)
            write(*,'(a,es12.4)') " ===== total run time [s]:          ", timer(4) - timer(1)
        else
            write(*,'(a,es12.4)') " ===== time for creating array [s]: ", timer(3) - timer(2)
            write(*,'(a,es12.4)') " ===== time for writing file [s]:   ", timer(4) - timer(3)
            write(*,'(a,es12.4)') " ===== total run time [s]:          ", timer(4) - timer(1)
        end if
    end if

    !!! cleanup pepc and MPI
    call pepc_finalize()

end program pepc

