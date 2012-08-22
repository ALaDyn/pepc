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

!!!!!!!!!!!!!!!!!!!!
!! helper module
!!!!!!!!!!!!!!!!!!!!

module helper

  use module_pepc_types
  use variables

  contains


<<<<<<< .mine
=======

  
     ! set default parameter values
     
    tnpp             = 1337
    nt              = 20
    dt              = 1e-3
    Bx              = 0.
    By              = 0.
    Bz              = 0.
    delx            = 80.
    dely            = 3.
    delz            = 3.
    ni              = 1e17
    ne              = 1e17
    te_ev           = 80.
    ti_ev           = 80.
    quelltyp        = 0
    tnwpy           = 0
    tnwpz           = 0
    dx=0.
    dy=0.
    dz=0.
    diag_interval   =0
    checkp_interval =0
    open_sides =.false.
    guiding_centre_electrons=.false.
    
    ! read in namelist file
    call pepc_read_parameters_from_first_argument(read_para_file, para_file)

    if (read_para_file) then
      if(root) write(*,'(a)') " == reading parameter file, section pepc-f: ", para_file
      open(fid,file=para_file)
      read(fid,NML=pepcf)
      read(fid,NML=walk_para_smpss)
      close(fid)
    else
      if(root) write(*,*) " == no param file, using default parameter "
    end if    

    tnpp = tnpp

    !constants
    pi = 2.0_8*acos(0.0_8)
    fc = 0.25_8/eps0/pi 

    !calculate Plasmaparameters
    l_debye=sqrt(eps0*te_ev/e/ne)
    omega_p=sqrt(ne*e*e/eps0/me)
    B=sqrt(Bx**2+By**2+Bz**2)
    
    IF (B.ne.0.) THEN
        r_lamor=mp*sqrt(ti_ev*e/mp)/e/B
    ELSE
        r_lamor=0.
    END IF
    IF (dx.eq.0. .or. dy.eq.0. .or. dz.eq.0.) THEN  !one of dx, dy, dz not set
        IF (B.ne.0.) THEN
            dy=dely*l_debye
            dz=delz*l_debye
        ELSE
            dy=dely*l_debye
            dz=delz*l_debye
        END IF
        dx=delx*l_debye
    ELSE                                             ! dx, dy, dz set in input file
        IF (B.ne.0.) THEN
            dely=dy/l_debye
            delz=dz/l_debye
        ELSE
            dely=dy/l_debye
            delz=dz/l_debye
        END IF                           
        delx=dx/l_debye
    END IF

    IF (tnwpy.eq.0 .or. tnwpz.eq.0) THEN
        tnwpy=int(10.*dy/l_debye)
        tnwpz=int(10.*dz/l_debye)
    END IF

    tnwp=tnwpy*tnwpz

    fsup = (ne+ni)*dx*dy*dz/tnpp

!========== used if domain does not start in 0,0,0; origin can be set in input later
!right now it is still 0,0,0, but the subroutines are already adapted
    xmin=0.0_8
    xmax=dx+xmin
    ymin=0.0_8
    ymax=dy+ymin
    zmin=0.0_8
    zmax=dz+zmin

    call init_periodicity()

    call pepc_prepare(3)

    if(root) then
      write(*,'(a,i12)')    " == total number of plasma particles : ", tnpp
      write(*,'(a,i12)')    " == number of time steps             : ", nt
      write(*,'(a,es12.4)') " == time step                        : ", dt
      write(*,'(a,es12.4)') " == timestep * omega_p               : ", dt*omega_p
      write(*,'(a,es12.4)') " == simulation volume                : ", dx*dy*dz
      write(*,'(a,es12.4)') " == superparticle factor             : ", fsup
      write(*,'(a,es12.4)') " == particles / debye spehere        : ", 0.5*tnpp/(dx*dy*dz)*l_debye**3
      write(*,'(a,i12)')    " == type of source distribution      : ", quelltyp
      write(*,'(a,l)')      " == open side walls                  : ", open_sides
      write(*,*)
      write(*,*) "========== Magnetic Field ========="
      write(*,'(a,f12.4)')  " == Bx                               : ", Bx
      write(*,'(a,f12.4)')  " == By                               : ", By
      write(*,'(a,f12.4)')  " == Bz                               : ", Bz
      write(*,*)
      write(*,*) "========== Simulation Domain ========="
      if (B.ne.0.) then
          write(*,'(a,12X,es10.2,a,es10.2)')  " == dx (debye_lengths,m)             : ", delx,", ",dx
          write(*,'(a,es10.2,a,es10.2,a,es10.2)')  " == dy (gyro_radii, debye_legths, m) : ", dy/r_lamor,", ",dy/l_debye,", ",dy
          write(*,'(a,es10.2,a,es10.2,a,es10.2)')  " == dz (gyro_radii, debye_legths, m) : ", dz/r_lamor,", ",dz/l_debye,", ",dz
      else
          write(*,'(a,es10.2,a,es10.2)')  " == dx (debye_lengths,m)             : ", delx,", ",dx
          write(*,'(a,es10.2,a,es10.2)')  " == dy (debye_lengths,m)             : ", dely,", ",dy
          write(*,'(a,es10.2,a,es10.2)')  " == dz (debye_lengths,m)             : ", delz,", ",dz
      end if
      write(*,*)
      write(*,*) "========== Plasmaparameters ========="
      write(*,'(a,f12.4)')  " == TE [eV] (sourcetemperature)      : ", te_ev
      write(*,'(a,f12.4)')  " == TI [eV] (sourcetemperature)      : ", ti_ev
      write(*,'(a,se12.4)') " == Ion density [m^-3]               : ", ni
      write(*,'(a,es12.4)') " == Electron density [m^-3]          : ", ne
      write(*,'(a,es12.4)') " == Gyro radius [m]                  : ", r_lamor
      write(*,'(a,es12.4)') " == Debye length [m]                 : ", l_debye
      write(*,'(a,es12.4)') " == Plasmafrequency [s^-1]           : ", omega_p
      write(*,*)
      write(*,*) "========== Wall Particles ========="
      write(*,'(a,i12)')    " == number of wall partiles          : ", tnwp
      write(*,'(a,es12.4)') " == wall particles / l_debye in y    : ", tnwpy/dy*l_debye
      write(*,'(a,es12.4)') " == wall particles / l_debye in z    : ", tnwpz/dz*l_debye
      

    end if
    startstep=0

 
  end subroutine


!====================================================================================== 

>>>>>>> .r3385
  subroutine init_after_resume()
      
    use module_pepc
    use module_interaction_specific
    use module_checkpoint

    implicit none
    include 'mpif.h'
      
<<<<<<< .mine
    vtk=.false.
    cmd_args = COMMAND_ARGUMENT_COUNT()
    if (cmd_args > 1) then
        call GET_COMMAND_ARGUMENT(1, file_in)
        call GET_COMMAND_ARGUMENT(2, file_out)
    end if
    if (file_out=="vtk") vtk=.true.
    call read_particles_mpiio_from_filename(MPI_COMM_WORLD,my_rank,n_ranks,step,np,npart,particles,file_in)
=======
    integer, parameter :: fid = 666
    integer :: global_max_label,local_max_label,i
>>>>>>> .r3385

<<<<<<< .mine
    if (vtk) then
        filenameh = trim(file_in)//".h"
        open(123, file=trim(filenameh),action='read')
        read(123,NML=pepcf)
        close(123)
    end if
=======
    ! set default parameter values
    tnpp             = 1337
    nt              = 20
    dt              = 1e-3
    Bx              = 0.
    By              = 0.
    Bz              = 0.
    delx            = 80.
    dely            = 3.
    delz            = 3.
    ni              = 1e17
    ne              = 1e17
    te_ev           = 80.
    ti_ev           = 80.
    quelltyp        = 0
    tnwpy           = 0
    tnwpz           = 0
    dx=0.
    dy=0.
    dz=0.
    diag_interval   =0
    checkp_interval =0
    open_sides =.false.
    guiding_centre_electrons=.false.


    if (root) CALL SYSTEM ('python last_succ_ts.py > resume_int.dat')
    if (root) open(unit = 666, file = "resume_int.dat")
    if (root) read(666,*)startstep
    if (root) close(666)
    if (root) write(*,*)
    if (root) write(*,*)"=================================================="
    if (root) write(*,*)"Restarting from checkpoint at timestep",startstep
    if (root) write(*,*)"=================================================="
    if (root) write(*,*)
    if (root) CALL SYSTEM ('rm resume_int.dat')
    call MPI_BCAST(startstep, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      
    call read_particles_mpiio(startstep,MPI_COMM_WORLD,my_rank,n_ranks,step,np,npart,particles,filename)    
    startstep=step

    cmd_args = COMMAND_ARGUMENT_COUNT()
    if (cmd_args > 1) then
        call GET_COMMAND_ARGUMENT(2, argument2)
        !if (root) write (*,*) "cmd_args > 1",cmd_args
        filename=argument2
        call pepc_read_parameters_from_file_name(filename)
    end if

    if(root) write(*,'(a)') " == reading parameter file, section pepc-f: ", filename
    open(fid,file=filename)
    read(fid,NML=pepcf)
    close(fid)

    !constants
    pi = 2.0_8*acos(0.0_8)
    fc = 0.25_8/eps0/pi 

    !calculate Plasmaparameters
    l_debye=sqrt(eps0*te_ev/e/ne)
    omega_p=sqrt(ne*e*e/eps0/me)
    B=sqrt(Bx**2+By**2+Bz**2)
    
    IF (B.ne.0.) THEN
        r_lamor=mp*sqrt(ti_ev*e/mp)/e/B
    ELSE
        r_lamor=0.
    END IF
    IF (dx.eq.0. .or. dy.eq.0. .or. dz.eq.0.) THEN  !one of dx, dy, dz not set
        IF (B.ne.0.) THEN
            dy=dely*l_debye
            dz=delz*l_debye
        ELSE
            dy=dely*l_debye
            dz=delz*l_debye
        END IF
        dx=delx*l_debye
    ELSE                                             ! dx, dy, dz set in input file
        IF (B.ne.0.) THEN
            dely=dy/l_debye
            delz=dz/l_debye
        ELSE
            dely=dy/l_debye
            delz=dz/l_debye
        END IF                           
        delx=dx/l_debye
    END IF

    IF (tnwpy.eq.0 .or. tnwpz.eq.0) THEN
        tnwpy=int(10.*dy/l_debye)
        tnwpz=int(10.*dz/l_debye)
    END IF

    tnwp=tnwpy*tnwpz
    tnp=tnwp+tnpp

    fsup = (ne+ni)*dx*dy*dz/tnpp

    call init_rng()
    !call random_seed(size = rsize)
    !allocate(rseed(rsize))
    !rseed = my_rank + [(i*144,i=1,rsize)]
    !call random_seed(put = rseed)
    !deallocate(rseed)

!========== used if domain does not start in 0,0,0; origin can be set in input later
!right now it is still 0,0,0, but the subroutines are already adapted
    xmin=0.0_8
    xmax=dx+xmin
    ymin=0.0_8
    ymax=dy+ymin
    zmin=0.0_8
    zmax=dz+zmin

    call init_periodicity()

>>>>>>> .r3385
    call pepc_prepare(3)

<<<<<<< .mine
=======
    global_max_label = 0
    local_max_label  = 0
    DO i=1,np
        if (particles(i)%label > local_max_label) local_max_label = particles(i)%label
    END DO
    call MPI_ALLREDUCE(local_max_label, global_max_label, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
    next_label = global_max_label + 1

    if(root) then
      write(*,'(a,i12)')    " == total number of plasma particles : ", tnpp
      write(*,'(a,i12)')    " == number of time steps             : ", nt
      write(*,'(a,es12.4)') " == time step                        : ", dt
      write(*,'(a,es12.4)') " == timestep * omega_p               : ", dt*omega_p
      write(*,'(a,es12.4)') " == simulation volume                : ", dx*dy*dz
      write(*,'(a,es12.4)') " == superparticle factor             : ", fsup
      write(*,'(a,es12.4)') " == particles / debye spehere        : ", 0.5*tnpp/(dx*dy*dz)*l_debye**3
      write(*,'(a,i12)')    " == type of source distribution      : ", quelltyp
      write(*,'(a,l)')      " == open side walls                  : ", open_sides
      write(*,*)
      write(*,*) "========== Magnetic Field ========="
      write(*,'(a,f12.4)')  " == Bx                               : ", Bx
      write(*,'(a,f12.4)')  " == By                               : ", By
      write(*,'(a,f12.4)')  " == Bz                               : ", Bz
      write(*,*)
      write(*,*) "========== Simulation Domain ========="
      if (B.ne.0.) then
          write(*,'(a,12X,es10.2,a,es10.2)')  " == dx (debye_lengths,m)             : ", delx,", ",dx
          write(*,'(a,es10.2,a,es10.2,a,es10.2)')  " == dy (gyro_radii, debye_legths, m) : ", dy/r_lamor,", ",dy/l_debye,", ",dy
          write(*,'(a,es10.2,a,es10.2,a,es10.2)')  " == dz (gyro_radii, debye_legths, m) : ", dz/r_lamor,", ",dz/l_debye,", ",dz
      else
          write(*,'(a,es10.2,a,es10.2)')  " == dx (debye_lengths,m)             : ", delx,", ",dx
          write(*,'(a,es10.2,a,es10.2)')  " == dy (debye_lengths,m)             : ", dely,", ",dy
          write(*,'(a,es10.2,a,es10.2)')  " == dz (debye_lengths,m)             : ", delz,", ",dz
      end if
      write(*,*)
      write(*,*) "========== Plasmaparameters ========="
      write(*,'(a,f12.4)')  " == TE [eV] (sourcetemperature)      : ", te_ev
      write(*,'(a,f12.4)')  " == TI [eV] (sourcetemperature)      : ", ti_ev
      write(*,'(a,se12.4)') " == Ion density [m^-3]               : ", ni
      write(*,'(a,es12.4)') " == Electron density [m^-3]          : ", ne
      write(*,'(a,es12.4)') " == Gyro radius [m]                  : ", r_lamor
      write(*,'(a,es12.4)') " == Debye length [m]                 : ", l_debye
      write(*,'(a,es12.4)') " == Plasmafrequency [s^-1]           : ", omega_p
      write(*,*)
      write(*,*) "========== Wall Particles ========="
      write(*,'(a,i12)')    " == number of wall partiles          : ", tnwp
      write(*,'(a,es12.4)') " == wall particles / l_debye in y    : ", tnwpy/dy*l_debye
      write(*,'(a,es12.4)') " == wall particles / l_debye in z    : ", tnwpz/dz*l_debye
      write(*,*) "========== Random Number Generator ========="
      write(*,'(a,i12)')  " == Random Number Generator          : ", rng
      
    end if

>>>>>>> .r3385
  end subroutine

!===============================================================================

    subroutine write_particles_vtk(p)
        use module_vtk
        implicit none

        type(t_particle), allocatable, intent(in) :: p(:)

        integer :: i
        type(vtkfile_unstructured_grid) :: vtk
        integer :: vtk_step
        real*8 :: time
        real*8 :: ta, tb

        ta = get_time()

        time = dt * step

        vtk_step = VTK_STEP_NORMAL

        call vtk%create_parallel("particles", step, my_rank, n_ranks, time, vtk_step)
        call vtk%write_headers(np, 0)
        call vtk%startpoints()
        call vtk%write_data_array("xyz", np, p(:)%x(1), p(:)%x(2), p(:)%x(3))
        call vtk%finishpoints()
        call vtk%startpointdata()
        call vtk%write_data_array("velocity", np, p(:)%data%v(1), p(:)%data%v(2), p(:)%data%v(3))
        call vtk%write_data_array("el_field", np, p(:)%results%e(1),p(:)%results%e(2), p(:)%results%e(3))
        call vtk%write_data_array("el_pot", np, p(:)%results%pot)
        call vtk%write_data_array("charge", np, p(:)%data%q)
        call vtk%write_data_array("mass", np, p(:)%data%m)
        call vtk%write_data_array("work", np, p(:)%work)
        call vtk%write_data_array("pelabel", np, p(:)%label)
        call vtk%write_data_array("local index", np, [(i,i=1,np)])
        call vtk%write_data_array("processor", np, p(:)%pid)
        call vtk%finishpointdata()
        call vtk%dont_write_cells()
        call vtk%write_final()
        call vtk%close()

        tb = get_time()

        !if(root) write(*,'(a,es12.4)') " == [write particles] time in vtk output [s]      : ", tb - ta

    end subroutine write_particles_vtk

!======================================================================================
  real*8 function get_time()
    implicit none
    include 'mpif.h'
    
    get_time = MPI_WTIME()
    
  end function get_time

<<<<<<< .mine
=======
!======================================================================================
  logical function hit_wall(p)
    implicit none
    type(t_particle), intent(in) :: p

    if (p%x(1) > xmax .and. p%label > 0) then
        hit_wall = .true.
    else
        hit_wall = .false.
    end if
     
  end function hit_wall

!======================================================================================
  logical function hit_src(p)
    implicit none
    type(t_particle), intent(in) :: p

    if (p%x(1) < xmin .and. p%label > 0) then
        hit_src = .true.
    else
        hit_src = .false.
    end if
     
  end function hit_src


!======================================================================================
  logical function hit_side(p)
    implicit none
    type(t_particle), intent(in) :: p

    if (((p%x(2) < ymin).or.(p%x(2) > ymax).or.(p%x(3) < zmin).or.(p%x(3) > zmax)) .and. p%label > 0) then
        hit_side = .true.
    else
        hit_side = .false.
    end if

  end function hit_side

!======================================================================================
  subroutine reallocate_particles(list, oldn, newn)
    implicit none
  
    type(t_particle), allocatable, intent(inout) :: list(:)
    integer, intent(in) :: oldn, newn
    
    type(t_particle) :: tmp_p(oldn)
    
    tmp_p(1:oldn) = list(1:oldn)
    deallocate(list)
    allocate(list(newn))
    list(1:oldn) = tmp_p
  
  end subroutine

>>>>>>> .r3385
end module
