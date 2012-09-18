! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2012 Juelich Supercomputing Centre, 
!                         Forschungszentrum Juelich GmbH,
!                         Germany
! 
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
  use zufall
  use module_cmdline

  contains

!======================================================================================
  subroutine set_parameter()
      
    use module_pepc
    use module_interaction_specific
    use module_mirror_boxes
    use module_interaction_specific
    implicit none
      
    integer, parameter :: fid = 12
  
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
    mirror_layers=1
    
    ! read in namelist file
    !call pepc_read_parameters_from_first_argument(read_para_file, para_file)
    call pepc_read_parameters_from_file_name(trim(input_file))

    if(root) write(*,'(a)') " == reading parameter file, section pepc-f: ", trim(input_file)
    open(fid,file=trim(input_file))
    read(fid,NML=pepcf)
    read(fid,NML=walk_para_smpss)
    close(fid)

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
      write(*,'(a,l)')      " == far_field_if_periodic            : ", include_far_field_if_periodic
      write(*,'(a,i12)')    " == mirror_layers                    : ", mirror_layers
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
    step=0

 
  end subroutine


!====================================================================================== 

  subroutine init_after_resume()
      
    use module_pepc
    use module_interaction_specific
    use module_mirror_boxes
    use module_checkpoint
    use module_interaction_specific
    implicit none
    include 'mpif.h'
      
    integer, parameter :: fid = 666
    integer :: global_max_label,local_max_label,i

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
    mirror_layers=1

    startstep=resume_step

    if (root) write(*,*)"=================================================="
    if (root) write(*,*)"Restarting from checkpoint at timestep",startstep
    if (root) write(*,*)"=================================================="

    call read_particles_mpiio(startstep,MPI_COMM_WORLD,my_rank,n_ranks,step,np,npart,particles,filename)    
    startstep=step

    call pepc_read_parameters_from_file_name(trim(input_file))

    if(root) write(*,'(a)') " == reading parameter file, section pepc-f: ", trim(input_file)
    open(fid,file=trim(input_file))
    read(fid,NML=pepcf)
    read(fid,NML=walk_para_smpss)
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
      write(*,'(a,l)')      " == far_field_if_periodic            : ", include_far_field_if_periodic
      write(*,'(a,i12)')    " == mirror_layers                    : ", mirror_layers
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

  end subroutine


!======================================================================================

    subroutine set_checkpoint()
        use module_checkpoint

        implicit none
        include 'mpif.h'

        integer, parameter :: fid = 666 
     
        npart=tnp
        !call write_particles_ascii(my_rank,step,np,particles,filename)
        call write_particles_mpiio(MPI_COMM_WORLD,my_rank,step,np,npart,particles,filename)
        
        if (root) then
            open(fid,file=trim(filename),STATUS='UNKNOWN', POSITION = 'APPEND')
            write(fid,NML=pepcf)
            write(fid,NML=walk_para_smpss)
            close(fid)

        endif

        if (root) write(*,*)
        if (root) write(*,*)"=================================================="
        if (root) write(*,'(a,i6)')"Saving particle data: checkpoint at timestep",step
        if (root) write(*,*)"=================================================="
        if (root) write(*,*)
    end subroutine

!====================================================================================== 
 
   subroutine init()
        implicit none
        integer :: rc,i,iy,iz
        !integer :: rsize
        !integer, allocatable :: rseed(:)
        real*8  :: dist_wpz,dist_wpy

        dist_wpz=dz/tnwpz
        dist_wpy=dy/tnwpy

        if(my_rank.eq.0) write(*,'(a)') " == [init] init particles "
    
        ! set initially number of local wall particles
        nwp = tnwp / n_ranks
        if(my_rank.eq.(n_ranks-1)) nwp = nwp + MOD(tnwp, n_ranks)
        allocate(wall_particles(nwp), stat=rc)
        if(rc.ne.0) write(*,*) " === wall particle allocation error!"
        ! set possible wall_particle positions
        allocate(wall_pos(tnwp,2), stat=rc)
        if(rc.ne.0) write(*,*) " === wallpos allocation error!"
 
        iz=1
        iy=1
        DO i=1,tnwp
            wall_pos(i,1) = iy*dist_wpy-dist_wpy/2.
            wall_pos(i,2) = iz*dist_wpz-dist_wpz/2.

            IF (iy==tnwpy) THEN
                iy=0
                iz=iz+1
            END IF        
            iy=iy+1
        END DO

        ! set initially number of local plasma particles
        npp = tnpp / n_ranks
        if(my_rank.eq.(n_ranks-1)) npp = npp + MOD(tnpp, n_ranks)
        allocate(plasma_particles(npp), stat=rc)
        if(rc.ne.0) write(*,*) " === plasma particle allocation error!"
  
        ! set initially number of local particles (wall+plasma)
        tnp=tnwp+tnpp
        np = npp + nwp
        allocate(particles(np), stat=rc)
        if(rc.ne.0) write(*,*) " === particle allocation error!"

        call init_rng()
        if (root) write(*,*) "========== Random Number Generator ========="
        if (root) write(*,'(a,i12)') " == Random Number Generator          : ", rng
        !call random_seed(size = rsize)
        !allocate(rseed(rsize))
        !rseed = my_rank + [(i*144,i=1,rsize)]
        !call random_seed(put = rseed)
        !deallocate(rseed)

    end subroutine


!======================================================================================

    subroutine init_periodicity()
        use module_mirror_boxes

        implicit none
        mirror_box_layers=mirror_layers
        periodicity=[.false.,.true.,.true.]
         
        t_lattice_1=[dx,0.0_8,0.0_8]
        t_lattice_2=[0.0_8,dy,0.0_8]
        t_lattice_3=[0.0_8,0.0_8,dz]

        LatticeOrigin=[xmin,ymin,zmin]

        !used to check if problem with potential is due to asymmetries when using mirror_layers
        !basically switches from nearest box to nearest image method
        !spatial_interaction_cutoff(2)=0.5*dy
        !spatial_interaction_cutoff(3)=0.5*dz


    end subroutine
        

!======================================================================================
  real*8 function get_time()
    implicit none
    include 'mpif.h'
    
    get_time = MPI_WTIME()
    
  end function get_time

!======================================================================================
  logical function hit_r(p)
    implicit none
    type(t_particle), intent(in) :: p

    if (p%x(1) > xmax .and. p%data%species /= 0) then
        hit_r = .true.
    else
        hit_r = .false.
    end if
     
  end function hit_r

!======================================================================================
  logical function hit_l(p)
    implicit none
    type(t_particle), intent(in) :: p

    if (p%x(1) < xmin .and. p%data%species /= 0) then
        hit_l = .true.
    else
        hit_l = .false.
    end if
     
  end function hit_l


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

end module
