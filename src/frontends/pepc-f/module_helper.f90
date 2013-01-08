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
  subroutine init_files()
      implicit none

      !if(root) open(unit=timing_out,file='timing.out',status='UNKNOWN',position='APPEND')
      !if(root) open(unit=recycling_out,file='recycling.out',status='UNKNOWN',position='APPEND')
      if(root) open(unit=out,file='pepcf.out',status='UNKNOWN',position='APPEND')

  end subroutine init_files

!======================================================================================
  subroutine close_files()
      implicit none

      !if(root) close(timing_out)
      !if(root) close(recycling_out)
      if(root) close(out)

  end subroutine close_files

!======================================================================================
  subroutine flush_files()
      implicit none

      call close_files()
      call init_files()

  end subroutine flush_files

!======================================================================================
  subroutine set_default_parameters()
      implicit none

      !set default parameter values
      fixed_npp=.true.
      diag_interval   =0
      checkp_interval =0
      open_sides =.false.
      guiding_centre_electrons=.false.
      mirror_layers=1
      Bx              = 0.
      By              = 0.
      Bz              = 0.
      nt              = 20
      dt              = 1e-3
      te_ev           = 80.
      ti_ev           = 80.
      quelltyp        = 0
      tnwpy           = 0
      tnwpz           = 0
      fsup            = 1859.
      tnpp            = 1
      tfpp            = 1
      dx              = 0.
      dy              = 0.
      dz              = 0.

  end subroutine set_default_parameters

!======================================================================================

  subroutine set_parameters()
      
    use module_pepc
    use module_interaction_specific
    use module_mirror_boxes

    implicit none
      
    integer, parameter :: fid = 12

    !constants
    pi = 2.0_8*acos(0.0_8)
    fc = 0.25_8/eps0/pi
  
    call pepc_read_parameters_from_file_name(trim(input_file))

    if(root) write(*,'(a)') " == reading parameter file, section pepc-f: ", trim(input_file)
    open(fid,file=trim(input_file))
    read(fid,NML=pepcf)
    read(fid,NML=walk_para_smpss)
    close(fid)

    IF (dx.eq.0. .or. dy.eq.0. .or. dz.eq.0.) THEN  !one of dx, dy, dz not set
        IF (root) write(*,*) "Geometry not set correctly. dx, dy and dz cannot be 0"
        STOP
    END IF
    !========== used if domain does not start in 0,0,0; origin can be set in input later
    !right now it is still 0,0,0, but the subroutines are already adapted
    xmin=0.0_8
    xmax=dx+xmin
    ymin=0.0_8
    ymax=dy+ymin
    zmin=0.0_8
    zmax=dz+zmin

    IF (tnwpy<=0 .or. tnwpz<=0) THEN  !no wall particles
        IF (root) write(*,*) "No Wallparticles set. tnwpy and tnwpz have to be positive."
        STOP
    END IF
    tnwp=tnwpy*tnwpz

    B=sqrt(Bx**2+By**2+Bz**2)
    IF (B.ne.0.) THEN
        r_lamor=mp*sqrt(ti_ev*e/mp)/e/B
    ELSE
        r_lamor=0.
    END IF

    IF (fixed_npp) THEN
        !calculate Plasmaparameters
        ne=0.5*tnpp*fsup/(dx*dy*dz)
        l_debye=sqrt(eps0*te_ev/e/ne)
        omega_p=sqrt(ne*e*e/eps0/me)
        tfpp=0
    END IF

    startstep=0
    step=0

    call init_rng()
    call init_periodicity()
    call pepc_prepare(3)


    IF (root) call write_parameters()

  end subroutine set_parameters

!======================================================================================

  SUBROUTINE write_parameters()
      use module_interaction_specific

      implicit none

      if (fixed_npp) write(*,'(a,i12)')       " == total number of simulated plasma particles: ", tnpp
      if (.not. fixed_npp) write(*,'(a,i12)') " == simulated plasma-particleflux per timestep: ", tfpp
      write(*,'(a,i12)')    " == number of time steps             : ", nt
      write(*,'(a,es12.4)') " == time step                        : ", dt
      if (fixed_npp) write(*,'(a,es12.4)') " == timestep * omega_p               : ", dt*omega_p
      write(*,'(a,es12.4)') " == simulation volume                : ", dx*dy*dz
      write(*,'(a,es12.4)') " == superparticle factor             : ", fsup
      if (fixed_npp) write(*,'(a,es12.4)') " == particles / debye spehere        : ", 0.5*tnpp/(dx*dy*dz)*l_debye**3
      write(*,'(a,i12)')    " == type of source distribution      : ", quelltyp
      write(*,'(a,l6)')     " == open side walls                  : ", open_sides
      write(*,'(a,l6)')     " == far_field_if_periodic            : ", include_far_field_if_periodic
      write(*,'(a,l6)')     " == do_restore_particles             : ", do_restore_particles
      write(*,'(a,i12)')    " == mirror_layers                    : ", mirror_layers
      write(*,*)
      write(*,*) "========== Magnetic Field ========="
      write(*,'(a,f12.4)')  " == Bx                               : ", Bx
      write(*,'(a,f12.4)')  " == By                               : ", By
      write(*,'(a,f12.4)')  " == Bz                               : ", Bz
      write(*,*)
      write(*,*) "========== Simulation Domain ========="
      if (B.ne.0.) then
          if (fixed_npp) write(*,'(a,12X,es10.2,a,es10.2)')  " == dx (debye_lengths,m)             : ", dx/l_debye,", ",dx
          if (fixed_npp) write(*,'(a,es10.2,a,es10.2,a,es10.2)')  " == dy (gyro_radii, debye_legths, m) : ", dy/r_lamor,", ",dy/l_debye,", ",dy
          if (fixed_npp) write(*,'(a,es10.2,a,es10.2,a,es10.2)')  " == dz (gyro_radii, debye_legths, m) : ", dz/r_lamor,", ",dz/l_debye,", ",dz
          if (.not.fixed_npp) write(*,'(a,12X,es10.2)')  " == dx (m)             : ", dx
          if (.not.fixed_npp) write(*,'(a,es10.2,a,es10.2)')  " == dy (gyro_radii, m) : ", dy/r_lamor,", ",dy
          if (.not. fixed_npp) write(*,'(a,es10.2,a,es10.2)')  " == dz (gyro_radii, m) : ", dz/r_lamor,", ",dz
      else
          if (fixed_npp) write(*,'(a,es10.2)')  " == dx (m)             : ", dx
          if (fixed_npp) write(*,'(a,es10.2)')  " == dy (m)             : ", dy
          if (fixed_npp) write(*,'(a,es10.2)')  " == dz (m)             : ", dz
          if (.not. fixed_npp) write(*,'(a,es10.2)')  " == dx (m)             : ", dx
          if (.not. fixed_npp) write(*,'(a,es10.2)')  " == dy (m)             : ", dy
          if (.not. fixed_npp) write(*,'(a,es10.2)')  " == dz (m)             : ", dz
      end if
      write(*,*)
      write(*,*) "========== Plasmaparameters ========="
      write(*,'(a,f12.4)')  " == TE [eV] (sourcetemperature)      : ", te_ev
      write(*,'(a,f12.4)')  " == TI [eV] (sourcetemperature)      : ", ti_ev
      if (fixed_npp) write(*,'(a,es12.4)') " == Electron/Ion density [m^-3]      : ", ne
      write(*,'(a,es12.4)') " == Gyro radius [m]                  : ", r_lamor
      if (fixed_npp) write(*,'(a,es12.4)') " == Debye length [m]                 : ", l_debye
      if (fixed_npp) write(*,'(a,es12.4)') " == Plasmafrequency [s^-1]           : ", omega_p
      write(*,*)
      write(*,*) "========== Wall Particles ========="
      write(*,'(a,i12)')    " == number of wall partiles          : ", tnwp
      if (fixed_npp) write(*,'(a,es12.4)') " == wall particles / l_debye in y    : ", tnwpy/dy*l_debye
      if (fixed_npp) write(*,'(a,es12.4)') " == wall particles / l_debye in z    : ", tnwpz/dz*l_debye
      write(*,*) "========== Random Number Generator ========="
      write(*,'(a,i12)') " == Random Number Generator          : ", rng

  END SUBROUTINE write_parameters

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

    end subroutine init


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
    integer :: nwp=0,npp=0

    ! set default parameter values
    call set_default_parameters()

    startstep=resume_step

    if (root) write(*,*)"=================================================="
    if (root) write(*,*)"Restarting from checkpoint at timestep",startstep
    if (root) write(*,*)"=================================================="

    call read_particles_mpiio(startstep,MPI_COMM_WORLD,my_rank,n_ranks,step,np,npart,particles,filename)    

    call pepc_read_parameters_from_file_name(trim(input_file))

    if(root) write(*,'(a)') " == reading parameter file, section pepc-f: ", trim(input_file)
    open(fid,file=trim(input_file))
    read(fid,NML=pepcf)
    read(fid,NML=walk_para_smpss)
    close(fid)

    !constants
    pi = 2.0_8*acos(0.0_8)
    fc = 0.25_8/eps0/pi 

    tnp=npart
    DO i=1,np
        IF (particles(i)%data%species==0) THEN
            nwp=nwp+1
        ELSE
            npp=npp+1
        END IF
    END DO

    allocate(wall_particles(nwp))
    allocate(plasma_particles(npp))

    call MPI_ALLREDUCE(nwp, tnwp, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(npp, tnpp, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

    IF (tnwp/=(tnwpy*tnwpz)) THEN
        IF (root) write(*,*) "Number of Wallparticles in checkpoint does not match number in input."
        IF (root) write(*,*) "You can't change the number of Wallparticles when resuming"
        STOP
    END IF

    IF (dx.eq.0. .or. dy.eq.0. .or. dz.eq.0.) THEN  !one of dx, dy, dz not set
        IF (root) write(*,*) "Geometry not set correctly. dx, dy and dz cannot be 0"
        STOP
    END IF

    !========== used if domain does not start in 0,0,0; origin can be set in input later
    !right now it is still 0,0,0, but the subroutines are already adapted
    xmin=0.0_8
    xmax=dx+xmin
    ymin=0.0_8
    ymax=dy+ymin
    zmin=0.0_8
    zmax=dz+zmin

    B=sqrt(Bx**2+By**2+Bz**2)
    IF (B.ne.0.) THEN
        r_lamor=mp*sqrt(ti_ev*e/mp)/e/B
    ELSE
        r_lamor=0.
    END IF

    IF (fixed_npp) THEN
        !calculate Plasmaparameters
        ne=0.5*tnpp*fsup/(dx*dy*dz)
        l_debye=sqrt(eps0*te_ev/e/ne)
        omega_p=sqrt(ne*e*e/eps0/me)
        tfpp=0
    END IF

    call init_rng()
    call init_periodicity()
    call pepc_prepare(3)

    global_max_label = 0
    local_max_label  = 0
    DO i=1,np
        if (particles(i)%label > local_max_label) local_max_label = particles(i)%label
    END DO
    call MPI_ALLREDUCE(local_max_label, global_max_label, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
    next_label = global_max_label + 1

    if (root) call write_parameters()

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

        if(root) write(*,'(a,i6)') " == [set_checkpoint] checkpoint at timestep",step

    end subroutine



!======================================================================================

    subroutine init_periodicity()
        use module_mirror_boxes

        implicit none
        mirror_box_layers=mirror_layers
        periodicity=periodicity_in
         
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


!======================================================================================
!sortiert wallparticles nach hinten
  subroutine sort_particles(p)
    implicit none

    type(t_particle), allocatable, intent(inout) :: p(:)
    type(t_particle)                               :: p_help_vorn


    integer :: i,iback

    iback=np

    DO i=1,np
        if (p(i)%data%species == 0) then
            p_help_vorn=p(i)
            do while (p(iback)%data%species==0)
                iback=iback-1
            end do
            if (iback<=i) exit
            p(i)=p(iback)
            p(iback)=p_help_vorn
        end if
    END DO
  end subroutine

!=======================================================================================
  subroutine pepc_tree_diagnostics()
      use module_vtk
      use module_treediags
      use module_pepc

      implicit none

      integer :: vtk_step

      if (step == 0) then
          vtk_step = VTK_STEP_FIRST
      else if (step == nt + startstep) then
          vtk_step = VTK_STEP_LAST
      else
          vtk_step = VTK_STEP_NORMAL
      endif

      call pepc_statistics(step)
      call write_branches_to_vtk(step,   step*dt, vtk_step)
      call write_spacecurve_to_vtk(step, step*dt, vtk_step, particles)

  end subroutine pepc_tree_diagnostics

end module
