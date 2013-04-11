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

module module_initialization

  use module_pepc_types
  use module_geometry
  use variables
  use zufall
  use helper
  use module_cmdline
  use particlehandling

  contains

!======================================================================================

  subroutine set_default_parameters()
      implicit none

      !constants
      pi = 2.0_8*acos(0.0_8)
      fc = 0.25_8/eps0/pi

      !set default parameter values
      diag_interval   =0
      checkp_interval =0
      guiding_centre_electrons=.false.
      Bx              = 0.
      By              = 0.
      Bz              = 0.
      nt              = 20
      dt              = 1e-3
      quelltyp        = 0
      fsup            = 1859.
      dx              = 0.
      dy              = 0.
      dz              = 0.

  end subroutine set_default_parameters

!======================================================================================

    subroutine init_source()
        implicit none

        integer :: fid=12

        x0_src=0.
        e1_src=0.
        e2_src=0.
        e3_src=0.
        quelltyp=0
        src_boundary=0


        IF(root) write(*,'(a,a)') " == reading parameter file, section source: ", trim(input_file)
        open(fid,file=trim(input_file))
        read(fid,NML=source_nml)
        close(fid)

        IF ((quelltyp==0)) THEN !surface source
            x0_src=0.
            e1_src=0.
            e2_src=0.
            e3_src=0.
            IF ((src_boundary<=0).or.(src_boundary>nb)) THEN
                IF (root) write(*,'(a)') "You have to select one of the boundaries as surface source"
                STOP
            ELSE
                IF(boundaries(src_boundary)%type==2) THEN
                    IF (root) write(*,'(a)') "Periodic boundary cannot be used as surface source"
                    STOP
                ELSE IF (boundaries(src_boundary)%type==3) THEN
                    IF (root) write(*,'(a)') "Open boundary cannot be used as surface source"
                    STOP
                END IF
                IF (root) write(*,'(a,i3,a,i3)') "Boundary ",src_boundary," chosen as surface source of type ",quelltyp
            END IF
        ELSE IF ((quelltyp==1).or.(quelltyp==2).or.(quelltyp==3)) THEN
            src_boundary=0
            IF (root) write(*,'(a,i2,a)') " == Volume source of type ",quelltyp," set. Parameters:"
            IF (root) write(*,'(a,3(1pe14.5E3))') " == x0: ",x0_src
            IF (root) write(*,'(a,3(1pe14.5E3))') " == e1: ",e1_src
            IF (root) write(*,'(a,3(1pe14.5E3))') " == e2: ",e2_src
            IF (root) write(*,'(a,3(1pe14.5E3))') " == e3: ",e3_src
        ELSE
            IF (root) write(*,'(a,i3,a)') " Source cannot be set. Type ",quelltyp," not available."
            STOP
        END IF



    end subroutine init_source
!======================================================================================

   subroutine init()
        implicit none
        integer :: rc,ispecies

        ! set initially number of local wall particles
        DO ispecies=0,nspecies-1
            npps(ispecies) = tnpps(ispecies) / n_ranks
            if(my_rank.eq.(n_ranks-1)) npps(ispecies) = npps(ispecies) + MOD(tnpps(ispecies), n_ranks)
        END DO

        tnp= SUM(tnpps)
        np = SUM(npps)
        allocate(particles(np), stat=rc)
        if(rc.ne.0) write(*,*) " === particle allocation error!"

        call init_particles(particles)

    end subroutine init

!======================================================================================

  subroutine set_parameters()
      
    use module_pepc
    use module_interaction_specific
    use module_mirror_boxes
    use module_checkpoint

    implicit none
    include 'mpif.h'

    integer, parameter :: fid = 12

    IF (do_resume) THEN
        startstep=resume_step       !resume step is read in read_args (from filename, maybe not the best idea...)

        if (root) write(*,*)"=================================================="
        if (root) write(*,*)"Restarting from checkpoint at timestep",startstep
        if (root) write(*,*)"=================================================="

        call read_particles_mpiio(startstep,MPI_COMM_WORLD,my_rank,n_ranks,step,np,npart,particles,filename)
    ELSE
        startstep=0
        step=0
    END IF

    call pepc_read_parameters_from_file_name(trim(input_file))

    if(root) write(*,'(a,a)') " == reading parameter file, section pepc-f: ", trim(input_file)
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

    B=sqrt(Bx**2+By**2+Bz**2)
    IF (B.ne.0.) THEN
        r_lamor=1.
    ELSE
        r_lamor=0.
    END IF


    !ne=0.5*tnpp*fsup/(dx*dy*dz)
    !omega_p=sqrt(ne*e*e/eps0/me)



    call init_rng()
    call init_periodicity()
    call pepc_prepare(3)


  end subroutine set_parameters

!====================================================================================== 

  subroutine init_after_resume()
    use module_pepc
    use module_interaction_specific
    use module_mirror_boxes
    use module_checkpoint
    implicit none
    include 'mpif.h'
      
    integer, parameter :: fid = 666
    integer :: global_max_label,local_max_label,ip,ib
    real(KIND=8) :: q_loc(nb),q_glob(nb)
    logical :: hit


    tnp=npart
    npps=0
    tnpps=0

    DO ip=1,np
        npps(particles(ip)%data%species)=npps(particles(ip)%data%species)+1
    END DO

    call MPI_ALLREDUCE(npps, tnpps, nspecies, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

    write(*,*)my_rank,npps,tnpps

    IF (tnpps(0)/=count_wallparticles()) THEN
        IF (root) write(*,*) "Number of Wallparticles in checkpoint does not match number in input."
        IF (root) write(*,*) "You can't change the number of Wallparticles when resuming"
        STOP
    END IF



    q_loc=0.
    q_glob=0.
    DO ip=1, np                                            ! calculate charge by adding up charges of wallparticles after resume
        IF (particles(ip)%data%species/=0) CYCLE
        DO ib=1,nb
            call check_hit(particles(ip)%x(1),particles(ip)%x(2),particles(ip)%x(3),boundaries(ib),hit)
            IF(hit) THEN
                q_loc(ib) = q_loc(ib) + particles(ip)%data%q
            END IF
        END DO
    END DO

    DO ib=1,nb
        call MPI_ALLREDUCE(q_loc(ib), q_glob(ib), 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
        boundaries(ib)%q_tot=q_glob(ib)
    END DO


    global_max_label = 0
    local_max_label  = 0
    DO ip=1,np
        if (particles(ip)%label > local_max_label) local_max_label = particles(ip)%label
    END DO
    call MPI_ALLREDUCE(local_max_label, global_max_label, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
    next_label = global_max_label + 1

  end subroutine


!======================================================================================

    subroutine init_periodicity()
        use module_mirror_boxes

        implicit none
         
        t_lattice_1=[dx,0.0_8,0.0_8]
        t_lattice_2=[0.0_8,dy,0.0_8]
        t_lattice_3=[0.0_8,0.0_8,dz]

        LatticeOrigin=[xmin,ymin,zmin]


    end subroutine

end module
