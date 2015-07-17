! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2015 Juelich Supercomputing Centre, 
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

  use module_pepc_kinds
  use module_pepc_types
  use module_geometry
  use variables
  use zufall
  use helper
  use module_cmdline
  use particlehandling

  implicit none
  real(KIND=8)  :: fsup_at_checkpoint

  contains

!======================================================================================

  subroutine set_default_parameters()
      implicit none

      !constants
      pi = 2.0_8*acos(0.0_8)
      fc = 0.25_8/eps0/pi

      !set default parameter values
      rngseed = 0
      diag_interval   = 0
      checkp_interval = 0
      vtk_interval    = 0
      npy_interval    = 0
      reflux_interval = 1
      guiding_centre_electrons=.false.
      Bx              = 0.
      By              = 0.
      Bz              = 0.
      nt              = 20
      dt              = 1e-3
      fsup            = 1859.
      dx              = 0.
      dy              = 0.
      dz              = 0.
      xmax            = 0.
      xmin            = 0.
      ymax            = 0.
      ymin            = 0.
      zmax            = 0.
      zmin            = 0.

      bool_energy_resolved_hits = .true.
      bool_angle_resolved_hits = .true.
      bool_age_resolved_hits = .true.
      bool_space_resolved_hits = .true.
      nbins_energy_resolved_hits = 80
      nbins_angle_resolved_hits = 45
      nbins_age_resolved_hits = 100
      nbins_e1_space_resolved_hits = 10
      nbins_e2_space_resolved_hits = 10
      ehit_max_in_T = 10.
      agehit_max_in_t_trav_ion = 2.

      bool_diag_bins_cylinder = .false.
      bool_avg_btwn_diag_steps = .true.
      diag_bins_x=1
      diag_bins_y=1
      diag_bins_z=1

      diag_bins_vpar = 100
      v_grid_max = 4.0_8
      bool_velocity_diag = .false.

      bool_hockney_diag = .true.
      hockney_start_step = 100

  end subroutine set_default_parameters

!======================================================================================

   subroutine init()
        implicit none
        integer :: rc

        !read probe positions
        open(1234,file=trim(input_file))
        read(1234,NML=probe_positions)
        close(1234)

        allocate(particles(sum(npps)), stat=rc)
        if(rc.ne.0) write(*,*) " === particle allocation error!"

        call init_particles(particles)

        allocate(n_bins(0:nspecies-1,diag_bins_x, diag_bins_y, diag_bins_z))
        allocate(data_bins(0:nspecies-1,39,diag_bins_x, diag_bins_y, diag_bins_z))
        n_bins = 0
        data_bins = 0.0_8



        !if the following condition is met, particles can leave the system. Diagnostics for collisionality and heating
        !(as they are implemented right now) won't work. So bool_hockney_diag ist set to .false.
        IF ( ANY(boundaries(:)%type == 0) .OR. ANY(boundaries(:)%type == 1) .OR. &
             ANY(boundaries(:)%type == 2) .OR. ANY(boundaries(:)%type == 3) ) THEN
            bool_hockney_diag = .false.
        END IF

        !if the following condition is met, particles enter the system. Diagnostics for collisionality and heating
        !(as they are implemented right now) won't work. So bool_hockney_diag ist set to .false.
        IF ( ANY(species(:)%nfp /= 0) ) bool_hockney_diag = .false.

        IF (bool_hockney_diag) bool_velocity_diag = .true.
        IF (bool_velocity_diag) THEN
            allocate(n_bins_v(0:nspecies-1,3,0:diag_bins_vpar+1))
            allocate(data_bins_v(0:nspecies-1,6,0:diag_bins_vpar+1))
            n_bins_v = 0
            data_bins_v = 0.0_8
        END IF

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
    integer*8 :: npart
    real(KIND=8) :: eps=1.0e-10


    IF (do_resume) THEN
        startstep=resume_step       !resume step is read in read_args (from filename, maybe not the best idea...)

        if (root) write(*,*)"=================================================="
        if (root) write(*,*)"Restarting from checkpoint at timestep",startstep
        if (root) write(*,*)"=================================================="

        call read_particles_mpiio_from_filename(MPI_COMM_WORLD,step,npart,particles,resume_file,noparams=.true.)
        open(fid,file=trim(resume_file)//".h")
        read(fid,NML=pepcf)
        close(fid)
        fsup_at_checkpoint = fsup
    ELSE
        startstep=0
        step=0
    END IF
    last_diag_output = startstep

    call pepc_read_parameters_from_file_name(trim(input_file))

    if(root) write(*,'(a,a)') " == reading parameter file, section pepc-f: ", trim(input_file)
    open(fid,file=trim(input_file))
    read(fid,NML=pepcf)
    read(fid,NML=walk_para_smpss)
    close(fid)


    dx = xmax - xmin
    dy = ymax - ymin
    dz = zmax - zmin

    IF (dx < 0. .or. dy < 0. .or. dz < 0.) THEN  !one of dx, dy, dz < 0 (min > max)
        IF (root) write(*,*) "Geometry not set correctly. dx, dy and dz cannot be < 0"
        STOP
    END IF

    IF (real_equal_zero(dx,eps) .or. real_equal_zero(dy,eps) .or. real_equal_zero(dz,eps)) THEN  !one of dx, dy, dz not set
        IF (root) write(*,*) "Geometry not set correctly. dx, dy and dz cannot be 0"
        STOP
    END IF

    B=sqrt(Bx**2+By**2+Bz**2)

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
      integer(kind_particle) :: global_max_label,local_max_label
      integer :: ip,ib,ispecies,i,j,rc
      real(KIND=8) :: eps=1e-12

      !read probe positions
      open(1234,file=trim(input_file))
      read(1234,NML=probe_positions)
      close(1234)

      npps=0
      tnpps=0
      !hockney diag does not work when resuming (at least not with the current implementation)
      bool_hockney_diag = .false.

      !change charge and mass of particles if fsup was changed
      IF (real_unequal(fsup_at_checkpoint, fsup, eps)) THEN
          IF (root) write(*,*)
          IF (root) write(*,*) "fsup was changed after resume. Existing particles will be adjusted. "
          IF (root) write(*,*) "fsup_old:",fsup_at_checkpoint, "fsup:", fsup
          IF (root) write(*,*)
          DO ip=1,size(particles)
              particles(ip)%data%m = particles(ip)%data%m / fsup_at_checkpoint * fsup
              particles(ip)%data%q = particles(ip)%data%q / fsup_at_checkpoint * fsup
          END DO
      END IF

      !count particles per species (npps) in checkpoint data
      DO ip=1,size(particles)
          npps(particles(ip)%data%species)=npps(particles(ip)%data%species)+1
      END DO
      call MPI_ALLREDUCE(npps, tnpps, nspecies, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

      !Check wallparticles
      IF (tnpps(0)/=count_wallparticles()) THEN
          IF (root) write(*,*) "Number of Wallparticles in checkpoint does not match number in input.",tnpps(0),count_wallparticles()
          IF (root) write(*,*) "You can't change the number of Wallparticles when resuming"
          STOP
      END IF


      ! set charges of wallparticles after resume according to q_tot input
      DO ip=1, size(particles)
          IF (particles(ip)%data%species/=0) CYCLE
          ib = particles(ip)%data%mp_int1
          particles(ip)%data%q = boundaries(ib)%q_tot / boundaries(ib)%nwp
      END DO


      !find global maximum label
      global_max_label = 0
      local_max_label  = 0
      DO ip=1,size(particles)
          if (particles(ip)%label > local_max_label) local_max_label = particles(ip)%label
      END DO
      call MPI_ALLREDUCE(local_max_label, global_max_label, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
      next_label = global_max_label + 1

      !remove existing probes and re-add them (so the number of probes can be changed or new probe species can be added)
      DO ispecies=1,nspecies-1
          call remove_all_probes(ispecies,particles)
      END DO

      DO ispecies=1,nspecies-1
          IF (species(ispecies)%physical_particle == 2) THEN
              IF (root) THEN  !all probes are added on root and will be redistributed later automatically
                  npps(ispecies) = species(ispecies)%nip
                  call reallocate_particles(particles, size(particles), sum(npps))
                  j=0
                  DO i=1, npps(ispecies)
                      ip = size(particles) - npps(ispecies) + i
                      particles(ip)%data%B(1)=Bx
                      particles(ip)%data%B(2)=By
                      particles(ip)%data%B(3)=Bz
                      particles(ip)%label       = next_label
                      next_label = next_label + 1
                      particles(ip)%data%q      = species(ispecies)%q*fsup
                      particles(ip)%data%m      = species(ispecies)%m*fsup
                      particles(ip)%data%species= species(ispecies)%indx
                      particles(ip)%results%e   = 0.0_8
                      particles(ip)%results%pot = 0.0_8
                      particles(ip)%work        = 1.0_8
                      particles(ip)%x(1) = probe_start_x(ispecies) + (j+0.5) * (probe_end_x(ispecies) - probe_start_x(ispecies)) / species(ispecies)%nip
                      particles(ip)%x(2) = probe_start_y(ispecies) + (j+0.5) * (probe_end_y(ispecies) - probe_start_y(ispecies)) / species(ispecies)%nip
                      particles(ip)%x(3) = probe_start_z(ispecies) + (j+0.5) * (probe_end_z(ispecies) - probe_start_z(ispecies)) / species(ispecies)%nip
                      particles(ip)%data%v = 0.
                      j=j+1
                  END DO
              END IF
              tnpps(ispecies) = tnpps(ispecies) + species(ispecies)%nip
          END IF
      END DO

      allocate(n_bins(0:nspecies-1,diag_bins_x, diag_bins_y, diag_bins_z))
      allocate(data_bins(0:nspecies-1,39,diag_bins_x, diag_bins_y, diag_bins_z))
      n_bins = 0
      data_bins = 0.0_8

      call MPI_BCAST(next_label, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, rc)

      !if the following condition is met, particles can leave the system. Diagnostics for collisionality and heating
      !(as they are implemented right now) won't work. So bool_hockney_diag ist set to .false.
      IF ( ANY(boundaries(:)%type == 0) .OR. ANY(boundaries(:)%type == 1) .OR. &
          ANY(boundaries(:)%type == 2) .OR. ANY(boundaries(:)%type == 3) ) THEN
          bool_hockney_diag = .false.
      END IF

      !if the following condition is met, particles enter the system. Diagnostics for collisionality and heating
      !(as they are implemented right now) won't work. So bool_hockney_diag ist set to .false.
      IF ( ANY(species(:)%nfp /= 0) ) bool_hockney_diag = .false.

  end subroutine init_after_resume


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
