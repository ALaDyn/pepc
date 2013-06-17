!!!!!!!!!!!!!!!!!!!!
!! Output module
!!
!! Enthaelt Methoden fuer den Output
!!!!!!!!!!!!!!!!!!!!

MODULE output
    use variables
    use helper


    implicit none

    integer,allocatable      :: thits_out(:,:),treflux_out(:,:)
    real(KIND=8) :: energy_0(2,3)

    CONTAINS


!===============================================================================

    SUBROUTINE init_output_arrays()
        implicit none
        integer rc

        allocate(thits_out(0:nspecies-1,1:nb),stat=rc)
        allocate(treflux_out(0:nspecies-1,1:nb),stat=rc)
        thits_out=0
        treflux_out=0


    END SUBROUTINE init_output_arrays

!===============================================================================

    SUBROUTINE velocity_output(ispecies,filehandle)
        use diagnostics
        implicit none
        include 'mpif.h'

        integer,intent(in)      :: filehandle,ispecies
        real(KIND=8) :: v_mean(3)
        real(KIND=8) :: v_th


        v_mean=get_v_mean(particles,np,ispecies)
        v_th=sqrt(2*species(ispecies)%t_src*e/species(ispecies)%m)
        if(root) write(filehandle,'(a,3(1pe16.7E3))') "Average velocity (1,2,3) [m/s]     : ",v_mean
        if(root) write(filehandle,'(a,3(1pe16.7E3))') "Average velocity/v_th (1,2,3)      : ",v_mean/v_th




    END SUBROUTINE  velocity_output

!===============================================================================

    SUBROUTINE probe_output(ispecies,filehandle)
        use diagnostics
        use module_pepc_types
        implicit none
        include 'mpif.h'

        integer,intent(in)      :: filehandle,ispecies
        integer :: nprobes(0:n_ranks-1),rc, displs(0:n_ranks-1),i
        type(t_particle), allocatable :: probes(:)
        type(t_particle)  :: tprobes(tnpps(ispecies))

        displs=0

        call MPI_GATHER(npps(ispecies), 1, MPI_INTEGER, nprobes, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, rc)

        IF (root) THEN
            displs(0)=0
            DO i=1,n_ranks-1
                displs(i)=displs(i-1)+nprobes(i-1)
            END DO
        END IF

        call get_probe_particles(probes,ispecies)

        call MPI_GATHERV(probes, npps(ispecies), MPI_TYPE_PARTICLE, tprobes, nprobes, displs, MPI_TYPE_PARTICLE, 0, MPI_COMM_WORLD, rc)


        IF (root) THEN
            write(filehandle,'(a10,a10,7(a16))')"","Label","x","y","z","POT","Ex","Ey","Ez"
            DO i=1,tnpps(ispecies)
                write(filehandle,'(a10,i10.10,7(1pe16.7E3))')"Probe:    ",tprobes(i)%label,tprobes(i)%x,tprobes(i)%results%pot*fc,tprobes(i)%results%E*fc
            END DO
        END IF

    END SUBROUTINE  probe_output

!===============================================================================

    SUBROUTINE energy_output(ispecies,filehandle)
        use diagnostics
        implicit none
        include 'mpif.h'

        integer,intent(in)      :: filehandle,ispecies
        real(KIND=8) :: v2_mean(3),ekin(3),epot


        v2_mean=get_v2_mean(particles,np,ispecies)
        ekin=v2_mean*0.5*species(ispecies)%m/e
        if(root) write(filehandle,'(a,3(1pe16.7E3))') "Average kinetic energy (1,2,3) [eV]: ",ekin
        epot=get_epot(particles,np,ispecies)
        if(root) write(filehandle,'(a,(1pe16.7E3))') "Average potential energy [eV]      : ",epot




    END SUBROUTINE  energy_output

!===============================================================================

    SUBROUTINE avg_wallpotential_output(ib,filehandle)
        use diagnostics
        implicit none
        include 'mpif.h'

        integer, intent(in) :: filehandle,ib
        real(KIND=8) :: avg_potential

        avg_potential = get_avg_wallpotential(particles,ib)
        if (root)write(filehandle,'(a,i2,a,(1pe16.7E3))') "Average Potential [V] on boundary", ib,":",avg_potential

    END SUBROUTINE  avg_wallpotential_output

!===============================================================================

    SUBROUTINE timing_output(integrator,particlehandling,pepc_grow,pepc_traverse,pepc_diag,pepc_timber,output,filehandle)

        implicit none

        integer,intent(in)      :: filehandle
        real(kind=8),intent(in) :: integrator,particlehandling,pepc_grow,output,pepc_traverse,pepc_diag,pepc_timber
        real(kind=8)             :: timestep

        timestep=integrator+particlehandling+pepc_grow+output+pepc_traverse+pepc_diag+pepc_timber
        write(filehandle,'(a,es16.8,a,f8.5,a)') " == time in integrator [s], %            :", integrator,", ",100.*integrator/timestep," %"
        write(filehandle,'(a,es16.8,a,f8.5,a)') " == time in particlehandling [s], %      :", particlehandling,", ",100.*particlehandling/timestep," %"
        write(filehandle,'(a,es16.8,a,f8.5,a)') " == time in grow_tree routine [s], %     :", pepc_grow,", ",100.*pepc_grow/timestep," %"
        write(filehandle,'(a,es16.8,a,f8.5,a)') " == time in traverse_tree routine [s], % :", pepc_traverse,", ",100.*pepc_traverse/timestep," %"
        write(filehandle,'(a,es16.8,a,f8.5,a)') " == time in timber_tree routine [s], %   :", pepc_timber,", ",100.*pepc_timber/timestep," %"
        write(filehandle,'(a,es16.8,a,f8.5,a)') " == time in pepc diag routines [s], %    :", pepc_diag,", ",100.*pepc_diag/timestep," %"
        write(filehandle,'(a,es16.8,a,f8.5,a)') " == time in output routines [s], %       :", output,", ",100.*output/timestep," %"
        write(filehandle,'(a,es16.8)') " == total time in timestep [s]           :", timestep
        write(filehandle,*)

    END SUBROUTINE timing_output

!===============================================================================

    SUBROUTINE end_of_ts_output(timestep,filehandle)

        implicit none

        integer,intent(in)      :: timestep,filehandle

        write(filehandle,'(a,es16.8)') " == timestep                             :",dt
        write(filehandle,'(a,i16)') " == finished computing step              :",timestep
        write(filehandle,'(a)') " "
        write(filehandle,'(a)')"############################################################################################################"
        write(filehandle,'(a)')"    ====================================================================================================    "
        write(filehandle,'(a)')"############################################################################################################"
        write(filehandle,'(a)') " "

    END SUBROUTINE end_of_ts_output

!===============================================================================

    SUBROUTINE main_output(filehandle)

        implicit none

        integer,intent(in)      :: filehandle
        integer :: ib,ispecies

        IF(root) write(filehandle,'(a)')"================================================================================================"
        IF(root) write(filehandle,'(a)')"=================================== Info on particle-species ==================================="
        IF(root) write(filehandle,'(a)')"================================================================================================"
        DO ispecies=0,nspecies-1
            IF(root) THEN
                IF (species(ispecies)%physical_particle) THEN
                    write(filehandle,'(a,i2,a)')"----------------------------------- Species ",ispecies," -----------------"
                ELSE
                    write(filehandle,'(a,i2,a)')"----------------------------------- Species ",ispecies," --(unphysical)---"
                END IF
                write(filehandle,'(a,a)')"Name: ",TRIM(species(ispecies)%name)
                write(filehandle,'(a,i10)') "Number of particles:",tnpps(ispecies)
                write(filehandle,*)
            END IF
            IF (species(ispecies)%physical_particle) THEN
                call velocity_output(ispecies,filehandle)
                call energy_output(ispecies,filehandle)
                IF(root) write(filehandle,*)
                DO ib=1,nb
                    IF(root) write(filehandle,'(a,i2,a,i10)') "Hits on boundary ",ib,":",thits_out(ispecies,ib)
                END DO
                IF(root) write(filehandle,*)
                IF(root) write(filehandle,'(a,i10)') "Refluxed particles :",SUM(treflux_out(ispecies,1:nb))+species(ispecies)%nfp
                IF(root) write(filehandle,*)
            ELSE
                IF(diag_interval.ne.0) THEN
                    IF ((MOD(step,diag_interval)==0).or.(step==nt+startstep) .or. (step==1)) THEN
                        IF(species(ispecies)%indx /= 0) call probe_output(ispecies,filehandle)
                    END IF
                END IF
                IF(root) write(filehandle,*)
            END IF

        END DO

        IF(root) write(filehandle,'(a)')"================================================================================================"
        IF(root) write(filehandle,'(a)')"=================================== Info on boundaries ========================================="
        IF(root) write(filehandle,'(a)')"================================================================================================"
        DO ib=1,nb
            IF (boundaries(ib)%nwp/=0) THEN
                IF(root) write(filehandle,'(a,i2,a,es16.8)') "Total charge on boundary ",ib,":",boundaries(ib)%q_tot
                call avg_wallpotential_output(ib,filehandle)
            END IF
        END DO

        IF(root) write(filehandle,*)

        IF(root) write(filehandle,'(a)')"================================================================================================"
        IF(root) write(filehandle,'(a)')"=================================== General info ==============================================="
        IF(root) write(filehandle,'(a)')"================================================================================================"
        IF(root) write(filehandle,'(a,i16)')    " == Total number of particles            :", tnp

        thits_out=0
        treflux_out=0

    END SUBROUTINE main_output

!===============================================================================

    SUBROUTINE set_recycling_output_values(thits,treflux)

        implicit none

        integer,intent(in)      :: thits(0:,:),treflux(0:,:)

        thits_out=thits
        treflux_out=treflux


    END SUBROUTINE set_recycling_output_values

!===============================================================================

  SUBROUTINE write_parameters()
      use module_interaction_specific

      implicit none

      IF (root) THEN
          write(*,'(a,i12)')    " == total number of simulated plasma particles: ", SUM(tnpps(1:nspecies-1))
          write(*,'(a,i12)')    " == number of time steps             : ", nt
          write(*,'(a,i12)')    " == number of species                : ", nspecies
          write(*,'(a,i12)')    " == number of of boundaries          : ", nb
          write(*,'(a,es12.4)') " == time step                        : ", dt
          !write(*,'(a,es12.4)') " == timestep * omega_p               : ", dt*omega_p
          write(*,'(a,es12.4)') " == simulation volume                : ", dx*dy*dz
          write(*,'(a,es12.4)') " == superparticle factor             : ", fsup
          !write(*,'(a,es12.4)') " == particles / debye spehere        : ", 0.5*tnpp/(dx*dy*dz)*l_debye**3
          write(*,'(a,i12)')    " == type of source distribution      : ", quelltyp
          write(*,'(a,l6)')     " == far_field_if_periodic            : ", include_far_field_if_periodic
          write(*,'(a,l6)')     " == do_restore_particles             : ", do_restore_particles
          write(*,*)
          write(*,*) "========== Magnetic Field ========="
          write(*,'(a,f12.4)')  " == Bx                               : ", Bx
          write(*,'(a,f12.4)')  " == By                               : ", By
          write(*,'(a,f12.4)')  " == Bz                               : ", Bz
          write(*,*)
          write(*,*) "========== Simulation Domain ========="
          if (B.ne.0.) then
              write(*,'(a,12X,es10.2)')  " == dx (m)             : ", dx
              write(*,'(a,es10.2,a,es10.2)')  " == dy (gyro_radii, m) : ", dy/r_lamor,", ",dy
              write(*,'(a,es10.2,a,es10.2)')  " == dz (gyro_radii, m) : ", dz/r_lamor,", ",dz
          else
              write(*,'(a,es10.2)')  " == dx (m)             : ", dx
              write(*,'(a,es10.2)')  " == dy (m)             : ", dy
              write(*,'(a,es10.2)')  " == dz (m)             : ", dz
          end if
          write(*,*)
          write(*,*) "========== Plasmaparameters ========="
          write(*,'(a,es12.4)') " == Gyro radius [m]                  : ", r_lamor
          !write(*,'(a,es12.4)') " == Debye length [m]                 : ", l_debye
          !write(*,'(a,es12.4)') " == Plasmafrequency [s^-1]           : ", omega_p
          write(*,*)
          write(*,*) "========== Wall Particles ========="
          write(*,'(a,i12)')    " == total number of wall particles   : ", tnpps(0)
          write(*,*) "========== Random Number Generator ========="
          write(*,'(a,i12)') " == Random Number Generator          : ", rng
      END IF

  END SUBROUTINE write_parameters

!======================================================================================

    subroutine write_particles(p)
        use module_vtk
        implicit none
    
        type(t_particle), allocatable, intent(in) :: p(:)

        integer :: i,n
        type(vtkfile_unstructured_grid) :: vtk
        integer :: vtk_step
        real*8 :: time
        real*8 :: ta, tb
    
        ta = get_time()
    
        time = dt* step

        n=size(p)

        if (step .eq. 0) then
          vtk_step = VTK_STEP_FIRST
        else if (step .eq. nt-1) then
          vtk_step = VTK_STEP_LAST
        else
          vtk_step = VTK_STEP_NORMAL
        endif

        call vtk%create_parallel("particles", step, my_rank, n_ranks, time, vtk_step)
        call vtk%write_headers(n, 0)
        call vtk%startpoints()
        call vtk%write_data_array("xyz", p(:)%x(1), p(:)%x(2), p(:)%x(3))
        call vtk%finishpoints()
        call vtk%startpointdata()
        call vtk%write_data_array("velocity", p(:)%data%v(1), p(:)%data%v(2), p(:)%data%v(3))
        call vtk%write_data_array("el_field", p(:)%results%e(1),p(:)%results%e(2), p(:)%results%e(3))
        call vtk%write_data_array("el_pot", p(:)%results%pot)
        call vtk%write_data_array("charge", p(:)%data%q)
        call vtk%write_data_array("mass", p(:)%data%m)
        call vtk%write_data_array("work", p(:)%work)
        call vtk%write_data_array("pelabel", p(:)%label)
        call vtk%write_data_array("local index", [(i,i=1,n)])
        call vtk%write_data_array("processor", n, my_rank)
        call vtk%write_data_array("species", p(:)%data%species)
        call vtk%finishpointdata()
        call vtk%dont_write_cells()
        call vtk%write_final()
        call vtk%close()

        tb = get_time()

        if(root) write(*,'(a,i6)') " == [write_particles] vtk output at timestep",step

    end subroutine write_particles  

!======================================================================================

    subroutine set_checkpoint()
        use module_checkpoint
        use module_mirror_boxes, only: mirror_box_layers,periodicity

        implicit none
        include 'mpif.h'

        real(KIND=8) :: x0(nb,3)
        real(KIND=8) :: e1(nb,3)
        real(KIND=8) :: e2(nb,3)
        real(KIND=8) :: n(nb,3)
        integer :: type(nb)
        integer :: opposite_bnd(nb)
        logical :: reflux_particles(nb)
        integer :: nwp(nb)
        integer :: ib,nbnd

        integer :: nfp(0:nspecies-1)
        integer :: nip(0:nspecies-1)
        real(KIND=8) :: mass(0:nspecies-1)
        real(KIND=8) :: charge(0:nspecies-1)
        real(KIND=8) :: t_src(0:nspecies-1)
        logical :: physical_particle(0:nspecies-1)
        character(255) :: name(0:nspecies-1)
        integer :: ispecies,ns

        namelist /geometry/ x0,e1,e2,n,type,opposite_bnd,reflux_particles,nwp,nbnd,mirror_box_layers,periodicity
        namelist /species_nml/ ns,nip,nfp,mass,charge,physical_particle,name,t_src

        integer, parameter :: fid = 666

        npart=tnp
        call write_particles_mpiio(MPI_COMM_WORLD,my_rank,step,npart,particles,filename)

        if (root) then
            open(fid,file=trim(filename),STATUS='UNKNOWN', POSITION = 'APPEND')
            write(fid,NML=pepcf)
            write(fid,NML=probe_positions)

            nbnd=nb
            DO ib=1,nbnd
                x0(ib,1:3)=boundaries(ib)%x0
                e1(ib,1:3)=boundaries(ib)%e1
                e2(ib,1:3)=boundaries(ib)%e2
                n(ib,1:3)=boundaries(ib)%n
                type(ib)=boundaries(ib)%type
                opposite_bnd(ib)=boundaries(ib)%opp_bnd
                nwp(ib)=boundaries(ib)%nwp
                reflux_particles(ib)=boundaries(ib)%reflux_particles
            END DO
            write(fid,NML=geometry)


            ns=nspecies
            DO ispecies=0,ns-1
                name(ispecies)=trim(species(ispecies)%name)
                mass(ispecies)=species(ispecies)%m
                charge(ispecies)=species(ispecies)%q
                physical_particle(ispecies)=species(ispecies)%physical_particle
                nfp(ispecies)=species(ispecies)%nfp
                nip(ispecies)=species(ispecies)%nip
                t_src(ispecies)=species(ispecies)%t_src
            END DO
            write(fid,NML=species_nml)
            write(fid,NML=source_nml)
            write(fid,NML=walk_para_smpss)
            close(fid)

        endif

        if(root) write(*,'(a,i6)') " == [set_checkpoint] checkpoint at timestep",step

    end subroutine

!======================================================================================


END MODULE
