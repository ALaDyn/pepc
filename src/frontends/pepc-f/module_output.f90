!!!!!!!!!!!!!!!!!!!!
!! Output module
!!
!! Enthaelt Methoden fuer den Output
!!!!!!!!!!!!!!!!!!!!

MODULE output
    use variables
    use helper


    implicit none

    CONTAINS

!===============================================================================

    SUBROUTINE timing_output(integrator,particlehandling,pepc_grow,pepc_traverse,pepc_rest,output,filehandle)

        implicit none

        integer,intent(in)      :: filehandle
        real(kind=8),intent(in) :: integrator,particlehandling,pepc_grow,output,pepc_traverse,pepc_rest
        real(kind=8)             :: timestep

        timestep=integrator+particlehandling+pepc_grow+output+pepc_traverse+pepc_rest
        if(root) write(filehandle,'(a,es16.8,a,f8.5,a)') " == time in integrator [s], %            : ", integrator,", ",100.*integrator/timestep," %"
        if(root) write(filehandle,'(a,es16.8,a,f8.5,a)') " == time in particlehandling [s], %      : ", particlehandling,", ",100.*particlehandling/timestep," %"
        if(root) write(filehandle,'(a,es16.8,a,f8.5,a)') " == time in grow_tree routine [s], %     : ", pepc_grow,", ",100.*pepc_grow/timestep," %"
        if(root) write(filehandle,'(a,es16.8,a,f8.5,a)') " == time in traverse_tree routine [s], % : ", pepc_traverse,", ",100.*pepc_traverse/timestep," %"
        if(root) write(filehandle,'(a,es16.8,a,f8.5,a)') " == time in other pepc routines [s], %   : ", pepc_rest,", ",100.*pepc_rest/timestep," %"
        if(root) write(filehandle,'(a,es16.8,a,f8.5,a)') " == time in output routines [s], %       : ", output,", ",100.*output/timestep," %"
        if(root) write(filehandle,'(a,es16.8)') " == total time in timestep [s]           : ", timestep

    END SUBROUTINE timing_output

!===============================================================================

    SUBROUTINE end_of_ts_output(timestep,filehandle)

        implicit none

        integer,intent(in)      :: timestep,filehandle

        if(root) write(filehandle,'(a,i6)') " == finished computing step              : ",timestep
        if(root) write(filehandle,'(a)') " "
        if(root) write(filehandle,'(a)') "#############################################################################"
        if(root) write(filehandle,'(a)') " "

    END SUBROUTINE end_of_ts_output

!===============================================================================

    SUBROUTINE recycling_output(thits,treflux,filehandle)

        implicit none

        integer,intent(in)      :: thits(0:,:),treflux(0:,:),filehandle
        integer :: ib,ispecies

        IF(root) THEN
            write(filehandle,'(a)')"========================= Info on particle-species ========================="
            DO ispecies=0,nspecies-1
                write(filehandle,'(a,i2,a)')"========================= Species ",ispecies," ========================="
                write(filehandle,'(a,a)')"Name: ",TRIM(species(ispecies)%name)
                IF (.not. species(ispecies)%physical_particle) write(filehandle,'(a)')"(no physical species)"
                write(filehandle,*)
                IF (species(ispecies)%physical_particle) THEN
                    DO ib=1,nb
                        write(filehandle,'(a,i2,a,i10)') "Hits on boundary ",ib," :",thits(ispecies,ib)
                    END DO
                    write(filehandle,*)
                    write(filehandle,'(a,i10)') "Refluxed particles  :",SUM(treflux(ispecies,1:nb))+species(ispecies)%nfp
                    write(filehandle,*)
                    write(filehandle,'(a,i10)') "Number of particles  :",tnpps(ispecies)
                    write(filehandle,*)
                ELSE
                    write(filehandle,*)
                    write(filehandle,'(a,i10)') "Number of particles  :",tnpps(ispecies)
                    write(filehandle,*)
                END IF
                write(filehandle,'(a)')"=============================================================="
                write(filehandle,*)
            END DO

            write(filehandle,'(a)')"========================= Info on boundaries ========================="
            DO ib=1,nb
                IF (boundaries(ib)%nwp/=0) THEN
                    write(filehandle,'(a,i2,a,es16.8)') "Total charge on boundary ",ib," :",boundaries(ib)%q_tot
                END IF
            END DO

            write(filehandle,*)
            write(filehandle,'(a,i16)')    " == Total number of particles            : ", tnp
        END IF

    END SUBROUTINE recycling_output

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
          write(*,'(a,i12)')    " == mirror_layers                    : ", mirror_layers
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
          write(*,'(a,f12.4)')  " == TE [eV] (sourcetemperature)      : ", te_ev
          write(*,'(a,f12.4)')  " == TI [eV] (sourcetemperature)      : ", ti_ev
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

        integer :: i
        type(vtkfile_unstructured_grid) :: vtk
        integer :: vtk_step
        real*8 :: time
        real*8 :: ta, tb
    
        ta = get_time()
    
        time = dt* step

        if (step .eq. 0) then
          vtk_step = VTK_STEP_FIRST
        else if (step .eq. nt-1) then
          vtk_step = VTK_STEP_LAST
        else
          vtk_step = VTK_STEP_NORMAL
        endif

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
        call vtk%write_data_array("species", np, p(:)%data%species)
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
        logical :: physical_particle(0:nspecies-1)
        character(255) :: name(0:nspecies-1)
        integer :: ispecies,ns

        namelist /geometry/ x0,e1,e2,n,type,opposite_bnd,reflux_particles,nwp,nbnd
        namelist /species_nml/ ns,nip,nfp,mass,charge,physical_particle,name

        integer, parameter :: fid = 666

        npart=tnp
        !call write_particles_ascii(my_rank,step,np,particles,filename)
        call write_particles_mpiio(MPI_COMM_WORLD,my_rank,step,np,npart,particles,filename)

        if (root) then
            open(fid,file=trim(filename),STATUS='UNKNOWN', POSITION = 'APPEND')
            write(fid,NML=pepcf)

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
            END DO
            write(fid,NML=species_nml)

            write(fid,NML=walk_para_smpss)
            close(fid)

        endif

        if(root) write(*,'(a,i6)') " == [set_checkpoint] checkpoint at timestep",step

    end subroutine

!======================================================================================


END MODULE
