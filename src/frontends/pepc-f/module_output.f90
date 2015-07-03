!!!!!!!!!!!!!!!!!!!!
!! Output module
!!
!! Enthaelt Methoden fuer den Output
!!!!!!!!!!!!!!!!!!!!

MODULE output
    use module_pepc_kinds
    use variables
    use helper


    implicit none

    integer,allocatable      :: thits_out(:,:),treflux_out(:,:),trethermalized_out(:)
    real(KIND=8) :: energy_0(2,3)

    CONTAINS


!===============================================================================

    SUBROUTINE init_output_arrays()
        implicit none
        integer rc

        if (.not. allocated(thits_out)) allocate(thits_out(0:nspecies-1,1:nb),stat=rc)
        if (.not. allocated(treflux_out))allocate(treflux_out(0:nspecies-1,1:nb),stat=rc)
        thits_out=0
        treflux_out=0


    END SUBROUTINE init_output_arrays

!===============================================================================

    SUBROUTINE check_output_intervals()
        implicit none

        IF ( (checkp_interval /=0) .AND. ( ( MOD(step,checkp_interval) == 0 ) .OR. ( step == nt+startstep ) ) ) THEN
            checkpoint_now = .TRUE.
        ELSE
            checkpoint_now = .FALSE.
        END IF

        IF ( (npy_interval /= 0) .AND. ( ( MOD(step,npy_interval) == 0 ) .OR. ( step == nt+startstep ) ) ) THEN
            npy_now = .TRUE.
        ELSE
            npy_now = .FALSE.
        END IF

        IF ( (vtk_interval /= 0) .AND. ( ( MOD(step,vtk_interval) == 0 ) .OR. ( step == nt+startstep ) ) ) THEN
            vtk_now = .TRUE.
        ELSE
            vtk_now = .FALSE.
        END IF

        IF ( (diag_interval /= 0) .AND. ( ( MOD(step,diag_interval) == 0 ) .OR. ( step == nt+startstep ) .OR. ( step == 1 ) ) ) THEN
            diag_now = .TRUE.
        ELSE
            diag_now = .FALSE.
        END IF

    END SUBROUTINE check_output_intervals

!===============================================================================

    SUBROUTINE velocity_output(ispecies,filehandle)
        use diagnostics
        implicit none
        include 'mpif.h'

        integer,intent(in)      :: filehandle,ispecies
        real(KIND=8) :: v_mean(3)
        real(KIND=8) :: v_th


        v_mean=get_v_mean(particles,ispecies)
        v_th=sqrt(2*species(ispecies)%src_t*e/species(ispecies)%m)
        if(root) write(filehandle,'(a,3(1pe16.7E3))') "Average velocity (1,2,3) [m/s]     : ",v_mean
        if(root) write(filehandle,'(a,3(1pe16.7E3))') "Average velocity/v_th (1,2,3)      : ",v_mean/v_th


    END SUBROUTINE  velocity_output

    !===============================================================================

    SUBROUTINE plasma_props_output(ispecies,filehandle,nbs,dbs,write_data)
        use diagnostics
        use module_pepc_types
        implicit none
        include 'mpif.h'

        integer                    :: rc,ix,iy,iz

        integer,intent(in)         :: filehandle,ispecies
        integer                    :: npoints
        real(KIND=8),intent(inout) :: dbs(:,:,:,:)
        integer,intent(inout)      :: nbs(:,:,:)
        logical,intent(in)         :: write_data
        real(KIND=8)               :: tdata_bins(39,diag_bins_x,diag_bins_y,diag_bins_z)
        integer                    :: tn_bins(diag_bins_x,diag_bins_y,diag_bins_z)
        real(KIND=8)               :: tn_bins_dble(diag_bins_x,diag_bins_y,diag_bins_z)


        npoints = diag_bins_x * diag_bins_y * diag_bins_z

        IF (bool_diag_bins_cylinder) THEN
            call fill_data_bins_cylindrical(ispecies,nbs,dbs)
        ELSE
            call fill_data_bins(ispecies,nbs,dbs)
        END IF

        IF (.not. write_data) THEN
            RETURN
        ELSE
            call MPI_ALLREDUCE(nbs, tn_bins, npoints, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, rc)
            call MPI_ALLREDUCE(dbs, tdata_bins, 39*npoints, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, rc)

            IF (bool_avg_btwn_diag_steps) THEN
                tdata_bins = tdata_bins / (step - last_diag_output)
                tn_bins_dble = DBLE(tn_bins) / (step - last_diag_output)
            ELSE
                tn_bins_dble = DBLE(tn_bins)
            END IF

            IF (root) THEN
                IF (bool_diag_bins_cylinder) THEN
                    write(filehandle,'(a,3(i6.5))')"Average values for particles at equidistant points (nx,nr,ntheta):", diag_bins_x,diag_bins_y,diag_bins_z
                    write(filehandle,'(a12,3(a7),a12,39(a16))')"","ix","ir","itheta","n","vx","vy","vz","vpar", &
                                                               "vperp1","vperp2","vxvx","vyvy","vzvz","vxvy","vyvz", &
                                                               "vzvx","vparvpar","vperp1vperp1","vperp2vperp2", &
                                                               "vparvperp1","vperp1vperp2","vperp2vpar", &
                                                               "phi","Ex","Ey","Ez","Epar","Eperp1","Eperp2", &
                                                               "phiphi","ExEx","EyEy","EzEz","ExEy","EyEz","EzEx", &
                                                               "EparEpar","Eperp1Eperp1","Eperp2Eperp","EparEperp1", &
                                                               "Eperp1Eperp2","Eperp2Epar","age"
                ELSE
                    write(filehandle,'(a,3(i6.5))')"Average values for particles at equidistant points (nx,ny,nz):", diag_bins_x,diag_bins_y,diag_bins_z
                    write(filehandle,'(a12,3(a7),a12,39(a16))')"","ix","iy","iz","n","vx","vy","vz","vpar", &
                                                               "vperp1","vperp2","vxvx","vyvy","vzvz","vxvy","vyvz", &
                                                               "vzvx","vparvpar","vperp1vperp1","vperp2vperp2", &
                                                               "vparvperp1","vperp1vperp2","vperp2vpar", &
                                                               "phi","Ex","Ey","Ez","Epar","Eperp1","Eperp2", &
                                                               "phiphi","ExEx","EyEy","EzEz","ExEy","EyEz","EzEx", &
                                                               "EparEpar","Eperp1Eperp1","Eperp2Eperp","EparEperp1", &
                                                               "Eperp1Eperp2","Eperp2Epar","age"
                END IF
                DO iz=1,diag_bins_z
                    DO iy=1,diag_bins_y
                        DO ix=1,diag_bins_x
                            write(filehandle,'(a12,3(i7.5),F14.3,39(1pe16.7E3))')"Bins:       ",ix,iy,iz,tn_bins_dble(ix,iy,iz), &
                                                        tdata_bins(:,ix,iy,iz)
                        END DO
                    END DO
                END DO
            END IF
        END IF

    END SUBROUTINE  plasma_props_output

!===============================================================================

    SUBROUTINE write_velocity_bins(ispecies,filehandle,nbs,dbs,write_data)
        use diagnostics, only: fill_velocity_bins
        use module_pepc_types
        implicit none
        include 'mpif.h'

        integer                    :: rc,ivx,iv2

        integer,intent(in)         :: filehandle,ispecies
        integer                    :: npoints
        real(KIND=8),intent(inout) :: dbs(:,0:,:)
        integer,intent(inout)      :: nbs(0:,:)
        logical,intent(in)         :: write_data
        real(KIND=8)               :: tdata_bins(3,0:diag_bins_vx+1,diag_bins_v2+1)
        integer                    :: tn_bins(0:diag_bins_vx+1,diag_bins_v2+1)
        real(KIND=8)               :: tn_bins_dble(0:diag_bins_vx+1,diag_bins_v2+1)
        real(KIND=8)               :: vxmin, vxmax, v2min, v2max


        npoints = (diag_bins_vx+2) * (diag_bins_v2+1)

        call fill_velocity_bins(ispecies,nbs,dbs)

        IF (.not. write_data) THEN
            RETURN
        ELSE
            call MPI_ALLREDUCE(nbs, tn_bins, npoints, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, rc)
            call MPI_ALLREDUCE(dbs, tdata_bins, 3*npoints, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, rc)

            IF (bool_avg_btwn_diag_steps) THEN
                tdata_bins = tdata_bins / (step - last_diag_output)
                tn_bins_dble = DBLE(tn_bins) / (step - last_diag_output)
            ELSE
                tn_bins_dble = DBLE(tn_bins)
            END IF

            IF (root) THEN
                write(filehandle,*)
                write(filehandle,'(a,2(i6.5))')"collision analysis (nvx,nv2):", diag_bins_vx,diag_bins_v2
                write(filehandle,'(a12,2(a7),4(a16),a12,3(a16))')"","ivx","iv2","vx_min","vx_max","v2_min","v2_max","N","dvx","dvx2","dv2"
                DO ivx=0,diag_bins_vx+1
                    vxmin = (-v_grid_max + ((ivx-1) * (v_grid_max*2)/diag_bins_vx))! * species(ispecies)%v_th*sqrt(2.)
                    vxmax = (-v_grid_max + (ivx * (v_grid_max*2)/diag_bins_vx))! * species(ispecies)%v_th*sqrt(2.)
                    if (ivx == 0) vxmin = -huge(vxmin)
                    if (ivx == diag_bins_vx+1) vxmax = huge(vxmax)
                    DO iv2=1,diag_bins_v2+1
                        v2min = ((iv2-1) * (v_grid_max*3)/diag_bins_v2)! * species(ispecies)%v_th*sqrt(2.)
                        v2max = (iv2 * (v_grid_max*3)/diag_bins_v2)! * species(ispecies)%v_th*sqrt(2.)
                        if (iv2 == diag_bins_v2+1) v2max = huge(v2max)
                        write(filehandle,'(a12,2(i7.5),4(1pe16.7E3),0pF14.3,3(1pe16.7E3))')"Bins_v:     ",ivx,iv2,vxmin,vxmax,v2min,v2max,&
                                                                                          tn_bins_dble(ivx,iv2), tdata_bins(:,ivx,iv2)
                    END DO
                END DO
            END IF
        END IF

    END SUBROUTINE  write_velocity_bins
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


        v2_mean=get_v2_mean(particles,ispecies)
        ekin=v2_mean*0.5*species(ispecies)%m/e
        if(root) write(filehandle,'(a,3(1pe16.7E3))') "Average kinetic energy (1,2,3) [eV]: ",ekin
        epot=get_epot(particles,ispecies)
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

    SUBROUTINE timing_output(integrator,particlehandling,pepc_grow,pepc_traverse,pepc_diag,pepc_timber,boundary_field,output,filehandle)

        implicit none

        integer,intent(in)      :: filehandle
        real(kind=8),intent(in) :: integrator,particlehandling,pepc_grow,output,pepc_traverse,&
                                   pepc_diag,pepc_timber,boundary_field
        real(kind=8)            :: timestep

        timestep=integrator+particlehandling+pepc_grow+output+pepc_traverse+pepc_diag+pepc_timber+boundary_field
        write(filehandle,'(a,es16.8,a,f8.5,a)') " == time in integrator [s], %            :", integrator,", ",100.*integrator/timestep," %"
        write(filehandle,'(a,es16.8,a,f8.5,a)') " == time in particlehandling [s], %      :", particlehandling,", ",100.*particlehandling/timestep," %"
        write(filehandle,'(a,es16.8,a,f8.5,a)') " == time in grow_tree routine [s], %     :", pepc_grow,", ",100.*pepc_grow/timestep," %"
        write(filehandle,'(a,es16.8,a,f8.5,a)') " == time in traverse_tree routine [s], % :", pepc_traverse,", ",100.*pepc_traverse/timestep," %"
        write(filehandle,'(a,es16.8,a,f8.5,a)') " == time in timber_tree routine [s], %   :", pepc_timber,", ",100.*pepc_timber/timestep," %"
        write(filehandle,'(a,es16.8,a,f8.5,a)') " == time in pepc diag routines [s], %    :", pepc_diag,", ",100.*pepc_diag/timestep," %"
        write(filehandle,'(a,es16.8,a,f8.5,a)') " == time in bnd field routines [s], %    :", boundary_field,", ",100.*boundary_field/timestep," %"
        write(filehandle,'(a,es16.8,a,f8.5,a)') " == time in output routines [s], %       :", output,", ",100.*output/timestep," %"
        write(filehandle,'(a,es16.8)') " == total time in timestep [s]           :", timestep
        write(filehandle,*)

    END SUBROUTINE timing_output

!===============================================================================

    SUBROUTINE energy_resolved_hits_output(ispecies)
        use module_utils
        implicit none
        include 'mpif.h'

        integer, intent(in) :: ispecies
        integer :: rc,ibins,tmp_filehandle=2345
        real(KIND=8) :: binwidth,emin,emax
        integer :: hits(nb, nbins_energy_resolved_hits+1)
        character(100) :: tmp_file,format,dir

        hits=0

        call MPI_ALLREDUCE(energy_resolved_hits(ispecies,:,:), hits, (nbins_energy_resolved_hits+1)*nb, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, rc)

        dir = "./energy_resolved_hits/"
        write(tmp_file,'(a,"erh_species_",i3.3,".dat")') trim(dir), ispecies
        write(format,'(a,i3,a)') "(2es13.5,",  nb  ,"i9)"

        binwidth = ehit_max(ispecies)/nbins_energy_resolved_hits

        IF(root) THEN
            IF (last_diag_output == startstep) call create_directory(trim(dir))
            open(unit=tmp_filehandle,file=trim(tmp_file),status='UNKNOWN',position='APPEND')
            write(tmp_filehandle,'(a,i6,a,i6,a)')"---------------------- TIMESTEPS: ",last_diag_output+1," - ",step," -----------------"
            DO ibins=1,nbins_energy_resolved_hits+1
                emin = (ibins-1)*binwidth
                emax = ibins*binwidth
                IF (ibins == nbins_energy_resolved_hits+1) emax = 1.e32
                write(tmp_filehandle,format) emin,emax,hits(:,ibins)
            END DO
            write(tmp_filehandle,*)""
            write(tmp_filehandle,*)"############################################################################################################"
            write(tmp_filehandle,*)"    ====================================================================================================    "
            write(tmp_filehandle,*)"############################################################################################################"
            write(tmp_filehandle,*)""
            close(tmp_filehandle)
        END IF

    END SUBROUTINE energy_resolved_hits_output

!===============================================================================


    SUBROUTINE age_resolved_hits_output(ispecies)
        use module_utils
        implicit none
        include 'mpif.h'

        integer, intent(in) :: ispecies
        integer :: rc,ibins,tmp_filehandle=2345
        real(KIND=8) :: binwidth,agemin,agemax
        integer :: hits(nb, nbins_age_resolved_hits+1)
        character(100) :: tmp_file,format,dir

        hits=0

        call MPI_ALLREDUCE(age_resolved_hits(ispecies,:,:), hits, (nbins_age_resolved_hits+1)*nb, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, rc)

        dir = "./age_resolved_hits/"
        write(tmp_file,'(a,"agerh_species_",i3.3,".dat")') trim(dir), ispecies
        write(format,'(a,i3,a)') "(2es13.5,",  nb  ,"i9)"

        binwidth = agehit_max(ispecies)/nbins_age_resolved_hits

        IF(root) THEN
            IF (last_diag_output == startstep) call create_directory(trim(dir))
            open(unit=tmp_filehandle,file=trim(tmp_file),status='UNKNOWN',position='APPEND')
            write(tmp_filehandle,'(a,i6,a,i6,a)')"---------------------- TIMESTEPS: ",last_diag_output+1," - ",step," -----------------"
            DO ibins=1,nbins_age_resolved_hits+1
                agemin = (ibins-1)*binwidth
                agemax = ibins*binwidth
                IF (ibins == nbins_age_resolved_hits+1) agemax = 1.e32
                write(tmp_filehandle,format) agemin,agemax,hits(:,ibins)
            END DO
            write(tmp_filehandle,*)""
            write(tmp_filehandle,*)"############################################################################################################"
            write(tmp_filehandle,*)"    ====================================================================================================    "
            write(tmp_filehandle,*)"############################################################################################################"
            write(tmp_filehandle,*)""
            close(tmp_filehandle)
        END IF

    END SUBROUTINE age_resolved_hits_output


!===============================================================================
    SUBROUTINE angle_resolved_hits_output(ispecies)
        use module_utils
        implicit none
        include 'mpif.h'

        integer, intent(in) :: ispecies
        integer :: rc,ibins,tmp_filehandle=2345
        real(KIND=8) :: binwidth,betamin,betamax
        integer :: hits(nb, nbins_angle_resolved_hits)
        character(100) :: tmp_file,format,dir

        hits=0

        call MPI_ALLREDUCE(angle_resolved_hits(ispecies,:,:), hits, (nbins_angle_resolved_hits)*nb, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, rc)

        dir = "./angle_resolved_hits/"
        write(tmp_file,'(a,"arh_species_",i3.3,".dat")') trim(dir), ispecies
        write(format,'(a,i3,a)') "(2es13.5,",  nb  ,"i9)"

        binwidth = 90./nbins_angle_resolved_hits

        IF(root) THEN
            IF (last_diag_output == startstep) call create_directory(trim(dir))
            open(unit=tmp_filehandle,file=trim(tmp_file),status='UNKNOWN',position='APPEND')
            write(tmp_filehandle,'(a,i6,a,i6,a)')"---------------------- TIMESTEPS: ",last_diag_output+1," - ",step," -----------------"
            DO ibins=1,nbins_angle_resolved_hits
                betamin = (ibins-1)*binwidth
                betamax = ibins*binwidth
                write(tmp_filehandle,format) betamin,betamax,hits(:,ibins)
            END DO
            write(tmp_filehandle,*)""
            write(tmp_filehandle,*)"############################################################################################################"
            write(tmp_filehandle,*)"    ====================================================================================================    "
            write(tmp_filehandle,*)"############################################################################################################"
            write(tmp_filehandle,*)""
            close(tmp_filehandle)
        END IF

    END SUBROUTINE angle_resolved_hits_output

!===============================================================================

    SUBROUTINE space_resolved_hits_output(ispecies)
        use module_utils
        implicit none
        include 'mpif.h'

        integer, intent(in) :: ispecies
        integer :: rc,ibins,ibnd,tmp_filehandle=2345
        integer :: hits(nb, nbins_e1_space_resolved_hits, nbins_e2_space_resolved_hits)
        character(100) :: tmp_file,format,dir

        hits=0

        call MPI_ALLREDUCE(space_resolved_hits(ispecies,:,:,:), hits, nbins_e1_space_resolved_hits*nbins_e2_space_resolved_hits*nb, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, rc)

        dir = "./space_resolved_hits/"
        write(format,'(a,i3,a)') "(",nbins_e2_space_resolved_hits  ,"i9)"


        IF(root) THEN
            IF ((last_diag_output == startstep) ) call create_directory(trim(dir))
            DO ibnd=1,nb
                IF (boundaries(ibnd)%type > 0) CYCLE
                write(tmp_file,'(a,"srh_species_",i3.3,"_bnd_",i3.3,".dat")') trim(dir), ispecies, ibnd
                open(unit=tmp_filehandle,file=trim(tmp_file),status='UNKNOWN',position='APPEND')
                write(tmp_filehandle,'(a,i6,a,i6,a)')"---------------------- TIMESTEPS: ",last_diag_output+1," - ",step," -----------------"
                DO ibins=1,nbins_e1_space_resolved_hits
                    write(tmp_filehandle,format) hits(ibnd,ibins,:)
                END DO
                write(tmp_filehandle,*)""
                write(tmp_filehandle,*)"############################################################################################################"
                write(tmp_filehandle,*)"    ====================================================================================================    "
                write(tmp_filehandle,*)"############################################################################################################"
                write(tmp_filehandle,*)""
                close(tmp_filehandle)
            END DO
        END IF

    END SUBROUTINE space_resolved_hits_output

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
        use diagnostics, only: hockney_diag
        implicit none

        integer,intent(in)      :: filehandle
        integer :: ib,ispecies
        real(KIND=8) :: avg_1(nspecies-1), avg_2(nspecies-1), avg_3(nspecies-1), avg_4(nspecies-1), avg_5(nspecies-1)
        real(KIND=8) :: avg_6(nspecies-1), avg_7(nspecies-1), avg_8(nspecies-1), avg_9(nspecies-1), avg_10(nspecies-1)
        real(KIND=8) :: avg_11(nspecies-1), avg_12(nspecies-1), avg_13(nspecies-1), avg_14(nspecies-1), avg_15(nspecies-1)
        real(KIND=8) :: avg_16(nspecies-1),avg_17(nspecies-1),avg_18(nspecies-1),avg_19(nspecies-1),avg_20(nspecies-1)
        real(KIND=8) :: avg_21(nspecies-1), avg_fields(12, nspecies-1)

        IF(root) write(filehandle,'(a)')"================================================================================================"
        IF(root) write(filehandle,'(a)')"=================================== Info on particle-species ==================================="
        IF(root) write(filehandle,'(a)')"================================================================================================"
        IF ((step >= hockney_start_step) .AND. (bool_hockney_diag)) THEN
            call hockney_diag(particles, avg_1, avg_2, avg_3, avg_4, avg_5, avg_6, avg_7, avg_8, &
                                         avg_9, avg_10, avg_11, avg_12, avg_13, avg_14, avg_15, avg_16, &
                                         avg_17, avg_18, avg_19, avg_20, avg_21, avg_fields)
        END IF
        DO ispecies=0,nspecies-1
            IF(root) THEN
                IF (species(ispecies)%physical_particle == 1) THEN
                    write(filehandle,'(a,i2,a)')"----------------------------------- Species ",ispecies," -----------------"
                ELSE IF (species(ispecies)%physical_particle == 0) THEN
                    write(filehandle,'(a,i2,a)')"-------------(wallparticle)-------- Species ",ispecies," --(unphysical)---"
                ELSE IF (species(ispecies)%physical_particle == 2) THEN
                    write(filehandle,'(a,i2,a)')"-------------(probe particle)------ Species ",ispecies," --(unphysical)---"
                ELSE IF (species(ispecies)%physical_particle == 3) THEN
                    write(filehandle,'(a,i2,a)')"-------------(test particle)------- Species ",ispecies," --(unphysical)---"
                END IF
                write(filehandle,'(a,a)')"Name: ",TRIM(species(ispecies)%name)
                write(filehandle,'(a,i10)') "Number of particles:",tnpps(ispecies)
                write(filehandle,*)
            END IF
            IF (species(ispecies)%moving_particle) THEN
                call velocity_output(ispecies,filehandle)
                call energy_output(ispecies,filehandle)
                IF ((step >= hockney_start_step) .AND. (bool_hockney_diag)) THEN
                    IF(root) write(filehandle,*)
                    IF(root) write(filehandle,'(a,(1pe16.7E3))') "Hockney <beta(t)>: ",avg_1(ispecies)
                    IF(root) write(filehandle,'(a,(1pe16.7E3))') "Hockney <|beta(t)|>: ",avg_2(ispecies)
                    IF(root) write(filehandle,'(a,(1pe16.7E3))') "Hockney <beta(t)^2>^0.5: ",avg_3(ispecies)
                    IF(root) write(filehandle,'(a,(1pe16.7E3))') "Hockney <h(t)>: ",avg_4(ispecies)
                    IF(root) write(filehandle,'(a,(1pe16.7E3))') "Hockney <|h(t)|>: ",avg_5(ispecies)
                    IF(root) write(filehandle,'(a,(1pe16.7E3))') "Hockney <h(t)^2>^0.5: ",avg_6(ispecies)
                    IF(root) write(filehandle,'(a,(1pe16.7E3))') "Hockney <vpar(t)>: ",avg_7(ispecies)
                    IF(root) write(filehandle,'(a,(1pe16.7E3))') "Hockney <|vpar(t)|>: ",avg_8(ispecies)
                    IF(root) write(filehandle,'(a,(1pe16.7E3))') "Hockney <vpar(t)^2>^0.5: ",avg_9(ispecies)
                    IF(root) write(filehandle,'(a,(1pe16.7E3))') "Hockney <|vperp(t)|>: ",avg_10(ispecies)
                    IF(root) write(filehandle,'(a,(1pe16.7E3))') "Hockney <vperp(t)^2>^0.5: ",avg_11(ispecies)
                    IF(root) write(filehandle,'(a,(1pe16.7E3))') "Hockney <|v(t)|>: ",avg_12(ispecies)
                    IF(root) write(filehandle,'(a,(1pe16.7E3))') "Hockney <v(t)^2>^0.5: ",avg_13(ispecies)
                    IF(root) write(filehandle,'(a,(1pe16.7E3))') "Hockney 1/2 * m * <vpar(t)^2>: ",avg_14(ispecies)
                    IF(root) write(filehandle,'(a,(1pe16.7E3))') "Hockney 1/2 * m * <vperp(t)^2>: ",avg_15(ispecies)
                    IF(root) write(filehandle,'(a,(1pe16.7E3))') "Hockney 1/2 * m * <v^2>: ",avg_16(ispecies)
                    IF(root) write(filehandle,'(a,(1pe16.7E3))') "Hockney <Ex>: ",avg_fields(1,ispecies)
                    IF(root) write(filehandle,'(a,(1pe16.7E3))') "Hockney <Ey>: ",avg_fields(2,ispecies)
                    IF(root) write(filehandle,'(a,(1pe16.7E3))') "Hockney <Ez>: ",avg_fields(3,ispecies)
                    IF(root) write(filehandle,'(a,(1pe16.7E3))') "Hockney <ExEx>: ",avg_fields(4,ispecies)
                    IF(root) write(filehandle,'(a,(1pe16.7E3))') "Hockney <EyEy>: ",avg_fields(5,ispecies)
                    IF(root) write(filehandle,'(a,(1pe16.7E3))') "Hockney <EzEz>: ",avg_fields(6,ispecies)
                    IF(root) write(filehandle,'(a,(1pe16.7E3))') "Hockney <ExEy>: ",avg_fields(7,ispecies)
                    IF(root) write(filehandle,'(a,(1pe16.7E3))') "Hockney <EyEz>: ",avg_fields(8,ispecies)
                    IF(root) write(filehandle,'(a,(1pe16.7E3))') "Hockney <EzEx>: ",avg_fields(9,ispecies)
                    IF(root) write(filehandle,'(a,(1pe16.7E3))') "Hockney <|E|>: ",avg_fields(10,ispecies)
                    IF(root) write(filehandle,'(a,(1pe16.7E3))') "Hockney <E^2>: ",sum(avg_fields(4:6,ispecies))
                    IF(root) write(filehandle,'(a,(1pe16.7E3))') "Hockney <pot>: ",avg_fields(11,ispecies)
                    IF(root) write(filehandle,'(a,(1pe16.7E3))') "Hockney <pot^2>: ",avg_fields(12,ispecies)
                    IF(root) write(filehandle,'(a,(1pe16.7E3))') "Hockney <|v(t)|^3>: ",avg_17(ispecies)
                    IF(root) write(filehandle,'(a,(1pe16.7E3))') "Hockney <|v(t)|^4>: ",avg_18(ispecies)
                    IF(root) write(filehandle,'(a,(1pe16.7E3))') "Hockney <|v(t)|^5>: ",avg_19(ispecies)
                    IF(root) write(filehandle,'(a,(1pe16.7E3))') "Hockney <|v(t)|^6>: ",avg_20(ispecies)
                    IF(root) write(filehandle,'(a,(1pe16.7E3))') "Hockney <|v(t)|^2 |v(0)|^2>: ",avg_21(ispecies)
                END IF
                IF(root) write(filehandle,*)

                IF (diag_now) THEN
                    call plasma_props_output(ispecies,filehandle,n_bins(ispecies,:,:,:), &
                                             data_bins(ispecies,:,:,:,:),.True.)
                    IF (bool_velocity_diag) THEN
                        call write_velocity_bins(ispecies,filehandle,n_bins_v(ispecies,:,:), &
                                                 data_bins_v(ispecies,:,:,:),.True.)
                    END IF
                    IF (bool_energy_resolved_hits) call energy_resolved_hits_output(ispecies)
                    IF (bool_angle_resolved_hits) call angle_resolved_hits_output(ispecies)
                    IF (bool_space_resolved_hits) call space_resolved_hits_output(ispecies)
                    IF (bool_age_resolved_hits) call age_resolved_hits_output(ispecies)
                ELSE
                    IF (bool_avg_btwn_diag_steps) THEN
                        call plasma_props_output(ispecies,filehandle,n_bins(ispecies,:,:,:),  &
                                                 data_bins(ispecies,:,:,:,:),.False.)
                        IF (bool_velocity_diag) THEN
                            call write_velocity_bins(ispecies,filehandle,n_bins_v(ispecies,:,:), &
                                                     data_bins_v(ispecies,:,:,:),.False.)
                        END IF
                    END IF
                END IF

                IF(root) write(filehandle,*)
                DO ib=1,nb
                    IF(root) write(filehandle,'(a,i2,a,i10)') "Hits on boundary ",ib,":",thits_out(ispecies,ib)
                END DO
                IF(root) write(filehandle,*)
                IF(root) write(filehandle,'(a,i10)') "Refluxed particles :",SUM(treflux_out(ispecies,1:nb))+species(ispecies)%nfp
                IF(root) write(filehandle,'(a,i10)') "Rethermalized particles :",trethermalized_out(ispecies)
                IF(root) write(filehandle,*)
            ELSE
                IF (diag_now) THEN
                    IF(species(ispecies)%indx /= 0) call probe_output(ispecies,filehandle)
                END IF
                IF(root) write(filehandle,*)
            END IF

        END DO

        IF (diag_now) THEN
            last_diag_output = step
            IF (bool_energy_resolved_hits) energy_resolved_hits = 0
            IF (bool_angle_resolved_hits) angle_resolved_hits = 0
            IF (bool_space_resolved_hits) space_resolved_hits = 0
            IF (bool_age_resolved_hits) age_resolved_hits = 0
            n_bins = 0
            data_bins = 0.0_8
            IF (bool_velocity_diag) n_bins_v = 0
            IF (bool_velocity_diag) data_bins_v = 0.0_8
        END IF

        IF(root) write(filehandle,'(a)')"================================================================================================"
        IF(root) write(filehandle,'(a)')"=================================== Info on boundaries ========================================="
        IF(root) write(filehandle,'(a)')"================================================================================================"
        DO ib=1,nb
            IF ((boundaries(ib)%type == 0) .OR. (boundaries(ib)%type == 1)) THEN
                IF(root) write(filehandle,'(a,i2,a,es16.8)') "Total charge on boundary ",ib,":",boundaries(ib)%q_tot
                IF (boundaries(ib)%type == 0) call avg_wallpotential_output(ib,filehandle)
            END IF
        END DO

        IF(root) write(filehandle,*)

        IF(root) write(filehandle,'(a)')"================================================================================================"
        IF(root) write(filehandle,'(a)')"=================================== General info ==============================================="
        IF(root) write(filehandle,'(a)')"================================================================================================"
        IF(root) write(filehandle,'(a,i16)')    " == Total number of particles            :", sum(tnpps)

        thits_out = 0
        treflux_out = 0

    END SUBROUTINE main_output

!===============================================================================

    SUBROUTINE set_recycling_output_values(thits,treflux,trethermalized)

        implicit none
        integer :: rc
        integer,intent(in) :: thits(0:,:),treflux(0:,:),trethermalized(0:)

        if (.not. allocated(thits_out)) allocate(thits_out(0:nspecies-1,1:nb),stat=rc)
        if (.not. allocated(treflux_out)) allocate(treflux_out(0:nspecies-1,1:nb),stat=rc)
        if (.not. allocated(trethermalized_out)) allocate(trethermalized_out(0:nspecies-1),stat=rc)

        thits_out(:,:) = thits
        treflux_out(:,:) = treflux
        trethermalized_out(:) = trethermalized(:)


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
          write(*,'(a,es12.4)') " == simulation volume                : ", dx*dy*dz
          write(*,'(a,es12.4)') " == superparticle factor             : ", fsup
          write(*,'(a,l6)')     " == far_field_if_periodic            : ", include_far_field_if_periodic
          write(*,*)
          write(*,*) "========== Magnetic Field ========="
          write(*,'(a,f12.4)')  " == Bx                               : ", Bx
          write(*,'(a,f12.4)')  " == By                               : ", By
          write(*,'(a,f12.4)')  " == Bz                               : ", Bz
          write(*,*)
          write(*,*) "========== Simulation Domain ========="
          write(*,'(a,es10.2)')  " == dx (m)             : ", dx
          write(*,'(a,es10.2)')  " == dy (m)             : ", dy
          write(*,'(a,es10.2)')  " == dz (m)             : ", dz
          write(*,*)
          write(*,*) "========== Wall Particles ========="
          write(*,'(a,i12)')    " == total number of wall particles   : ", tnpps(0)
          write(*,*) "========== Random Number Generator ========="
          write(*,'(a,i12)') " == Random Number Generator          : ", rng
      END IF

  END SUBROUTINE write_parameters

  !======================================================================================

    subroutine write_particles_npy(p, my_rank, itime)
        use module_pepc_types
        use module_utils
        use module_debug, only : pepc_status

        implicit none
        integer(kind_pe), intent(in) :: my_rank
        integer(kind_default), intent(in) :: itime
        type(t_particle), allocatable, intent(in) :: p(:)

        character(12), parameter :: directory = './particles'
        integer :: filehandle = 91

        character(255) :: filename
        logical :: firstcall  = .true.
        character(50) :: dir, format
        integer :: n, i

        integer(1) :: MAGIC1 = -109 !eqv to X'93'
        character(5) :: MAGIC2 = "NUMPY"
        integer(1) :: MAGIC3 = X'1'
        integer(1) :: MAGIC4 = X'0'
        integer(1) :: HEAD_LEN0
        integer(1) :: HEAD_LEN1
        integer(1) :: SPACE = X'20'
        integer(1) :: NEWLINE = X'0A'
        character(512) :: HEADER
        integer :: header_length, total_length_is, total_length_needed, spaces_needed

        dir = trim(directory)//"/npy/"
        write(filename,'(a,"particle_",i12.12,"_",i6.6,".npy")') trim(dir), itime, my_rank
        call pepc_status("DUMP PARTICLES NUMPY: "//trim(filename))

        if (firstcall) then
            call create_directory(trim(directory))
            call create_directory(trim(dir))
            firstcall = .false.
        endif

        if(root) write(*,'(a,i6)') " == [npy particle output] npy output at timestep",step

        open(filehandle, file=trim(filename), STATUS='UNKNOWN', ACCESS="STREAM", POSITION="APPEND", CONVERT="BIG_ENDIAN")

        n = size(p)

        !create npy file header
        !https://github.com/numpy/numpy/blob/master/doc/neps/npy-format.rst for more info on the file format
        write(format,'(a,i2,a)') "(3a,i",9,",a)"
        write(HEADER,format) "{'descr': [('label', '>i8'), ('x', '>f8'), ('y', '>f8'), ('z', '>f8'), ('vx', '>f8'), ('vy', '>f8'), ", &
                           "('vz', '>f8'), ('q', '>f8'), ('m', '>f8'), ('age', '>f8'), ('species', '>i4'), ('mp_int1', '>i4'), ", &
                           "('Ex', '>f8'), ('Ey', '>f8'), ('Ez', '>f8'), ('phi', '>f8')], 'fortran_order': False, 'shape': (", &
                           n, &
                           ",), }"

        header_length = len_trim(HEADER)
        total_length_is = header_length + 6 + 4
        total_length_needed = ((total_length_is / 16) + 1) * 16
        spaces_needed = total_length_needed - total_length_is - 1

        HEAD_LEN1 = INT((total_length_needed - 10) / 256,kind=1)
        HEAD_LEN0 = INT(MOD((total_length_needed - 10), 256),kind=1)

        write(filehandle) MAGIC1,MAGIC2,MAGIC3,MAGIC4,HEAD_LEN0, HEAD_LEN1
        write(filehandle) trim(HEADER)
        DO i=1, spaces_needed
            write(filehandle) SPACE
        END DO
        write(filehandle) NEWLINE


        !write data
        DO i=1, n
            write(filehandle) p(i)%label, p(i)%x, p(i)%data%v, p(i)%data%q, p(i)%data%m, p(i)%data%age, &
                              p(i)%data%species, p(i)%data%mp_int1, p(i)%results%e, p(i)%results%pot
        END DO
        close(filehandle)

    end subroutine write_particles_npy


!======================================================================================

    subroutine write_particles_vtk(p, vtk_mask)
        use module_vtk
        implicit none
    
        type(t_particle), allocatable, intent(in) :: p(:)
        integer(kind_particle), intent(in) :: vtk_mask

        type(t_particle), allocatable :: p_out(:)
        integer :: n_out
        integer :: i,n,rc,j
        type(vtkfile_unstructured_grid) :: vtk
        integer :: vtk_step
        real*8 :: time
        real*8 :: ta, tb
    
        ta = get_time()
        time = dt* step

        if(root) write(*,'(a,i6)') " == [write_particles] vtk output at timestep",step

        n=size(p)
        allocate(p_out(1:n),stat=rc)
        n_out = 0
        j=1
        do i=1,n
            if (MOD(p(i)%label,vtk_mask)==0) then
                n_out = n_out+1
                p_out(j) = p(i)
                j = j+1
            end if
        end do

        call reallocate_particles(p_out, n_out, n_out)

        if (step .eq. 0) then
          vtk_step = VTK_STEP_FIRST
        else if (step .eq. nt-1) then
          vtk_step = VTK_STEP_LAST
        else
          vtk_step = VTK_STEP_NORMAL
        endif

        call vtk%create_parallel("particles", step, my_rank, n_ranks, time, vtk_step)
        call vtk%write_headers(n_out, 0)
        call vtk%startpoints()
        call vtk%write_data_array("xyz", p_out(:)%x(1), p_out(:)%x(2), p_out(:)%x(3))
        call vtk%finishpoints()
        call vtk%startpointdata()
        call vtk%write_data_array("velocity", p_out(:)%data%v(1), p_out(:)%data%v(2), p_out(:)%data%v(3))
        call vtk%write_data_array("el_field", p_out(:)%results%e(1),p_out(:)%results%e(2), p_out(:)%results%e(3))
        call vtk%write_data_array("el_pot", p_out(:)%results%pot)
        call vtk%write_data_array("charge", p_out(:)%data%q)
        call vtk%write_data_array("mass", p_out(:)%data%m)
        call vtk%write_data_array("work", p_out(:)%work)
        call vtk%write_data_array("pelabel", p_out(:)%label)
        call vtk%write_data_array("local index", [(i,i=1,n_out)])
        call vtk%write_data_array("processor", n_out, my_rank)
        call vtk%write_data_array("species", p_out(:)%data%species)
        call vtk%write_data_array("mp_int1", p_out(:)%data%mp_int1)
        call vtk%write_data_array("age", p_out(:)%data%age)
        call vtk%finishpointdata()
        call vtk%dont_write_cells()
        call vtk%write_final()
        call vtk%close()

        tb = get_time()


    end subroutine write_particles_vtk

!======================================================================================

    subroutine set_checkpoint()
        use module_checkpoint

        implicit none
        include 'mpif.h'

        integer :: ispecies,ns
        integer*8 :: npart
        integer, parameter :: fid = 666

        integer :: ib,nbnd
        real(KIND=8), allocatable :: x0(:,:)
        real(KIND=8), allocatable :: e1(:,:)
        real(KIND=8), allocatable :: e2(:,:)
        real(KIND=8), allocatable :: n(:,:)
        integer, allocatable :: type(:)
        integer, allocatable :: opposite_bnd(:)
        logical, allocatable :: reflux_particles(:)
        integer, allocatable :: nwp(:)
        real(KIND=8), allocatable :: q_tot(:)

        integer, allocatable :: nfp(:)
        integer, allocatable :: nip(:)
        real(KIND=8), allocatable :: mass(:)
        real(KIND=8), allocatable :: charge(:)
        real(KIND=8), allocatable :: src_t(:)
        integer, allocatable :: physical_particle(:)
        character(255), allocatable :: name(:)
        real(KIND=8), allocatable :: src_x0(:,:)
        real(KIND=8), allocatable :: src_e1(:,:)
        real(KIND=8), allocatable :: src_e2(:,:)
        real(KIND=8), allocatable :: src_e3(:,:)
        integer, allocatable :: src_type_x(:), src_type_v(:)
        integer, allocatable :: src_bnd(:)


        namelist /geometry/ x0,e1,e2,n,type,opposite_bnd,reflux_particles,nwp,nbnd,q_tot
        namelist /species_nml/ ns,nip,nfp,mass,charge,physical_particle,name,src_t,src_x0,src_e1,src_e2,src_e3,src_bnd,src_type_x, src_type_v


        if(root) write(*,'(a,i6)') " == [set_checkpoint] checkpoint at timestep",step

        npart=sum(tnpps)
        call write_particles_mpiio(MPI_COMM_WORLD,step,npart,particles,filename)

        if (root) then
            allocate(x0(nb,3))
            allocate(e1(nb,3))
            allocate(e2(nb,3))
            allocate(n(nb,3))
            allocate(type(nb))
            allocate(opposite_bnd(nb))
            allocate(reflux_particles(nb))
            allocate(nwp(nb))
            allocate(q_tot(nb))

            allocate(nfp(0:nspecies-1))
            allocate(nip(0:nspecies-1))
            allocate(mass(0:nspecies-1))
            allocate(charge(0:nspecies-1))
            allocate(src_t(0:nspecies-1))
            allocate(physical_particle(0:nspecies-1))
            allocate(name(0:nspecies-1))
            allocate(src_x0(0:nspecies-1,3))
            allocate(src_e1(0:nspecies-1,3))
            allocate(src_e2(0:nspecies-1,3))
            allocate(src_e3(0:nspecies-1,3))
            allocate(src_type_x(0:nspecies-1))
            allocate(src_type_v(0:nspecies-1))
            allocate(src_bnd(0:nspecies-1))

            open(fid,file=trim(filename),STATUS='UNKNOWN', POSITION = 'APPEND')
            write(fid,NML=pepcf,DELIM="QUOTE")
            write(fid,NML=probe_positions,DELIM="QUOTE")

            nbnd=nb
            DO ib=1,nbnd
                x0(ib,1:3)=boundaries(ib)%x0
                e1(ib,1:3)=boundaries(ib)%e1
                e2(ib,1:3)=boundaries(ib)%e2
                n(ib,1:3)=boundaries(ib)%n
                type(ib)=boundaries(ib)%type
                opposite_bnd(ib)=boundaries(ib)%opp_bnd
                nwp(ib)=boundaries(ib)%nwp
                q_tot(ib)=boundaries(ib)%q_tot
                reflux_particles(ib)=boundaries(ib)%reflux_particles
            END DO
            write(fid,NML=geometry,DELIM="QUOTE")

            ns=nspecies
            DO ispecies=0,ns-1
                name(ispecies)=trim(species(ispecies)%name)
                mass(ispecies)=species(ispecies)%m
                charge(ispecies)=species(ispecies)%q
                physical_particle(ispecies)=species(ispecies)%physical_particle
                nfp(ispecies)=species(ispecies)%nfp
                nip(ispecies)=species(ispecies)%nip
                src_t(ispecies)=species(ispecies)%src_t
                src_x0(ispecies,1:3)=species(ispecies)%src_x0
                src_e1(ispecies,1:3)=species(ispecies)%src_e1
                src_e2(ispecies,1:3)=species(ispecies)%src_e2
                src_e3(ispecies,1:3)=species(ispecies)%src_e3
                src_type_x(ispecies)=species(ispecies)%src_type_x
                src_type_v(ispecies)=species(ispecies)%src_type_v
                src_bnd(ispecies)=species(ispecies)%src_bnd
            END DO
            write(fid,NML=species_nml,DELIM="QUOTE")
            write(fid,NML=walk_para_smpss,DELIM="QUOTE")
            close(fid)

            deallocate(x0)
            deallocate(e1)
            deallocate(e2)
            deallocate(n)
            deallocate(type)
            deallocate(opposite_bnd)
            deallocate(reflux_particles)
            deallocate(nwp)
            deallocate(q_tot)

            deallocate(nfp)
            deallocate(nip)
            deallocate(mass)
            deallocate(charge)
            deallocate(src_t)
            deallocate(physical_particle)
            deallocate(name)
            deallocate(src_x0)
            deallocate(src_e1)
            deallocate(src_e2)
            deallocate(src_e3)
            deallocate(src_type_x)
            deallocate(src_type_v)
            deallocate(src_bnd)

        endif


    end subroutine

!======================================================================================


END MODULE
