
        case(20)  ! disc + tube: inset target
            !    ====================================

            write(ipefile,'(/a/)') "Setting up flange"

            target_geometry=2
            velocity_config=1
            plasma_centre =  (/ xl/2., yl/2., zl/2. /) 
            offset_e = me*nep + ne_rest
            offset_i = ne + me*nip + ni_rest

            ! Electrons 
            call plasma_start( 1, nep, ne, offset_e, target_geometry, velocity_config, idim, &
                -rho0, -1.0, 1.0, vte, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre, &
                number_faces, Vplas, Aplas, Qplas, qe, mass_e, a_ee )
            ! Ions
            call plasma_start( nep+1, nip, ni, offset_i, target_geometry, velocity_config, idim, &
                rho0, 1.0, mass_ratio, vti, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre, &
                number_faces, Vplas, Aplas, Qplas, qi, mass_i, a_ii )

            ! Tube = cylinder with inset

            write(ipefile,'(/a/)') "Setting up tube"

            ! Adjust local numbers if total non-multiple of # PEs
            if (my_rank==0) then
     		np_rest = mod(n_layer(1),n_cpu)
            else
     		np_rest = 0
            endif

            nlayp = n_layer(1)/n_cpu + np_rest  ! total # ions on this CPU
            nep0 = nlayp
            !  Make particle numbers on root known (in case of unequal particle #s - need for label offset)
            call MPI_BCAST( nep0, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
            ne_rest = nep0-nlayp
            ni_rest = nep0-nlayp

            ipstart = nep+nip+1
            label_offset = ne+ni+me*nlayp +n_layer(1) + ni_rest

            ! Place on rear of front disc
            displace = (/ x_plasma/2.+x_layer(1)/2.,0.,0. /)
            layer_geometry = 12  ! tube

            call plasma_start( ipstart, nlayp, n_layer(1), label_offset, layer_geometry, velocity_config, idim, &
                rho_layer(1), 1.0, mratio_layer(1), vti, x_layer(1), y_layer(1), z_layer(1), r_layer(1), plasma_centre+displace, &
                faces(1), V_layer(1), A_layer(1), Q_layer(1), qpart_layer(1), mass_layer(1), ai_layer(1) )

            ! Equal number of neutralising electrons 
            label_offset = ne+ni+me*nlayp + ne_rest 
            call plasma_start( ipstart+nlayp, nlayp, n_layer(1), label_offset, layer_geometry, velocity_config, idim, &
                -rho_layer(1), -1.0, 1.0, vte, x_layer(1), y_layer(1), z_layer(1), r_layer(1), plasma_centre+displace, &
                faces(1), V_layer(1), A_layer(1), Q_layer(1), qpart_layer(1), mass_layer(1), ai_layer(1) )

            npp=npp + 2*nlayp  ! Total # local particles
            ne = ne + n_layer(1)  ! Global # particles
            ni = ni + n_layer(1)  ! Global # particles
            npart = ni+ne

            if (scheme /= 5 .and. ramp) then
                call add_ramp(x_plasma)     ! add exponential ramp to target (stretch container)
            endif


        case(21)  ! disc + tube: tincan target, 50:50 ion, proton mix
            !    ====================================

            if (my_rank==0) write(*,'(/a/)') "Setting up flange"

            target_geometry=2
            velocity_config=1
            plasma_centre =  (/ xl/2.+x_offset, yl/2., zl/2. /) 
            offset_e = me*nep + ne_rest
            offset_i = ne + me*nip + ni_rest

            ! Electrons 
            call plasma_start( 1, nep, ne, offset_e, target_geometry, velocity_config, idim, &
                -rho0, -1.0, 1.0, vte, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre, &
                number_faces, Vplas, Aplas, Qplas, qe, mass_e, a_ee )
            ! Ions
            call plasma_start( nep+1, nip, ni, offset_i, target_geometry, velocity_config, idim, &
                rho0, 1.0, mass_ratio, vti, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre, &
                number_faces, Vplas, Aplas, Qplas, qi, mass_i, a_ii )

!  Turn 1/2 of ions into protons - assumes 50:50 mix @ same density
            do i=nep+1,nep+nip,2
                m(i) = m(i)/mass_ratio*1836
            end do

            if (debug_level.ge.1 .and. me==0) then
                write(*,*) "cap charge ",qi
                write(*,*) "cap mass ",mass_i
                write(*,*) "cap spacing",a_ii
            endif

            if (scheme /= 5 .and. ramp) then
                call add_ramp(x_plasma)     ! add exponential ramp to target (stretch container)
            endif

            ! Tube = cylinder with inset

            if (my_rank==0) write(*,'(/a/)') "Setting up tube"

            ! Adjust local numbers if total non-multiple of # PEs
            if (my_rank==0) then
     		np_rest = mod(n_layer(1),n_cpu)
            else
     		np_rest = 0
            endif

            nlayp = n_layer(1)/n_cpu + np_rest  ! total # ions on this CPU
            nep0 = nlayp
            !  Make particle numbers on root known (in case of unequal particle #s - need for label offset)
            call MPI_BCAST( nep0, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
            ne_rest = nep0-nlayp
            ni_rest = nep0-nlayp

            ipstart = nep+nip+1
            label_offset = ne+ni+me*nlayp +n_layer(1) + ni_rest

            ! Place on rear of front disc
            displace = (/ x_plasma/2.+x_layer(1)/2.,0.,0. /)
            layer_geometry = 12  ! tube

            call plasma_start( ipstart, nlayp, n_layer(1), label_offset, layer_geometry, velocity_config, idim, &
                rho_layer(1), 1.0, mratio_layer(1), vti, x_layer(1), y_layer(1), z_layer(1), r_layer(1), plasma_centre+displace, &
                faces(1), V_layer(1), A_layer(1), Q_layer(1), qpart_layer(1), mass_layer(1), ai_layer(1) )

!  Turn 1/2 of ions into protons - assumes 50:50 mix @ same density
            do i=ipstart,ipstart+nlayp,2
                m(i) = m(i)/mratio_layer(1)*1836
            end do

            ! Equal number of neutralising electrons 
            label_offset = ne+ni+me*nlayp + ne_rest 
            call plasma_start( ipstart+nlayp, nlayp, n_layer(1), label_offset, layer_geometry, velocity_config, idim, &
                -rho_layer(1), -1.0, 1.0, vte, x_layer(1), y_layer(1), z_layer(1), r_layer(1), plasma_centre+displace, &
                faces(1), V_layer(1), A_layer(1), Q_layer(1), qpart_layer(1), mass_layer(1), ai_layer(1) )

            npp=npp + 2*nlayp  ! Total # local particles
            ne = ne + n_layer(1)  ! Global # particles
            ni = ni + n_layer(1)  ! Global # particles
            npart = ni+ne

            if (debug_level.ge.1 .and. me==0) then
                write(*,*) "tube charge ",qpart_layer(1)
                write(*,*) "tube mass ",mass_layer(1)
                write(*,*) "tube spacing",ai_layer(1)
            endif




            !####################################################################################################
        case(12)  ! A.P.L.R's set-up (8th March 2006)
            !=====================================================

            target_geometry=0
            velocity_config=1
            plasma_centre =  (/ xl/2., yl/2., zl/2. /) 
            offset_e = me*nep + ne_rest
            offset_i = ne + me*nip + ni_rest

            ! Electrons 
            call plasma_start( 1, nep, ne, offset_e, target_geometry, velocity_config, idim, &
                -rho0, -1.0, 1.0, vte, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre, &
                number_faces, Vplas, Aplas, Qplas, qe, mass_e, a_ee )
            ! Ions
            call plasma_start( nep+1, nip, ni, offset_i, target_geometry, velocity_config, idim, &
                rho0, 1.0, mass_ratio, vti, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre, &
                number_faces, Vplas, Aplas, Qplas, qi, mass_i, a_ii )

            ! Protons 
            ! Adjust local numbers if total non-multiple of # PEs
            if (my_rank == 0) then
     		!np_rest = mod(n_layer(1),n_cpu)
                np_rest = mod(n_layer(1),n_cpu) !2 layers = twice the no. of protons

            else
     		np_rest = 0
            endif

            !nlayp = n_layer(1)/n_cpu + np_rest  ! total # protons on this CPU
            nlayp = 2*(n_layer(1)/n_cpu + np_rest)  !2 layers = twice the number of protons
            !nlayp = total # of protons in both layers
            !assume that np_rest << n_layer(1)
            nlaypfront = n_layer(1)/n_cpu + np_rest
            nlaypback = n_layer(1)/n_cpu + np_rest
            nep0 = nlayp
            !  Make particle numbers on root known (in case of unequal particle #s - need for label offset)
            call MPI_BCAST( nep0, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
            ne_rest = nep0-nlayp
            ni_rest = nep0-nlayp

            ipstart = nep+nip+1
            ! Place on rear of main slab 
            !-------------------------------------------------
            displace = (/ x_plasma/2.+x_layer(1)/2,0.,0. /)
            label_offset = ne+ni+me*nlaypback+ni_rest 

            ! Equal number of neutralising electrons 
            call plasma_start( ipstart+nlaypback, nlaypback, n_layer(1), label_offset, target_geometry, velocity_config, idim, &
                -rho_layer(1), -1.0, 1.0, vte, x_layer(1), y_layer(1), z_layer(1), r_layer(1), plasma_centre+displace, &
                faces(1), V_layer(1), A_layer(1), Q_layer(1), qpart_layer(1), mass_layer(1), ai_layer(1) )

            !And now for the protons:
            label_offset = ne+ni+n_layer(1) + me*nlaypback + ne_rest 
            call plasma_start( ipstart, nlaypback, n_layer(1), label_offset, target_geometry, velocity_config, idim, &
                rho_layer(1), 1.0, mratio_layer(1), vti, x_layer(1), y_layer(1), z_layer(1), r_layer(1), plasma_centre+displace, &
                faces(1), V_layer(1), A_layer(1), Q_layer(1), qpart_layer(1), mass_layer(1), ai_layer(1) )

            !We also need a layer on the front of the main slab:
            !--------------------------------------------------------
            displace = (/ -x_plasma/2.-x_layer(1)/2,0.,0. /)
            label_offset = ne+ni+2*n_layer(1)+ni_rest+me*nlaypfront
            ! Equal number of neutralising electrons 
            call plasma_start( ipstart+2*nlaypback+nlaypfront, nlaypfront, n_layer(1), label_offset, target_geometry, velocity_config, idim, &
                -rho_layer(1), -1.0, 1.0, vte, x_layer(1), y_layer(1), z_layer(1), r_layer(1), plasma_centre+displace, &
                faces(1), V_layer(1), A_layer(1), Q_layer(1), qpart_layer(1), mass_layer(1), ai_layer(1) )

            !And now for the protons:
            label_offset = ne+ni+3*n_layer(1) +me*nlaypfront + ne_rest 
            call plasma_start( ipstart+2*nlaypback, nlaypfront, n_layer(1), label_offset, target_geometry, velocity_config, idim, &
                rho_layer(1), 1.0, mratio_layer(1), vti, x_layer(1), y_layer(1), z_layer(1), r_layer(1), plasma_centre+displace, &
                faces(1), V_layer(1), A_layer(1), Q_layer(1), qpart_layer(1), mass_layer(1), ai_layer(1) )


            if (debug_level==2 .and. me==0) then
                write(*,*) "proton charge ",qpart_layer(1)
                write(*,*) "proton mass ",mass_layer(1)
                write(*,*) "spacing",ai_layer(1)
            endif


            npp=npp + 2*nlayp  ! Total # local particles
            ne = ne + 2*n_layer(1)  ! Global # particles
            ni = ni + 2*n_layer(1)  ! Global # particles
            npart = ni+ne

            if (scheme /= 5 .and. ramp) then
                call add_ramp(x_plasma)     ! add exponential ramp to target (stretch container)
            endif





        case(13)  ! Add proton layer to slab INSIDE the target
            !    ====================================

            target_geometry=0
            velocity_config=1
            plasma_centre =  (/ xl/2., yl/2., zl/2. /) 
            offset_e = me*nep + ne_rest
            offset_i = ne + me*nip + ni_rest

            ! Electrons 
            call plasma_start( 1, nep, ne, offset_e, target_geometry, velocity_config, idim, &
                -rho0, -1.0, 1.0, vte, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre, &
                number_faces, Vplas, Aplas, Qplas, qe, mass_e, a_ee )
            ! Ions
            call plasma_start( nep+1, nip, ni, offset_i, target_geometry, velocity_config, idim, &
                rho0, 1.0, mass_ratio, vti, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre, &
                number_faces, Vplas, Aplas, Qplas, qi, mass_i, a_ii )

            ! Protons 
            ! Adjust local numbers if total non-multiple of # PEs
            if (my_rank==0) then
     		np_rest = mod(n_layer(1),n_cpu)
            else
     		np_rest = 0
            endif

            nlayp = n_layer(1)/n_cpu + np_rest  ! total # protons on this CPU
            nep0 = nlayp
            !  Make particle numbers on root known (in case of unequal particle #s - need for label offset)
            call MPI_BCAST( nep0, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
            ne_rest = nep0-nlayp
            ni_rest = nep0-nlayp

            ipstart = nep+nip+1
            ! Place on rear of main slab 
            displace = (/ x_plasma/2.-x_layer(1)/2,0.,0. /)
            label_offset = ne+ni+me*nlayp+ni_rest 

            ! Equal number of neutralising electrons 
            call plasma_start( ipstart+nlayp, nlayp, n_layer(1), label_offset, target_geometry, velocity_config, idim, &
                -rho_layer(1), -1.0, 1.0, vte, x_layer(1), y_layer(1), z_layer(1), r_layer(1), plasma_centre+displace, &
                faces(1), V_layer(1), A_layer(1), Q_layer(1), qpart_layer(1), mass_layer(1), ai_layer(1) )

            label_offset = ne+ni+n_layer(1) +me*nlayp + ne_rest 
            call plasma_start( ipstart, nlayp, n_layer(1), label_offset, target_geometry, velocity_config, idim, &
                rho_layer(1), 1.0, mratio_layer(1), vti, x_layer(1), y_layer(1), z_layer(1), r_layer(1), plasma_centre+displace, &
                faces(1), V_layer(1), A_layer(1), Q_layer(1), qpart_layer(1), mass_layer(1), ai_layer(1) )

            if (debug_level==2 .and. me==0) then
                write(*,*) "proton charge ",qpart_layer(1)
                write(*,*) "proton mass ",mass_layer(1)
                write(*,*) "spacing",ai_layer(1)
            endif


            npp=npp + 2*nlayp  ! Total # local particles
            ne = ne + n_layer(1)  ! Global # particles
            ni = ni + n_layer(1)  ! Global # particles
            npart = ni+ne

            if (scheme /= 5 .and. ramp) then
                call add_ramp(x_plasma)     ! add exponential ramp to target (stretch container)
            endif






        case(32)  ! A.P.L.R's SECOND set-up (8th March 2006)
            !=====================================================

            target_geometry=0
            velocity_config=1
            plasma_centre =  (/ xl/2., yl/2., zl/2. /) 
            offset_e = me*nep + ne_rest
            offset_i = ne + me*nip + ni_rest

            ! Electrons 
            call plasma_start( 1, nep, ne, offset_e, target_geometry, velocity_config, idim, &
                -rho0, -1.0, 1.0, vte, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre, &
                number_faces, Vplas, Aplas, Qplas, qe, mass_e, a_ee )
            ! Ions
            call plasma_start( nep+1, nip, ni, offset_i, target_geometry, velocity_config, idim, &
                rho0, 1.0, mass_ratio, vti, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre, &
                number_faces, Vplas, Aplas, Qplas, qi, mass_i, a_ii )

            ! Protons 
            ! Adjust local numbers if total non-multiple of # PEs
            if (my_rank == 0) then
     		!np_rest = mod(n_layer(1),n_cpu)
                np_rest = mod(n_layer(1),n_cpu) !2 layers = twice the no. of protons

            else
     		np_rest = 0
            endif

            !nlayp = n_layer(1)/n_cpu + np_rest  ! total # protons on this CPU
            nlayp = 2*(n_layer(1)/n_cpu + np_rest)  !2 layers = twice the number of protons
            !nlayp = total # of protons in both layers
            !assume that np_rest << n_layer(1)
            nlaypfront = n_layer(1)/n_cpu + np_rest
            nlaypback = n_layer(1)/n_cpu + np_rest
            nep0 = nlayp
            !  Make particle numbers on root known (in case of unequal particle #s - need for label offset)
            call MPI_BCAST( nep0, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
            ne_rest = nep0-nlayp
            ni_rest = nep0-nlayp

            ipstart = nep+nip+1
            ! Place on rear of main slab INSIDE THE TARGET
            !-------------------------------------------------
            displace = (/ x_plasma/2.-x_layer(1)/2,0.,0. /)
            label_offset = ne+ni+me*nlaypback+ni_rest 

            ! Equal number of neutralising electrons 
            call plasma_start( ipstart+nlaypback, nlaypback, n_layer(1), label_offset, target_geometry, velocity_config, idim, &
                -rho_layer(1), -1.0, 1.0, vte, x_layer(1), y_layer(1), z_layer(1), r_layer(1), plasma_centre+displace, &
                faces(1), V_layer(1), A_layer(1), Q_layer(1), qpart_layer(1), mass_layer(1), ai_layer(1) )

            !And now for the protons:
            label_offset = ne+ni+n_layer(1) + me*nlaypback + ne_rest 
            call plasma_start( ipstart, nlaypback, n_layer(1), label_offset, target_geometry, velocity_config, idim, &
                rho_layer(1), 1.0, mratio_layer(1), vti, x_layer(1), y_layer(1), z_layer(1), r_layer(1), plasma_centre+displace, &
                faces(1), V_layer(1), A_layer(1), Q_layer(1), qpart_layer(1), mass_layer(1), ai_layer(1) )

            !We also need a layer on the front of the main slab:
            !INSIDE THE TARGET
            !--------------------------------------------------------
            displace = (/ -x_plasma/2.+x_layer(1)/2,0.,0. /)
            label_offset = ne+ni+2*n_layer(1)+ni_rest+me*nlaypfront
            ! Equal number of neutralising electrons 
            call plasma_start( ipstart+2*nlaypback+nlaypfront, nlaypfront, n_layer(1), label_offset, target_geometry, velocity_config, idim, &
                -rho_layer(1), -1.0, 1.0, vte, x_layer(1), y_layer(1), z_layer(1), r_layer(1), plasma_centre+displace, &
                faces(1), V_layer(1), A_layer(1), Q_layer(1), qpart_layer(1), mass_layer(1), ai_layer(1) )

            !And now for the protons:
            label_offset = ne+ni+3*n_layer(1) +me*nlaypfront + ne_rest 
            call plasma_start( ipstart+2*nlaypback, nlaypfront, n_layer(1), label_offset, target_geometry, velocity_config, idim, &
                rho_layer(1), 1.0, mratio_layer(1), vti, x_layer(1), y_layer(1), z_layer(1), r_layer(1), plasma_centre+displace, &
                faces(1), V_layer(1), A_layer(1), Q_layer(1), qpart_layer(1), mass_layer(1), ai_layer(1) )


            if (debug_level==2 .and. me==0) then
                write(*,*) "proton charge ",qpart_layer(1)
                write(*,*) "proton mass ",mass_layer(1)
                write(*,*) "spacing",ai_layer(1)
            endif


            npp=npp + 2*nlayp  ! Total # local particles
            ne = ne + 2*n_layer(1)  ! Global # particles
            ni = ni + 2*n_layer(1)  ! Global # particles
            npart = ni+ne

            if (scheme /= 5 .and. ramp) then
                call add_ramp(x_plasma)     ! add exponential ramp to target (stretch container)
            endif




        case(35)  ! A.P.L.R's third set-up (19th May 2006)
            !=====================================================

            target_geometry=0
            velocity_config=1
            plasma_centre =  (/ xl/2., yl/2., zl/2. /) 
            offset_e = me*nep + ne_rest
            offset_i = ne + me*nip + ni_rest

            ! Electrons 
            call plasma_start( 1, nep, ne, offset_e, target_geometry, velocity_config, idim, &
                -rho0, -1.0, 1.0, vte, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre, &
                number_faces, Vplas, Aplas, Qplas, qe, mass_e, a_ee )
            ! Ions
            call plasma_start( nep+1, nip, ni, offset_i, target_geometry, velocity_config, idim, &
                rho0, 1.0, mass_ratio, vti, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre, &
                number_faces, Vplas, Aplas, Qplas, qi, mass_i, a_ii )

            ! Protons 
            ! Adjust local numbers if total non-multiple of # PEs
            if (my_rank == 0) then
     		!np_rest = mod(n_layer(1),n_cpu)
                np_rest = mod(n_layer(1),n_cpu) !2 layers = twice the no. of protons

            else
     		np_rest = 0
            endif

            !nlayp = n_layer(1)/n_cpu + np_rest  ! total # protons on this CPU
            nlayp = 2*(n_layer(1)/n_cpu + np_rest)  !2 layers = twice the number of protons
            !nlayp = total # of protons in both layers
            !assume that np_rest << n_layer(1)
            nlaypfront = n_layer(1)/n_cpu + np_rest
            nlaypback = n_layer(1)/n_cpu + np_rest
            nep0 = nlayp
            !  Make particle numbers on root known (in case of unequal particle #s - need for label offset)
            call MPI_BCAST( nep0, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
            ne_rest = nep0-nlayp
            ni_rest = nep0-nlayp

            ipstart = nep+nip+1
            ! Place on rear of main slab 
            !-------------------------------------------------
            displace = (/ x_plasma/2.+x_layer(1)/2,0.,0. /)
            label_offset = ne+ni+me*nlaypback+ni_rest 

            ! Equal number of neutralising electrons 
            call plasma_start( ipstart+nlaypback, nlaypback, n_layer(1), label_offset, target_geometry, velocity_config, idim, &
                -rho_layer(1), -1.0, 1.0, vte, x_layer(1), y_layer(1), z_layer(1), r_layer(1), plasma_centre+displace, &
                faces(1), V_layer(1), A_layer(1), Q_layer(1), qpart_layer(1), mass_layer(1), ai_layer(1) )

            !And now for the protons:
            label_offset = ne+ni+n_layer(1) + me*nlaypback + ne_rest 
            call plasma_start( ipstart, nlaypback, n_layer(1), label_offset, target_geometry, velocity_config, idim, &
                rho_layer(1), 1.0, mratio_layer(1), vti, x_layer(1), y_layer(1), z_layer(1), r_layer(1), plasma_centre+displace, &
                faces(1), V_layer(1), A_layer(1), Q_layer(1), qpart_layer(1), mass_layer(1), ai_layer(1) )

            !We also need a another layer in the same place
            !--------------------------------------------------------
            displace = (/ x_plasma/2.+x_layer(2)/2,0.,0. /)
            label_offset = ne+ni+2*n_layer(2)+ni_rest+me*nlaypfront
            ! Equal number of neutralising electrons 
            call plasma_start( ipstart+2*nlaypback+nlaypfront, nlaypfront, n_layer(2), label_offset, target_geometry, velocity_config, idim, &
                -rho_layer(2), -1.0, 1.0, vte, x_layer(2), y_layer(2), z_layer(2), r_layer(2), plasma_centre+displace, &
                faces(2), V_layer(2), A_layer(2), Q_layer(2), qpart_layer(2), mass_layer(2), ai_layer(2) )

            !And now for the protons:
            label_offset = ne+ni+3*n_layer(2) +me*nlaypfront + ne_rest 
            call plasma_start( ipstart+2*nlaypback, nlaypfront, n_layer(2), label_offset, target_geometry, velocity_config, idim, &
                rho_layer(2), 1.0, mratio_layer(2), vti, x_layer(2), y_layer(2), z_layer(2), r_layer(2), plasma_centre+displace, &
                faces(2), V_layer(2), A_layer(2), Q_layer(2), qpart_layer(2), mass_layer(2), ai_layer(2) )


            if (debug_level==2 .and. me==0) then
                write(*,*) "proton charge ",qpart_layer(1)
                write(*,*) "proton mass ",mass_layer(1)
                write(*,*) "spacing",ai_layer(1)
            endif


            npp=npp + 2*nlayp  ! Total # local particles
            ne = ne + 2*n_layer(1)  ! Global # particles
            ni = ni + 2*n_layer(1)  ! Global # particles
            npart = ni+ne

            if (scheme /= 5 .and. ramp) then
                call add_ramp(x_plasma)     ! add exponential ramp to target (stretch container)
            endif


            !###########################################################################################################



        case(15)  ! 2-layer slab with foam
            !    ==================================

            target_geometry=0
            velocity_config=1
            plasma_centre =  (/ xl/2., yl/2., zl/2. /) 
            offset_e = me*nep + ne_rest
            offset_i = ne + me*nip + ni_rest

            !  Main block
            write(ipefile,'(/a/)') "Setting up main slab"

            ! Electrons 
            call plasma_start( 1, nep, ne, offset_e, target_geometry, velocity_config, idim, &
                -rho0, -1.0, 1.0, vte, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre, &
                number_faces, Vplas, Aplas, Qplas, qe, mass_e, a_ee )
            ! Ions
            call plasma_start( nep+1, nip, ni, offset_i, target_geometry, velocity_config, idim, &
                rho0, 1.0, mass_ratio, vti, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre, &
                number_faces, Vplas, Aplas, Qplas, qi, mass_i, a_ii )

            ! Protons
            write(ipefile,'(/a/)') "Setting up proton layer"

            ! Adjust local numbers if total non-multiple of # PEs
            if (my_rank==0) then
     		np_rest = mod(n_layer(1),n_cpu)
            else
     		np_rest = 0
            endif

            nlayp = n_layer(1)/n_cpu + np_rest  ! total # protons on this CPU

            ipstart = nep+nip+1  ! local index
            label_offset = ne+ni+me*nlayp ! global label

            ! Place on front of main slab 
            !	displace = plasma_centre - (/ x_plasma/2.+x_layer(1)/2,0.,0. /)

            call plasma_start( ipstart, nlayp, n_layer(1), label_offset, target_geometry, velocity_config, idim, &
                rho_layer(1), 1.0, mratio_layer(1), vti, x_layer(1), y_layer(1), z_layer(1), r_layer(1), plasma_centre+displace, &
                faces(1), V_layer(1), A_layer(1), Q_layer(1), qpart_layer(1), mass_layer(1), a_layer(1) )

            ! Equal number of neutralising electrons 
            label_offset = ne+ni+n_layer(1) +me*nlayp 
            call plasma_start( ipstart+nlayp, nlayp, n_layer(1), label_offset, target_geometry, velocity_config, idim, &
                -rho_layer(1), -1.0, 1.0, vte, x_layer(1), y_layer(1), z_layer(1), r_layer(1), plasma_centre+displace, &
                faces(1), V_layer(1), A_layer(1), Q_layer(1), qpart_layer(1), mass_layer(1), ai_layer(1) )

            npp=npp + 2*nlayp  ! Total # local particles
            ne = ne + n_layer(1)  ! Global # particles
            ni = ni + n_layer(1)  ! Global # particles
            npart = ni+ne

            if (scheme /= 5 .and. ramp) then
                call add_ramp(x_plasma+x_layer(1))     ! add exponential ramp to target (stretch container)
            endif


            !  Add foam layer (2) to rear side

            write(ipefile,'(/a/)') "Setting up foam layer"

            target_geometry=7
            velocity_config=1   ! Ions, electrons thermal
            plasma_centre =  (/ xl/2., yl/2., zl/2. /) ! Centre of plasma
            plasma_centre=plasma_centre+displace+r_layer(2) ! centre of 1st shell 


            nshell_z=1
            nshell_y=1
            nshell_x=2
            nshell = nshell_x * nshell_y * nshell_z
            ne_shell = n_layer(2)/nshell  ! # particles allocated to foam array must be exact multiple
            ni_shell = n_layer(2)/nshell  ! of # electrons or ions per cell

            ! Adjust local numbers if total non-multiple of # PEs
            if (my_rank==0) then
     		ne_rest = mod(ne_shell,n_cpu)
     		ni_rest = mod(ni_shell,n_cpu)
            else
     		ne_rest = 0
     		ni_rest = 0
            endif
            nep_shell = ne_shell/n_cpu + ne_rest  ! # electrons per shell allocated to this CPU
            nip_shell = ni_shell/n_cpu + ni_rest 
            nlayp = (nep_shell+nip_shell)*nshell  ! total # particles in foam on this CPU 

            nep0 = nep_shell
            nip0 = nip_shell
            !  Make particle numbers on root known (in case of unequal particle #s - need for label offset)
            call MPI_BCAST( nep0, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
            call MPI_BCAST( nip0, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
            ne_rest = nep0-nep_shell
            ni_rest = nip0-nip_shell

            ipstart_e = npp + 1  ! local index, electrons
            ipstart_i = ipstart_e + nep_shell*nshell  ! local index, ions
            offset_e = npart + me*ne_shell + ne_rest  ! Global label 
            ! TODO: - may need to change this later to have sequentially numbered species
            offset_i = npart + n_layer(2) + me*ni_shell + ni_rest

            do ishell = 0,nshell_z-1
                do jshell = 0,nshell_y-1
                    do kshell = 0,nshell_x-1
                        displace = (/ 2*ishell*r_sphere, 2*jshell*r_sphere, 2*kshell*r_sphere /)

                        ! Electrons 

                        call plasma_start( ipstart_e, nep_shell, ne_shell, offset_e, target_geometry, velocity_config, idim, &
                            -rho_layer(2), -1.0, 1.0, vte, x_layer(2), y_layer(2), z_layer(2), r_layer(2), plasma_centre+displace, &
                            faces(2), V_layer(2), A_layer(2), Q_layer(2), qpart_layer(2), mass_layer(2), a_layer(2) )

                        ipstart_e = ipstart_e + nep_shell  ! index start
                        offset_e = offset_e + nep_shell  ! label offset

                        ! Ions
                        call plasma_start( ipstart_i, nip_shell, ni_shell, offset_i, target_geometry, 0, idim, &
                            rho_layer(2), 1.0, mratio_layer(2), vti, x_layer(2), y_layer(2), z_layer(2), r_layer(2), plasma_centre+displace, &
                            faces(2), V_layer(2), A_layer(2), Q_layer(2), qpart_layer(2), mass_layer(2), a_layer(2) )

                        ipstart_i = ipstart_i + nip_shell  ! index start
                        offset_i = offset_i + nip_shell  ! label offset

                    end do
                end do
                ! shift for FCC
                if (mod(ishell,2)==0) then
                    plasma_centre = plasma_centre +  r_layer(2)/sqrt(2.)*(/ -1., 1., 1. /)
                else
                    plasma_centre = plasma_centre +  r_layer(2)/sqrt(2.)*(/ -1., -1., -1. /)
                endif
            end do


            npp=npp + nlayp  ! Total # local particles
            ne = ne + n_layer(2)  ! Global # particles
            ni = ni + n_layer(2)  ! Global # particles
            npart = ni+ne
            if (debug_level==2) then
                write(*,*)  "# particles increased to (cpu, global):",npp,npart
            endif
