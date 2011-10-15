! =============================================
!
!                CONFIGURE
!
!  Sets up physical system: particle positions, velocities
!
! ==============================================

subroutine configure

    use module_particle_props
    use module_physvars
    use module_geometry
    use module_velocity_setup
    use module_laser
    use module_io
    use module_particle_beam
 
    implicit none
    include 'mpif.h'


    integer :: i, ierr

    integer :: label_offset
    integer :: faces(maxlayers)
    real ::  V_layer(maxlayers), A_layer(maxlayers), Q_layer(maxlayers)
    real ::  qpart_layer(maxlayers), mass_layer(maxlayers), ai_layer(maxlayers)
    integer :: ipstart, ipstart_e, offset_e, ipstart_i, offset_i, nlayp, np_rest, ne_rest, ni_rest, nep0, nip0

    ! foam stuff
    integer :: nshell, ne_shell, ni_shell, ishell, jshell, kshell, nshell_x, nshell_y, nshell_z
    integer :: nep_shell, nip_shell, jion
    real :: Vshell, Ashell, Qshell, qe_shell, mass_e_shell, a_ee_shell, qi_shell, mass_i_shell, a_ii_shell


    np_local = nep+nip  ! initial # particles on this CPU
    nep0 = nep
    nip0 = nip
    !  Make particle numbers on root known (in case of unequal particle #s - need for label offset)
    call MPI_BCAST( nep0, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
    call MPI_BCAST( nip0, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
    ne_rest = nep0-nep
    ni_rest = nip0-nip


    if (restart) then

        call predef_parts    ! Predefined particle properties read from peXX/parts_dump.NNNNNN

        !  Find critical/shelf densities if laser switched on

        if (beam_config ==4 ) call track_nc 
        if (beam_config ==5) focus(1)=focus(1)+x_offset

    else
        ! Setup new plasma or special configuration

        ! Thermal velocities - Maxwellian 
	if (vte.eq.0) then
          vte = sqrt(Te_keV/511.)  ! convert from keV to /c
          vti = sqrt(Ti_keV/511./mass_ratio)
	else
	  ! Debye units: velocities normalised to vte
	  vti = vte/sqrt(mass_ratio)
	endif 

        new_config: select case(plasma_config)




        case(1)        ! Set up single neutral plasma target according to geometry
            !     =========================================================================

            if (debug_level==2 .and. my_rank==0) then
                write(*,*) "Setting up single plasma target"
            endif
            plasma_centre =  (/ xl/2., yl/2., zl/2. /) ! Centre of plasma

            offset_e = my_rank*nep + ne_rest
            offset_i = ne + my_rank*nip + ni_rest

            ! Electrons 
            call plasma_start( my_rank, 1, nep, ne, offset_e, target_geometry, idim, &
                -rho0, -1.0, 1.0, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre+displace, &
                number_faces, Vplas, Aplas, Qplas, qe, mass_e, a_ee )

	    call set_velocities(1, nep, vte, velocity_config, plasma_centre)

            !write (*,*) 'Electrons Vplas, Qplas:',Vplas, Qplas
            ! Ions
            call plasma_start( my_rank, nep+1, nip, ni, offset_i, target_geometry, idim, &
                rho0, 1.0, mass_ratio, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre+displace, &
                number_faces, Vplas, Aplas, Qplas, qi, mass_i, a_ii )

	    call set_velocities(nep+1, nip, vti, velocity_config, plasma_centre)

            !write (*,*) 'ions Vplas, Qplas:',Vplas, Qplas

            if (scheme /= 5 .and. ramp) then
                call add_ramp(x_plasma)     ! add exponential ramp to target (stretch container)
            endif

        case(2)        ! Special test config
            call special_start(ispecial)

        case(3)        ! Spatially separated ion and electron slabs

            target_geometry=0
            velocity_config=2   ! Ions cold, electrons with vx=vosc
            plasma_centre =  (/ xl/4., yl/2., zl/2. /) ! Centre of plasma
            offset_i = my_rank*nip + ni_rest
            offset_e = ni + my_rank*nep + ne_rest

            ! Ions
            call plasma_start( my_rank, 1, nip, ni, offset_i, target_geometry, idim, &
                rho0, 1.0, mass_ratio, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre, &
                number_faces, Vplas, Aplas, Qplas, qi, mass_i, a_ii )

	    call set_velocities(1, nip, vti, velocity_config, plasma_centre)

            ! Electrons shifted by displace vector 
            call plasma_start( my_rank, nip+1, nep, ne, offset_e, target_geometry, idim, &
                -rho0, -1.0, 1.0, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre+displace, &
                number_faces, Vplas, Aplas, Qplas, qe, mass_e, a_ee )

	    call set_velocities(nip+1, nep, vosc, velocity_config, plasma_centre)

        case(7)        ! Spherical Coulomb implosion


            if (debug_level==2 .and. my_rank==0) then
                write(*,*) "Setting up Coulomb implosion"
            endif
!            target_geometry=1
            velocity_config=3   ! Ions only with v_r=-v0 
            plasma_centre =  (/ xl/2., yl/2., zl/2. /) ! Centre of plasma

            ! Ions
            call plasma_start( my_rank, 1, nip, ni, 0, target_geometry, idim, &
                rho0, 1.0, mass_ratio, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre, &
                number_faces, Vplas, Aplas, Qplas, qi, mass_i, a_ii )

	    call set_velocities(1, nip, vosc, velocity_config, plasma_centre)


        case(4)        ! Spherical Coulomb explosion
            !    ============================================


            if (debug_level==2 .and. my_rank==0) then
                write(*,*) "Setting up ion sphere"
            endif
            target_geometry=1
            velocity_config=0   ! Ions cold, electrons with radial v0=vosc
            plasma_centre =  (/ xl/2., yl/2., zl/2. /) ! Centre of plasma
            !        plasma_centre =  (/ 0., 0., 0. /) ! Centre of plasma
            offset_e = my_rank*nep + ne_rest
            offset_i = ne + my_rank*nip + ni_rest

            ! Electrons can be shifted by displace vector 
            call plasma_start( my_rank, 1, nep, ne, offset_e, target_geometry, idim, &
                -rho0, -1.0, 1.0, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre+displace, &
                number_faces, Vplas, Aplas, Qplas, qe, mass_e, a_ee )

	    call set_velocities(1, nep, vosc, velocity_config, plasma_centre)

            ! Ions
            call plasma_start( my_rank, nep+1, nip, ni, offset_i, target_geometry, idim, &
                rho0, 1.0, mass_ratio, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre, &
                number_faces, Vplas, Aplas, Qplas, qi, mass_i, a_ii )

	    call set_velocities(nep+1, nip, vti, velocity_config, plasma_centre)

            if (scheme /= 5 .and. ramp) then
                call stretch_sphere(r_sphere)     ! create spherically symmetric density profile 
            endif




        case(40)        ! Foam: array of spherical shells 
            !    ================================================

            velocity_config=1   ! Ions, electrons thermal
            nshell_z=foam_geom(3)
            nshell_y=foam_geom(2)
            nshell_x=foam_geom(1)
            nshell = nshell_x * nshell_y * nshell_z
            ne_shell = ne/nshell
            ni_shell = ni/nshell
            nep_shell = nep/nshell
            nip_shell = nip/nshell
	    if (my_rank==0) then
	      write(*,*) "Particles/sphere/cpu: ",nep_shell
	      write(*,*) "Particles/sphere: ",ne_shell
	    endif
            ipstart_e = 1
            ipstart_i = nep+1
            offset_e = my_rank*nep
            offset_i = ne + my_rank*nip

! Physical transverse extent of target
	    
	    y_plasma = 2*r_sphere*nshell_y
	    z_plasma = 2*r_sphere*nshell_z
            plasma_centre =  (/ xl/2., yl/2.-y_plasma/2., zl/2.-z_plasma/2. /) ! Centre of first sphere 

            do ishell = 0,nshell_x-1
                do jshell = 0,nshell_y-1
                    do kshell = 0,nshell_z-1
                        ! offset for next sphere
                        displace = (/ 2*ishell*r_sphere, 2*jshell*r_sphere, 2*kshell*r_sphere /)

                        ! Electrons 

                        call plasma_start( my_rank, ipstart_e, nep_shell, ne_shell, offset_e, target_geometry, idim, &
                            -rho0, -1.0, 1.0, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre+displace, &
                            number_faces, Vshell, Ashell, Qshell, qe_shell, mass_e_shell, a_ee_shell )

	    		call set_velocities(ipstart_e, nep_shell, vte, velocity_config, plasma_centre)

                        ipstart_e = ipstart_e + nep_shell  ! index start
                        offset_e = offset_e + nep_shell  ! label offset

                        ! Ions
                        call plasma_start( my_rank, ipstart_i, nip_shell, ni_shell, offset_i, target_geometry, idim, &
                            rho0, 1.0, mass_ratio, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre+displace, &
                            number_faces, Vshell, Ashell, Qshell, qi_shell, mass_i_shell, a_ii_shell )

	    		call set_velocities(ipstart_i, nip_shell, vti, velocity_config, plasma_centre)

                        ipstart_i = ipstart_i + nip_shell  ! index start
                        offset_i = offset_i + nip_shell  ! label offset

                    end do
                end do
                ! shift for FCC
                if (mod(ishell,2)==0) then
                    plasma_centre = plasma_centre +  r_sphere/sqrt(2.)*(/ -1., 1., 1. /)
                else
                    plasma_centre = plasma_centre +  r_sphere/sqrt(2.)*(/ -1., -1., -1. /)
                endif
            end do
            Vplas = Vshell
            Aplas = Ashell
            Qplas = Qshell
            qe=qe_shell
            mass_e = mass_e_shell
            a_ee = a_ee_shell
            qi=qi_shell
            mass_i = mass_i_shell
            a_ii = a_ii_shell




        case(6)        ! Spherical cluster with Andreev profile
            !    ============================================
            !  r_sphere is radius of equivalent uniform sphere
            !  - used to define particle charges before stretching outer portion

            if (debug_level>=1 .and. my_rank==0) then
                write(*,*) "Setting up Andreev cluster"
            endif

            target_geometry=1
            velocity_config=1   ! Ions cold, electrons with vte
            plasma_centre =  (/ xl/2., yl/2., zl/2. /) ! Centre of plasma
            !        plasma_centre =  (/ 0., 0., 0. /) ! Centre of plasma
            offset_e = my_rank*nep + ne_rest
            offset_i = ne + my_rank*nip + ni_rest
	    vti = Ti_keV  ! Te represents vte/vmax in CE units
	    vte = Te_keV

! Ion mass set to unity to define CE timescale (omega_pi^-1)
            call plasma_start( my_rank, nep+1, nip, ni, offset_i, target_geometry, idim, &
                rho0, 1.0, 1.0, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre, &
                number_faces, Vplas, Aplas, Qplas, qi, mass_i, a_ii )

	    call set_velocities(nep+1, nip, vti, 0, plasma_centre)

            ! create spherically symmetric cluster with Andreev profile
            ! r_layer(1) is characteristic radius r0

            call cluster_sa(my_rank,n_cpu, nep+1,nip,r_layer(1),rho_layer(1),r_sphere,qi,Qplas,plasma_centre,mass_ratio)

            ! Electrons: use same geometry but reduced charge density
            ! - should get qe=-qi; masses reduced by miome

            call plasma_start( my_rank, 1, nep, ne, offset_e, target_geometry, idim, &
                -rho0*ne/ni, -1.0, 1.0/mass_ratio, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre, &
                number_faces, Vplas, Aplas, Q_layer(1), qe, mass_e, a_ee )

	    call set_velocities(1, nep, vte, velocity_config, plasma_centre)

            ! Stretch electron positions to match ions
            do i=1,nep
                jion=min(nep+nip/nep*i,nep+nip)
                x(i) = x(jion)+eps
                y(i) = y(jion)+ eps
                z(i) = z(jion) + eps
            end do



        case(10)  ! 2 overlapping ion species
            !    ===============================================

!            target_geometry=0
            velocity_config=1
            plasma_centre =  (/ xl/2., yl/2., zl/2. /) 
            offset_e = my_rank*nep + ne_rest
            offset_i = ne + my_rank*nip + ni_rest

            ! Electrons 
            call plasma_start( my_rank, 1, nep, ne, offset_e, target_geometry, idim, &
                -rho0, -1.0, 1.0, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre, &
                number_faces, Vplas, Aplas, Qplas, qe, mass_e, a_ee )

	    call set_velocities(1, nep, vte, velocity_config, plasma_centre)

            ! Ions
            call plasma_start( my_rank, nep+1, nip, ni, offset_i, target_geometry, idim, &
                rho0, 1.0, mass_ratio, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre, &
                number_faces, Vplas, Aplas, Qplas, qi, mass_i, a_ii )

	    call set_velocities(nep+1, nip, vosc, velocity_config, plasma_centre)

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

            ! Place on top of first slab 
            displace = (/0.,0.,0./)
            !       displace = (/ x_plasma/2.+x_layer(1)/2,0.,0. /)
            label_offset = ne+ni+my_rank*nlayp+ni_rest 

            ! Equal number of neutralising electrons 
            call plasma_start( my_rank, ipstart+nlayp, nlayp, n_layer(1), label_offset, target_geometry, idim, &
                -rho_layer(1), -1.0, 1.0, x_layer(1), y_layer(1), z_layer(1), r_layer(1), plasma_centre+displace, &
                faces(1), V_layer(1), A_layer(1), Q_layer(1), qpart_layer(1), mass_layer(1), ai_layer(1) )

	    call set_velocities(ipstart+nlayp, nlayp, vte, velocity_config, plasma_centre)

            label_offset = ne+ni+n_layer(1) +my_rank*nlayp + ne_rest 
            call plasma_start( my_rank, ipstart, nlayp, n_layer(1), label_offset, target_geometry, idim, &
                rho_layer(1), 1.0, mratio_layer(1), x_layer(1), y_layer(1), z_layer(1), r_layer(1), plasma_centre+displace, &
                faces(1), V_layer(1), A_layer(1), Q_layer(1), qpart_layer(1), mass_layer(1), ai_layer(1) )

	    call set_velocities(ipstart, nlayp, vti, velocity_config, plasma_centre)

            if (debug_level==2 .and. my_rank==0) then
                write(*,*) "proton charge ",qpart_layer(1)
                write(*,*) "proton mass ",mass_layer(1)
                write(*,*) "spacing",ai_layer(1)
            endif


            np_local=np_local + 2*nlayp  ! Total # local particles
            ne = ne + n_layer(1)  ! Global # particles
            ni = ni + n_layer(1)  ! Global # particles
            npart_total = ni+ne

            if (scheme /= 5 .and. ramp) then
                call add_ramp(x_plasma)     ! add exponential ramp to target (stretch container)
            endif



            !###########################################################################################################

        case(11)  ! Add proton disc to slab
            !    ====================================

            write(ipefile,'(/a/)') "Setting up main slab"

            target_geometry=0
            velocity_config=1
            plasma_centre =  (/ xl/2., yl/2., zl/2. /) 
            offset_e = my_rank*nep + ne_rest
            offset_i = ne + my_rank*nip + ni_rest

            ! Electrons 
            call plasma_start( my_rank, 1, nep, ne, offset_e, target_geometry, idim, &
                -rho0, -1.0, 1.0, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre, &
                number_faces, Vplas, Aplas, Qplas, qe, mass_e, a_ee )

	    call set_velocities(1, nep, vte, velocity_config, plasma_centre)

            ! Ions
            call plasma_start( my_rank, nep+1, nip, ni, offset_i, target_geometry, idim, &
                rho0, 1.0, mass_ratio, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre, &
                number_faces, Vplas, Aplas, Qplas, qi, mass_i, a_ii )

	    call set_velocities(nep+1, nip, vti, velocity_config, plasma_centre)

            ! Proton disc 

            write(ipefile,'(/a/)') "Setting up proton disc"

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
            label_offset = ne+ni+my_rank*nlayp +n_layer(1) + ni_rest

            ! Place on rear of main slab 
            displace = (/ x_plasma/2.+x_layer(1)/2,0.,0. /)
            !	layer_geometry = 2  ! disc

            call plasma_start( my_rank, ipstart, nlayp, n_layer(1), label_offset, layer_geometry, idim, &
                rho_layer(1), 1.0, mratio_layer(1), x_layer(1), y_layer(1), z_layer(1), r_layer(1), plasma_centre+displace, &
                faces(1), V_layer(1), A_layer(1), Q_layer(1), qpart_layer(1), mass_layer(1), ai_layer(1) )
	    call set_velocities(ipstart, nlayp, vti, velocity_config, plasma_centre)

            ! Equal number of neutralising electrons 
            label_offset = ne+ni+my_rank*nlayp + ne_rest 
            call plasma_start( my_rank, ipstart+nlayp, nlayp, n_layer(1), label_offset, layer_geometry, idim, &
                -rho_layer(1), -1.0, 1.0, x_layer(1), y_layer(1), z_layer(1), r_layer(1), plasma_centre+displace, &
                faces(1), V_layer(1), A_layer(1), Q_layer(1), qpart_layer(1), mass_layer(1), ai_layer(1) )
	    call set_velocities(ipstart+nlayp, nlayp, vte, velocity_config, plasma_centre)

            np_local=np_local + 2*nlayp  ! Total # local particles
            ne = ne + n_layer(1)  ! Global # particles
            ni = ni + n_layer(1)  ! Global # particles
            npart_total = ni+ne

            if (scheme /= 5 .and. ramp) then
                call add_ramp(x_plasma)     ! add exponential ramp to target (stretch container)
            endif

            !###########################################################################################################

        case(16)  ! Add proton disc to bowl
            !    ====================================

            write(ipefile,'(/a/)') "Setting up bowl"

            target_geometry=36
            velocity_config=1
            plasma_centre =  (/ xl/2., yl/2., zl/2. /) 
            offset_e = my_rank*nep + ne_rest
            offset_i = ne + my_rank*nip + ni_rest

            ! Electrons 
            call plasma_start( my_rank, 1, nep, ne, offset_e, target_geometry, idim, &
                -rho0, -1.0, 1.0, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre, &
                number_faces, Vplas, Aplas, Qplas, qe, mass_e, a_ee )

	    call set_velocities(1, nep, vte, velocity_config, plasma_centre)

            ! Ions
            call plasma_start( my_rank, nep+1, nip, ni, offset_i, target_geometry, idim, &
                rho0, 1.0, mass_ratio, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre, &
                number_faces, Vplas, Aplas, Qplas, qi, mass_i, a_ii )

	    call set_velocities(nep+1, nip, vti, velocity_config, plasma_centre)

            ! Proton disc 

            write(ipefile,'(/a/)') "Setting up proton disc"

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
            label_offset = ne+ni+my_rank*nlayp +n_layer(1) + ni_rest

            ! Place on rear (inside) surface of shell  
            displace = (/ x_plasma+x_layer(1)/2-r_sphere/2.,0.,0. /)
            !	layer_geometry = 2  ! disc

            call plasma_start( my_rank, ipstart, nlayp, n_layer(1), label_offset, layer_geometry, idim, &
                rho_layer(1), 1.0, mratio_layer(1), x_layer(1), y_layer(1), z_layer(1), r_layer(1), plasma_centre+displace, &
                faces(1), V_layer(1), A_layer(1), Q_layer(1), qpart_layer(1), mass_layer(1), ai_layer(1) )

	    call set_velocities(ipstart, nlayp, vti, velocity_config, plasma_centre)

            ! Equal number of neutralising electrons 
            label_offset = ne+ni+my_rank*nlayp + ne_rest 
            call plasma_start( my_rank, ipstart+nlayp, nlayp, n_layer(1), label_offset, layer_geometry, idim, &
                -rho_layer(1), -1.0, 1.0, x_layer(1), y_layer(1), z_layer(1), r_layer(1), plasma_centre+displace, &
                faces(1), V_layer(1), A_layer(1), Q_layer(1), qpart_layer(1), mass_layer(1), ai_layer(1) )

	    call set_velocities(ipstart+nlayp, nlayp, vte, velocity_config, plasma_centre)

            np_local=np_local + 2*nlayp  ! Total # local particles
            ne = ne + n_layer(1)  ! Global # particles
            ni = ni + n_layer(1)  ! Global # particles
            npart_total = ni+ne

            if (scheme /= 5 .and. ramp) then
                call add_ramp(x_plasma)     ! add exponential ramp to target (stretch container)
            endif

#ifdef HAVE_PRIVATE_CONFIGS
	include 'config_private.h'   ! Non-public configuration cases
#endif

        case default     ! Default = 0 - no plasma target
            if (my_rank==0) then
		write (6,*) 'Warning: no plasma set up'
		write (6,*) 'Chosen config ',plasma_config,' not found'
		write (6,*) '- check you have built pepcb with OPT="-DHAVE_PRIVATE_CONFIGS"'
		write (6,*) '  and that config exists in config_private.h'
		call closefiles
  	  	call MPI_FINALIZE(ierr)
		stop
	    endif

            npart_total=0
            np_local = 0
            Vplas = x_plasma * y_plasma * z_plasma
            Aplas = x_plasma * y_plasma
            plasma_centre =  (/ xl/2., yl/2., zl/2. /) ! Centre of plasma
            Qplas = 1.
        end select new_config

 ! zero fields
  if (idim==2) uz(1:np_local)=0.
!  Ex(1:np_local) = 0
!  Ey(1:np_local) = 0
!  Ez(1:np_local) = 0
!  Bx(1:np_local) = 0
!  By(1:np_local) = 0
!  Bz(1:np_local) = 0
!  Ax(1:np_local) = 0
!  Ay(1:np_local) = 0
!  Az(1:np_local) = 0
!  Axo(1:np_local) = 0
! Ayo(1:np_local) = 0
!  Azo(1:np_local) = 0
!  pot(1:np_local) = 0
  work(1:np_local) = 1.   ! set work load balanced initially

        laser_config: select case(plasma_config)
        case(0)
            focus = (/ xl/4., yl/2., zl/2. /) ! Centre of laser focal spot

        case(1,3,4)
            ! Centre of laser focal spot
            laser_focus: select case(target_geometry)

            case(0) ! slab
                focus = (/ xl/2. + x_offset, yl/2., zl/2. /) 
            case(1,11) ! sphere
                focus = (/ xl/2.-r_sphere, yl/2., zl/2./) 
            case(2,12) ! disc, tube
                focus = (/ xl/2.-x_plasma/2., yl/2., zl/2. /) 
            case(3,13) ! wire
                focus = (/ xl/2.-r_sphere+x_offset, yl/2., zl/2.+z_offset /) 
            case(4) ! ellipsoid
                focus = (/ xl/2.-x_plasma*r_sphere, yl/2., zl/2. /) 
            case(5) ! wedge
                focus = (/ xl/2.-x_plasma/2., yl/2., zl/2. /)
            case(6,16,26,36) ! hemisphere/hollow hemisphere
                focus = (/ xl/2.-r_sphere/2., yl/2., zl/2. /)

            end select laser_focus
        case(2)
	case(10)
            focus = plasma_centre -  (/ x_plasma/2.+x_layer(1),0.,0. /)
	case(11)
            focus = plasma_centre -  (/ x_plasma/2., 0., 0. /)
        case default
            focus = (/ xl/4., yl/2., zl/2. /) ! Centre of laser focal spot
        end select laser_config
    endif


    ! Parameters & conversion factors derived from target/particle config.

    window_min = plasma_centre(1) - x_plasma/2.
    propag_laser=focus(1)
    convert_fs = 10.*omega*lambda/(6*pi)     ! convert from wp^-1 to fs
    convert_mu = omega/2./pi*lambda          ! convert from c/wp to microns
    r_neighbour = fnn*a_ii  ! Nearest neighbour search radius
    navcycle = 2*pi/dt/omega  ! # timesteps in a laser cycle
    nu_ei = 1./40./pi*a_ii**3/max(vte,1.e-8)/eps**2  ! collision frequency (fit to Okuda & Birdsall)
    sigma_e = 1./nu_ei   ! Spitzer conductivity
    intensity = 0.2*vosc**2*omega**2  ! normalised laser intensity

    energy_units: select case(plasma_config)  ! TODO: Should have extra switch for unit system

    case(6) ! Andreev profile - allow for x 1/2 factor in energies with CE norms 
	convert_erg = 2*(ni*4.8e-10)**2/r_layer(1)  ! Radius in cm
	convert_kev = convert_erg/1.6e-9

    case default
      if (ne>0) then
        convert_keV = 511./abs(qe)     ! convert from code energy units to keV
        convert_erg = 8.16e-7/abs(qe)
      else
        convert_keV = 511./abs(qi)
        convert_erg = 8.16e-7/qi
      endif

    end select energy_units


    if (te_perturb) then
        call perturb_temp    ! Impose perturbation on Te for transport test
#ifdef VISIT_NBODY
        if ( vis_on ) then
            call densities
            call sum_fields
            call vis_fields_nbody(itime+itime_start)
        endif
#endif

    endif

    ! Ion config - make sure potential repulsive
    if (scheme==5) then
        eps = a_ii*2.
    endif


    if (launch) then ! Ignore if still in config loop

        beamconf: select case(beam_config)  ! Configure laser or particle beam

        case(1)
            call beam           ! Fixed beam
#ifdef VISIT_NBODY
            if (steering) call vis_control   ! Display default parameters
#endif
        case(2)
#ifdef VISIT_NBODY
            if (steering) call vis_control   ! Constant particle source
#endif
            if (my_rank==0) write(*,'(//a)') '===> Particle beam switched on' 

        case(8)
            call beam_dust   ! Dust particle

        case(3:6) ! laser on

            if (my_rank==0) write(*,'(//a)') '===> Laser switched on' 
#ifdef VISIT_NBODY
            if (steering) call vis_control 
#endif

        end select beamconf


      endif
   End subroutine configure
   
