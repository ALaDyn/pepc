! =============================================
!
!                CONFIGURE
!
!  Sets up physical system: particle positions, velocities
!
! ==============================================

subroutine configure

  use physvars   ! Use internal particle arrays from lpepc
  use treevars
  implicit none
  include 'mpif.h'

  integer :: i, ipe, jion, idummy=0, ierr, ifile
  real :: t_walk, t_walkc, t_force, t_domain,t_build,t_prefetch
  integer :: label_offset
  integer :: faces(maxlayers)
  real ::  V_layer(maxlayers), A_layer(maxlayers), Q_layer(maxlayers)
  real ::  qpart_layer(maxlayers), mass_layer(maxlayers), ai_layer(maxlayers)
  integer :: ipstart, ipstart_e, offset_e, ipstart_i, offset_i, nlayp, np_rest, ne_rest, ni_rest, nep0, nip0
  !APLR stuff:
  integer :: nlaypfront,nlaypback
! foam stuff
  integer :: nshell, ne_shell, ni_shell, ishell, jshell, kshell, nshell_x, nshell_y, nshell_z
  integer :: nep_shell, nip_shell
  real :: Vshell, Ashell, Qshell, qe_shell, mass_e_shell, a_ee_shell, qi_shell, mass_i_shell, a_ii_shell

  npp = nep+nip  ! initial # particles on this CPU
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

    vte = sqrt(Te_keV/511.)  ! convert from keV to /c
    vti = sqrt(Ti_keV/511./mass_ratio)

     new_config: select case(plasma_config)

     case(1)        ! Set up single neutral plasma target according to geometry
!     =========================================================================

     if (debug_level==2 .and. me==0) then
	write(*,*) "Setting up single plasma target"
     endif
     plasma_centre =  (/ xl/2., yl/2., zl/2. /) ! Centre of plasma

     velocity_config=1
     offset_e = me*nep + ne_rest
     offset_i = ne + me*nip + ni_rest

 ! Electrons 
        call plasma_start( 1, nep, ne, offset_e, target_geometry, velocity_config, idim, &
               -rho0, -1.0, 1.0, vte, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre+displace, &
               number_faces, Vplas, Aplas, Qplas, qe, mass_e, a_ee )
!write (*,*) 'Electrons Vplas, Qplas:',Vplas, Qplas
 ! Ions
        call plasma_start( nep+1, nip, ni, offset_i, target_geometry, velocity_config, idim, &
               rho0, 1.0, mass_ratio, vti, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre+displace, &
               number_faces, Vplas, Aplas, Qplas, qi, mass_i, a_ii )
     
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
	offset_i = me*nip + ni_rest
	offset_e = ni + me*nep + ne_rest

 ! Ions
        call plasma_start( 1, nip, ni, offset_i, target_geometry, 0, idim, &
               rho0, 1.0, mass_ratio, vti, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre, &
               number_faces, Vplas, Aplas, Qplas, qi, mass_i, a_ii )

 ! Electrons shifted by displace vector 
        call plasma_start( nip+1, nep, ne, offset_e, target_geometry, velocity_config, idim, &
               -rho0, -1.0, 1.0, vosc, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre+displace, &
               number_faces, Vplas, Aplas, Qplas, qe, mass_e, a_ee )


     case(4)        ! Spherical Coulomb explosion
!    ============================================


     if (debug_level==2 .and. me==0) then
	write(*,*) "Setting up ion sphere"
     endif
	target_geometry=1
        velocity_config=0   ! Ions cold, electrons with vte
        plasma_centre =  (/ xl/2., yl/2., zl/2. /) ! Centre of plasma
!        plasma_centre =  (/ 0., 0., 0. /) ! Centre of plasma
	offset_e = me*nep + ne_rest
	offset_i = ne + me*nip + ni_rest

 ! Electrons can be shifted by displace vector 
        call plasma_start( 1, nep, ne, offset_e, target_geometry, 1, idim, &
               -rho0, -1.0, 1.0, vte, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre+displace, &
               number_faces, Vplas, Aplas, Qplas, qe, mass_e, a_ee )
 
 ! Ions
        call plasma_start( nep+1, nip, ni, offset_i, target_geometry, 0, idim, &
               rho0, 1.0, mass_ratio, vti, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre, &
               number_faces, Vplas, Aplas, Qplas, qi, mass_i, a_ii )


      if (scheme /= 5 .and. ramp) then
           call stretch_sphere(r_sphere)     ! create spherically symmetric density profile 
      endif

     case(5)        ! Foam: array of spherical shells 
!    ================================================

	target_geometry=7
        velocity_config=1   ! Ions, electrons thermal
        plasma_centre =  (/ xl/2., yl/2., zl/2. /) ! Centre of plasma
	nshell_z=5
	nshell_y=3
	nshell_x=2
	nshell = nshell_x * nshell_y * nshell_z
	ne_shell = ne/nshell
	ni_shell = ni/nshell
	nep_shell = nep/nshell
	nip_shell = nip/nshell
        ipstart_e = 1
        ipstart_i = nep+1
	offset_e = me*nep
	offset_i = ne + me*nip

	do ishell = 0,nshell_z-1
	do jshell = 0,nshell_y-1
	do kshell = 0,nshell_x-1
	  displace = (/ 2*ishell*r_sphere, 2*jshell*r_sphere, 2*kshell*r_sphere /)

 ! Electrons 

        call plasma_start( ipstart_e, nep_shell, ne_shell, offset_e, target_geometry, velocity_config, idim, &
               -rho0, -1.0, 1.0, vte, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre+displace, &
               number_faces, Vshell, Ashell, Qshell, qe_shell, mass_e_shell, a_ee_shell )

	ipstart_e = ipstart_e + nep_shell  ! index start
	offset_e = offset_e + nep_shell  ! label offset
 
 ! Ions
        call plasma_start( ipstart_i, nip_shell, ni_shell, offset_i, target_geometry, 0, idim, &
               rho0, 1.0, mass_ratio, vti, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre+displace, &
               number_faces, Vshell, Ashell, Qshell, qi_shell, mass_i_shell, a_ii_shell )

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


     case(6)        ! Spherical cluster with Andreev profile
!    ============================================
!  r_sphere is radius of equivalent uniform sphere
!  - used to define particle charges before stretching outer portion

     if (debug_level>=1 .and. me==0) then
	write(*,*) "Setting up Andreev cluster"
     endif

	target_geometry=1
        velocity_config=1   ! Ions cold, electrons with vte
        plasma_centre =  (/ xl/2., yl/2., zl/2. /) ! Centre of plasma
!        plasma_centre =  (/ 0., 0., 0. /) ! Centre of plasma
	offset_e = me*nep + ne_rest
	offset_i = ne + me*nip + ni_rest

 
 ! Ions
        call plasma_start( nep+1, nip, ni, offset_i, target_geometry, 0, idim, &
               rho0, 1.0, mass_ratio, vti, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre, &
               number_faces, Vplas, Aplas, Qplas, qi, mass_i, a_ii )

! create spherically symmetric cluster with Andreev profile
! r_layer(1) is characteristic radius r0
      
        call cluster_sa(nep+1,nip,r_layer(1),r_sphere,qi,Qplas,plasma_centre)
   

! Electrons: use same geometry but reduced charge density
! - should get qe=-qi
 
        call plasma_start( 1, nep, ne, offset_e, target_geometry, velocity_config, idim, &
               -rho0*ne/ni, -1.0, 1.0, vte, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre, &
               number_faces, Vplas, Aplas, Q_layer(1), qe, mass_e, a_ee )

! Stretch electron positions to match ions
	do i=1,nep
	   jion=min(nep+nip/nep*i,nep+nip)
	   x(i) = x(jion)+eps
	   y(i) = y(jion)+ eps
	   z(i) = z(jion) + eps
        end do

     case(10)  ! Add proton layer to primary target - arbitrary geometry for both layers
!    =====================================================================

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
!	displace = (/ x_plasma/2.+x_layer(1)/2,0.,0. /)
	label_offset = ne+ni+me*nlayp+ni_rest 

 ! Equal number of neutralising electrons 
        call plasma_start( ipstart+nlayp, nlayp, n_layer(1), label_offset, layer_geometry, velocity_config, idim, &
               -rho_layer(1), -1.0, 1.0, vte, x_layer(1), y_layer(1), z_layer(1), r_layer(1), plasma_centre+displace, &
               faces(1), V_layer(1), A_layer(1), Q_layer(1), qpart_layer(1), mass_layer(1), ai_layer(1) )

	label_offset = ne+ni+n_layer(1) +me*nlayp + ne_rest 
        call plasma_start( ipstart, nlayp, n_layer(1), label_offset, layer_geometry, velocity_config, idim, &
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


!###########################################################################################################
     case(11)  ! Add proton disc to slab
!    ====================================

   if (debug_level>1) write(ipefile,'(/a/)') "Setting up main slab"

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

 ! Proton disc 

     if(debug_level>1) write(ipefile,'(/a/)') "Setting up proton disc"

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
	label_offset = ne+ni+me*nlayp +n_layer(1) + ni_rest
 
! Place on rear of main slab 
	displace = (/ x_plasma/2.+x_layer(1)/2,0.,0. /)
!	layer_geometry = 2  ! disc

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

!###########################################################################################################
     case(21)  ! Add proton disc to slab
!    ====================================

   if (debug_level>1) write(ipefile,'(/a/)') "Setting up main slab"

!     target_geometry=0
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

 ! Proton disc 

     if(debug_level>1) write(ipefile,'(/a/)') "Setting up proton disc"

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
	label_offset = ne+ni+me*nlayp +n_layer(1) + ni_rest
 
! Place on rear of main slab 
	displace = displace + (/ x_plasma/2.,0.,0. /)
!	layer_geometry = 2  ! disc

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




     case(15)  ! 2-layer slab with foam
!    ==================================

     target_geometry=0
     velocity_config=1
     plasma_centre =  (/ xl/2., yl/2., zl/2. /) 
     offset_e = me*nep + ne_rest
     offset_i = ne + me*nip + ni_rest

!  Main block
     if (debug_level>1) write(ipefile,'(/a/)') "Setting up main slab"

 ! Electrons 
        call plasma_start( 1, nep, ne, offset_e, target_geometry, velocity_config, idim, &
               -rho0, -1.0, 1.0, vte, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre, &
               number_faces, Vplas, Aplas, Qplas, qe, mass_e, a_ee )
 ! Ions
        call plasma_start( nep+1, nip, ni, offset_i, target_geometry, velocity_config, idim, &
               rho0, 1.0, mass_ratio, vti, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre, &
               number_faces, Vplas, Aplas, Qplas, qi, mass_i, a_ii )

 ! Protons
      if (debug_level>1) write(ipefile,'(/a/)') "Setting up proton layer"

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

     if (debug_level>1) write(ipefile,'(/a/)') "Setting up foam layer"

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

     case default     ! Default = 0 - no plasma target
        if (me==0) write (6,*) 'Warning: no plasma set up'
        npart=0
        npp = 0
        Vplas = x_plasma * y_plasma * z_plasma
        Aplas = x_plasma * y_plasma
        plasma_centre =  (/ xl/2., yl/2., zl/2. /) ! Centre of plasma
        Qplas = 1.
     end select new_config



     laser_config: select case(plasma_config)
        case(0)
        focus = (/ xl/4., yl/2., zl/2. /) ! Centre of laser focal spot

        case(1,3,4)
! Centre of laser focal spot
        laser_focus: select case(target_geometry)

          case(0) ! slab
           focus = (/ xl/2. + x_offset, yl/2., zl/2. /) 
          case(1) ! sphere
           focus = (/ xl/2.-r_sphere, yl/2., zl/2./) 
          case(2) ! disc
           focus = (/ xl/2.-x_plasma/2., yl/2., zl/2. /) 
          case(3) ! wire
           focus = (/ xl/2.-r_sphere+x_offset, yl/2., zl/2.+z_offset /) 
          case(4) ! ellipsoid
           focus = (/ xl/2.-x_plasma*r_sphere, yl/2., zl/2. /) 
          case(5) ! wedge
           focus = (/ xl/2.-x_plasma/2., yl/2., zl/2. /)
          case(6) ! hemisphere
           focus = (/ xl/2.-r_sphere/2., yl/2., zl/2. /)
          case(7) ! hollow sphere
           focus = (/ xl/2.-r_sphere/2., yl/2., zl/2. /)
          case(8) ! hollow hemisphere
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
  lolam = lolam*2.*pi/omega  ! normalise scale-length
  if (ne>0) then
	convert_keV = 511./abs(qe)     ! convert from code energy units to keV
  	convert_erg = 8.16e-7/abs(qe)
  else
	convert_keV = 511./abs(qi)
  	convert_erg = 8.16e-7/qi
  endif 
  r_neighbour = fnn*a_ii  ! Nearest neighbour search radius
  navcycle = 2*pi/dt/omega  ! # timesteps in a laser cycle
  nu_ei = 1./40./pi*a_ii**3/max(vte,1.e-8)/eps**2  ! collision frequency (fit to Okuda & Birdsall)
  sigma_e = 1./nu_ei   ! Spitzer conductivity
  intensity = 0.2*vosc**2*omega**2  ! normalised laser intensity


  if (mc_init) call mc_config  ! Do MC min-PE initialisation depending on config

  if (te_perturb) then
     call perturb_temp    ! Impose perturbation on Te for transport test
     if ( vis_on ) then
        call densities
        call sum_fields
        call vis_fields_nbody(itime+itime_start)
     endif

  endif

! Ion config - make sure potential repulsive
  if (scheme==5) then
     eps = a_ii*2.
  endif

  beamconf: select case(beam_config)  ! Configure laser or particle beam

  case(1)
     call beam           ! Fixed beam
     if (steering) call vis_control   ! Display default parameters

  case(2)
     if (steering) call vis_control   ! Constant particle source
     if (me==0) write(*,'(//a)') '===> Particle beam switched on' 

  case(8)
     call beam_dust   ! Dust particle

  case(3:6) ! laser on

     if (me==0) write(*,'(//a)') '===> Laser switched on' 
     if (steering) call vis_control 

  end select beamconf


! local and global # particles for  main routine
  np_local = npp
  npart_total = npart

  call param_dump      ! Dump initial data

! Compute initial field values - need these to get vec. pots consistent with velocities

  if (me==0) write(*,*) 'Computing initial fields'
  call pepc_fields_p(np_local,  walk_scheme, mac, theta, ifreeze, eps, force_tolerance, balance,&
          force_const, bond_const, &
          dt, xl, yl, zl, 0, &
          coulomb, bfields, bonds, lenjones, &
          t_domain,t_build,t_prefetch,t_walk,t_walkc,t_force, iprot, work_tot)   

!  Initialise vec. pots. to avoid jump in induced E-field
  Axo(1:npp) = Ax(1:npp)
  Ayo(1:npp) = Ay(1:npp)
  Azo(1:npp) = Az(1:npp)

end subroutine configure

