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

  integer :: i, ipe, idummy=0, ierr, ifile
  real :: t_walk, t_walkc, t_force, t_domain,t_build,t_prefetch
  integer :: label_offset
  integer :: faces(maxlayers)
  real ::  V_layer(maxlayers), A_layer(maxlayers), Q_layer(maxlayers)
  real ::  qpart_layer(maxlayers), mass_layer(maxlayers), ai_layer(maxlayers)
  integer :: ipstart, nlayp, np_rest

  npp = nep+nip  ! initial # particles on this CPU


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

     if (debug_level==2) then
	write(*,*) "Setting up single plasma target"
     endif
     plasma_centre =  (/ xl/2., yl/2., zl/2. /) ! Centre of plasma

     velocity_config=1
 ! Electrons 
        call plasma_start( 1, nep, ne, 0, target_geometry, velocity_config, idim, &
               -rho0, -1.0, 1.0, vte, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre, &
               number_faces, Vplas, Aplas, Qplas, qe, mass_e, a_ee )
write (*,*) 'Electrons Vplas, Qplas:',Vplas, Qplas
 ! Ions
        call plasma_start( nep+1, nip, ni, ne, target_geometry, velocity_config, idim, &
               rho0, 1.0, mass_ratio, vti, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre, &
               number_faces, Vplas, Aplas, Qplas, qi, mass_i, a_ii )
     
write (*,*) 'ions Vplas, Qplas:',Vplas, Qplas
        if (scheme /= 5 .and. ramp) then
           call add_ramp     ! add exponential ramp to target (stretch container)
        endif

     case(2)        ! Special test config
        call special_start(ispecial)

     case(3)        ! Spatially separated ion and electron slabs

	target_geometry=0
        velocity_config=2   ! Ions cold, electrons with vx=vosc
        plasma_centre =  (/ xl/4., yl/2., zl/2. /) ! Centre of plasma
 ! Ions
        call plasma_start( 1, nip, ni, 0, target_geometry, 0, idim, &
               rho0, 1.0, mass_ratio, vti, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre, &
               number_faces, Vplas, Aplas, Qplas, qi, mass_i, a_ii )

 ! Electrons shifted by displace vector 
        call plasma_start( nip+1, nep, ne, ni, target_geometry, velocity_config, idim, &
               -rho0, -1.0, 1.0, vosc, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre+displace, &
               number_faces, Vplas, Aplas, Qplas, qe, mass_e, a_ee )

     case(4)        ! Spherical Coulomb explosion

	target_geometry=1
        velocity_config=3   ! Ions cold, electrons with radial v0=vosc
!        plasma_centre =  (/ xl/2., yl/2., zl/2. /) ! Centre of plasma
        plasma_centre =  (/ 0., 0., 0. /) ! Centre of plasma

 ! Electrons can be shifted by displace vector 
        call plasma_start( 1, nep, ne, 0, target_geometry, velocity_config, idim, &
               -rho0, -1.0, 1.0, vosc, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre+displace, &
               number_faces, Vplas, Aplas, Qplas, qe, mass_e, a_ee )
 
 ! Ions
        call plasma_start( nep+1, nip, ni, ne, target_geometry, 0, idim, &
               rho0, 1.0, mass_ratio, vti, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre, &
               number_faces, Vplas, Aplas, Qplas, qi, mass_i, a_ii )


     case(10)  ! Add proton layer to slab

     target_geometry=0
     velocity_config=1
     plasma_centre =  (/ xl/2., yl/2., zl/2. /) 

 ! Electrons 
        call plasma_start( 1, nep, ne, 0, target_geometry, velocity_config, idim, &
               -rho0, -1.0, 1.0, vte, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre, &
               number_faces, Vplas, Aplas, Qplas, qe, mass_e, a_ee )
 ! Ions
        call plasma_start( nep+1, nip, ni, ne, target_geometry, velocity_config, idim, &
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

  	ipstart = nep+nip+1
	label_offset = ne+ni 
! Place on front of main slab 
	displace = plasma_centre - (/ x_plasma/2.+x_layer(1)/2,0.,0. /)

        call plasma_start( ipstart, nlayp, n_layer(1), label_offset, target_geometry, velocity_config, idim, &
               rho_layer(1), 1.0, mratio_layer(1), vti, x_layer(1), y_layer(1), z_layer(1), r_layer(1), displace, &
               faces(1), V_layer(1), A_layer(1), Q_layer(1), qpart_layer(1), mass_layer(1), a_layer(1) )

 ! Equal number of neutralising electrons 
	label_offset = ne+ni+n_layer(1)  
        call plasma_start( ipstart+nlayp, nlayp, n_layer(1), label_offset, target_geometry, velocity_config, idim, &
               -rho_layer(1), -1.0, 1.0, vte, x_layer(1), y_layer(1), z_layer(1), r_layer(1), displace, &
               faces(1), V_layer(1), A_layer(1), Q_layer(1), qpart_layer(1), mass_layer(1), ai_layer(1) )

       npp=npp + 2*nlayp  ! Total # local particles
       npart = npart + 2 * n_layer(1)  ! Global # particles

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

        case(1,3:20)
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
     end select laser_config
  endif


! Parameters & conversion factors derived from target/particle config.

  window_min = plasma_centre(1) - x_plasma/2.
  propag_laser=focus(1)
  convert_fs = 10.*omega*lambda/(6*pi)     ! convert from wp^-1 to fs
  convert_mu = omega/2./pi*lambda          ! convert from c/wp to microns
  lolam = lolam*2.*pi/omega  ! normalise scale-length
  convert_keV = 2./3./abs(Qplas)*511     ! convert from code energy units to keV/particle (Temperature)
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
     eps = a_ii*2
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
  call pepc_fields_p(np_local, mac, theta, ifreeze, eps, force_tolerance, balance,&
          force_const, bond_const, &
          dt, xl, yl, zl, 0, &
          coulomb, bfields, bonds, lenjones, &
          t_domain,t_build,t_prefetch,t_walk,t_walkc,t_force, iprot)   

!  Initialise vec. pots. to avoid jump in induced E-field
  Axo(1:npp) = Ax(1:npp)
  Ayo(1:npp) = Ay(1:npp)
  Azo(1:npp) = Az(1:npp)

end subroutine configure

