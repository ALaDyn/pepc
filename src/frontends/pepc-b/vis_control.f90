! ==============================================
!
!                BEAM_CONTROL
!
!  Interactive control of particle source
!
! ==============================================

subroutine vis_control

  use physvars
  use treevars
  use utils
  implicit none
  include 'mpif.h'

  integer :: i, p, iseed1, iseed2, ierr
  real :: Volb, dpx, yt, zt, vosc_old, sigma_old, tpulse_old, u_old, theta_old, phi_old
  integer :: lvisit_active
  real :: ct, st, cp, sp, vx_beam, vy_beam, vz_beam, th_beam, xb, yb, zb
  logical :: beam_on = .true.
  logical :: beam_debug = .true.
  integer :: isteer1 = 0, isteer2 =0, isteer3=0, isteer4=0
  real*8 :: dsteer1, dsteer2, dsteer3, dsteer4

 
  integer, save :: np_beam_dt  ! current # beam particles

  ! First check for VISIT connection

!  if (me==0)   call flvisit_spk_check_connection(lvisit_active)
#ifdef VISIT_NBODY
  if (me==0)   call flvisit_nbody2_check_connection(lvisit_active)
  call MPI_BCAST( lvisit_active, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
#endif

  if (lvisit_active==0 )then
     if (me==0) write(*,*) 'VISNB | No connection to visualization'
     return
  endif

  if (itime == 0) then
     np_beam_dt = np_beam  ! initial # beam particles to introduce per timestep
  endif

  !  if (.not. beam_on) return

  !  if (npart + np_beam_dt > npartm .or. max_list_length > nintmax-20) then
  !     if (me==0) write(*,*) 'Array limit reached: switching off beam'
  !     beam_on = .false.
  !     return
  !  endif


!  if (beam_config==7) 
!  Constant Bz
!     dsteer1 = vosc
!     dsteer2 = sigma
!     dsteer3 = tpulse

  if ( beam_config >=5 .and. beam_config <=7) then

!  else if ( beam_config == 5 .or. beam_config==6) then
      ! Define beam from laser parameters
     vosc_old = vosc
     sigma_old = sigma
     dsteer1 = vosc
     dsteer2 = sigma
     dsteer3 = tpulse  ! incidence angle instead of pulse duration

  else if ( beam_config >=3 .and. beam_config <=4) then
      ! Define beam from laser parameters
     vosc_old = vosc
     sigma_old = sigma
     dsteer1 = vosc
     dsteer2 = sigma
     dsteer3 = theta_beam  ! incidence angle instead of pulse duration
     theta_old = th_beam

  else if (scheme==4 .and. .not. (beam_config<=3 .and. beam_config>0)) then
     ! Temperature clamp mode - laser should be off
     dsteer1 = Te_kev

  else if (scheme == 5) then
     ! ion crystal eqm  mode:
     !  r_beam is mean ion spacing
     !  u_beam is ion temperature (eV)
     !  rho_beam is potential constant
     dsteer2 = Ti_keV
     dsteer3 = dt 
     dsteer4 = eps

  else if (beam_config==8) then ! dust particle
     u_old = u_beam
     theta_old = theta_beam
     phi_old = phi_beam
  endif

#ifdef VISIT_NBODY
  if (itime == 0 .and. me==0 )  then
     call flvisit_nbody2_check_connection(lvisit_active)
  ! Specify default parameters at beginning of run

!     call flvisit_spk_beam_paraminit_send(th_beam,phi_beam,r_beam,rho_beam,u_beam)
  endif


  if (me==0) then
     call flvisit_nbody2_check_connection(lvisit_active)

     ! Fetch real-time, user-specified control parameters

     if (lvisit_active /= 0) then 
        call flvisit_nbody2_steering_recv( dsteer1,dsteer2,dsteer3, dsteer4,isteer1,isteer2,isteer3,isteer4)
	th_beam = dsteer1
	u_beam = dsteer2
	r_beam = dsteer4
	rho_beam = dsteer3
	ivis=max(isteer2,2)
	ivis_fields=max(isteer3,2)
        if (beam_debug) then
	   write(*,*) 'VISNB | steer1 = ',dsteer1
	   write(*,*) 'VISNB | steer2 = ',dsteer2
	   write(*,*) 'VISNB | steer3 = ',dsteer3
	   write(*,*) 'VISNB | steer4 = ',dsteer4
	   write(*,*) 'VISNB | switches = ',isteer1,isteer2,isteer3,isteer4
	endif
     else
        write(*,*) 'VISNB | No connection to visualization'
        return
     endif
     
  endif

#endif

  ! Broadcast beam parameters to all other PEs
  call MPI_BARRIER( MPI_COMM_WORLD, ierr)   ! Synchronize first
  call MPI_BCAST( lvisit_active, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
  if (lvisit_active /= 0) then
     call MPI_BCAST( th_beam, 1, MPI_REAL, 0, MPI_COMM_WORLD,ierr)
     call MPI_BCAST( phi_beam, 1, MPI_REAL, 0, MPI_COMM_WORLD,ierr)
     call MPI_BCAST( r_beam, 1, MPI_REAL, 0, MPI_COMM_WORLD,ierr)
     call MPI_BCAST( rho_beam, 1, MPI_REAL, 0, MPI_COMM_WORLD,ierr)
     call MPI_BCAST( u_beam, 1, MPI_REAL, 0, MPI_COMM_WORLD,ierr)
     call MPI_BCAST( ivis, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
     call MPI_BCAST( ivis_fields, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
  else
     if (me==0) write(*,*) ' No Connection to Visualization'
     return
  endif

!  if (rho_beam==0 )then
!     if (me==0) write(*,*) ' Switching off beam'
!
!     return
!  endif

  if (beam_config ==7) then
	vosc = u_beam

  else if (beam_config >= 3 .and. beam_config <=6 ) then
     ! laser standing wave or pond bullet

     !     u_beam = max(abs(u_beam),0.1)
     !     rho_beam = max(abs(rho_beam),0.5)
     !     r_beam = max(abs(r_beam),0.1)
     
     if (u_beam>0 ) then
      if ((u_beam-vosc)/u_beam <5.0 ) then
        vosc=u_beam ! limit amplitude change but allow switch-off
        if (me==0 .and. vosc /= vosc_old) write(*,*) 'VISNB | Laser amplitude changed'
      else
        if (me==0) write(*,*) 'VISNB | Amplitude change too big - readjust!'
      endif
     else
	vosc=0.
        if (me==0) write(*,*) 'VISNB | Laser switched off'
     endif

     if (sigma/rho_beam > 0.1 .and. sigma/rho_beam < 10.0) then
        sigma=rho_beam
        if (me==0 .and. sigma /= sigma_old) write(*,*) 'VISNB | Laser spot size changed'
     else
        if (me==0) write(*,*) 'VISNB | Spot size change too big - readjust!'
     endif
     if (tpulse/r_beam > 0.1 .and. tpulse/r_beam < 10.0) then
        tpulse=r_beam
     else
        if (me==0) write(*,*) 'VISNB | Pulse length change too big - readjust!'
        if (me==0 .and. tpulse /= tpulse_old) write(*,*) 'VISNB | Laser pulse length changed'
     endif

     if (th_beam /= theta_beam) then
        theta_beam=th_beam
        if (me==0 .and. theta_beam /= theta_old) write(*,*) 'VISNB | Incidence angle changed'
     endif

  else if (beam_config ==2) then
     nb_pe = np_beam_dt/num_pe  ! # beam particles to load per PE

     if (me ==0 .and. beam_debug) then
        write(*,*) 'VISNB | beam parts ',nb_pe,' theta ',theta_beam,' phi ', phi_beam
        write(*,*) 'VISNB | r ',r_beam,' rho ',rho_beam,' u ',u_beam, ' dt ',np_beam_dt
     endif

     Volb = pi*r_beam**2*x_beam       !  beam cylinder volume:  r_beam is radius
     qeb = Volb*rho_beam/np_beam    ! charge
     dpx=x_beam/nb_pe           ! x-axis spacing
     iseed2 = -131 - np_beam -3*me ! seed
     iseed1 = -333 - np_beam -3*me

     if (qeb <0 ) then
        mass_beam = 1.  ! electrons
     else
        mass_beam = 1836.  ! protons
     endif



     ! particle beam initialised along x-axis and rotated by theta, phi
     ct = cos(theta_beam)
     st = sin(theta_beam)
     cp = cos(phi_beam)
     sp = sin(phi_beam)
     vz_beam = -u_beam*st
     vx_beam = u_beam*ct*cp
     vy_beam = u_beam*ct*sp


     i = 0

     do while (i < nb_pe)
        yt = r_beam*(2*rano(iseed2)-1.)          
        zt = r_beam*(2*rano(iseed1)-1.)  
        if (yt**2 + zt**2 <= r_beam**2 ) then
           i = i+1
           p = npp+i  ! put them after last particle on this PE
           xb=  dpx*i + dpx/num_pe*me
           yb = yt 
           zb = zt 

           ! Now rotate disc about Euler angles

           x(p) =  xb*ct*cp - yb*ct*sp + zb*st 
           y(p) =  xb*sp    + yb*cp             
           z(p) =  -xb*st*cp + yb*st*sp + zb*ct 

           ! Add starting point
           x(p) = x(p) + start_beam
           y(p) = y(p) + yl/2.
           z(p) = z(p) + zl/2.

           q(p) = qeb
           m(p) = abs(qeb)*mass_beam
           ux(p)=vx_beam
           uy(p)=vy_beam
           uz(p)=vz_beam
           pepid(p) = me                ! processor ID
           pelabel(p) = npart+me*nb_pe+i  ! labels
           Ex(p) = 0.
           Ey(p) = 0.
           Ez(p) = 0.
           Bx(p) = 0.
           By(p) = 0.
           Bz(p) = 0.
           Ax(p) = 0.
           Ay(p) = 0.
           Az(p) = 0.
           pot(p) =0.
           work(p) =0.
        endif
     end do

     ! Augment total # particles - have to limit increase to max array size
     npp = npp+nb_pe
     np_beam = np_beam + np_beam_dt
     npart = npart + np_beam_dt


  else if (beam_config ==8) then
     ! Dust particle - # beam particles constant; infinite mass

     if (me==0 .and. u_beam /= u_old) write(*,*) 'VISNB | Beam velocity changed'

     if (me ==0 .and. beam_debug) then
        write(*,*) 'VISNB |  theta ',theta_beam,' phi ', phi_beam,' u ',u_beam
     endif
     ! dust particle velocity rotated by theta, phi
     ! beam particles could be sitting anywhere now
     Volb = 4*pi/3.*r_beam**3
     qeb = Volb*rho_beam/np_beam    ! new charge
     ct = cos(theta_beam)
     st = sin(theta_beam)
     cp = cos(phi_beam)
     sp = sin(phi_beam)
     vy_beam = u_beam*st*cp*vte*10   ! Scale by thermal velocity
     vx_beam = u_beam*st*sp*vte*10
     vz_beam = u_beam*ct*vte*10
     do i=1,npp
        if (pelabel(i)>ne+ni) then
           ux(i) = vx_beam
           uy(i) = vy_beam
           uz(i) = vz_beam
	   q(i) = qeb
        endif

     end do

  else if (scheme == 5) then
     ! ion crystal eqm  mode:
     !  r_beam is mean ion spacing
     !  u_beam is ion temperature (eV)
     !  rho_beam is dt 
     Ti_kev = u_beam
     eps = r_beam
     dt = rho_beam 
     if (me==0) write(*,*) 'VISNB | Steering pars: eps=',eps,' Ti=',Ti_kev,' dt=',dt
  else if (scheme==4) then
     ! Electron temp clamp
     Te_kev = u_beam
  endif


end subroutine vis_control