! ==============================================
!
!                BEAM_CONTROL
!
!  Interactive control of particle source
!
! ==============================================

subroutine beam_control

  use treevars
  use utils
  implicit none
  integer :: i, p, iseed1, iseed2
  real :: Volb, dpx, yt, zt
  integer :: lvisit_active
  real :: ct, st, cp, sp, vx_beam, vy_beam, vz_beam, xb, yb, zb
  logical :: beam_on = .true.
  logical :: beam_debug = .true.

  integer, save :: np_beam_dt  ! current # beam particles


  if (itime == 0) then
     np_beam_dt = np_beam  ! initial # beam particles to introduce per timestep
  endif

  if (.not. beam_on) return

  if (npart + np_beam_dt > npartm .or. max_list_length > nintmax-20) then
     if (me==0) write(*,*) 'Array limit reached: switching off beam'
     beam_on = .false.
     return
  endif

  if (itime == 0 .and. me==0 )  then

     ! Specify default parameters at beginning of run
     call flvisit_spk_check_connection(lvisit_active)
     call flvisit_spk_beam_paraminit_send(theta_beam,phi_beam,r_beam,rho_beam,u_beam)
  endif



  if (me==0) then
     call flvisit_spk_check_connection(lvisit_active)

     ! Fetch real-time, user-specified beam parameters
     if (lvisit_active /= 0) call flvisit_spk_beam_param_recv( theta_beam,phi_beam,r_beam,rho_beam,u_beam)
  endif


  ! Broadcast beam parameters to all other PEs
  call MPI_BARRIER( MPI_COMM_WORLD, ierr)   ! Synchronize first

  call MPI_BCAST( theta_beam, one, MPI_REAL8, root, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( phi_beam, one, MPI_REAL8, root, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( r_beam, one, MPI_REAL8, root, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( rho_beam, one, MPI_REAL8, root, MPI_COMM_WORLD,ierr)
  call MPI_BCAST( u_beam, one, MPI_REAL8, root, MPI_COMM_WORLD,ierr)

  call MPI_BCAST( lvisit_active, one, MPI_INTEGER8, root, MPI_COMM_WORLD,ierr)

  if (lvisit_active==0 )then
     if (me==0) write(*,*) ' No Connection to Visualization'
     !     return
  endif

  if (rho_beam==0 )then
     if (me==0) write(*,*) ' Switching off beam'

     return
  endif

  if (beam_config == 4) then
 ! laser standing wave
!     vosc = rho_beam
 !    sigma = r_beam

  else if (beam_config ==2) then
     nb_pe = np_beam_dt/num_pe  ! # beam particles to load per PE

     if (me ==0 .and. beam_debug) then
        write(*,*) 'beam parts ',nb_pe,' theta ',theta_beam,' phi ', phi_beam
        write(*,*) 'r ',r_beam,' rho ',rho_beam,' u ',u_beam, ' dt ',np_beam_dt
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
           ax(p) = 0.
           ay(p) = 0.
           az(p) = 0.
        endif
     end do

     ! Augment total # particles - have to limit increase to max array size
     npp = npp+nb_pe
     np_beam = np_beam + np_beam_dt
     npart = npart + np_beam_dt
  else
  endif
end subroutine beam_control
