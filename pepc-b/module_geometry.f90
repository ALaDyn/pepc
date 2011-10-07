!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>  Encapsulates all particle configuration routines
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_geometry

      implicit none
      save
      private

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	public plasma_start	!< Plasma target geometries
	public special_start	!< User-defined configs
	public add_electrons	!< Adds electrons to ion lattice
	public add_ramp		!< Add density ramp to target
	public cluster_sa	!< Stretch homogeneous sphere into 'Andreev' profile
	public stretch_sphere	!< Convert uniform sphere to tapered density profile
	public double_target	!< Add second target shifted by displacement vector displace(1:3)
!	public mc_config	!< OBSOLETE:  Do Monte-Carlo initialisation of particle positions
	public cutvector	!< Cuts a vector on the surface face_nr
	public face	 	!< Function to set up bounded geometry targets in conjunction with plasma_start


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      contains

! ==============================================
!
!                PLASMA_START
!
!   $Revision: 1.5  (Previously: randion)
!
!  Sets up 3D random distribution of plasma particles
!  according to target_geometry
! 
! ==============================================

subroutine plasma_start(my_rank, i1, n, nglobal, label_off, target_geometry, idim, &
     rho0, Zion, mass_ratio, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre, &
     number_faces, Vplas, Aplas, Qplas, q_part, m_part, a_ave )

  use module_particle_props
  use module_utilities
  use module_io

  implicit none

  integer, intent(in) :: my_rank ! rank
  integer, intent(in) :: i1 ! particle index start
  integer, intent(in) :: n  ! # particles to set up
  integer, intent(in) :: nglobal  ! global # particles of this species (determining charge weight)
  integer, intent(in) :: label_off ! offset for pelabel
  integer, intent(in) :: target_geometry ! geometry

  integer, intent(in) :: idim ! # dimensions
  real, intent(in) :: rho0 ! plasma density
  real, intent(in) :: Zion ! particle atomic charge (-1, Zion)

  real, intent(in) :: mass_ratio ! mass ratio (m_q/m_e)
  real, intent(in) :: x_plasma, y_plasma, z_plasma, r_sphere ! container params
  real, intent(in) :: plasma_centre(3) ! target offset vector
  real, intent(out) :: Vplas, Aplas, Qplas ! Plasma volume, area, charge
  real, intent(out) :: q_part, m_part, a_ave ! Derived particle charge, mass, spacing
  integer, intent(out) :: number_faces   ! # faces of target container

  integer              :: i, face_nr, n1
  integer              :: iseed1, iseed3, p, k

  real*8               :: xt, yt, zt
  real*8 :: c_status
  real :: pi=3.141592654
  logical :: start_debug=.false.

  iseed1 =  -99 - label_off      ! Select seed depending on PE and offset

! TODO  - need to pass ipefile for debug

  if (start_debug) then
    write (file_ipefile, '(a,3i8)') 'Seed: ', iseed1
    write(file_ipefile,'(a,i5/6(a20,i8/),6(a20,f12.3/),a20,3f12.3)') "Call parameters on PE ",my_rank, &
	"start index:",i1, &
	"# parts:",n, &
	"label offset:",label_off, &
	"geometry:",target_geometry, &
	"dimensions:",idim, &
	"density:",rho0, &
	"Z:",zion, &
	"mass ratio:",mass_ratio, &
	"xplas:",x_plasma,&
	"yplas:",y_plasma,&
	"zplas:",z_plasma,&
	"radius:",r_sphere,&
	"centre:",plasma_centre
 endif

  !  Predefined geometries:
  !  ---------------------
  !       0 = random slab
  !       1 = random sphere
  !       2 = random disc
  !       3 = random wire
  !       4 = ellipsoid
  !       5 = wedge
  !       6, 26 = hemisphere
  !       11 = hollow sphere
  !       16, 36 = hollow hemisphere

  !  2D geometries allowed in x-y plane: 
  !    0 = slab
  !    3 = disc

  !  Define target container parameters: volume, area, centre, # faces

  geom_container: select case(target_geometry)

  case(0) ! slab
     Vplas = x_plasma * y_plasma * z_plasma  ! plasma volume
     Aplas = x_plasma * y_plasma ! plasma area
     if (idim.eq.2) then
	number_faces = 4 ! In 2D, just check particles inside region xplasma*yplasma
     else 
     	number_faces = 6
     endif

  case(1) ! sphere
     Vplas = 4 * pi * r_sphere**3 / 3.
     Aplas = pi*r_sphere**2
     number_faces = 1

  case(11) ! hollow sphere
     Vplas = (4 * pi / 3.) * (r_sphere**3 - (r_sphere - x_plasma)**3)
     Aplas = pi*(r_sphere**2-(r_sphere-x_plasma)**2)
     number_faces = 2

  case(2) ! disc
     Vplas = pi * r_sphere**2 * x_plasma
     Aplas = pi*r_sphere**2
     number_faces = 3

  case(12) ! hollow tube
     Vplas = pi * (r_sphere**2 - (r_sphere-z_plasma)**2) * x_plasma
     Aplas = pi*(r_sphere**2 - (r_sphere-z_plasma)**2)
     number_faces = 4

  case(3) ! wire
     Vplas = pi * r_sphere**2 * z_plasma
     Aplas = pi*r_sphere**2
     if (idim.eq.2) then
	number_faces = 1 ! In 2D, just check particles inside circle radius r_sphere
     else 
	number_faces = 3 ! Otherwise check z-faces too
     endif

  case(13) ! hollow wire
     Vplas = pi * (r_sphere**2 - (r_sphere-x_plasma)**2) * z_plasma
     Aplas = pi*(r_sphere**2 - (r_sphere-x_plasma)**2)
     number_faces = 4

  case(4) ! ellipsoid
     Vplas = 4 * pi * x_plasma * y_plasma * z_plasma / 3.
     Aplas = pi*x_plasma*y_plasma*2
     number_faces = 1

  case(5) ! wedge
     Vplas = .5 * x_plasma * y_plasma * z_plasma
     Aplas = .5*x_plasma*y_plasma
     number_faces = 5

  case(6,26) ! hemisphere
     Vplas = 4 * pi * r_sphere**3 / 6.
     Aplas = pi*r_sphere**2/2.
     number_faces = 2


  case(16,36) ! hollow hemisphere
     Vplas = (4 * pi / 6.) * (r_sphere**3 - (r_sphere - x_plasma)**3)
     Aplas = pi/2.*(r_sphere**2-(r_sphere-x_plasma)**2)
     number_faces = 3

  case default ! slab
     Vplas = x_plasma * y_plasma * z_plasma  ! plasma volume
     Aplas = x_plasma * y_plasma ! plasma area
     number_faces = 6

  end select geom_container

! ensure soft exit if particle # = zero
  if (n==0) then
     q_part=1.
     m_part=1.
     a_ave=1.
     return
  endif

  !  Initialise particles according to target geometry


 i=0 

  do while (i < n)

     geom_parts: select case (target_geometry)
     case(0, 5) ! slab or wedge 
        xt = .5 * x_plasma * (2 * rano(iseed1) - 1.) + plasma_centre(1)
        yt = .5 * y_plasma * (2 * rano(iseed1) - 1.) + plasma_centre(2)         
        zt = plasma_centre(3)
        if (idim.eq.3) zt = zt + .5 * z_plasma * (2 * rano(iseed1) - 1.) 

     case(1, 11) ! sphere and hollow sphere
        xt = r_sphere * (2 * rano(iseed1) - 1.) + plasma_centre(1)
        yt = r_sphere * (2 * rano(iseed1) - 1.) + plasma_centre(2)        
        zt = r_sphere * (2 * rano(iseed1) - 1.) + plasma_centre(3)

     case(2,12) ! disc, tube
        xt = .5 * x_plasma * (2 * rano(iseed1) - 1.) + plasma_centre(1)
        yt = r_sphere * (2 * rano(iseed1) - 1.) + plasma_centre(2)         
        zt = r_sphere * (2 * rano(iseed1) - 1.) + plasma_centre(3)

     case(3,13) ! wire
        xt = r_sphere * (2 * rano(iseed1) - 1.) + plasma_centre(1)         
        yt = r_sphere * (2 * rano(iseed1) - 1.) + plasma_centre(2)
        zt =  plasma_centre(3)
	if (idim.eq.3) zt=zt + .5 * z_plasma * (2 * rano(iseed1) - 1.) ! Add z-component in 3D

     case(4) ! ellipsoid
        xt = r_sphere * (2 * rano(iseed1) - 1.) * x_plasma + plasma_centre(1)
        yt = r_sphere * (2 * rano(iseed1) - 1.) * y_plasma + plasma_centre(2)
        zt = r_sphere * (2 * rano(iseed1) - 1.) * z_plasma + plasma_centre(3)

     case(6, 16, 26, 36) ! hemisphere and hollow hemisphere
        xt = .5 * r_sphere * (2 * rano(iseed1) - 1.) + plasma_centre(1)
        yt = r_sphere * (2 * rano(iseed1) - 1.) + plasma_centre(2)
        zt = r_sphere * (2 * rano(iseed1) - 1.) + plasma_centre(3)

     end select geom_parts

     ! check if the new particle is inside container
     ! - if so then add it

     do face_nr = 1, number_faces
        call face(  (/xt, yt, zt/), c_status, face_nr, &
	            target_geometry, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre )

        if (c_status .ge. 0.) exit ! particle outside
     end do

     if (c_status .lt. 0) then
	i=i+1
        p=i+i1-1
        z(p) = zt 
        y(p) = yt
        x(p) = xt
     end if
  end do


  ! scramble indices to remove correlations
  iseed3 = -17 - 4 * my_rank -i1
  n1 = n

  !  exclude odd one out
  if (mod(n,2) .ne. 0) then
     n1 = n1 - 1
  end if

  do i = 1, n1
     p = i + i1-1
     k = n1 * rano(iseed3) + i1
     xt = x(p)
     yt = y(p)
     zt = z(p)
     x(p) = x(k)
     y(p) = y(k)
     z(p) = z(k)
     x(k) = xt
     y(k) = yt
     z(k) = zt
  end do


  ! Derive electron/ion charges and masses from target parameters

  if (idim==3 .or. idim==1) then    ! 3D

     Qplas = Vplas*rho0      ! total target charge
     q_part = Qplas/nglobal   ! particle charge
     m_part = mass_ratio*abs(q_part) ! particle mass
     a_ave = (Vplas/max(1,nglobal))**(1./3.)  ! ave. spacing

  else        ! 2D - use areal density instead of volume

     Qplas = Aplas*rho0
     q_part = Qplas/nglobal   ! particle charge
     m_part = mass_ratio*abs(q_part) ! particle mass
     a_ave = (Aplas/max(1,nglobal))**(1./2.) ! ave spacing
  endif


  q(i1:i1+n-1)       = q_part     ! charge
  m(i1:i1+n-1)       = m_part     ! mass
  pepid(i1:i1+n-1)   = my_rank    ! processor ID
  pelabel(i1:i1+n-1) = label_off + (/(i, i = 1, n)/) ! unique labels



!  if (start_debug) then
!	write(ipefile,'(a/(i6,8f15.5,i6))') "Initial particle positions, velocities:", &
!	  (i,x(i),y(i),z(i),ux(i),uy(i),uz(i),q(i),m(i),pelabel(i),i=i1,i1+n-1)
!  endif
end subroutine plasma_start


! ==============================================
!
!                SPECIAL_START
!
!  Initialise set of particles for special test configs
!
! ==============================================

subroutine special_start(iconf)

  use module_particle_props
  use module_physvars
  use module_utilities
  use module_velocity_setup
  implicit none

  integer, intent(in) :: iconf
  integer :: i,iseed
  real :: rs,v0,gamma0,vt,xt,yt,zt,thetvel,phivel
  real :: dpx, vx_beam, vy_beam, vz_beam, st, ct, sp, cp, volb
  integer :: npb
  real :: gamma, Ndebye

  iseed = -17-my_rank

  ! Define default container parameters in absence of plasma target
  Vplas = x_plasma * y_plasma * z_plasma
  Aplas = x_plasma * y_plasma
  focus = (/ xl/4., yl/2., zl/2. /) ! Centre of laser focal spot
  plasma_centre =  (/ xl/2., yl/2., zl/2. /) ! Centre of plasma
  Qplas = ne

  config: select case(iconf)

  case(1)        ! electron disc at x=0 with 2pi phase spread
     qe=-1.0
     mass_e=1.0
     qi=1.0
     mass_i=mass_ratio

     r_beam=sigma/2.
     gamma0=sqrt(1+vosc**2/2.)

     do while (i < nep)
        yt = r_beam*(2*rano(iseed)-1.)

        zt = r_beam*(2*rano(iseed)-1.)
        if (yt**2 + zt**2 <= r_beam**2 ) then
           i = i+1
           x(i) = 2*pi*rano(iseed)
           y(i) = yt + focus(2)
           z(i) = zt + focus(3)
           ux(i) = 0.
           uy(i) = 0.
           uz(i) = 0.

        endif
     end do
     q(1:nep) = qe           ! plasma electrons
     m(1:nep) = mass_e            ! electron mass

  case(2)        ! electron disc at laser focus with 2pi phase spread

     qe=-1.0
     mass_e=1.0
     qi=1.0
     mass_i=mass_ratio

     r_beam=sigma/2.
     gamma0=sqrt(1+vosc**2/2.)


     x_crit = focus(1)
     do while (i < nep)
        yt = r_beam*(2*rano(iseed)-1.)

        zt = r_beam*(2*rano(iseed)-1.)
        if (yt**2 + zt**2 <= r_beam**2 ) then
           i = i+1
           x(i) = focus(1)
           y(i) = yt + focus(2)
           z(i) = zt + focus(3)
           ux(i) = 0.
           uy(i) = 0.
           uz(i) = 0.

        endif
     end do
     q(1:nep) = qe           ! plasma electrons
     m(1:nep) = mass_e            ! electron mass


  case(3)       ! random placement in box with vte in x,y plane for B test 

     qe=-1.0
     mass_e=1.0
     qi=1.0
     mass_i=mass_ratio

     do i=1,np_local
        if (i<=nep) then
           vt=vte
           q(i) = qe        ! plasma electrons
           m(i) = mass_e 
        else
           vt = vti
           q(i) = qi        ! plasma ions (need Z* here)
           m(i) = mass_i    ! ion mass

        endif
        xt =  xl * rano(iseed) 
        yt =  yl * rano(iseed)
!        zt =  zl * rano(iseed)
        zt = zl/2.-z_plasma/2.+z_plasma*i/1./np_local
        x(i) = xt
        y(i) = yt
        z(i) = zt
        rs = rano(iseed)
        v0 = vt*sqrt(-2.0*log(rs))
        thetvel = 2*pi*rano(iseed)
        phivel = pi*rano(iseed)
        ux(i) = v0*cos(thetvel)
        uy(i) = v0*sin(thetvel)
!        uz(i) = v0/10.
        uz(i) = v0*cos(phivel)
     end do

  case(4)      !  Special config for FMM comparison
     open(30,file='fmm_config.data')
     do i=1,np_local
        read(30,*) x(i), y(i), z(i), q(i)
        m(i) = 1.
        ux(i)=0.
        uy(i)=0.
        uz(i)=0.
     end do
     force_const=1.
     dt=0.
     eps=0.
     close(30)

  case(5)     ! beam initialised along x-axis and rotated by theta, phi


     write (*,*) 'Preparing particle beam'
     npb = max(ni,ne)  ! Take whichever species switched on
     np_beam=0 ! Discard normal beam mode (plasma+beam)
     np_local = max(nep,nip)

     ct = cos(theta_beam)
     st = sin(theta_beam)
     cp = cos(phi_beam)
     sp = sin(phi_beam)
     vz_beam = u_beam*st
     vx_beam = u_beam*ct*cp
     vy_beam = u_beam*ct*sp
     Volb = pi*r_beam**2*x_beam
     qe = Volb*rho_beam/npb    ! charge
     mass_e = mass_beam
     mass_i = mass_beam
     qi=-qe
     dpx=x_beam/npb           ! x-axis spacing
     i=0
     do while (i < np_local)
        yt = r_beam*(2*rano(iseed)-1.)
        zt = r_beam*(2*rano(iseed)-1.)

        if (yt**2 + zt**2 <= r_beam**2 ) then
           i = i+1
           x(i)= start_beam + x_beam*rano(iseed)
           y(i) = yt + focus(2)
           z(i) = zt + focus(3)
           ux(i) = vx_beam
           uy(i) = vy_beam
           uz(i) = vz_beam

        endif
     end do
     q(1:np_local) = qe           ! plasma electrons
     m(1:np_local) = abs(qe*mass_beam)  ! electron mass

     !     write (*,'(a30/(5f13.5))') 'Beam config:', &
     !          (x(i),y(i),z(i),q(i),m(i),i=1,npp)


  case(6)       ! random placement in box xl*xl*xl for periodic system 

     ! Use Debye units 
     ! - lengths normalized to lambda_De
     ! - velocities to v_te
     ! - time to wpe^-1
     ! Effective density adjusted via box length
     ! TODO:   gamma, Ndebye need to be made global (physvars.f90)

     yl=xl
     zl=xl
     qe=-1.0
     mass_e=1.0
     qi=1.0
     mass_i=mass_ratio
     !  Inter-electron spacing
     a_ee = xl/(4*pi/3.*ne)**(1./3.)
     !  Inter-ion spacing
     a_ii = a_ee*Zion**(1./3.) 

     !  Renormalise softening parameter
     eps = eps * a_ee

     !  Equivalent Gamma
     gamma= a_ee**2/3.
     !  # electrons in Debye sphere
     Ndebye = (3*gamma)**(-1.5)
     !  Adjust force constant
     force_const = 1./3./Ndebye

     write (*,*) "Ewald Setup:"
     write (*,'(5(a30,f15.5/))') "Particle spacing ",a_ee, &
          "Softening parameter:",eps, &
          "Gamma:",gamma, &
          "# electrons in Debye sphere:",Ndebye, &
          "force const:",force_const

     do i=1,np_local
        if (i<=nep) then
           q(i) = qe        ! plasma electrons
           m(i) = mass_e 
        else
           q(i) = qi        ! plasma ions (need Z* here)
           m(i) = mass_i    ! ion mass

        endif
        xt =  xl * rano(iseed) 
        yt =  xl * rano(iseed)
        zt =  xl * rano(iseed)
        x(i) = xt
        y(i) = yt
        z(i) = zt
     end do


     ! Setup 3v Maxwellian electrons
     ! TODO: max random seed processor-dependent

     call maxwell1(ux(1:nep),nep,vte)
     call maxwell1(uy(1:nep),nep,vte)
     call maxwell1(uz(1:nep),nep,vte)
     call scramble_v(ux(1:nep),uy(1:nep),uz(1:nep),nep)   ! remove x,y,z correlations
     ! Ions cold    
     call cold_start(nep+1,nip)


  end select config



  pelabel(1:nep)       = my_rank * nep + (/(i, i = 1, nep)/)      ! Electron labels
  pelabel(nep + 1:np_local) = ne + my_rank * nip + (/(i, i = 1, nip)/) ! Ion labels
  pepid(1:np_local) = my_rank                ! processor ID
  work(1:np_local) = 1.


end subroutine special_start


!  ===============================================================
!
!                        ADD_ELECTRONS
!
!
!   Add neutralising electrons after ions-only const-temp eqm phase
!
!  ===============================================================

subroutine add_electrons

  use module_physvars
  use module_particle_props
  use module_velocity_setup
  use module_utilities

  integer :: iseed1, npp
  real :: xt, yt, zt




  ! Double # particles/PE

  ne = ni
  nip = np_local
  nep = nip
  npart_total = ni+ne
  np_local = nep+nip
  npp = np_local

  !  now have nep=nip;  ions (1:nip); electrons  (nip+1:npp) 
  !  reverse of usual order, but get mixed up by sorting later

  q(nip+1:npp) = qe        ! plasma electrons

  m(nip+1:npp) = mass_e      ! electron mass

  pepid(nip+1:npp) = my_rank                ! processor ID


  pelabel(nip+1:npp) =  pelabel(1:nip)  ! Electron labels: 1->ne copied from ions
  pelabel(1:nip) = pelabel(1:nip) + ni  ! Augment ion labels: ne+1 -> npart

  ! zero accelerations - should really compute these for electrons
  Ex(nip+1:npp) = 0.
  Ey(nip+1:npp) = 0.
  Ez(nip+1:npp) = 0.
  pot(nip+1:npp) = 0.

  work(1:npp) = 1.   ! set work load balanced initially
  iseed1 = -7901-my_rank

  !  Place electrons within 1/10 of ave. ion spacing in vicinity of ions
  xt = .1*a_ii*(2*rano(iseed1)-1.)
  yt = .1*a_ii*(2*rano(iseed1)-1.)          
  zt = .1*a_ii*(2*rano(iseed1)-1.)  
  x(nip+1:npp) = x(1:nip)+xt
  y(nip+1:npp) = y(1:nip)+yt
  z(nip+1:npp) = z(1:nip)+zt

  !  Set up thermal distribution
  if (vte > 0) then
     call maxwell1(ux(nip+1:nip+nep),nep,vte)
     call maxwell1(uy(nip+1:nip+nep),nep,vte)
     call maxwell1(uz(nip+1:nip+nep),nep,vte)
     call scramble_v(ux(nip+1:nip+nep),uy(nip+1:nip+nep),uz(nip+1:nip+nep),nep)   ! remove x,y,z correlations
  else
     call cold_start(nip+1,nep)
  endif

  !  Set ion velocities to zero
  ux(1:nip) = 0.
  uy(1:nip) = 0.
  uz(1:nip) = 0.

end subroutine add_electrons


! ==============================================
!
!                ADD_RAMP
!
!   $Revision: 1.0
!
!  Add ramp region to front of profile
! 
! ==============================================

subroutine add_ramp(d_layer)

    use module_physvars
    use module_particle_props
    use module_utilities

    implicit none

    real, intent(in) :: d_layer
    integer              :: i
    integer              :: nramp
    real                 :: qtot, qramp, s, xi, xedge, rind

    
! Transform particles in layer at front of target
! to form exponential density ramp

    qtot = rho0*d_layer ! total line density
    s = lolam*(1-rho_min)  ! required layer thickness for stretching
    qramp = rho0*s
    nramp = npart_total*s/d_layer
    xedge = plasma_centre(1)-d_layer/2.  ! leading edge (slab & disc)

    do i=1,np_local
       if (x(i).le.xedge+s) then
          rind = (x(i)-xedge)/s  ! fractional index ~ i/nramp
          xi = -lolam*log(1. - rind*(1.-rho_min))
          x(i) = xedge+s - xi   ! transform coordinate out into vacuum
       endif
    end do



end subroutine add_ramp

! ==============================================
!
!               CLUSTER_SA 
!
!   $Revision: 1.0
!
!  Convert uniform sphere to tapered density profile
!  according to Sasha Andreev parametric form
! 
! ==============================================

subroutine cluster_sa(my_rank,num_pe, i1,nip,r0,n0,r_sphere,qi,Qplas,plasma_centre,miome)

    use module_particle_props
    use module_utilities

    implicit none

    integer, intent(in) :: my_rank, num_pe
    real, intent(in) :: r0  ! characteristic radius (cm)
    real, intent(in) :: n0  ! characteristic cluster density (cm**-3)
    real, intent(in) :: plasma_centre(3) !  sphere centre
    real, intent(in) :: r_sphere ! equivalent sphere radius for qi 
    real, intent(in) :: qi ! particle charge - already assigned in configure (plasma_start)
    real, intent(in) :: Qplas ! total plasma charge = 1
    real, intent(in) :: miome ! mass ratio

    integer, intent(in) :: i1  ! index offset
    integer, intent(in) :: nip ! # particles to set up
    real :: n0n ! normalised central density
    real :: n_max ! central density normalised to n0
    integer              :: i, j
    integer              :: p
    real*8 :: rt, xt, yt, zt, pi, xi
    real*8 :: r_c, Q_c, Q_s, dens, Qnorm, dens0, v_max
    integer :: nitot, N_c, N_s, npi_c, npi_s, nbin, np_rest, offset
    real*8 :: zeta_max, zeta, dzeta, t_start, Qt
    real*8 :: S, S_c, T, ch, sh2, th, ch2, rp, integ
    
! Andreev cluster:
! uniform up to r=0.1125 r0
! tapered from r=0.1125 r0 to 2.5 r0
    pi = 2.*asin(1.d0)
    t_start = 0.05  
    zeta_max = 0.23  ! start point r(0.23)=0.1125

! # ions in central region
!    r_c = 0.1125  ! centre radius normalised to r0 
    r_c = t_start*cosh(zeta_max)**2/(sinh(2*zeta_max)/2. + zeta_max)  ! Radius associated with zeta_max
    dens0 = 1./3./t_start**2*(sinh(2*zeta_max)/2. + zeta_max)**2/(cosh(zeta_max)**4*(1.-zeta_max*tanh(zeta_max)))
 
    v_max = sqrt(2./3.)/sqrt(dens0)/sqrt(miome)*r0

  ! # ions in central, self-sim regions

    n_max = n0*dens0
!    N_c = 4*pi/3.*n_max*(r_c*r0)**3 ! # ions in central portion
    Q_c = 0.041  ! Charge in centre (fixed to give correct field - 4pi/3 rho_m r_c^3; rho_m=6.917)
    N_c = Q_c/qi
    nitot = Qplas/qi  ! total ions (check)
    n0n = Q_c*3/4/pi/r_c**3 ! normalised central density

! # charge in outer region
    Q_s = Qplas - Q_c
    N_s = Q_s/qi
    Qnorm = 1.0 ! normalisation factor for charge integral

! # put remainder central particles on last cpu 
    if (my_rank==num_pe-1) then
      np_rest = mod(N_c,num_pe)
    else
      np_rest = 0
    endif
    npi_c = N_c/num_pe + np_rest
    npi_s = nip-npi_c
    offset = my_rank*(nip-N_c/num_pe)  ! shell offset same for all

    if (my_rank==0) then
      write(*,'(a30,3i12)') "Total ions, in centre/outside:",nitot,N_c,N_s
      write(*,'(a42,3i12)') "Local ions in centre/outside:",npi_c,npi_s
       write(*,'(a30,4(1pe12.4))') "Densities:", n0, n_max, dens0, n0n, Qnorm
     write(*,'(a30,3f12.4)') "Radii r_eff, r_c, r0:", r_sphere, r_c, r0
      write(*,'(a30,3f12.4)') "Charge Q_tot, Q_c, Q_s",Qplas, Q_c, Q_s
    endif

! # rescale radii of inner ions
    do i=1,npi_c
	p = i + i1-1  ! local index (including offset)
        j = my_rank*(npi_c-np_rest) + i   ! Unique global particle #
       	xt = x(p)-plasma_centre(1)
	yt = y(p)-plasma_centre(2)
	zt = z(p)-plasma_centre(3)
	xi = (1.*j/N_c)**(1./3.) !  inversion for uniform density
	rt = sqrt(xt**2+yt**2+zt**2)
        x(p) = plasma_centre(1) + xt/rt*r_c*xi  ! keep direction vector; scale by inner sphere radius 
	y(p) = plasma_centre(2) + yt/rt*r_c*xi
	z(p) = plasma_centre(3) + zt/rt*r_c*xi
! Scale velocities according to ss value at radius r_c 
!          ux(i) = tanh(zeta_max)*xt/rt*xi
!          uy(i) = tanh(zeta_max)*yt/rt*xi
!          uz(i) = tanh(zeta_max)*zt/rt*xi
	ux(i)=0.
	uy(i)=0.
	uz(i)=0.
!        write (ipefile,'(i6,a8,f15.5)') j,'r/r0, ',r_c*xi
     end do


! # Integrate outer zone numerically, according to Andreev parametric form
! # Each CPU will place particles in different segments of integrated shell
     nbin = 200*N_s
!     nbin = 50*
     dzeta = -zeta_max/nbin
     zeta=zeta_max
     Qt = 0.
     j = offset + 1  ! Global starting index of ions
     i = npi_c+i1 ! Local index
        if(my_rank==0) write (*,*) 'zeta_max, dzeta, nbin',zeta_max,dzeta,nbin

     S_c =  sinh(2*zeta_max)/2. + zeta_max  ! Sinh function for integral

!    do while (zeta>0.004 .and. j<=offset + npi_s)
    do while (j<=offset + npi_s)

! Integrand
        sh2 = sinh(2*zeta)
        th = tanh(zeta)
        ch2 = cosh(2*zeta)
        ch = cosh(zeta)
        S = sh2/2. + zeta
        T = 1.-zeta*th

	dens = 1./3./t_start**2*S**2/(ch**4*(1.-zeta*th))
        integ = -2*ch**2/S**2  ! Simplified integrand
!        integ =  (sh2*sz - ch**2*(ch2 + 1)) / (T*sz**2) 
!        Qt = Qt + integ*t_start*dzeta/Qnorm
        Qt = t_start*(1./S - 1./S_c)
        rp = t_start*ch**2/S ! Radius associated with zeta
!        if(me==0) write (*,'(8f12.5)') rp,zeta,integ,Qt

     if (Qt > j*qi) then
! place particle
          xt = x(i)-plasma_centre(1)  ! Get direction vector
          yt = y(i)-plasma_centre(2)
          zt = z(i)-plasma_centre(3)
          rt = sqrt(xt**2+yt**2+zt**2)
!         write (ipefile,'(i6,a23,4f15.5)') i,'r/r0, zeta, integ, Qt',rp,zeta,integ,Qt
!         if (me==0) write (*,'(2i6,a23,5(1pe15.7))') i,j,' r/r0, dens, integ, Qt/qi',rp,dens,Qt/qi

          x(i) = plasma_centre(1) + xt/rt*rp  ! scale by sphere radius 
          y(i) = plasma_centre(2) + yt/rt*rp
          z(i) = plasma_centre(3) + zt/rt*rp
! Scale velocities according to ss value at radius rp 
!          ux(i) = th*xt/rt
!          uy(i) = th*yt/rt
!          uz(i) = th*zt/rt
	ux(i)=0.
	uy(i)=0.
	uz(i)=0.
          j=j+1
          i=i+1
       endif
       zeta = zeta+dzeta
    end do



end subroutine cluster_sa 

! ==============================================
!
!               STRETCH_SPHERE 
!
!   $Revision: 1.0
!
!  Convert uniform sphere to tapered density profile
! 
! ==============================================

subroutine stretch_sphere(r0)

    use module_physvars
    use module_particle_props

    implicit none

    real, intent(in) :: r0 
    integer              :: i
    real                 :: xi
    real*8 :: rt, xt, yt, zt, s
    
! Transform particles in layer at front of target
! to form exponential density ramp


    do i=1,np_local
	xt = x(i)-plasma_centre(1)
	yt = y(i)-plasma_centre(2)
	zt = z(i)-plasma_centre(3)
	xi = pelabel(i)/1.001/npart_total ! impose cutoff
	rt = sqrt(xt**2+yt**2+zt**2)
! inverted  coord for 1/(1+s^3)^2 distribution
	s = (1./(1.-xi) - 1.)**(1./3.)
	x(i) = plasma_centre(1) + xt/rt*r0*s  ! keep direction vector; scale by sphere radius 
	y(i) = plasma_centre(2) + yt/rt*r0*s
	z(i) = plasma_centre(3) + zt/rt*r0*s
    end do

end subroutine stretch_sphere 

!  ===============================================================
!
!                       DOUBLE_TARGET 
!
!   $Revision: 206 $
!
!   Add secondary target shifted by displacement vector displace(1:3) 
!
!  ===============================================================


subroutine double_target

  use module_physvars
  use module_particle_props

  integer ::  iseed1, npp

  npp=np_local

! Make copy of target displaced by vector displace(1:3)

  x(npp+1:2*npp) = x(1:npp) + displace(1)
  y(npp+1:2*npp) = y(1:npp) + displace(2)
  z(npp+1:2*npp) = z(1:npp) + displace(3)

  ux(npp+1:2*npp) = ux(1:npp) ! same velocities      
  uy(npp+1:2*npp) = uy(1:npp) ! 
  uz(npp+1:2*npp) = uz(1:npp) !
 
  q(npp+1:2*npp) = q(1:npp) ! same charge      
  m(npp+1:2*npp) = m(1:npp)  ! same mass

  pepid(npp+1:2*npp) = my_rank                ! processor ID


  pelabel(npp+1:2*npp) =  pelabel(1:npp) + npart_total  ! labels: electrons_1, ions_1, electrons_2, ions_2 

  ! zero accelerations - should really compute these for electrons
  Ex(npp+1:2*npp) = 0.
  Ey(npp+1:2*npp) = 0.
  Ez(npp+1:2*npp) = 0.
  pot(npp+1:2*npp) = 0.


  work(1:2*npp) = 1.   ! set work load balanced initially
  iseed1 = -7901-my_rank



  ! Double # particles/PE
  nip = 2*nip
  nep = 2*nep
  ni = 2*ni
  ne = 2*ne
  npart_total = ni+ne
  np_local = 2*np_local

end subroutine double_target




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   FACE
!
!  Function to set up bounded geometry targets in conjunction with 
!  PLASMA_START. 
!  M. Hammes, U. Wuppertal (Guest student programme, 2002)
!
!
! Generalised description of convex faces described 
! by a function f(x, y, z) = 0. 
! Saves f(r) in c_status
! and means c_status > 0.: particle is out,
! c_status < 0.: particle is in.

subroutine face(r, c_status, face_nr, &
    target_geometry, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre )

    implicit none 

    real*8, intent(in), dimension(1:3)	    :: r
    real*8, intent(out)			    :: c_status
    integer, intent(in)                     :: face_nr
    integer, intent(in)	                    :: target_geometry
    real, intent(in) :: x_plasma, y_plasma, z_plasma, r_sphere
    real, intent(in) :: plasma_centre(3)

    ! vectors and other values for the slab
    real*8, dimension(1:3)	    :: r_diff
    real*8, dimension(1:3)            :: normal_vector, offset_vector
    real*8                            :: gamma

    ! matrix for the ellipsoid
    real, dimension(1:3, 1:3)       :: A

    select case (target_geometry)
    case (0) ! slab
        select case (face_nr)
        case(1) ! x-y-plane, negative z
            normal_vector = (/ 0., 0., +1. /)
            offset_vector = (/ 0., 0., -z_plasma/2. /) + plasma_centre
        case(2) ! x-y-plane, positive z
            normal_vector = (/ 0., 0., -1./)
            offset_vector = (/ 0., 0., +z_plasma/2. /) + plasma_centre
        case(3) ! x-z-plane, negative y
            normal_vector = (/ 0., +1., 0. /)
            offset_vector = (/ 0., -y_plasma/2., 0. /) + plasma_centre
        case(4) ! x-z-plane, positive y
            normal_vector = (/ 0., -1., 0. /)
            offset_vector = (/ 0., +y_plasma/2., 0. /) + plasma_centre
        case(5) ! y-z-plane, negative x
            normal_vector = (/ +1., 0., 0. /)
            offset_vector = (/ -x_plasma/2., 0., 0. /) + plasma_centre
        case(6) ! y-z-plane, positive x
            normal_vector = (/ -1., 0., 0. /)
            offset_vector = (/ +x_plasma/2., 0., 0. /) + plasma_centre
        end select

        c_status = dot_product(normal_vector, offset_vector - r)

    case(1) ! sphere
        r_diff = r - plasma_centre
        c_status = dot_product(r_diff, r_diff) - r_sphere**2

    case(11) ! hollow sphere
        r_diff = r - plasma_centre
        select case (face_nr)
        case(1) ! outer_sphere
            c_status = dot_product(r_diff, r_diff) - r_sphere**2
        case(2) ! inner_sphere
            c_status = (r_sphere - x_plasma)**2 - dot_product(r_diff, r_diff)
        end select

    case(2) ! disc
        select case (face_nr)
        case(1) ! the tube in x-direction
            r_diff = r - plasma_centre
            c_status = r_diff(2)**2 + r_diff(3)**2 - r_sphere**2
        case(2) ! y-z-plane, negative x
            normal_vector = (/ +1., 0., 0. /)
            offset_vector = (/ -x_plasma/2., 0., 0. /) + plasma_centre
            c_status = dot_product(normal_vector, offset_vector - r)
        case(3) ! y-z-plane, positive x
            normal_vector = (/ -1., 0., 0. /)
            offset_vector = (/ +x_plasma/2., 0., 0. /) + plasma_centre
            c_status = dot_product(normal_vector, offset_vector - r)
        end select

    case(12)  ! hollow disc/tube, thickness z_plasma
        select case (face_nr)
        case(1) ! the tube in x-direction
            r_diff = r - plasma_centre
            c_status = r_diff(2)**2 + r_diff(3)**2 - r_sphere**2
        case(2) ! the inner tube in z-direction
            r_diff = r - plasma_centre
            c_status =  (r_sphere-z_plasma)**2 - (r_diff(2)**2 + r_diff(3)**2)
        case(3) ! y-z-plane, negative x
            normal_vector = (/ +1., 0., 0. /)
            offset_vector = (/ -x_plasma/2., 0., 0. /) + plasma_centre
            c_status = dot_product(normal_vector, offset_vector - r)
        case(4) ! y-z-plane, positive x
            normal_vector = (/ -1., 0., 0. /)
            offset_vector = (/ +x_plasma/2., 0., 0. /) + plasma_centre
            c_status = dot_product(normal_vector, offset_vector - r)
        end select

    case(3)  ! wire
        select case (face_nr)
        case(1) ! the tube in z-direction
            r_diff = r - plasma_centre
            c_status = r_diff(1)**2 + r_diff(2)**2 - r_sphere**2
        case(2) ! x-y-plane, negative z
            normal_vector = (/ 0., 0., +1. /)
            offset_vector = (/ 0., 0., -z_plasma/2. /) + plasma_centre
            c_status = dot_product(normal_vector, offset_vector - r)
        case(3) ! x-y-plane, positive z
            normal_vector = (/ 0., 0., -1. /)
            offset_vector = (/ 0., 0., +z_plasma/2. /) + plasma_centre
            c_status = dot_product(normal_vector, offset_vector - r)
        end select

    case(13)  ! hollow wire
        select case (face_nr)
        case(1) ! the tube in z-direction
            r_diff = r - plasma_centre
            c_status = r_diff(1)**2 + r_diff(2)**2 - r_sphere**2
        case(2) ! the inner tube in z-direction
            r_diff = r - plasma_centre
            c_status =  (r_sphere-x_plasma)**2 - (r_diff(1)**2 + r_diff(2)**2)
        case(3) ! x-y-plane, negative z
            normal_vector = (/ 0., 0., +1. /)
            offset_vector = (/ 0., 0., -z_plasma/2. /) + plasma_centre
            c_status = dot_product(normal_vector, offset_vector - r)
        case(4) ! x-y-plane, positive z
            normal_vector = (/ 0., 0., -1. /)
            offset_vector = (/ 0., 0., +z_plasma/2. /) + plasma_centre
            c_status = dot_product(normal_vector, offset_vector - r)
        end select

    case(4) ! ellipsoid
        A = reshape((/ 1./x_plasma, 0., 0., 0., 1./y_plasma, 0., 0., 0., 1./z_plasma /), (/3, 3/))
        r_diff = matmul(A, r - plasma_centre)
        c_status = dot_product(r_diff, r_diff) - r_sphere**2

    case(5) ! prism in x-y plane
        gamma = atan(y_plasma / (2 * x_plasma))
        select case (face_nr)
        case(1) ! x-y-plane, positive z
            normal_vector = (/ 0., 0., -1. /)
            offset_vector = (/ 0., 0., +z_plasma/2. /) + plasma_centre
        case(2) ! x-y-plane, negative z
            normal_vector = (/ 0., 0., +1. /)
            offset_vector = (/ 0., 0., -z_plasma/2. /) + plasma_centre
        case(3) ! y-z-plane, negative x
            normal_vector = (/ +1., 0., 0. /)
            offset_vector = (/ -x_plasma/2., 0., 0. /) + plasma_centre
        case(4) ! negative y
            normal_vector = (/ -sin(gamma), +cos(gamma), 0.d0 /)
            offset_vector = (/ -x_plasma/2., -y_plasma/2., 0.0 /) + plasma_centre
        case(5) ! positive y
            normal_vector = (/ -sin(gamma), -cos(gamma), 0.d0 /)
            offset_vector = (/ -x_plasma/2., +y_plasma/2.,0. /) + plasma_centre
        end select
        c_status = dot_product(normal_vector, offset_vector - r)


    case(6) ! hemisphere: flat side facing laser
        select case (face_nr)
        case(1) ! hemisphere
            r_diff = r + (/ r_sphere/2., 0., 0. /) - plasma_centre
            c_status = dot_product(r_diff, r_diff) - r_sphere**2
        case(2) ! y-z-plane
            normal_vector = (/ +1., 0., 0. /)
            offset_vector = (/ -r_sphere/2., 0., 0. /) + plasma_centre
            c_status = dot_product(normal_vector, offset_vector - r)
        end select

    case(16) ! hollow hemisphere
        select case (face_nr)
        case(1) ! outer hemisphere
            r_diff = r + (/ r_sphere/2., 0., 0. /) - plasma_centre
            c_status = dot_product(r_diff, r_diff) - r_sphere**2
        case(2) ! inner hemisphere
            r_diff = r + (/ r_sphere/2., 0., 0. /) - plasma_centre
            c_status = (r_sphere - x_plasma)**2 - dot_product(r_diff, r_diff)
        case(3) ! y-z-plane
            normal_vector = (/ 1., 0., 0. /)
            offset_vector = (/-r_sphere/2., 0., 0./) + plasma_centre
            c_status = dot_product(normal_vector, offset_vector - r)
        end select

    case(26) ! flipped hemisphere: curved side facing laser
        select case (face_nr)
        case(1) ! hemisphere
            r_diff = r - (/ r_sphere/2., 0., 0. /) - plasma_centre
            c_status = dot_product(r_diff, r_diff) - r_sphere**2
        case(2) ! y-z-plane
            normal_vector = (/-1., 0., 0./)
            offset_vector = (/ +r_sphere/2., 0., 0./) + plasma_centre
            c_status = dot_product(normal_vector, offset_vector - r)
        end select


    case(36) ! flipped hollow hemisphere: curved surface facing laser
        select case (face_nr)
        case(1) ! outer hemisphere
            r_diff = r - (/r_sphere/2., 0., 0./) - plasma_centre
            c_status = dot_product(r_diff, r_diff) - r_sphere**2
        case(2) ! inner hemisphere
            r_diff = r - (/r_sphere/2., 0., 0./) - plasma_centre
            c_status = (r_sphere - x_plasma)**2 - dot_product(r_diff, r_diff)
        case(3) ! y-z-plane
            normal_vector = (/-1., 0., 0./)
            offset_vector = (/+r_sphere/2., 0., 0./) + plasma_centre
            c_status = dot_product(normal_vector, offset_vector - r)
        end select

    case default  ! revert to slab
        select case (face_nr)
        case(1) ! x-y-plane, negative z
            normal_vector = (/ 0., 0., +1. /)
            offset_vector = (/ 0., 0., -z_plasma/2. /) + plasma_centre
        case(2) ! x-y-plane, positive z
            normal_vector = (/ 0., 0., -1./)
            offset_vector = (/ 0., 0., +z_plasma/2. /) + plasma_centre
        case(3) ! x-z-plane, negative y
            normal_vector = (/ 0., +1., 0./)
            offset_vector = (/ 0., -y_plasma/2., 0. /) + plasma_centre
        case(4) ! x-z-plane, positive y
            normal_vector = (/ 0., -1., 0./)
            offset_vector = (/ 0., +y_plasma/2., 0. /) + plasma_centre
        case(5) ! y-z-plane, negative x
            normal_vector = (/ +1., 0., 0. /)
            offset_vector = (/ -x_plasma/2., 0., 0. /) + plasma_centre
        case(6) ! y-z-plane, positive x
            normal_vector = (/ -1., 0., 0. /)
            offset_vector = (/ +x_plasma/2., 0., 0. /) + plasma_centre
        end select

        c_status = dot_product(normal_vector, offset_vector - r)

    end select

end subroutine face

subroutine cutvector(r_in, face_nr, n, r_out, &
        geometry, x_plasma, y_plasma, z_plasma, r_sphere, centre )
    !
    ! cuts a vector on the face face_nr
    ! and optionally puts out thr normal vector
    !

    implicit none
    
    real*8, dimension(1:3), intent(in)  :: r_in
    real*8, dimension(1:3), intent(out) :: n, r_out
    real*8, dimension(1:3)              :: r_affin, temp
    integer, intent(in)               :: face_nr
    integer, intent(in)               :: geometry
    real, intent(in) :: x_plasma, y_plasma, z_plasma, r_sphere
    real, intent(in) ::  centre(3)

    real*8, dimension(1:3, 1:3)	      :: T
    real*8                              :: phi, psi, gamma

    r_affin = r_in - centre
    select case (geometry) 
    case (0) ! slab
        select case (face_nr)
        case (1, 2) ! x-y-planes
            n = (/0., 0., 1./)
            if (face_nr == 2) n = -n
            r_affin = .5 * z_plasma * r_affin / abs(dot_product(r_affin, n))
        case (3, 4) ! x-z-planes
            n = (/0., 1., 0./)
            if (face_nr == 4) n = -n
            r_affin = .5 * y_plasma * r_affin / abs(dot_product(r_affin, n))
        case (5, 6) ! y-z-planes
            n = (/1., 0., 0./)
            if (face_nr == 6) n = -n
            r_affin = .5 * x_plasma * r_affin / abs(dot_product(r_affin, n))
        end select
        r_out = r_affin + centre
        
    case (1) ! sphere
        n = - r_affin / sqrt(dot_product(r_affin, r_affin))
        r_out = r_sphere * r_affin / abs(dot_product(r_affin, n)) + centre
        
    case (2) ! disc
        select case (face_nr)
        case (1) ! zylinder
            n = (/0.d0, -r_affin(2), -r_affin(3)/)
            n = n / sqrt(dot_product(n, n))
            r_affin = r_sphere * r_affin / abs(dot_product(r_affin, n))
        case (2, 3) ! y-z-planes
            n = (/1., 0., 0./)
            if (face_nr == 3) n = -n
            r_affin = .5 * x_plasma * r_affin / abs(dot_product(r_affin, n))
        end select
        r_out = r_affin + centre
        
    case (3) ! wire
        select case (face_nr)
        case (1) ! zylinder
            n = (/-r_affin(1), -r_affin(2), 0.d0/)
            n = n / sqrt(dot_product(n, n))
            r_affin = r_sphere * r_affin / abs(dot_product(r_affin, n))
        case (2, 3) ! x-y-planes
            n = (/0., 0., 1./)
            if (face_nr == 3) n = -n
            r_affin = .5 * z_plasma * r_affin / abs(dot_product(r_affin, n))
        end select
        r_out = r_affin + centre
        
    case (4) ! ellipsoid
        T = reshape((/1., 0., 0., 0., x_plasma / y_plasma, 0., 0., 0., x_plasma / z_plasma/), (/3, 3/))
        temp = matmul(T, r_affin)
        r_out = r_sphere * x_plasma * r_affin / sqrt(dot_product(temp, temp))
        temp = r_sphere * x_plasma * temp / sqrt(dot_product(temp, temp))
        phi = atan(temp(1) / sqrt(dot_product(temp, temp) - temp(1)**2))
        psi = atan(temp(1) * sqrt(dot_product(r_out, r_out) - temp(1)**2) / (dot_product(temp, temp) - temp(1)**2))
        phi = phi - psi
        call rotate(phi, temp, n)
        n = -n / sqrt(dot_product(n, n))
        r_out = r_out + centre
        
    case (5) ! prism
        gamma = atan(y_plasma / (2 * x_plasma))
        select case (face_nr)
        case (1, 2) ! x-y-planes
            n = (/0., 0.,-1./)
            if (face_nr == 2) n = -n
            r_affin = .5 * z_plasma * r_affin / abs(dot_product(r_affin, n))
        case (3) ! y-z-plane
            n = (/1., 0., 0./)
            r_affin = .5 * x_plasma * r_affin / abs(dot_product(r_affin, n))
        case (4, 5) ! the other planes
            n = (/-sin(gamma), +cos(gamma),0.d0/)
            r_affin(2) = tan(gamma) * r_affin(1) - y_plasma / 4.
            if (face_nr == 5) then
                n(2) = -n(2)
                r_affin(2) = -r_affin(2)
            end if
        end select
        r_out = r_affin + centre
        
    case (6) ! semisphere
        temp = r_affin - (/-r_sphere / 2.d0, 0.d0, 0.d0/)
        select case (face_nr)
        case (1) ! semisphere
            n = -temp
            n = n / sqrt(dot_product(n, n))
            r_affin = r_sphere * temp / abs(dot_product(temp, n)) + (/-r_sphere / 2.d0, 0.d0, 0.d0/)
        case (2) ! plane
            n = (/-1., 0., 0./)
            r_affin = .5 * r_sphere * r_affin / abs(dot_product(r_affin, n))
        end select
        r_out = r_affin + centre
        
    case (7) ! hollow ball
        select case (face_nr)
        case (1) ! outer sphere
            n = -r_affin
            n = n / sqrt(dot_product(n, n))
            r_affin = r_sphere * r_affin / abs(dot_product(r_affin, n))
        case (2) ! inner sphere
            n = r_affin
            n = n / sqrt(dot_product(n, n))
            r_affin = (r_sphere - x_plasma) * r_affin / abs(dot_product(r_affin, n))
        end select
        r_out = r_affin + centre
        
    case (8) ! hollow semisphere
        temp = r_affin - (/-r_sphere / 2.d0, 0.d0, 0.d0/)
        select case (face_nr)
        case (1) ! outer semisphere
            n = -temp
            n = n / sqrt(dot_product(n, n))
            r_affin = r_sphere * temp / abs(dot_product(temp, n)) + (/-r_sphere / 2., 0., 0./)
        case (2) ! inner semisphere
            n = temp
            n = n / sqrt(dot_product(n, n))
            r_affin = (r_sphere - x_plasma) * temp / abs(dot_product(temp, n)) + (/-r_sphere / 2., 0., 0./)
        case (3) ! y-z-plane
            n = (/-1., 0., 0./)
            r_affin = .5 * r_sphere * r_affin / abs(dot_product(r_affin, n))
        end select
        r_out = r_affin + centre
        
    end select
    
contains

    !
    ! Rotates a vector around the y-axis
    ! with an angle phi (needed for the 
    ! ellipsoid)
    !
    subroutine rotate(phi, v_in, v_out)
        real*8				:: phi
        real*8, dimension(1:3)		:: v_in, v_out
        real*8				:: c, s
        real*8, dimension(1:3, 1:3)	:: R

        c = cos(phi)
        s = sin(phi)

        R = reshape((/1.d0, 0.d0, 0.d0, 0.d0, c, s, 0.d0, -s, c/), (/3, 3/))
        v_out = matmul(R, v_in)
    end subroutine rotate

end subroutine cutvector
end module module_geometry
