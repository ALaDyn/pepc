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

subroutine plasma_start(i1, n, nglobal, label_off, target_geometry, velocity_config, idim, &
     rho0, Zion, mass_ratio, vt, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre, &
     number_faces, Vplas, Aplas, Qplas, q_part, m_part, a_ave )

  use treevars
  use utils

  implicit none

  integer, intent(in) :: i1 ! particle index start
  integer, intent(in) :: n  ! # particles to set up
  integer, intent(in) :: nglobal  ! global # particles of this species (determining charge weight)
  integer, intent(in) :: label_off ! offset for pelabel
  integer, intent(in) :: target_geometry ! geometry
  integer, intent(in) :: velocity_config ! Maxwellian/cold start switch
  integer, intent(in) :: idim ! # dimensions
  real, intent(in) :: rho0 ! plasma density
  real, intent(in) :: Zion ! particle atomic charge (-1, Zion)
  real, intent(in) :: vt ! thermal or drift velocity
  real, intent(in) :: mass_ratio ! mass ratio (m_q/m_e)
  real, intent(in) :: x_plasma, y_plasma, z_plasma, r_sphere ! container params
  real, intent(in) :: plasma_centre(3) ! target offset vector
  real, intent(out) :: Vplas, Aplas, Qplas ! Plasma volume, area, charge
  real, intent(out) :: q_part, m_part, a_ave ! Derived particle charge, mass, spacing
  integer, intent(out) :: number_faces   ! # faces of target container

  integer              :: i, face_nr, n1
  integer              :: iseed1, iseed3, p, k

  real*8               :: xt, yt, zt, rt
  real*8 :: c_status
  real :: pi=3.141592654
  logical :: start_debug=.false.

  iseed1 =  -99 - label_off      ! Select seed depending on PE and offset

  if (start_debug) then
    write (ipefile, '(a,3i8)') 'Seed: ', iseed1
    write(ipefile,'(a,i5/6(a20,i8/),8(a20,f12.3/),a20,3f12.3)') "Call parameters on PE ",me, &
	"start index:",i1, &
	"# parts:",n, &
	"label offset:",label_off, &
	"geometry:",target_geometry, &
	"vel. config:",velocity_config,&
	"dimensions:",idim, &
	"density:",rho0, &
	"Z:",zion, &
	"vthermal:",vt,&
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


  !  Define target container parameters: volume, area, centre, # faces

  geom_container: select case(target_geometry)

  case(0) ! slab
     Vplas = x_plasma * y_plasma * z_plasma  ! plasma volume
     Aplas = x_plasma * y_plasma ! plasma area
     number_faces = 6

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
     Aplas = x_plasma*y_plasma
     number_faces = 3

  case(12) ! hollow tube
     Vplas = pi * (r_sphere**2 - (r_sphere-z_plasma)**2) * x_plasma
     Aplas = pi*(r_sphere**2 - (r_sphere-z_plasma)**2)
     number_faces = 4

  case(3) ! wire
     Vplas = pi * r_sphere**2 * z_plasma
     Aplas = pi*r_sphere**2
     number_faces = 3

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
        zt = .5 * z_plasma * (2 * rano(iseed1) - 1.) + plasma_centre(3)

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
        zt = .5 * z_plasma * (2 * rano(iseed1) - 1.) + plasma_centre(3)

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
  iseed3 = -17 - 4 * me -i1
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
  pepid(i1:i1+n-1)   = me    ! processor ID
  pelabel(i1:i1+n-1) = label_off + (/(i, i = 1, n)/) ! unique labels



  velocities: select case (velocity_config)

  case(0) ! Default is cold start
     call cold_start(i1,n)

  case(1) ! isotropic Maxwellian

     if (vt > 0) then
        call maxwell1(ux,nppm,i1,n,vt)
        call maxwell1(uy,nppm,i1,n,vt)
        call maxwell1(uz,nppm,i1,n,vt)
        call scramble_v(i1,n)   ! remove x,y,z correlations
     else
        call cold_start(i1,n)
     endif

  case(2)  ! drift with vt
	ux(i1:i1+n-1) = vt
	uy(i1:i1+n-1) = 0.
	uz(i1:i1+n-1) = 0.

  case(3)  ! radial velocity spread
     do i = 1, n
        p = i + i1-1
        xt = x(p)-plasma_centre(1)
        yt = y(p)-plasma_centre(2)
        zt = z(p)-plasma_centre(3)

        rt = sqrt(xt**2 + yt**2 + zt**2)
! Set velocity direction according to position vector from plasma centre
        ux(p) = vt*xt/rt
        uy(p) = vt*yt/rt
        uz(p) = vt*zt/rt
     end do

  end select velocities

 

  if (start_debug) then
	write(ipefile,'(a/(i6,8f15.5,i6))') "Initial particle positions, velocities:", &
	  (i,x(i),y(i),z(i),ux(i),uy(i),uz(i),q(i),m(i),pelabel(i),i=i1,i1+n-1)
  endif
end subroutine plasma_start

