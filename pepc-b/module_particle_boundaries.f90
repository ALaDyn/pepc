!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>  Encapsulates anything concerning the particle bcs
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_particle_boundaries
 
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

      public constrain
      public earth_plate
      public reinject


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!  CONSTRAIN
!
!!!!!!!!!!!!!!!!

subroutine constrain
    use module_physvars
    use module_particle_props
    use module_geometry

    implicit none

    real*8, dimension(1:3)		:: r_new, r_old, r_test, v
    real*8 				:: c_status		! particle in or out?
    integer 			        :: p, face_nr

    ! bisections
    real*8, dimension(1:3)		:: r_d, temp	        ! crossing point of particle
    real*8				:: p_new, p_old, diff

    ! reflection plane
    real*8, dimension(1:3)		:: n

    ! new direction of v
    real*8, dimension(1:3)                :: c

    ! file id
    integer                             ::  nr_out
    logical :: constrain_debug=.false.

    if (constrain_debug) write(*, *) 'in constrain'

    ! walk through all particles
    nr_out = 0
    do p = 1, np_local
        r_new = (/x(p), y(p), z(p)/)
        v = (/ux(p), uy(p), uz(p)/)

        do face_nr = 1, number_faces ! walk through all faces
            call face(r_new, c_status, face_nr, & 
              target_geometry, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre )

            if (c_status .ge. 0.) then
                nr_out = nr_out + 1

                ! get a good r_old
                r_test = r_new - dt * v
                call face(r_test, c_status, face_nr, & 
                  target_geometry, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre )
                if (c_status .gt. 0.) then
                    call cutvector(r_test, face_nr, n, r_old, &
                        target_geometry, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre )
                   temp = r_new - r_old
                    v = sqrt(dot_product(v, v)) * temp / sqrt(dot_product(temp, temp))
                else
                    r_old = r_test
                end if
        
                ! bisection for the particle
                p_new = 1.
                p_old = 0.
                r_test = r_old + p_new * (r_new - r_old)
                call face(r_test, c_status, face_nr, &
                  target_geometry, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre )
                do while (abs(c_status) .gt. constrain_proof)
                    diff = abs(p_new - p_old) / 2.
                    p_old = p_new
                    if (c_status .gt. 0.) then
                        p_new = p_new - diff
                    else
                        p_new = p_new + diff
                    end if
                    r_test = r_old + p_new * (r_new - r_old)
                    call face(r_test, c_status, face_nr, &
                       target_geometry, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre )
                end do
!                if (constrain_debug) write(*, *) 'Bisection for particle done.'
                call cutvector(r_test, face_nr, n, r_d, & 
                  target_geometry, x_plasma, y_plasma, z_plasma, r_sphere, plasma_centre )
                ! Reflect
                v = v - 2. * dot_product(v, n) * n
                if (constrain_debug) write(*, *) 'Reflecting .. new v: v     = ', v(1:3)
                c = v / sqrt(dot_product(v, v))
                r_new = r_d + sqrt(dot_product(r_d - r_new, r_d - r_new)) * c
                x(p) = r_new(1)
                y(p) = r_new(2)
                z(p) = r_new(3)
                ux(p) = v(1)
                uy(p) = v(2)
                uz(p) = v(3)
            end if
        end do
    end do
 if (constrain_debug)    write(*, *) 'Number of reflections: ', nr_out
 
end subroutine constrain


!  ===============================================================
!
!                          EARTH_PLATE 
!
!    Emulate grounded earth plate at z=0 
!
!  ===============================================================

subroutine earth_plate

  use module_physvars
  use module_particle_props
  use module_utilities
  implicit none

  integer :: p
  integer, save :: iseed1=5
  integer, save :: iseed2=7
  real :: x_limit, y_limit, r_limit, xt, yt, zt, xt2, yt2

  do p=1,np_local
     xt = x(p) - plasma_centre(1)  ! shift origin to 1st target centre
     yt = y(p) - plasma_centre(2)
     zt = z(p) - plasma_centre(3)

     xt2 = xt - displace(1)
     yt2 = yt - displace(2)

     r_limit = r_sphere
     x_limit = x_plasma/2.
     y_limit = y_plasma/2.

     if (target_geometry==1) then
        ! sphere
        !	      constrained = (xt**2 + yt**2 + zt**2 <= r_limit**2) 


     else if (target_geometry==2) then

     else if (target_geometry==3) then
        ! wire

        if ( q(p)>0 .and. zt < -x_limit ) then
 ! Reflect ions that hit earth plate
           z(p) = z(p) + 2*(-x_limit - zt)
           uz(p) = -uz(p)

        else if ( q(p)<0 .and. zt <  -x_limit ) then 

 ! Absorb electrons and reinject with thermal velocity
 	   call reinject(vte,1,ux(p),uy(p),uz(p))
!write(*,*) 'new velocities:', ux(p),uy(p),uz(p)
           z(p)= plasma_centre(3) -x_plasma/2 + uz(p)*dt/2. 

           if (pelabel(p) <= npart_total/2 .and. xt**2+yt**2 > r_sphere**2) then
 ! Wire 1 electrons falling outside wire1 get reinjected inside wire1
             do while (xt**2+yt**2 > r_sphere**2)
                xt = r_sphere*(2*rano(iseed2)-1.)
                yt = r_sphere*(2*rano(iseed1)-1.)
             end do

             y(p) = yt + plasma_centre(2)
             x(p) = xt + plasma_centre(1)

           else if (pelabel(p) <= npart_total/2 .and. xt**2+yt**2 <= r_sphere**2) then
  ! Wire 1 electrons falling inside get put back in  wire 2
	     y(p) = y(p) + displace(2)
   	     x(p) = x(p) + displace(1)
!             pelabel(p) = pelabel(p) + npart/2  ! Relabel as wire2 electron
! write(*,*) 'wire 1 reflect:',x(p),y(p),z(p),ux(p),uy(p),uz(p),pelabel(p)
 
           else if (pelabel(p) > npart_total/2 .and. xt2**2+yt2**2 > r_sphere**2) then
 ! Wire 2 electrons falling outside wire2 get reinjected inside wire2
             do while (xt2**2+yt2**2 > r_sphere**2)
                xt2 = r_sphere*(2*rano(iseed2)-1.)
                yt2 = r_sphere*(2*rano(iseed1)-1.)
             end do
             y(p) = yt2 + plasma_centre(2)+displace(2)
             x(p) = xt2 + plasma_centre(1)+displace(1)

           else if (pelabel(p) > npart_total/2 .and. xt2**2+yt2**2 <= r_sphere**2) then
  ! Wire 2 electrons falling inside get put back in  wire 1
	     y(p) = y(p) - displace(2)
   	     x(p) = x(p) - displace(1)
 !            pelabel(p) = pelabel(p) - npart/2  ! Relabel as wire1 electron
 ! write(*,*) 'wire 2 reflect:',x(p),y(p),z(p),ux(p),uy(p),uz(p),pelabel(p)
           endif

        endif


     else
 
    endif
 
end do

end subroutine earth_plate


!  ========================================
!
!     Thermal particle reinjection
!    - conserves  background temp. 
!    - pairwise injection for uy,uz to give zero transverse current
!
!  ========================================

 subroutine reinject(vt,idir,uxi,uyi,uzi)

 use module_utilities

 implicit none
 real, parameter :: pi=3.14159265
 integer, parameter :: nsamp=300000
 integer :: idir, ntrial
 real*8, save :: usamp(nsamp+1)
 real*8  :: uxi, uyi
 real*8, save :: theta
 real*8 :: df,f0,A,u,du,g,uzi,g0
 real :: rs
 real :: ute,  pio2, umax
 real, intent(in) :: vt
 integer, save :: idum=-13
 integer, save :: isamp=1
 integer, save :: i2
 integer :: i,k, nsteps

!  conserve flux in x dirn: maintain const. Te

!  Tabulate momenta on first call
    ute = vt**2  ! Temperature normalised to mc^2
    pio2 = pi/2
   if (isamp.eq.1) then

        nsteps=100000
!  dimensionless electron/ion temp
        umax=6.*sqrt(ute + 4*ute**2)
        du=umax/nsteps
!  distn norm factor
        A = 1.0001*nsamp/2./pi/ute/(1+ute)
        f0=0.
        k=1

!  load loop
        do i=1,nsteps

!  cumulative integral of rel. distn function f(u)
	  u = i*du
          g = sqrt(u**2 + 1.0)
	  df = 2*pi*A*du*u*exp((1.-g)/ute)

	  f0 = f0+df

	  if(f0.ge.k) then
!  store u corresponding to integer values of f0
            usamp(k) = u
	    k=k+1
	  endif
        end do
!        write(40,*) (usamp(k),k=1,nsamp)
      endif

!  Flux in x-dirn: invert Int (vx f(ux) dux) directly
     
      rs=rano(idum)
      g0 = max(1.0,1.0-ute*alog(rs))
      uzi = idir*sqrt(g0**2-1.d0)

!  pick random  u from sample:  
!  usamp contains integrated momentum flux Int(u f(u) du)
      if (mod(isamp,2).eq.1) then
	ntrial = nsamp*rano(idum)+1
       i2=min(nsamp,ntrial)
!  components for uy,uz
       theta=2*pi*rano(idum)
       uxi = usamp(i2)*cos(theta)
       uyi = usamp(i2)*sin(theta)
!  write(40,*) 'isamp=',isamp,' i2=',i2,' idum=',idum,'ux,uy=',uxi, uyi
      else
!  use -ve of previous sample
       uxi = -usamp(i2)*cos(theta)
       uyi = -usamp(i2)*sin(theta)
!  write(40,*) 'isamp=',isamp,' i2=',i2,' idum=',idum, 'ux,uy=',uxi, uyi
      endif
  
!  increment sample index
      isamp = isamp + 1
      if (isamp.gt.nsamp) isamp = 1

  end subroutine reinject


!  ===============================================================
!
!                           RESET_IONS
!
!   Reset ion mass and velocity after const-temp eqm phase
!
!  ===============================================================

subroutine reset_ions

  use module_physvars
  use module_particle_props

  integer :: p
  real :: ratio_clamp  ! mass ratio used for NVT dynamics

  ratio_clamp = mass_i/mass_e  ! (typically 10)
  do p=1,np_local

     if (q(p) > 0) then   ! only constrain ions to target boundaries
        m(p) = mass_ratio*mass_e  ! mass_ratio preserved from initial inputs
        ux(p) = ux(p)/sqrt(mass_ratio/ratio_clamp)  ! scale velocities back so that T_i = mv^2 conserved
        uy(p) = uy(p)/sqrt(mass_ratio/ratio_clamp) 
        uz(p) = uz(p)/sqrt(mass_ratio/ratio_clamp) 
     endif
  end do
end subroutine reset_ions
end module module_particle_boundaries
