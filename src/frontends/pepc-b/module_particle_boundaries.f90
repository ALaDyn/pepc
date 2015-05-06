! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2015 Juelich Supercomputing Centre, 
!                         Forschungszentrum Juelich GmbH,
!                         Germany
! 
! PEPC is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! PEPC is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
! 
! You should have received a copy of the GNU Lesser General Public License
! along with PEPC.  If not, see <http://www.gnu.org/licenses/>.
!

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

      public check_particle_bounds
      public constrain
      public wrap
      public reinject_rel
      public reinject_2v

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      contains

!  Entry routine

subroutine check_particle_bounds
    use module_physvars
    use module_particle_props
    implicit none

  boundaries: select case(particle_bcs)

  case(2)
     call constrain   ! relective particle bcs for temperature-clamp mode

  case(3)
     call wrap  ! periodic particles

  case(4)
 ! TODO need to generalize boundary directions
     call wrap ! periodic in y
     call reflect_thermal  ! thermal reinjection on RHB of x-axis

  case default
     ! open (vacuum) boundaries - do nothing 

  end select boundaries

end subroutine check_particle_bounds


!> Periodic wrapper

subroutine wrap
    use module_physvars
    use module_particle_props
    use module_fmm_framework
    implicit none

    integer p

! TODO: need to include periodicity vector

    do p = 1, np_local
	if (particle_wrap(1)) then
	  if (x(p).lt.0) x(p)=x(p)+xl
	  if (x(p).gt.xl) x(p)=x(p)-xl
	endif
	if (particle_wrap(2)) then
	  if (y(p).lt.0) y(p)=y(p)+yl
	  if (y(p).gt.yl) y(p)=y(p)-yl
 	endif
	if (particle_wrap(3)) then
	  if (z(p).lt.0) z(p)=z(p)+zl
	  if (z(p).gt.zl) z(p)=z(p)-zl
	endif

    end do
	
end subroutine wrap

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

! =========================================
!> Reflecting RHB with thermal reinjection
! =========================================

subroutine reflect_thermal
    use module_physvars
    use module_particle_props
    use module_fmm_framework
    implicit none

    integer p

    do p = 1, np_local

      if (x(p).gt.xl) then
	if (q(p).lt.0) then
	  call reinject_2v(vte,-1,ux(p),uy(p),my_rank)
        else 
	  call reinject_2v(vti,-1,ux(p),uy(p),my_rank)
	endif
	x(p)=xl+dt*ux(p)
      end if
    end do
	
end subroutine reflect_thermal


!  ========================================
!
!     2V thermal particle reinjection via direct inversion
!    - conserves  background temp. 
!
!  ========================================

 subroutine reinject_2v(vt,idir,vxi,vyi,my_rank)

 use module_utilities

 implicit none
 real, parameter :: pi=3.14159265
 integer :: my_rank, idir
 real :: vt
 real*8 :: theta, v0, rs, vxi, vyi
 integer, save :: idum=-3193191

!  Flux in x-dirn: invert Int (vx f(ux) dux) directly
     
  rs=rano(idum)
  v0 = vt*sqrt(-2.*log(rs))
  theta=pi*rano(idum)
  vxi = idir*v0*sin(theta)
  vyi = v0*cos(theta)
		
 end subroutine reinject_2v



!  ========================================
!
!     Relativistic thermal particle reinjection
!    - conserves  background temp. 
!    - pairwise injection for uy,uz to give zero transverse current
!
!  ========================================

 subroutine reinject_rel(vt,idir,uxi,uyi,uzi)

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

  end subroutine reinject_rel


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
