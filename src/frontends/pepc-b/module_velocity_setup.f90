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
!>  Encapsulates functions for setting up particle velocities with different models
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_velocity_setup
      use module_particle_props
      use module_physvars
      use module_utilities
      implicit none
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

      public set_velocities
      public cold_start
      public maxwell1
      public maxwell2
      public maxwell3
      public scramble_v
      public scramble_v2
      public perturb_temp



      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      contains

      subroutine set_velocities(i1,n,vt,velocity_config, plasma_centre)

      use module_particle_props
      implicit none

  integer, intent(in) :: i1 ! particle index start
  integer, intent(in) :: n  ! # particles to set up
  integer, intent(in) :: velocity_config ! Maxwellian/cold start switch
  real, intent(in) :: vt ! thermal or drift velocity
  real, intent(in) :: plasma_centre(3) ! target offset vector
	
  integer :: i, p
  real*8 :: rt, xt, yt, zt

  velocities: select case (velocity_config)

  case(0) ! Default is cold start
     call cold_start(i1,n)

  case(1) ! isotropic Maxwellian

     if (vt > 0) then
        call maxwell1(ux(i1:i1+n),n,vt)
        call maxwell1(uy(i1:i1+n),n,vt)
        call maxwell1(uz(i1:i1+n),n,vt)
        call scramble_v(ux(i1:i1+n),uy(i1:i1+n),uz(i1:i1+n),n)   ! remove x,y,z correlations
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

  case(4) ! 2D isotropic Maxwellian

     if (vt > 0) then
        call maxwell2(ux(i1:i1+n),uy(i1:i1+n),n,vt)
        call scramble_v2(ux(i1:i1+n),uy(i1:i1+n),n)   ! remove x,y correlations
	uz(i1:i1+n-1) = 0.
     else
        call cold_start(i1,n)
     endif


  end select velocities

end subroutine set_velocities

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !>   COLD_START
        !>   initialises velocity array slice to 0.
        !>   @param u array of velocities to be initialized
        !>   @param i1 minimal index in u to be used
        !>   @param n maximum index in u to be used
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		subroutine cold_start(i1,n)

		  use module_particle_props
		  implicit none
		  integer, intent(in) :: i1,n

		  ux(i1:i1+n) = 0.
		  uy(i1:i1+n) = 0.
		  uz(i1:i1+n) = 0.

		end subroutine cold_start


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !>   MAXWELL1
        !>   initialises 1D Maxwellian velocity distribution
        !>   @param u array of velocities to be initialized
        !>   @param n maximum index in u to be used
        !>   @param vt desired average velocity of particles
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		subroutine maxwell1(u,n,vt)

		  implicit none
		  integer, intent(in) :: n
		  real, intent(in) :: vt
		  real*8 :: u(n)
		  real, parameter :: pi=3.141592654

		  integer :: ip1, ip2, i, cntr, nv
		  real :: f0 ,df, v, dv, vmax, finf, deltv
		  real*8 :: vip

		  if (n.eq.0) return
		  nv = 30*n
		  vmax = 4.0
		  dv = vmax/nv
		  f0 = 0.5
		  finf = sqrt(pi/2.0)
		  cntr = 1
		  ip1 = n/2
		  ip2 = n/2+1

		  do  i=1,nv
		     v = vmax - (i-0.5)*dv
		     deltv = vmax/10.
		     df = exp( max(-30.0,-0.5*v**2) )*dv/finf*n/2.
		     f0 = f0 + df       ! integrate dist. fn.
		     if(f0.ge.cntr) then
		        vip = vt*(v-dv*(f0-cntr)/df)
		        u(ip1) = vip
		        u(ip2) = -vip
		        cntr = cntr + 1
		        ip1 = max(ip1-1,1)
		        ip2 = min(ip2+1,n)
		     endif
		  end do

		  !  odd one out
		  if (mod(n,2).eq.1) then
		     u(n)=0.
		  endif
		end subroutine maxwell1


       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !>   MAXWELL2
        !>   initialises 2D Maxwellian velocity distribution by direct inversion
	!>   assumes thermal distribution exp{ -(vx**2+vy**2)/2v_t**2 }
        !>   @param ux,uy array of velocities to be initialized
        !>   @param n maximum index in ux to be used
        !>   @param vt desired thermal velocity of particles
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		subroutine maxwell2(u1,u2,n,vt)

		  implicit none
		  integer, intent(in) :: n
		  real, intent(in) :: vt
		  real*8 :: u1(n),u2(n)
		  real, parameter :: pi=3.141592654

		  integer :: i
		  real*8 :: theta, u0
		  integer :: dum1=-319

		  if (n.eq.0) return

		  do  i=1,n
		     u0 = vt*sqrt(-2.*log((i-0.5)/n))
		     theta=2*pi*rano(dum1)
		     u1(i) = u0*cos(theta)
		     u2(i) = u0*sin(theta)
		  end do

		end subroutine maxwell2


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !>   MAXWELL3
        !>   initialises 3D Maxwellian velocity distribution
        !>   @param ux,uy,uz arrays of velocities to be initialized
        !>   @param n maximum index in u to be used
        !>   @param vt desired average velocity of particles
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine maxwell3(ux,uy,uz,n,vt)
          implicit none
          integer, intent(in) ::  n
          real, intent(in) :: vt
          real*8 :: ux(n),uy(n),uz(n)
 	  real :: vc

          vc = vt / sqrt(3.) ! homogeneous: vx=vy=vz = sqrt(|v|^2/3)

          call maxwell1(ux,n,vc)
          call maxwell1(uy,n,vc)
          call maxwell1(uz,n,vc)
          call scramble_v(ux,uy,uz,n)   ! remove x,y,z correlations

        end subroutine maxwell3


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !>   SCRAMBLE_V
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		subroutine scramble_v(ux,uy,uz,n)

		  implicit none

		  integer :: dum1, dum2, dum3, n
		  real*8 :: uxt, uyt, uzt
          	  real*8 :: ux(n), uy(n), uz(n)
		  integer :: j, k, kk, p, n1

		  dum1 = -71 - 10*my_rank
		  dum2 = -113301 - 10*my_rank
		  dum3 = -8651 - 10*my_rank
		  !  exclude odd one out
		  if (mod(n,2).ne.0) then
		     n1=n-1
		  else
		    n1 = n
		  endif

		  !  scramble indices to remove correlation between ux,uy,uz
		  do p=1,n1
		     j=int(n1*rano(dum1)+1)
		     k=int(n1*rano(dum2)+1)
		     kk=int(n1*rano(dum3)+1)
		     uxt=ux(p)
		     uyt=uy(p)
		     uzt=uz(p)
		     ux(p)=ux(kk)
		     uy(p)=uy(j)
		     uz(p)=uz(k)
		     ux(kk)=uxt
		     uy(j)=uyt
		     uz(k)=uzt
		  end do

		end subroutine scramble_v

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>
        !>   SCRAMBLE_V2
        !>
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		subroutine scramble_v2(ux,uy,n)

		  implicit none

		  integer :: dum1, dum2, n
		  real*8 :: uxt, uyt
          	  real*8 :: ux(n), uy(n)
		  integer :: j, k, p, n1

		  dum1 = -71 - 10*my_rank
		  dum2 = -113301 - 10*my_rank

		  !  exclude odd one out
		  if (mod(n,2).ne.0) then
		     n1=n-1
		  else
		    n1 = n
		  endif

		  !  scramble indices to remove correlation between ux,uy
		  do p=1,n1
		     j=int(n1*rano(dum1)+1)
		     k=int(n1*rano(dum2)+1)
		     uxt=ux(p)
		     uyt=uy(p)
		     ux(p)=ux(k)
		     uy(p)=uy(j)
		     ux(k)=uxt
		     uy(j)=uyt
		  end do

		end subroutine scramble_v2
!  ===============================================================
!
!                       PERTURB_TEMP 
!
!   Add temperature perturbation to isotropic thermal distrib. 
!
!  ===============================================================


subroutine perturb_temp 

  implicit none
  integer :: i
  real :: Te0, deltaT, k_therm, s
  real :: lambdaT, xpert, xpe, deltaT0


! Scale velocities by space-dep. temperature perturbation 

Te0=Te_keV
xpert=xl*0.8
lambdaT = xpert/kpert
k_therm = 2*pi/lambdaT   ! Perturbation wavenumber - leave 10% buffer at either end 
deltaT0 = tpert*Te0  ! 50% temperature variation 
if (my_rank==0) then 
  write(*,*) 'PEPC-B | Doing electron temperature perturbation'
  write(*,*) 'PEPC-B | Wavelength: ',lambdaT
endif

do i=1,np_local
  if (q(i)<0 .and. x(i) >= 0.1*xl .and. x(i) < 0.9*xl) then
      xpe = x(i)-xl/10. 
     deltaT = deltaT0 * sin(k_therm*xpe)
     s = sqrt((deltaT+Te0)/Te0)
     ux(i) = s*ux(i)
     uy(i) = s*uy(i)
     uz(i) = s*uz(i)
  else
! buffer region at Te0
  endif
enddo
 

end subroutine perturb_temp



end module module_velocity_setup
