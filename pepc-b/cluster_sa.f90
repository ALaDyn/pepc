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

subroutine cluster_sa(i1,nip,r0,r_sphere,qi,Qplas,plasma_centre)

    use treevars
    use utils

    implicit none

    real, intent(in) :: r0  ! characteristic radius
    real, intent(in) :: plasma_centre(3) !  sphere centre
    real, intent(in) :: r_sphere ! equivalent sphere radius for qi 
    real, intent(in) :: qi ! particle charge
    real, intent(in) :: Qplas ! total plasma charge

    integer, intent(in) :: i1  ! index offset
    integer, intent(in) :: nip ! # particles to set up
    integer              :: i, j, ierr
    integer              :: p, k, nramp, nx, ny
    real*8 :: rt, xt, yt, zt, s, pi, xi
    real*8 :: r_c, Q_c, Q_s, dens
    integer :: nitot, N_c, N_s, npi_c, npi_s, nbin, np_rest, offset
    real*8 :: zeta_max, zeta, dzeta, t_start, a_const, b_const, Qt
    real*8 :: SZ, T, ch, sh, sh2, th, ch2, rp, integ, nom, dom
    
! Andreev cluster:
! uniform up to r=0.1125 r0
! tapered from r=0.1125 r0 to 2.5 r0
    pi = 2.*asin(1.d0)
    zeta_max = 0.23  ! start point r(0.23)=0.1125
! # ions in central region
    r_c = 0.1125*r0
    Q_c = 4*pi/3.*r_c**3 ! (rho=1)
    N_c = Q_c/qi
    nitot = Qplas/qi

! # charge in outer region
    Q_s = Qplas - Q_c
    N_s = Q_s/qi

! # put remainder central particles on last cpu 
    if (me==num_pe-1) then
      np_rest = mod(N_c,num_pe)
    else
      np_rest = 0
    endif
    npi_c = N_c/num_pe + np_rest
    npi_s = nip-npi_c
    offset = me*(nip-N_c/num_pe)  ! shell offset same for all

    if (me==0) then
      write(*,'(a30,3i12)') "Total ions, in centre/outside:",nitot,N_c,N_s
      write(*,'(a30,3i12)') "Local ions in centre/outside:",npi_c,npi_s
      write(*,'(a30,3f12.4)') "Radii r_eff, r_c, r0:", r_sphere, r_c, r0
      write(*,'(a30,3f12.4)') "Charge Q_tot, Q_c, Q_s",Qplas, Q_c, Q_s
    endif

! # rescale radii of inner ions
    do i=1,npi_c
	p = i + i1-1  ! local index (including offset)
        j = me*(npi_c-np_rest) + i   ! Unique global particle #
       	xt = x(p)-plasma_centre(1)
	yt = y(p)-plasma_centre(2)
	zt = z(p)-plasma_centre(3)
	xi = (1.*j/N_c)**(1./3.) !  inversion for uniform density
	rt = sqrt(xt**2+yt**2+zt**2)
        x(p) = plasma_centre(1) + xt/rt*r_c*xi  ! keep direction vector; scale by inner sphere radius 
	y(p) = plasma_centre(2) + yt/rt*r_c*xi
	z(p) = plasma_centre(3) + zt/rt*r_c*xi
        write (ipefile,'(i6,a8,f15.5)') j,'r/r0, ',r_c*xi/r0
     end do


! # Integrate outer zone numerically, according to Andreev parametric form
! # Each CPU will place particles in different segments of integrated shell
!     nbin = 30*N_s
     nbin = 100000
     dzeta = -zeta_max/nbin
     zeta=zeta_max
     Qt = 0.
     t_start = 0.05  
     j = offset + 1  ! Global starting index of ions
     i = npi_c+i1 ! Local index
        if(me==0) write (*,*) 'zeta_max, dzeta, nbin',zeta_max,dzeta,nbin

    do while (zeta>0.005 .and. j<=offset + npi_s)

! Integrand
        sh2 = sinh(2*zeta)
        th = tanh(zeta)
        ch2 = cosh(2*zeta)
        ch = cosh(zeta)
        SZ = sh2/2. + zeta
        T = 1.-zeta*th

	dens = 1./3./t_start**2*SZ**2/(ch**4*(1.-zeta*th))
        integ =  (sh2*S - ch**2*(ch2 + 1)) / (T*SZ**2) 
        Qt = Qt + integ*t_start*dzeta*N_s
        rp = t_start*ch**2/SZ ! Radius associated with zeta
!        if(me==0) write (*,'(8f12.5)') rp,zeta,integ,Qt

        if (Qt > j) then      
! place particle
          xt = x(i)-plasma_centre(1)  ! Get direction vector
          yt = y(i)-plasma_centre(2)
          zt = z(i)-plasma_centre(3)
          rt = sqrt(xt**2+yt**2+zt**2)
!         write (ipefile,'(i6,a23,4f15.5)') i,'r/r0, zeta, integ, Qt',rp,zeta,integ,Qt
         write (*,'(i6,a23,5f15.5)') i,'r/r0, zeta, dens, integ, Qt',rp,zeta,dens,integ,Qt

          x(i) = plasma_centre(1) + xt/rt*rp*r0  ! scale by sphere radius 
          y(i) = plasma_centre(2) + yt/rt*rp*r0
          z(i) = plasma_centre(3) + zt/rt*rp*r0
          j=j+1
          i=i+1
       endif
       zeta = zeta+dzeta
    end do



end subroutine cluster_sa 


