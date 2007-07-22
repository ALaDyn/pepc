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

subroutine cluster_sa(i1,nip,r0,n0,r_sphere,qi,Qplas,plasma_center,miome)

    use treevars
    use utils

    implicit none

    real, intent(in) :: r0  ! characteristic radius (cm)
    real, intent(in) :: n0  ! characteristic cluster density (cm**-3)
    real, intent(in) :: plasma_center(3) !  sphere center
    real, intent(in) :: r_sphere ! equivalent sphere radius for qi 
    real, intent(in) :: qi ! particle charge - already assigned in configure (plasma_start)
    real, intent(in) :: Qplas ! total plasma charge = 1
    real, intent(in) :: miome ! mass ratio

    integer, intent(in) :: i1  ! index offset
    integer, intent(in) :: nip ! # particles to set up
    real :: n0n ! normalised central density
    real :: n_max ! central density normalised to n0
    integer              :: i, j, ierr
    integer              :: p, k, nramp, nx, ny
    real*8 :: rt, xt, yt, zt, pi, xi
    real*8 :: r_c, Q_c, Q_s, dens, Qnorm, dens0, v_max
    integer :: nitot, N_c, N_s, npi_c, npi_s, nbin, np_rest, offset
    real*8 :: zeta_max, zeta, dzeta, t_start, a_const, b_const, Qt
    real*8 :: S, S_c, T, ch, sh, sh2, th, ch2, rp, integ, nom, dom
    
! Andreev cluster:
! uniform up to r=0.1125 r0
! tapered from r=0.1125 r0 to 2.5 r0
    pi = 2.*asin(1.d0)
    t_start = 0.05  
    zeta_max = 0.23  ! start point r(0.23)=0.1125

! # ions in central region
!    r_c = 0.1125  ! center radius normalised to r0 
    r_c = t_start*cosh(zeta_max)**2/(sinh(2*zeta_max)/2. + zeta_max)  ! Radius associated with zeta_max
    dens0 = 1./3./t_start**2*(sinh(2*zeta_max)/2. + zeta_max)**2/(cosh(zeta_max)**4*(1.-zeta_max*tanh(zeta_max)))
 
    v_max = sqrt(2./3.)/sqrt(dens0)/sqrt(miome)*r0

  ! # ions in central, self-sim regions

    n_max = n0*dens0
!    N_c = 4*pi/3.*n_max*(r_c*r0)**3 ! # ions in central portion
    Q_c = 0.041  ! Charge in center (fixed to give correct field - 4pi/3 rho_m r_c^3; rho_m=6.917)
    N_c = Q_c/qi
    nitot = Qplas/qi  ! total ions (check)
    n0n = Q_c*3/4/pi/r_c**3 ! normalised central density

! # charge in outer region
    Q_s = Qplas - Q_c
    N_s = Q_s/qi
    Qnorm = 1.0 ! normalisation factor for charge integral

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
      write(*,'(a30,3i12)') "Total ions, in center/outside:",nitot,N_c,N_s
      write(*,'(a42,3i12)') "Local ions in center/outside:",npi_c,npi_s
       write(*,'(a30,4(1pe12.4))') "Densities:", n0, n_max, dens0, n0n, Qnorm
     write(*,'(a30,3f12.4)') "Radii r_eff, r_c, r0:", r_sphere, r_c, r0
      write(*,'(a30,3f12.4)') "Charge Q_tot, Q_c, Q_s",Qplas, Q_c, Q_s
    endif

! # rescale radii of inner ions
    do i=1,npi_c
	p = i + i1-1  ! local index (including offset)
        j = me*(npi_c-np_rest) + i   ! Unique global particle #
       	xt = x(p)-plasma_center(1)
	yt = y(p)-plasma_center(2)
	zt = z(p)-plasma_center(3)
	xi = (1.*j/N_c)**(1./3.) !  inversion for uniform density
	rt = sqrt(xt**2+yt**2+zt**2)
        x(p) = plasma_center(1) + xt/rt*r_c*xi  ! keep direction vector; scale by inner sphere radius 
	y(p) = plasma_center(2) + yt/rt*r_c*xi
	z(p) = plasma_center(3) + zt/rt*r_c*xi
! Scale velocities according to ss value at radius r_c 
!          ux(i) = tanh(zeta_max)*xt/rt*xi
!          uy(i) = tanh(zeta_max)*yt/rt*xi
!          uz(i) = tanh(zeta_max)*zt/rt*xi
	ux(i)=0.
	uy(i)=0.
	uz(i)=0.
        write (ipefile,'(i6,a8,f15.5)') j,'r/r0, ',r_c*xi
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
        if(me==0) write (*,*) 'zeta_max, dzeta, nbin',zeta_max,dzeta,nbin

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
          xt = x(i)-plasma_center(1)  ! Get direction vector
          yt = y(i)-plasma_center(2)
          zt = z(i)-plasma_center(3)
          rt = sqrt(xt**2+yt**2+zt**2)
!         write (ipefile,'(i6,a23,4f15.5)') i,'r/r0, zeta, integ, Qt',rp,zeta,integ,Qt
!         if (me==0) write (*,'(2i6,a23,5(1pe15.7))') i,j,' r/r0, dens, integ, Qt/qi',rp,dens,Qt/qi

          x(i) = plasma_center(1) + xt/rt*rp  ! scale by sphere radius 
          y(i) = plasma_center(2) + yt/rt*rp
          z(i) = plasma_center(3) + zt/rt*rp
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


