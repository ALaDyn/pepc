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

subroutine cluster_sa(r0)

    use physvars
    use treevars
    use utils

    implicit none

    real, intent(in) :: r0 
    integer              :: i, j, ierr
    integer              :: idum, iseed1, iseed2, iseed3, i1, n1,p, k, nramp, nx, ny
    real                 :: qtot, qramp, xi, xedge, rcut 
    real*8 :: rt, xt, yt, zt, s
    real :: r_c, Q_c, Q_s
    integer :: N_c, N_s, npi_c, npi_s, nbin
    real*8 :: zeta, dzeta, a_const, b_const, Qt
    real*8 :: ch, sh, sh2, th, ch2, rp, integ, nom, dom
    
! Andreev cluster:
! uniform up to r=0.1 r0
! tapered from r=0.1 r0 to 2.5 r0

! # ions in central region
    r_c = 0.1125*r_layer(1)
    Q_c = 4*pi/3.*r_c**3 ! (rho=1)
    N_c = Q_c/qi

! # charge in outer region
    Q_s = Qplas - Q_c
    N_s = Q_s/qi

! # particles on this cpu in central, outer regions
! must be exactly divisible
    npi_c = N_c/num_pe
    npi_s = nip-npi_c

    if (me==0) then
       write(*,*) "Total ions in centre/outside:",N_c,N_s
       write(*,*) "Local ions in centre/outside:",npi_c,npi_s
       write(*,*) "Radii r_c, r_eff, r0:",r_c, r_sphere, r_layer(1)
       write(*,*) "Charge Q_c, Q_s",Q_c, Q_s
    endif

! # rescale radii of inner ions
    do i=1,npi_c
        j = me*npi_c + i   ! Unique global particle #
       	xt = x(i)-plasma_centre(1)
	yt = y(i)-plasma_centre(2)
	zt = z(i)-plasma_centre(3)
	xi = (1.*j/N_c)**(1./3.) !  inversion for uniform density
	rt = sqrt(xt**2+yt**2+zt**2)
        x(i) = plasma_centre(1) + xt/rt*r_c*xi  ! keep direction vector; scale by inner sphere radius 
	y(i) = plasma_centre(2) + yt/rt*r_c*xi
	z(i) = plasma_centre(3) + zt/rt*r_c*xi
     end do


! # Integrate outer zone numerically, according to Andreev parametric form
! # Each CPU will place particles in different segments of integrated shell
     nbin = 100*N_s
     zeta = 0.23
     dzeta = -zeta/nbin
     Qt = 0.
     a_const = 4.805*4*pi/3./Q_s
     b_const = 0.05
     j = me*npi_s + 1  ! Global starting index of ions
     i = npi_c+1 ! Local index

    do while (zeta > 0.01 .and. j<=(me+1)*npi_s)

! Integrand
        sh2 = sinh(2*zeta)
        th = tanh(zeta)
        ch2 = cosh(2*zeta)
        ch = cosh(zeta)
        integ = a_const*b_const**3/(1.-zeta*th) * ( ( sh2*(.5*sh2 + zeta) ) - ch**2*(ch2 + 1) ) &
             / (.5*sh2 + zeta)**2
        Qt = Qt + integ*N_s*dzeta
        nom = b_const*ch**2
        dom = (sh2/2. + zeta)
        rp = nom/dom
!        rp = b_const*ch**2/(sh2/2. + zeta)  ! Radius associated with zeta
!         write (*,'(9f12.5)') rp,zeta,ch,sh2,ch2,nom,dom,integ,Qt

        if (Qt > j) then      
! place particle
          xt = x(i)-plasma_centre(1)  ! Get direction vector
          yt = y(i)-plasma_centre(2)
          zt = z(i)-plasma_centre(3)
          rt = sqrt(xt**2+yt**2+zt**2)
         write (ipefile,'(i6,a23,4f12.5)') i,'r/r0, zeta, integ, Qt',rp,zeta,integ,Qt

          x(i) = plasma_centre(1) + xt/rt*rp*r_layer(1)  ! scale by sphere radius 
          y(i) = plasma_centre(2) + yt/rt*rp*r_layer(1)
          z(i) = plasma_centre(3) + zt/rt*rp*r_layer(1)
          j=j+1
          i=i+1
       endif
       zeta = zeta+dzeta
    end do



end subroutine cluster_sa 


