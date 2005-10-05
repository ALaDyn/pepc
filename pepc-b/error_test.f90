!  =========================
!
!  Estimate force errors
!
!  =========================

      subroutine error_test(ntest)

      use physvars
      use treevars
      use utils
      implicit none
!      integer, parameter :: ntest = 3
      integer, intent(in) ::  ntest
      integer :: listerr(ntest)
      real, dimension(ntest) ::  potd, exd, eyd, ezd
      real :: dfx2, dfy2, dfz2, fxs, fys, fzs, dpot, spot
      real :: errfx, errfy, errfz, errf_ave, err_pot
      integer :: i, iseed = -317, isamp


! direct force evaluation
! TODO: make random list of particles for sample

      do i=1,ntest
         isamp = npp*rano(iseed)+1
         listerr(i) = isamp

      end do

! Get direct forces, potential
      call force_direct(npp,ntest,x(1:npp),y(1:npp),z(1:npp),q(1:npp),listerr(1:ntest),eps,force_const,exd,eyd,ezd,potd)

! find rms error

      dfx2 = 0.
      dfy2 = 0.
      dfz2 = 0.
      fxs = 0.
      fys = 0.
      fzs = 0.
      dpot=0.
      spot=0.

      do i=1,ntest
        dfx2 = dfx2 + (exd(i) - ex(listerr(i)))**2
        dfy2 = dfy2 + (eyd(i) - ey(listerr(i)))**2
        dfz2 = dfz2 + (ezd(i) - ez(listerr(i)))**2
        dpot = dpot + potd(i) - pot(listerr(i))
	spot = spot+ potd(i)
        fxs = fxs + exd(i)**2
        fys = fys + eyd(i)**2
        fzs = fzs + ezd(i)**2
      end do

      write (*,*) '  i   list pot_t    pot_d    ex_t    ex_d'
      write (*,'((2i8,4(1pe14.4)))') (i,listerr(i),pot(listerr(i)),potd(i),ex(listerr(i)),exd(i),i=1,ntest)
      err_pot = dpot/spot
      errfx = sqrt(dfx2/fxs)
      errfy = sqrt(dfy2/fys)
      errfz = sqrt(dfz2/fzs)
      errf_ave = (errfx+errfy+errfz)/3.
      write (6,'(a/a20,1pe13.6/a20,4(1pe13.6))') 'Relative errors:','Potential ',err_pot, &
	'Forces (ave,x,y,z) ',errf_ave,errfx,errfy,errfz
  
      end subroutine error_test


