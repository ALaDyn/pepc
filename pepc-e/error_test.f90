!  =========================
!
!  Estimate force errors
!
!  =========================

      subroutine error_test(ntest)

      use physvars
      implicit none
!      integer, parameter :: ntest = 3
      integer ntest
      real, dimension(ntest) ::  potd, exd, eyd, ezd
      real :: dfx2, dfy2, dfz2, fxs, fys, fzs, dpot
      real :: errfx, errfy, errfz
      integer :: i



! direct force evaluation

      call force_direct(ntest,x(1:ntest),y(1:ntest),z(1:ntest),q(1:ntest),eps,force_const,exd(1:ntest),eyd(1:ntest),ezd(1:ntest),potd(1:ntest))


! find rms error

      dfx2 = 0.
      dfy2 = 0.
      dfz2 = 0.
      fxs = 0.
      fys = 0.
      fzs = 0.
      do i=1,ntest
        dfx2 = dfx2 + (exd(i) - ex(i))**2
        dfy2 = dfy2 + (eyd(i) - ey(i))**2
        dfz2 = dfz2 + (ezd(i) - ez(i))**2
        dpot = dpot + (potd(i) - pot(i))/potd(i)
        fxs = fxs + exd(i)**2
        fys = fys + eyd(i)**2
        fzs = fzs + ezd(i)**2
      end do

      write (*,*) '  i   pot_t    pot_d    ex_t    ex_d'
      write (*,'((i5,4(1pe14.4)))') (i,pot(i),potd(i),ex(i),exd(i),i=1,ntest)

      errfx = sqrt(dfx2/fxs)
      errfy = sqrt(dfy2/fys)
      errfz = sqrt(dfz2/fzs)
      write (6,101) errfx,errfy,errfz,(errfx+errfy+errfz)/3.
 101  format (/'Force errors - x,y,z,ave: ',4(1pe15.5))
  
      end subroutine error_test


