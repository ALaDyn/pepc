


program rantest
!VAMPINST include
      INCLUDE 'VTcommon.h'
      INTEGER VTIERR
!
  interface
!VAMPINST program_start
      CALL VTdefs()
      CALL VTENTER(IF_rantest,VTNOSCL,VTIERR)
      write(*,*) 'VT: rantest P>',VTIERR, &
         IF_rantest,ICLASSH
!
!VAMPINST program_end
      CALL VTLEAVE(ICLASSH,VTIERR)
      write(*,*) 'VT: rantest P<',VTIERR,ICLASSH
!
  end interface

  parameter (n=10)
  integer*4 v(n), iflag

  real r(n),x(n)

  write (*,*) '# bits = ',bit_size(1)

  maxnumber = 2147483647
  write (*,'(z10)') maxnumber
    iseed = -74093875
  do while(iseed.ne.0)
     write (*,*) 'Give seed:'
     read (*,*) iseed
     iflag = 0
     do i = 1, n

        r(i) = rano ( iseed )      ! 0-> 1, double precision
     end do

     !  call times2(v)

     write(*,'(i5,f20.16)') (i,r(i),i=1,n)
  end do
     subroutine times2(k)
!VAMPINST include
      INCLUDE 'VTcommon.h'
      INTEGER VTIERR
!
       integer*4, intent(inout) :: k(:)
!VAMPINST subroutine_start
       CALL VTENTER(IF_times2,VTNOSCL,VTIERR)
!      write(*,*) 'VT: times2 S>',VTIERR,
!     *    IF_times2,ICLASSH
!
!VAMPINST subroutine_end
       CALL VTLEAVE(ICLASSH,VTIERR)
!      write(*,*) 'VT: times2 S<',VTIERR,ICLASSH
!
     end subroutine times2

!VAMPINST program_end
      CALL VTLEAVE(ICLASSH,VTIERR)
      write(*,*) 'VT: rantest P<',VTIERR,ICLASSH
!
end program rantest


subroutine times2(k)
!VAMPINST include
      INCLUDE 'VTcommon.h'
      INTEGER VTIERR
!

  integer*4, dimension(:) :: k
!VAMPINST subroutine_start
       CALL VTENTER(IF_times2,VTNOSCL,VTIERR)
!      write(*,*) 'VT: times2 S>',VTIERR,
!     *    IF_times2,ICLASSH
!
  maxnumber = 2147483647/100

  k = k/maxnumber        ! return number between 1 and 100

!VAMPINST subroutine_end
       CALL VTLEAVE(ICLASSH,VTIERR)
!      write(*,*) 'VT: times2 S<',VTIERR,ICLASSH
!
end subroutine times2


! Random number scrambler
! =======================
!
!  called with: 
!               x=rano(iseed)
!
!  returns real number in the interval (0.0 - 1.0)
!  set iseed = -ve integer before using call
!
!  Example of calling routine:
!
!      subroutine coords
!      include 'common.h'
!
!      iseed1 = -11
!      iseed2 = -7
!
!
!      do i=1,n
!	x(i)=xlen*rano(iseed1)
!       y(i)=ylen*rano(iseed2)
!      end do
!
!
!      end


! - routine taken from Numerical Recipies, p195
!
function rano(idum)
!VAMPINST include
      INCLUDE 'VTcommon.h'
      INTEGER VTIERR
!
  real dseed
  dimension v(97)
  data iff,icall/0,0/
!VAMPINST subroutine_start
       CALL VTENTER(IF_rano,VTNOSCL,VTIERR)
!      write(*,*) 'VT: rano S>',VTIERR,
!     *    IF_rano,ICLASSH
!
  if (idum.lt.0.or.iff.eq.0) then
     iff = 1
     dseed=abs(idum)*1.0
     idum=1
     do  j=1,97
        dum=genran(dseed)
     end do
     do j=1,97
        v(j)=genran(dseed)
     end do
     y=genran(dseed)
  endif

  !  next index - make sure we don't overstep array bounds if
  !  generator returns a 0.0 or 1.0

  j=max(mod(1+int(97.*y),98),1)
  if(j.gt.97.or.j.lt.1) then
     write (6,*) 'Call: ',icall
     write (6,*) 'idum = ',idum,'j = ',j,' y= ',y
     write (6,*) 'Random No. generator not initialised properly'
     write (6,*) 'dummy =',dum,' dseed=',dseed
     write (6,100) (i,v(i),i=1,97)
100  format (i4,f10.6)
!VAMPINST stop
!      CALL VTLEAVE(ICLASSH,VTIERR)
!      write(*,*) 'VT: rano S<',VTIERR,ICLASSH
!
     stop
  endif
  !  get next variate and generate new one to fill gap

  y=v(j)
  rano=y
  v(j)=genran(dseed)
  icall = icall + 1

!VAMPINST return
       CALL VTLEAVE(ICLASSH,VTIERR)
!      write(*,*) 'VT: rano S<',VTIERR,ICLASSH
!
  return
!VAMPINST subroutine_end
       CALL VTLEAVE(ICLASSH,VTIERR)
!      write(*,*) 'VT: rano S<',VTIERR,ICLASSH
!
end function rano

                                                                          
real function genran (dseed)                                      
!VAMPINST include
      INCLUDE 'VTcommon.h'
      INTEGER VTIERR
!
!                                  specifications for arguments         
  real   dseed                                          
!                                  specifications for local variables   
  real   d2p31m,d2p31                                   
!                                  d2p31m=(2**31) - 1                   
!                                  d2p31 =(2**31)(or an adjusted value) 
  data               d2p31m/2147483647.0/                          
  data               d2p31 /2147483648.0/                          
!                                  first executable statement           
!VAMPINST subroutine_start
       CALL VTENTER(IF_genran,VTNOSCL,VTIERR)
!      write(*,*) 'VT: genran S>',VTIERR,
!     *    IF_genran,ICLASSH
!
  dseed = mod(16807.0*dseed,d2p31m)                               
  genran = dseed / d2p31                                            
!VAMPINST return
       CALL VTLEAVE(ICLASSH,VTIERR)
!      write(*,*) 'VT: genran S<',VTIERR,ICLASSH
!
  return                                                            
!VAMPINST subroutine_end
       CALL VTLEAVE(ICLASSH,VTIERR)
!      write(*,*) 'VT: genran S<',VTIERR,ICLASSH
!
end function genran
