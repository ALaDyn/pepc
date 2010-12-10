!
!                          General utility module
!                       *   I/O routines
!                       *   Random number gen.
!

module utils


  interface blank
     module procedure blankn, blank6
  end interface

contains



  subroutine blankn(ichan)
    integer :: ichan
    write(ichan,'(/)')
  end subroutine blankn

  subroutine blank6
    write(6,'(/)')
  end subroutine blank6


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
  real function rano(idum)
    implicit none
    integer :: idum
    real, save :: dseed, dum
    real, save :: v(97), y
    integer, save :: iff, icall, i, j
    data iff,icall/0,0/
    if (idum.lt.0.or.iff.eq.0) then
       !  first call generates new sequence
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

    !  next index - make sure we do not overstep array bounds if
    !  generator returns a 0.0 or 1.0

    j=max(mod(1+int(97.*y),98),1)
    if(j.gt.97.or.j.lt.1) then
       write (6,*) 'Call: ',icall
       write (6,*) 'idum = ',idum,'j = ',j,' y= ',y
       write (6,*) 'Random No. generator not initialised properly'
       write (6,*) 'dummy =',dum,' dseed=',dseed
       write (6,100) (i,v(i),i=1,97)
100    format (i4,f10.6)
       stop
    endif
    !  get next variate and generate new one to fill gap

    y=v(j)
    rano=y
    v(j)=genran(dseed)
    icall = icall + 1

    return
  end function rano


  real function genran (dseed)                                      
    !                                  specifications for arguments         
    real ::  dseed                                          
    !                                  specifications for local variables   
    !    real ::  d2p31m,d2p31                                   
    real ::  d2p31m,d2p31    
    !                                  d2p31m=(2**31) - 1                   
    !                                  d2p31 =(2**31)(or an adjusted value) 
    data               d2p31m/2147483647.0/                          
    data               d2p31 /2147483648.0/                          
    !                                  first executable statement           
    dseed = mod(16807.0*dseed,d2p31m)                               
    genran = dseed / d2p31                                            
    return                                                            
  end function genran


end module utils
