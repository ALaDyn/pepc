! ==============================================
!
!                SCRAMBLE_V
!
!  Sets up physical system: particle positions, velocities
!
! ==============================================

subroutine scramble_v(i1,n)

  use physvars
  use utils
  implicit none

 integer :: dum1, dum2, dum3  
real :: uxt, uyt, uzt
  integer :: i, j, k, kk, p, i1, n, n1

  dum1 = -71 - 10*my_rank
  dum2 = -113301 - 10*my_rank
  dum3 = -8651 - 10*my_rank
  !  exclude odd one out
  if (mod(n,2).ne.0) then
     n1=n-1
  else
    n1 = n
  endif
 write(15,*) 'scrambling ',n1,' particles'
  !  scramble indices to remove correlation between ux,uy,uz
  do i=1,n1
     p=i+i1-1
     j=n1*rano(dum1)+i1
     k=n1*rano(dum2)+i1
     kk=n1*rano(dum3)+i1
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
