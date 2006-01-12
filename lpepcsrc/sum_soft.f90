!  ===================================================================
!
!                              SUM_SOFT
!
!   Soft potential for initialising ions in container
!
!  ===================================================================

subroutine sum_soft( p, n, inode, a_bond, sumfx, sumfy, sumfz, sumphi )

  use treevars
  implicit none
  integer, intent(in) :: p  ! particle number
  integer, intent(in) :: n  !  # terms on interaction list
  integer, dimension(1:n) ::  inode
  real, intent(in) :: a_bond
  integer :: jnode, i

  real*8, intent(out) ::  sumfx,sumfy,sumfz,sumphi 
  real*8 ::  d2, d, pots, fxs, fys, fzs, epsc, plj, dx, dy, dz

  sumfx = 0
  sumfy = 0
  sumfz = 0
  sumphi = 0

  do i=1,n

     !  preprocess distances
     jnode = inode(i)
     dx =  x(p) - xcoc( jnode ) 
     dy =  y(p) - ycoc( jnode ) 
     dz =  z(p) - zcoc( jnode ) 

     d2 = dx**2+dy**2+dz**2
     d = sqrt(d2) 

     pots = -d*exp(-d/a_bond)
     fxs = -dx/d*exp(-d/a_bond)*( d/a_bond -1 )
     fys = -dy/d*exp(-d/a_bond)*( d/a_bond -1 )
     fzs = -dz/d*exp(-d/a_bond)*( d/a_bond -1 )



     ! potential
     sumphi = sumphi + pots

     !  forces

     sumfx = sumfx + fxs
     sumfy = sumfy + fys
     sumfz = sumfz + fzs


  end do


end subroutine sum_soft



