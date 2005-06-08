!  ===================================================================
!
!                              SUMPOT
!
!   Calculate potential of particle from interaction list
!   Pseudoparticles are given by: intlist(j,p), j=1,nterm(p)
!
!  ===================================================================

subroutine sumpot(p, n, inode, pote)

  use treevars

  integer, intent(in) :: p  ! particle label
  integer, dimension(1:n) :: inode

  !  work arrays
  real, dimension(1:n) :: dx,dy,dz,d,dx2,dy2,dz2,d3,d5,eq1,eq2,eq3,eq4,eq5,eq6,petotal

  real pote, eps2

  eps2=eps**2

 !  preprocess distances - loops implicit from 1 -> n

  dx = x(p) - xcoc( inode )
  dy = y(p) - ycoc( inode )
  dz = z(p) - zcoc( inode ) 

  !  sum potential to quadrupole
  !       ep(p)=ep(p)+Epm+Epd+Epq

  dx2 = dx**2
  dy2 = dy**2
  dz2 = dz**2
  d = sqrt(dx2+dy2+dz2+eps2)
  d3 = d**3
  d5 = d**5

  eq1 = 0.5*( 3*dx2/d5 - 1/d3 )
  eq2 = 0.5*( 3*dy2/d5 - 1/d3 )
  eq3 = 0.5*( 3*dz2/d5 - 1/d3 )
  eq4 = 3*dx*dy/d5
  eq5 = 3*dy*dz/d5
  eq6 = 3*dx*dz/d5

  petotal =   charge( inode )/d    &                                     !  monopole term
       !
            + (dx*xdip( inode ) + dy*ydip( inode ) + dz*zdip( inode ))/d3  &    !  dipole
       !
            + eq1*xxquad( inode ) + eq2*yyquad( inode ) + eq3*zzquad( inode )  &  !  quadrupole
            + eq4*xyquad( inode ) + eq5*yzquad( inode ) + eq6*zxquad( inode )      

  !  add up contributions
  pote = SUM(petotal)

end subroutine sumpot
