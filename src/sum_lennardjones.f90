!  ===================================================================
!
!                              SUM_LENNARD-JONES
!
!   Calculate Lennard-Jones forces of particle from interaction list
!
!  ===================================================================

subroutine sum_lennardjones( p, n, inode, sumfx, sumfy, sumfz, sumphi )
  use treevars

  integer, intent(in) :: p  ! particle number
  integer, intent(in) :: n  !  # terms on interaction list
  integer, dimension(1:n) ::  inode
  integer :: jnode

  real, intent(out) ::  sumfx,sumfy,sumfz,sumphi 
  real :: a_bond, d2, dlj2, flj, eps2

!  a_bond = 1.1*(Vplas/ni)**(1./3.)  ! mean interparticle spacing - slightly bigger to push particles onto boundaries
  a_bond=a_ii*1.05
  eps2 = (a_bond/4)**2 ! cutoff
  sumfx = 0
  sumfy = 0
  sumfz = 0
  sumphi = 0

  do i=1,n

     !  preprocess distances
     jnode = inode(i)
     dx = x(p) - xcoc( jnode )
     dy = y(p) - ycoc( jnode )
     dz = z(p) - zcoc( jnode ) 

     d2 = dx**2+dy**2+dz**2
     dlj2 = d2/a_bond**2  ! reduced distance
     plj = (1./dlj2**6 - 1./dlj2**3)
     flj = (2./dlj2**6 - 1./dlj2**3)
     plj = min(plj,100.)
     flj = max(min(flj,100.),-100.)
     !   natoms = abs_charge( jnode) / qi    ! # particles in monopole cluster


     ! potential
     sumphi = sumphi + plj

     !  forces

     sumfx = sumfx + 6*dx/d2*flj

     sumfy = sumfy + 6*dy/d2*flj

     sumfz = sumfz + 6*dz/d2*flj


  end do


end subroutine sum_lennardjones
