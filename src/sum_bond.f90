!  ===================================================================
!
!                              SUM_BOND
!
!   Calculate SHO bonding forces of particle from interaction list
!
!  ===================================================================

subroutine sum_bond( p, n, inode, sumfx, sumfy, sumfz, sumphi )
  use treevars

  integer, intent(in) :: p  ! particle number
  integer, intent(in) :: n  !  # terms on interaction list
  integer, dimension(1:n) ::  inode
  integer :: jnode

  real, intent(out) ::  sumfx,sumfy,sumfz,sumphi 
  real :: a_bond

  a_bond = Vplas**(1./3.)  ! mean interparticle spacing
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

     d = sqrt(dx**2+dy**2+dz**2)

!   natoms = abs_charge( jnode) / qi    ! # particles in monopole cluster

     if ( d <= 2*a_bond) then
        ! potential
        sumphi = sumphi + 0.5*(d - a_bond)**2

        !  forces

        sumfx = sumfx - dx/d*(d - a_bond)

        sumfy = sumfy - dy/d*(d - a_bond)

        sumfz = sumfz - dz/d*(d - a_bond)
     else
        sumphi = sumphi + a_bond**2/2.
     endif

  end do


end subroutine sum_bond
