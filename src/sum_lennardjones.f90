!  ===================================================================
!
!                              SUM_LENNARD-JONES
!
!   Calculate Lennard-Jones forces of particle from interaction list
!
!  ===================================================================

subroutine sum_lennardjones( p, n, inode, sumfx, sumfy, sumfz, sumphi )
  use physvars
  use treevars
  implicit none
  integer, intent(in) :: p  ! particle number
  integer, intent(in) :: n  !  # terms on interaction list
  integer, dimension(1:n) ::  inode
  integer :: jnode, i

  real, intent(out) ::  sumfx,sumfy,sumfz,sumphi 
  real :: a_bond, d2, d, dlj2, dlj, flj, epsc, plj, dx, dy, dz

 ! mean interparticle spacing - slightly bigger to push particles onto boundaries
  a_bond=a_ii*sqrt(2.)
  epsc = 0.8 ! cutoff
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
     dlj = d/a_bond

     if (dlj >= epsc) then 
       plj = 4*(1./dlj**12 - 1./dlj**6)
       flj = 24./a_bond*(2./dlj**13 - 1./dlj**7)
!     else if (dlj >= 0.5) then
!       plj = 48./pi*cos(pi*dlj/2.) 
 !      flj = 24./a_bond*sin(pi*dlj/2.)
     else
       plj = 4*(1./epsc**12 - 1./epsc**6)
       flj = 24./a_bond*(2./epsc**13-1./epsc**7)
     endif
     !   natoms = abs_charge( jnode) / qi    ! # particles in monopole cluster


     ! potential
     sumphi = sumphi + plj

     !  forces

     sumfx = sumfx + dx/d*flj

     sumfy = sumfy + dy/d*flj

     sumfz = sumfz + dz/d*flj


  end do


end subroutine sum_lennardjones
