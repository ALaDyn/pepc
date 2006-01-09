!  ===================================================================
!
!                              SUM_BFIELD
!
!   Calculate vector potential and b-field of particle from interaction list
!
!  ===================================================================

subroutine sum_bfield( p, n, inode, eps, sumbx, sumby, sumbz, sumax, sumay, sumaz )
  use treevars
  implicit none
  integer, intent(in) :: p  ! particle label 
  integer, intent(in) :: n  !  # terms on interaction list
  integer, dimension(1:n) ::  inode
  real, intent(in) :: eps ! smoothing parameter
  integer :: jnode, i,j,k 

  real*8 :: rd,dx,dy,dz,d,dx2,dy2,dz2 
 real*8 :: dx3,dy3,dz3,rd3,rd5,rdotj,rdotm
 real*8 :: fsx,fsy,fsz,phi
 real, dimension(n*10) :: mult 
 real, dimension(n*3) :: coc
  real*8, intent(out) ::  sumax,sumay,sumaz,sumbx, sumby, sumbz 
  real :: eps2

  eps2=eps**2
  sumax = 0
  sumay = 0
  sumaz = 0
  sumbx = 0
  sumby = 0
  sumbz = 0

! copy multipole moments into stride 1 array
  do j=1,n
    jnode=inode(j)
    i = (j-1)*10+1
!    k = (j-1)*3+1
!    coc(k) = xcoc(jnode)
!    coc(k+1) = ycoc(jnode)
!    coc(k+2) = zcoc(jnode)
    mult(i) = jx(jnode)
    mult(i+1) = jy(jnode)
    mult(i+2) = jz(jnode)
    mult(i+3) = magmx(jnode)
    mult(i+4) = magmy(jnode)
    mult(i+5) = magmz(jnode)
  end do

  do j=1,n
     
  !  preprocess distances
     i = 10*(j-1) + 1  ! multipole index
!     k = 3*(j-1) + 1  ! coc index
     jnode=inode(j)
     dx = x(p) - xcoc(jnode)
     dy = y(p) - ycoc(jnode)
     dz = z(p) - zcoc(jnode) 

     d = sqrt(dx**2+dy**2+dz**2+eps2)
     rd = 1./d
     rd3 = rd**3
     rd5 = rd**5

     rdotm = dx*magmx(jnode) + dy*magmy(jnode) + dz*magmz(jnode)
     rdotj = dx*jx(jnode) + dy*jy(jnode) + dz*jz(jnode)

     ! vector potential (magnetoinductive)

     sumax = sumax + 0.5*jx(jnode)*rd + 0.5*dx*rdotj*rd3        !  monopole term
!     + ( dz*magmy(jnode) - dy*magmz(jnode) )*rd3   !  dipole 

     sumay = sumay + 0.5*jy(jnode)*rd + 0.5*dy*rdotj*rd3          !  monopole term
!     + ( dx*magmz(jnode) - dz*magmx(jnode) )*rd3  !  dipole 

     sumaz = sumaz + 0.5*jz(jnode)*rd + 0.5*dz*rdotj*rd3         !  monopole term
!     + ( dy*magmx(jnode) - dx*magmy(jnode) )*rd3  !  dipole 



     ! B-field

     sumbx = sumbx + ( jy(jnode)*dz - jz(jnode)*dy)*rd3       ! monopole term
!       + 3*dx*rd5*rdotm - magmx(jnode)*rd3

     sumby = sumby + ( jz(jnode)*dx - jx(jnode)*dz)*rd3       ! monopole term
!       + 3*dy*rd5*rdotm - magmy(jnode)*rd3

     sumbz = sumbz + ( jx(jnode)*dy - jy(jnode)*dx)*rd3       ! monopole term
!       + 3*dz*rd5*rdotm - magmz(jnode)*rd3




  end do



end subroutine sum_bfield
