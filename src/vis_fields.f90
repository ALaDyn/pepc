! ======================
!
!   VIS_FIELDS
!
!   Send field data to VISIT for surface visualisation
!
!
! ======================

subroutine vis_fields


  use physvars
  use treevars
  implicit none   


  real, dimension(ngx*ngy*ngz) :: qvis,mvis
  real, dimension(0:ngx+1,0:ngy+1,0:ngz+1) :: bzg
  real :: s, simtime, dummy, xd,yd,zd, dx, dz, dy, epond_max
  real :: box_max, az_em, ez_em, by_em, bx_em, phipond
  real :: epon_x, epon_y, epon_z, tpon, amp_las
  integer, parameter :: ngmax=100
  integer :: i, j, k, ioffset,ixd, iyd, izd, ilev, lcount, iskip,itlas
  integer :: lvisit_active, ierr
  integer :: npx, npy, npz, ng
  integer :: iskip_x, iskip_y, iskip_z

  simtime = dt*(itime+itime_start)
  amp_las = vosc*min(1.,simtime/tpulse)

  ng = (ngx+2)*(ngy+2)*(ngz+2)                         ! total # gridpoints
  ! Merge sums for Bz
  call MPI_ALLREDUCE(bz_loc, bzg, ng, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

  if (me==0 ) then

     box_max = xl/25.
     epond_max = sqrt(1.+vosc**2/2.)/2.
     lcount=0

     ! limit size of field data to 100^3
     !     iskip_x = ngx/ngmax + 1
     !     iskip_y = ngy/ngmax + 1
     !     iskip_z = ngz/ngmax + 1
     iskip_x = 1
     iskip_y = 1
     iskip_z = 1


     !     npx = ngx/iskip_x + mod(ngx,2)
     !     npy = ngy/iskip_y + mod(ngy,2)
     !     npz = ngz/iskip_z + mod(ngz,2)
     npx = ngx
     npy = ngy
     npz = ngz

     dx = xl/ngx
     dz = zl/ngz
     dy = yl/ngy
     do k=1,ngz,iskip_z
        do j=1,ngy,iskip_y
           do i=1,ngx,iskip_x
              lcount=lcount+1
              xd = (i-0.5)*dx - focus(1) ! position relative to laser focus
              yd = (j-0.5)*dy - focus(2)
              zd = (k-0.5)*dz - focus(3)
!              xd = (i-0.5)*dx - 50.

              laser: select case(beam_config)

              case(4)
                 call fpond( 1.57/omega, tpulse,sigma,vosc,omega, rho_upper, &
                      xd,yd,zd,epon_x,epon_y,epon_z,phipond)
                 Tpon = min(1.,tlaser/tpulse) * (sin(omega*tlaser))**2
                 !                 mvis(lcount) = epon_x/omega ! Pond field, EM norm
                 mvis(lcount) = phipond ! Pond potential

              case(5)  ! propagating fpond
                 call laser_bullet( tlaser, focus(1), tpulse,sigma,vosc,omega, &
                      xd,yd,zd,epon_x,epon_y,epon_z,phipond)
                 mvis(lcount) = phipond ! Pond potential

              case(7) ! Oblique incidence fpond
                 call emobliq(tlaser,tpulse,sigma,vosc,omega,theta_inc, &
                      xd,yd,zd,epon_x,epon_y,epon_z,phipond)
                 mvis(lcount) = phipond

              case(6) ! Plane wave
                 call emplane(tlaser,tpulse,sigma,vosc,omega,xd,yd,zd,ez_em,by_em,bx_em,az_em,phipond)
                 mvis(lcount) = by_em 

              case default ! Propagating fpond
                 phipond = 0  
                 mvis(lcount) = 0.
                 mvis(lcount) = rhoi(i,j,k)/omega**2   ! electron density /nc
              end select laser

              !              qvis(lcount) = ez_em 
              qvis(lcount) = 100*bzg(i,j,k)   ! Bz
           end do
        end do
     end do
     call flvisit_spk_check_connection(lvisit_active)
     !     call flvisit_spk_info_send(npart,xl,yl,zl,zl,ne,ni,np_beam,itime+itime_start)
     itlas=int(tlaser)

     !        call flvisit_spk_info_send(npart,xl,yl,zl, zl, &
     !             x_crit, amp_las, sigma, tpulse, &
     !             ne,ni,npart,itlas)

     !     if (beam_config>=4) call flvisit_spk_3dfieldA_send(mvis,npx,npy,npz)  ! laser potential
     call flvisit_spk_3dfieldB_send(mvis,npx,npy,npz)  ! ion density 
     call flvisit_spk_3dfieldA_send(qvis,npx,npy,npz)  ! ship selected field
  endif


end subroutine vis_fields

