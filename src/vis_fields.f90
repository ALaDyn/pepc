! ======================
!
!   VIS_FIELDS
!
!   Send field data to VISIT for surface visualisation
!
!
! ======================

subroutine vis_fields


  use treevars
  implicit none   


  real, dimension(ngx*ngy*ngz) :: qvis,mvis
  real :: s, simtime, dummy, xd,yd,zd, dx, dz, dy, epond_max,&
          box_max, epondx, epondy, epondz,phipond
  integer :: i, j, k, ioffset,ixd, iyd, izd, ilev, lcount, iskip
  integer :: lvisit_active
  integer :: npx, npy, npz
  if (me==0 ) then

     box_max = xl/25.
     epond_max = sqrt(1.+vosc**2/2.)/2.
     lcount=0
     iskip=2
     npx = ngx/iskip + mod(ngx,2)
     npy = ngy/iskip + mod(ngy,2)
     npz = ngz/iskip + mod(ngz,2)
     dx = xl/ngx
     dz = zl/ngz
     dy = yl/ngy
     do k=1,ngz,iskip
        do j=1,ngy,iskip
           do i=1,ngx,iskip
              lcount=lcount+1
              xd = (i-0.5)*dx - focus(1) ! position relative to laser focus
              yd = (j-0.5)*dy - focus(2)
              zd = (k-0.5)*dz - focus(3)
              call fpond(tlaser,tpulse,sigma,vosc,omega,-xd,yd,zd,epondx,epondy,epondz,phipond)  ! evaluate pond force
              mvis(lcount) = phipond  
              qvis(lcount) = abs(rhoe(i,j,k))/omega**2   ! electron density
           end do
        end do
     end do
     call flvisit_spk_check_connection(lvisit_active)
     call flvisit_spk_info_send(npart,xl,yl,zl,zl,ne,ni,np_beam,itime+itime_start)
!     if (beam_config==4) call flvisit_spk_3dfieldB_send(mvis,npx,npy,npz)  ! laser potential
!     call flvisit_spk_3dfieldA_send(qvis,npx,npy,npz)  ! ion density
  endif


end subroutine vis_fields





