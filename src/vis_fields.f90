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


  real, dimension(ngx*ngy*ngz) :: qvis,mvis,xvis,yvis,zvis,exvis,eyvis,ezvis, cvis
  integer, dimension(ngx*ngy*ngz) :: labvis, idvis
  real :: s, simtime, dummy, xd,yd,zd, dx, dz, dy, epond_max,&
       box_max, epondx, epondy, epondz,phipond
  integer :: i, j, k, ioffset,ixd, iyd, izd, ilev, lcount, iskip
  integer :: lvisit_active
  integer :: npx, npy, npz, ng
  if (me==0 ) then

     box_max = xl/25.
     epond_max = sqrt(1.+vosc**2/2.)/2.
     lcount=0
     iskip=1
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
              xvis(lcount) = convert_mu*(i-0.5)*dx
              yvis(lcount) = convert_mu*(j-0.5)*dy
              zvis(lcount) = convert_mu*(k-0.5)*dz
              labvis(lcount) = lcount
              xd = (i-0.5)*dx - focus(1) ! position relative to laser focus
              yd = (j-0.5)*dy - focus(2)
              zd = (k-0.5)*dz - focus(3)

              laser_model: select case(beam_config)
              case(4)  ! standing wave fpond
                 call fpond( tlaser, tpulse,sigma,vosc,omega,xd,yd,zd,epondx,epondy,epondz,phipond)
              case(5)  ! propagating fpond
                 call laser_bullet( tlaser, focus(1), tpulse,sigma,vosc,omega, & 
                      xd,yd,zd,epondx,epondy,epondz,phipond)
              case default
                 Epondx=0
                 Epondy=0
                 Epondz=0
                 phipond=0
              end select laser_model

              mvis(lcount) = 10*phipond/vosc**2  
              qvis(lcount) = 10*phi_g(i,j,k)   ! potential
              exvis(lcount) = 10*ex_g(i,j,k)
              eyvis(lcount)  = 10*ey_g(i,j,k)
              ezvis(lcount)  = 10*ez_g(i,j,k)
              cvis(lcount) = -1.
 !             qvis(lcount) = 10*(rhoi(i,j,k) + rhoe(i,j,k)) ! delta-n/n
           end do
        end do
     end do

     ng = ngx*ngy*ngz

     call flvisit_spk_check_connection(lvisit_active)
     call flvisit_spk_info_send(ng,xl*convert_mu,yl*convert_mu,zl*convert_mu, &
         zl*convert_mu,ng,0,np_beam,itime+itime_start)
     if (beam_config==4 .or. beam_config ==5) then
        call flvisit_spk_3dfieldB_send(mvis,ngx,ngy,ngz)  ! laser potential
     endif
     call flvisit_spk_3dfieldA_send(qvis,ngx,ngy,ngz)  ! electron density

! Ship electric field as particle data
!     call flvisit_spk_particles_send(tlaser*convert_fs,xvis,yvis,zvis, & 
!     exvis,eyvis,ezvis,cvis,idvis,labvis,ng)

  endif


end subroutine vis_fields





