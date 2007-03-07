! ======================
!
!   DUMP_FIELDS
!
!   Write out field data for 
!     postprocessing

!
!
! ======================

subroutine dump_fields(timestamp)


  use physvars
!  use treevars
  implicit none 
  include 'mpif.h'

  real, dimension(ngx) :: work1, work2
  real, dimension(0:nxh) :: phi_pond, ex_pond, ey_pond, ez_pond
  real, dimension(0:ngx+1) :: rhoi_slice, rhoe_slice, ex_slice, ey_slice, ez_slice
  real, dimension(0:ngx+1) :: jxe_slice, jye_slice, jze_slice 
  real, dimension(0:ngx+1,0:ngy+1,0:ngz+1) :: exg, eyg, ezg, jxeg, jyeg, jzeg

  real :: dx, dz, dy, xd, yd, zd, dummy, simtime, epon_x, epon_y, epon_z, phipond, epond_max, box_max
  real :: uxd
  real :: Qtot, Qbox, norm, rhonorm, tpon, bx_em, by_em, az_em,ez_em

  character(30) :: cfile
  character(5) :: cme
  character(6) :: cdump, cvis
  integer, intent(in) :: timestamp
  integer :: i, j, k, ioffset, idummy=0, ilev, ixd, iyd, izd, npx, npz, npy
  integer :: icall, lcount, ierr
  integer :: jfoc, kfoc, ng, nave

  icall = timestamp/ivis
  simtime = timestamp*dt

  ! Gather locally accumulated averages
  ng = (ngx+2)*(ngy+2)*(ngz+2)                         ! total # gridpoints

  call MPI_ALLREDUCE(rhoe_loc, rhoe, ng, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(rhoi_loc, rhoi, ng, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(jxe_loc, jxeg, ng, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(jye_loc, jyeg, ng, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(jze_loc, jzeg, ng, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(ex_loc, exg, ng, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(ey_loc, eyg, ng, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(ez_loc, ezg, ng, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)

  ! get filename suffix from dump counter
  do i=0,4
     cdump(6-i:6-i) =  achar(mod(timestamp/10**i,10) + 48)  
  end do
  cdump(1:1) = achar(timestamp/10**5 + 48)

  if (my_rank==0) then
     write(*,'(//a/a,f10.1,a1,f10.1,a1,f10.1)') 'Cycle-average field dump:', &
          'writing out densities, fields on grid ',xl,'*',yl,'*',zl
     ! dump electron and ion densities, electron currents, DC fields within xl*yl*zl
     cfile = "fields/"//cdump//".xyz"
     open (62,file=cfile)
     cfile = "fields/laser_"//cdump//".xyz"
     open (63,file=cfile)
     dx = xl/ngx
     dy = yl/ngy
     dz = zl/ngz
     Qbox = 0.
     do k=1,ngz
        do j=1,ngy
           do i=1,ngx
              Qbox = Qbox + rhoe(i,j,k)*dx*dy*dz
              write(62,'(8e13.3)') rhoe(i,j,k)/omega**2,rhoi(i,j,k)/omega**2, &
                   jxeg(i,j,k), jyeg(i,j,k), jzeg(i,j,k), &
                   exg(i,j,k),  eyg(i,j,k), ezg(i,j,k)
!              xd = (i-0.5)*dx - 20.-focus(1) ! position relative to laser focus
              xd = (i-0.5)*dx - focus(1) ! position relative to laser focus
              yd = (j-0.5)*dy - focus(2)
              zd = (k-0.5)*dz - focus(3)

              laser: select case(beam_config)

              case(4)
                 call fpond( 1.57/omega, tpulse,sigma,vosc,omega, rho_upper, &
                      xd,yd,zd,epon_x,epon_y,epon_z,phipond)
                 Tpon = min(1.,tlaser/tpulse) * (sin(omega*tlaser))**2
                 !                 mvis(lcount) = epon_x/omega ! Pond field, EM norm

              case(5)  ! propagating fpond
                 call laser_bullet( tlaser, focus(1), tpulse,sigma,vosc,omega, &
                      xd,yd,zd,epon_x,epon_y,epon_z,phipond)

              case(24) ! Standing wave fpond Ez, By, Az
                 call emobliq(tlaser,tpulse,sigma,vosc,omega,theta_inc, &
                      xd,yd,zd,epon_x,epon_y,epon_z,phipond)

              case(6) ! Plane wave
                 call emplane(tlaser,tpulse,sigma,vosc,omega,xd,yd,zd,ez_em,by_em,bx_em,az_em,phipond)

              case default ! Propagating fpond
                 phipond = 0  

              end select laser
              write(63,'(4e13.3)') phipond,epon_x,epon_y,epon_z
           end do
        end do
     end do
     Qtot = SUM(rhoe)*dx*dy*dz  ! including ghost cells
     write(ifile_cpu,'(4(a,f14.5/))') &
          'Total charge on grid:',Qbox, &
          '         ghost cells:',Qtot-Qbox, &
          '                 sum:',Qtot, &
          'Initial charge Q_s*Ne = rho0*Vplas = ',Vplas*rho0
     write(*,'(4(a,f14.5/))') &
          'Total charge on grid:',Qbox, &
          '         ghost cells:',Qtot-Qbox, &
          '                 sum:',Qtot, &
          'Initial charge Q_s*Ne = rho0*Vplas = ',Vplas*rho0

     close(62)
     close(63)


     ! x-slices along laser axis
     jfoc = focus(2)/dy
     kfoc = focus(3)/dz

     ! density average line-out along laser axis: nave*nave average, converted to n/nc

     rhoe_slice = 0.
     rhoi_slice = 0.
     ex_slice = 0.
     jxe_slice = 0.

     if (ngz<=5) then
        nave=0
     else
        nave = min(3,ngz/2)
     endif

     norm = (2*nave+1)**2
     rhonorm = norm*omega**2

     do k=kfoc-nave,kfoc+nave
        do j=jfoc-nave,jfoc+nave

           rhoe_slice(1:ngx) = rhoe_slice(1:ngx)+rhoe(1:ngx,j,k)/norm  ! density slice along laser axis: 5x5 average
           rhoi_slice(1:ngx) = rhoi_slice(1:ngx)+rhoi(1:ngx,j,k)/norm
           ex_slice(1:ngx) = ex_slice(1:ngx)+exg(1:ngx,j,k)/norm
           jxe_slice(1:ngx) = jxe_slice(1:ngx)+jxeg(1:ngx,j,k)/norm

        end do
     end do

     ! Laser fields on axis
     do i=1,ngx
        xd=i*dx-focus(1)
        !        yd=sigma/2.
        !        zd=sigma/2.
        yd = 0.
        zd = 0.
        call fpond( 1.57/omega,1.0,sigma,vosc,omega,rho_upper,xd,yd,zd,ex_pond(i),ey_pond(i), &
             ez_pond(i), phi_pond(i))
     end do
     ! Renormalise to EM units
     cfile = "fields/xslice."//cdump
     open (62,file=cfile)
     write(62,'(a)') '!   x      rho_e   rho_i  ex, jxe,  phi_p,  ex_p,  ey_p,   ez_p  '
     write(62,'((9(1pe12.4)))') &
          (i*dx+x_offset,rhoe_slice(i)/omega**2, rhoi_slice(i)/omega**2, ex_slice(i)/omega,&
          jxe_slice(i), phi_pond(i),ex_pond(i)/omega, ey_pond(i)/omega, ez_pond(i)/omega,i=1,ngx)
     close(62)

! Write out time-averaged E-field profiles to file
     write(*,'(a40,a10)') 'Cycle-averaged lineouts at',cdump
     dx = (xl-x_offset)/ngx ! spacing for time-ave grid
     cfile = "fields/linave."//cdump
     open (60,file=cfile)
     write(60,'(2(a12))') '!   x      ',' ex  '
     write(60,'((2(1pe12.4)))') &
          (i*dx+x_offset, ex_ave(i), i=0,ngx)
     close(60)

! Laser fields on Helmholtz grid
     do i=0,nxh
        xd=i*dxh+xh_start
        !        yd=sigma/2.
        !        zd=sigma/2.
        yd = 0.
        zd = 0.
        uxd = 0.
        call fpond_helm( tlaser, tpulse,sigma,vosc,omega, &
               xd,yd,zd,uxd,Az_helm,nxh,xh_start, xh_end, dxh, focus(1), &
	       ex_pond(i),ey_pond(i),ez_pond(i),phi_pond(i))


     end do
! Write out Helmholtz fields to file
     write(*,'(a40,a10)') 'Helmholtz lineouts at',cdump
     cfile = "fields/helmholtz."//cdump
     open (60,file=cfile)
     write(60,'(2(a12))') '!   x_helm  ',' rho         az^2  '
     dxh = (xh_end-xh_start)/nxh
     write(60,'((4(1pe12.3)))') &
!          (i*dxh+xh_start-xl/2.+x_plasma/2., rho_helm(i), abs(az_helm(i)**2), i=0,nxh)
          (i*dxh+xh_start, rho_helm(i), abs(az_helm(i)),ex_pond(i), i=0,nxh)
     close(60)
  endif

  icall = icall + 1

  ! Rezero local fields
  rhoe_loc = 0.

  ex_ave = 0.
  ex_loc = 0.
  ey_loc = 0.
  ez_loc = 0.
  jxe_loc = 0.
  jye_loc = 0.
  jze_loc = 0.


end subroutine dump_fields













