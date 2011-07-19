! ======================
!
!   SLICES
!
!   Write out field data for 
!   1D postprocessing

!
!
! ======================

subroutine slices(timestamp)


  use treevars
  implicit none   

  real, dimension(ngx) :: work1, work2
  real, dimension(ngx) :: phi_pond, ex_pond, ey_pond, ez_pond
  real :: dx, dz, dy, xd, yd, zd, dummy, simtime, epondx, epondy, epondz, phipond, epond_max, box_max
  real :: Qtot, Qbox

  character(30) :: cfile
  character(5) :: cme
  character(6) :: cdump, cvis
  integer, intent(in) :: timestamp
  integer :: i, j, k, ioffset, idummy=0, ilev, ixd, iyd, izd, npx, npz, npy
  integer :: icall, lcount
  integer :: jfoc, kfoc

  icall = timestamp/ivis
  simtime = timestamp*dt

  ! get filename suffix from dump counter
  do i=0,4
     cdump(6-i:6-i) =  achar(mod(timestamp/10**i,10) + 48)  
  end do
  cdump(1:1) = achar(timestamp/10**5 + 48)

  if (me==0) then
     !  Fields: electron and ion densities within xl*yl*zl
     cfile = "wf_field."//cdump
     open (62,file=cfile)

     dx = xl/ngx
     dy = yl/ngy
     dz = zl/ngz
     Qbox = 0.
     do k=1,ngz
        do j=1,ngy
           do i=1,ngx
              Qbox = Qbox + rhoe(i,j,k)*dx*dy*dz
              write(62,'(3f13.5,2e13.3)') i*dx,j*dy,k*dz,abs(rhoe(i,j,k)),rhoi(i,j,k)
           end do
        end do
     end do
     Qtot = SUM(rhoe)*dx*dy*dz  ! including ghost cells
     if (debug_level>1) write(ipefile,'(4(a,f14.5/))') &
          'Total charge on grid:',Qbox, &
          '         ghost cells:',Qtot-Qbox, &
          '                 sum:',Qtot, &
          'Initial charge Q_s*Ne = rho0*V = ',Vplas*rho0

     close(62)


     cfile = "field_slice."//cdump
     open (62,file=cfile)

     jfoc = focus(2)/dy
     kfoc = focus(3)/dz

     ! density average line-out along laser axis: 5x5 average, converted to n/nc

     work1 = 0.
     work2 = 0.
     do k=kfoc-2,kfoc+2
        do j=jfoc-2,jfoc+2
           work1(1:ngx) = work1(1:ngx)+rhoi(1:ngx,j,k)/25./omega**2  ! density slice along laser axis: 5x5 average
           work2(1:ngx) = work2(1:ngx)+rhoe(1:ngx,j,k)/25./omega**2
        end do
     end do

     do i=1,ngx
        xd=i*dx-focus(1)
        yd=sigma/2.
        zd=sigma/2.
        laser_model: select case(beam_config)

        case(4)  ! standing wave fpond
           call fpond(1.57/omega,1.0,sigma,vosc,omega,xd,yd,zd,ex_pond(i),ey_pond(i), &
                ez_pond(i), phi_pond(i))
        case(5)  ! propagating fpond
           call laser_bullet( tpulse,sigma,vosc, & 
                xd,yd,zd,ex_pond(i),ey_pond(i),ez_pond(i),phi_pond(i))
        case default
           ex_pond(i)=0
           ey_pond(i)=0
           ez_pond(i)=0
           phi_pond(i)=0
        end select laser_model

     end do

     write(62,'(7f13.5)') (i*dx,work1(i),work2(i), &
          phi_pond(i),ex_pond(i), ey_pond(i), ez_pond(i),i=1,ngx)
     close(62)

  endif
  icall = icall + 1

end subroutine slices



