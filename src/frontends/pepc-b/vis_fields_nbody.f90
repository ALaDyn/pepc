! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2012 Juelich Supercomputing Centre, 
!                         Forschungszentrum Juelich GmbH,
!                         Germany
! 
! PEPC is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! PEPC is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
! 
! You should have received a copy of the GNU Lesser General Public License
! along with PEPC.  If not, see <http://www.gnu.org/licenses/>.
!

! ======================
!
!   VIS_FIELDS
!
!   Send field data to VISIT for visualisation of scalar volumetric data
!
!
! ======================

subroutine vis_fields_nbody(timestamp)


  use physvars
  use treevars
  use module_utils
  implicit none   
  include 'mpif.h'

  real*4, dimension(ngx*ngy*ngz) :: field1, field2, field3, field4
  real*4, dimension(ngx,ngy,ngz) :: bzg, rhoeg, rhoig, Telec, Tion, ggelec,ggion
  real, dimension(0:ngx+1) :: te_slice 
  real :: s, simtime, dummy, xd,yd,zd, dx, dz, dy, epond_max
  real :: box_max, az_em, ez_em, by_em, bx_em, bz_em, ex_em, ey_em, phipond
  real :: epon_x, epon_y, epon_z, tpon, amp_las, field_laser
  integer, parameter :: ngmax=100
  integer :: i, j, k, ioffset,ixd, iyd, izd, ilev, lcount, iskip,itlas
  integer, intent(in) :: timestamp
  integer :: lvisit_active, ierr 
  integer :: npx, npy, npz, ng, jfoc, kfoc, nave
  real :: norm, uxd
  integer :: iskip_x, iskip_y, iskip_z
  integer :: fselect1=0,fselect2=0,fselect3=0,fselect4=0
  integer :: incdf
  real*4 :: grid_pars(24)  ! origins and mesh sizes of vis fields
  character(30) :: cfile
  character(5) :: cme
  character(6) :: cdump, cvis

  simtime = dt*(itime+itime_start)
  amp_las = vosc*min(1.,simtime/tpulse)

#ifdef VISIT_NBODY
  if (me==0)   call flvisit_nbody2_check_connection(lvisit_active)
  call MPI_BCAST( lvisit_active, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
#endif

  if (lvisit_active==0 )then
     if (me==0) write(*,'(a)') 'VIS_FIELDS  | No connection to visualization'
  endif

! Default field selection 
  fselect1=1
  fselect2=4
  fselect3=5
  fselect4=0

! Connected to vis, so proceed with field select & gather

#ifdef VISIT_NBODY
! Fetch user-selected config from vis
  if (me==0 .and. lvisit_active.ne.0) then
	call flvisit_nbody2_selectfields_recv(fselect1,fselect2,fselect3,fselect4)

  endif
#endif

  if (me==0) write(*,'(a,4i6)') 'VIS_FIELDS  | Selection:',fselect1,fselect2,fselect3,fselect4

  ! get filename suffix from dump counter
  do i=0,4
     cdump(6-i:6-i) =  achar(mod(timestamp/10**i,10) + 48)  
  end do
  cdump(1:1) = achar(timestamp/10**5 + 48)

  ng = ngx*ngy*ngz                         ! total # gridpoints
  ! Merge sums for gridded fields 

!  call MPI_ALLREDUCE(bz_loc(1:ngx,1:ngy,1:ngz), bzg, ng, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(rhoe_loc(1:ngx,1:ngy,1:ngz), rhoeg, ng, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(rhoi_loc(1:ngx,1:ngy,1:ngz), rhoig, ng, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(Te_loc(1:ngx,1:ngy,1:ngz), Telec, ng, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(Ti_loc(1:ngx,1:ngy,1:ngz), Tion, ng, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(g_ele(1:ngx,1:ngy,1:ngz), ggelec, ng, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(g_ion(1:ngx,1:ngy,1:ngz), ggion, ng, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)

! Normalise temperatures & convert to keV
  Telec = 2./3.*Telec/ggelec
!  Tion = 2./3.*Tion/ggion

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

  !   xl_laser = focus(1) + x_plasma/2. ! Reduced box size for laser
     dx = xl/ngx
     dz = zl/ngz
     dy = yl/ngy

     if (itime>0 .and. beam_config==4) focus(1) = x_crit  ! laser tracks n_c
     do k=1,ngz,iskip_z
        do j=1,ngy,iskip_y
           do i=1,ngx,iskip_x
              lcount=lcount+1
              xd = (i-0.5)*dx - focus(1) ! position relative to laser focus
              yd = (j-0.5)*dy - focus(2)
              zd = (k-0.5)*dz - focus(3)
!              xd = (i-0.5)*dx - 50.

              laser: select case(beam_config_in)

              case(4)
                 call fpond( 1.57/omega, tpulse,sigma,vosc,omega, rho_upper, &
                      xd,yd,zd,epon_x,epon_y,epon_z,phipond)
                 Tpon = min(1.,tlaser/tpulse) * (sin(omega*tlaser))**2
                 !                 mvis(lcount) = epon_x/omega ! Pond field, EM norm
                 field_laser = phipond ! Pond potential

              case(14)  ! linear build-up
                 call fpond_lin( tlaser, tpulse,sigma,vosc,omega, rho_upper, &
                      xd,yd,zd,epon_x,epon_y,epon_z,phipond)
!                 Tpon = min(1.,tlaser/tpulse) * (sin(omega*tlaser))**2
                 Tpon = min(1.,tlaser/tpulse) * sin(omega*tlaser-omega*xd)
                 field_laser = phipond*Tpon ! Pond potential

              case(94)  ! standing wave fpond with transverse fields artificially reduced
                 call fpond( tlaser, tpulse,sigma,vosc,omega,rho_upper, &
                 xd,yd,zd,epon_x,epon_y,epon_z,phipond)
                 Tpon = min(1.,tlaser/tpulse) * (sin(omega*tlaser))**2
                 field_laser = phipond*Tpon ! Pond potential

              case(44,54)  ! fpond derived from Az_helm; both linear & sin2 pulse forms
                 xd = (i-0.5)*dx ! absolute position
                 uxd = 0.
                 call fpond_helm( tlaser, tpulse,sigma,vosc,omega, &
                   xd,yd,zd,uxd,Az_helm,nxh,xh_start, xh_end, dxh, focus(1), &
                   epon_x,epon_y,epon_z,phipond)
                 field_laser = epon_x ! Pond potential

              case(5)  ! propagating fpond
                 call laser_bullet(tpulse,sigma,vosc,xd,yd,zd,epon_x,epon_y,epon_z,phipond)
                 field_laser = phipond ! Pond potential

              case(24) ! Oblique incidence fpond, s-pol
                 call emobliq(tlaser,tpulse,sigma,vosc,omega,theta_beam, rho_upper, &
                      xd,yd,zd,epon_x,epon_y,epon_z,phipond,ez_em,bx_em,by_em)
                 field_laser = ez_em**2

              case(6) ! Plane wave
                 call emplane(tlaser,tpulse,sigma,vosc,omega,xd,yd,zd,ez_em,by_em,bx_em,az_em,phipond)
                 field_laser = by_em 

              case(16) ! plane wave with Gaussian spot, linear rise-time
                 call emplane_lin(tlaser,tpulse,sigma,vosc,omega,xd,yd,zd,Epon_z,By_em,Bx_em,az_em,phipond)
	         Epon_x=0.
	         Epon_y=0.
	         field_laser = by_em

              case default ! No laser
                 phipond = 0  
                 field_laser = 0.
              end select laser

 ! Vis field selection
              f1: select case(fselect1)
	      case(1)  ! ion density / nc 
                field1(lcount) = rhoig(i,j,k)/omega**2   
	      case(2)  ! electron density / nc 
                field1(lcount) = rhoeg(i,j,k)/omega**2 
	      case(3)  ! Ion temperature  / MeV 
                field1(lcount) = Tion(i,j,k) 
	      case(4)  ! Electron temperature  / MeV 
                field1(lcount) = Telec(i,j,k) 
	      case(5)  ! laser 
                field1(lcount) = field_laser 
	      case(0) 
	        field1(lcount) = 0
	      end select f1

              f2: select case(fselect2)
	      case(1)  ! ion density / nc 
                field2(lcount) = rhoig(i,j,k)/omega**2   
	      case(2)  ! electron density / nc 
                field2(lcount) = rhoeg(i,j,k)/omega**2 
	      case(3)  ! Ion temperature  / MeV 
                field2(lcount) = Tion(i,j,k) 
	      case(4)  ! Electron temperature  / MeV 
                field2(lcount) = Telec(i,j,k) 
	      case(5)  ! laser 
                field2(lcount) = field_laser 
              case(0)
                field2(lcount) = 0
	      end select f2

              f3: select case(fselect3)
	      case(1)  ! ion density / nc 
                field3(lcount) = rhoig(i,j,k)/omega**2   
	      case(2)  ! electron density / nc 
                field3(lcount) = rhoeg(i,j,k)/omega**2 
	      case(3)  ! Ion temperature  / MeV 
                field3(lcount) = Tion(i,j,k) 
	      case(4)  ! Electron temperature  / MeV 
                field3(lcount) = Telec(i,j,k) 
	      case(5)  ! laser 
                field3(lcount) = field_laser 
              case(0)
                field3(lcount) = 0
	      end select f3

              f4: select case(fselect4)
	      case(1)  ! ion density / nc 
                field4(lcount) = rhoig(i,j,k)/omega**2   
	      case(2)  ! electron density / nc 
                field4(lcount) = rhoeg(i,j,k)/omega**2 
	      case(3)  ! Ion temperature  / MeV 
                field4(lcount) = Tion(i,j,k) 
	      case(4)  ! Electron temperature  / MeV 
                field4(lcount) = Telec(i,j,k) 
	      case(5)  ! laser 
                field4(lcount) = field_laser 
              case(0)
                field4(lcount) = 0
	      end select f4
 	   if (field1(lcount)<0) then
	    write(*,*) 'field -ve ',field1(lcount),' at i=',i,' j=',j,' k=',k
  	   endif
           end do
        end do
     end do


#ifdef VISIT_NBODY
      write(*,'(a)') 'VIS_FIELDS  | confirm field select'
      call flvisit_nbody2_check_connection(lvisit_active)

! Tell vis which fields are coming

  if (lvisit_active.ne.0) then
	 call flvisit_nbody2_selectedfields_send(fselect1,fselect2,fselect3,fselect4)
  endif
#endif

#ifdef NETCDFLIB
! Netcdf write
	write(*,'(a)') 'VIS_FIELDS  | NETCDF: confirm field select'
         if (netcdf) call ncnbody_putselfield( ncid, simtime, fselect1, fselect2, fselect3, fselect4, incdf )
#endif

!  Set up vis field grid
   do i=0,3
     grid_pars(6*i+1:6*i+3) = 0.
     grid_pars(6*i+4) = dx
     grid_pars(6*i+5) = dy
     grid_pars(6*i+6) = dz
   end do

#ifdef VISIT_NBODY
   call flvisit_nbody2_fielddesc_send(grid_pars,4,6)
#endif

#ifdef NETCDFLIB
	write(*,'(a)') 'VIS_FIELDS  | NETCDF: grid-pars'
   if (netcdf) call ncnbody_putfielddesc( ncid, simtime, grid_pars, incdf )
#endif

!   write(*,*) 'Grids: ',grid_pars
      if (fselect1>0) then
#ifdef VISIT_NBODY
       	 write (*,'(a)') "VIS_FIELDS  | Shipping field 1: min/max =", &
      minval(field1),maxval(field1)
         call flvisit_nbody2_field1_send(field1,npx,npy,npz)   
#endif
#ifdef NETCDFLIB
       	 write (*,'(a)') "VIS_FIELDS  | Writing field 1 to netcdf"
         if (netcdf) call ncnbody_putfield( ncid, simtime, 1, npx, npy, npz, field1, incdf )
#endif
      endif

      if (fselect2>0) then
#ifdef VISIT_NBODY
       	 write (*,'(a)') "VIS_FIELDS  | Shipping field 2"
         call flvisit_nbody2_field2_send(field2,npx,npy,npz)   
#endif
#ifdef NETCDFLIB
       	 write (*,'(a)') "VIS_FIELDS  | Writing field 2 to netcdf"
         if (netcdf) call ncnbody_putfield( ncid, simtime, 2, npx, npy, npz, field2, incdf )
#endif
      endif

      if (fselect3>0) then
#ifdef VISIT_NBODY
       	 write (*,'(a)') "VIS_FIELDS  | Shipping field 3"
         call flvisit_nbody2_field3_send(field3,npx,npy,npz)  
#endif
#ifdef NETCDFLIB
       	 write (*,'(a)') "VIS_FIELDS  | Writing field 3 to netcdf"
         if (netcdf) call ncnbody_putfield( ncid, simtime,3, npx, npy, npz, field3, incdf )
#endif
      endif

      if (fselect4>0) then
#ifdef VISIT_NBODY
       	 write (*,'(a)') "VIS_FIELDS  | Shipping field 4"
         call flvisit_nbody2_field4_send(field4,npx,npy,npz)
#endif
#ifdef NETCDFLIB
       	 write (*,'(a)') "VIS_FIELDS  | Writing field 4 to netcdf"
         if (netcdf) call ncnbody_putfield( ncid, simtime, 4, npx, npy, npz, field4, incdf )
#endif
      endif

  endif


! TODO : Put this in dump_fields
     return
     call create_directory("fields")
     cfile = "fields/tslice."//cdump
     open (62,file=cfile)

     ! x-slices along laser axis
     jfoc = focus(2)/dy
     kfoc = focus(3)/dz

     ! temperature average line-out along laser axis: nave*nave average, converted to n/nc

     te_slice = 0.


     if (ngz<=5) then
        nave=0
     else
        nave = ngz/3
     endif

     norm = (2*nave+1)**2  ! convert Temp to keV

     do k=kfoc-nave,kfoc+nave
        do j=jfoc-nave,jfoc+nave

           te_slice(1:ngx) = te_slice(1:ngx)+telec(1:ngx,j,k)/norm  ! slice along laser axis: 5x5 average

        end do
     end do
     write(62,'((2(1pe12.4)))') &
          (i*dx,te_slice(i),i=1,ngx)
     close(62)


end subroutine vis_fields_nbody

