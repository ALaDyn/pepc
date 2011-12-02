module diagnostics
contains


subroutine error_norms(itime)

   use physvars
   implicit none
   include 'mpif.h'

   real,intent(in) :: itime

   real*8,parameter :: U_norm = 0.1504581298
   integer :: i, ierr
   real*8 :: r1, u_err_part, u_err_all, u_err_rel_all, u_err_rel_part, u_err_part_max, u_err_all_max, u_max_part, u_max_all, &
             w_err_part, w_err_all, w_err_rel_all, w_err_rel_part, w_err_part_max, w_err_all_max, w_max_part, w_max_all
   real*8 :: ux_ori(1:np), uy_ori(1:np), uz_ori(1:np), wx_ori(1:np), wy_ori(1:np), wz_ori(1:np)

   u_err_part = 0.
   u_err_all = 0.
   u_err_part_max = 0.
   u_err_all_max = 0.
   u_max_part = 0.
   u_max_all = 0.
   u_err_rel_part = 0.
   u_err_rel_all = 0.

   w_err_part = 0.
   w_err_all = 0.
   w_err_part_max = 0.
   w_err_all_max = 0.
   w_max_part = 0.
   w_max_all = 0.
   w_err_rel_part = 0.
   w_err_rel_all = 0.

   do i = 1,np

      if (ispecial == 10) then
         r1 = sqrt(dot_product(vortex_particles(i)%x,vortex_particles(i)%x))
         if (r1 > 0) then
            if (r1 .lt. 1) then
               ux_ori(i) = -1.0/(8.0*r1)*(1.0-(1.0-r1)**4)*vortex_particles(i)%x(2)
               uy_ori(i) = 1.0/(8.0*r1)*(1.0-(1.0-r1)**4)*vortex_particles(i)%x(1)
               uz_ori(i) = 0.
            else
               ux_ori(i) = -1.0/(8.0*r1)*vortex_particles(i)%x(2)
               uy_ori(i) = 1.0/(8.0*r1)*vortex_particles(i)%x(1)
               uz_ori(i) = 0.
            end if
            u_err_part_max = max(u_err_part_max,sqrt((ux_ori(i)-vortex_particles(i)%results%u(1))**2 + (uy_ori(i)-vortex_particles(i)%results%u(2))**2))
            u_max_part = max(u_max_part,sqrt(ux_ori(i)**2 + uy_ori(i)**2))
            u_err_rel_part = u_err_rel_part + ux_ori(i)**2 + uy_ori(i)**2
            u_err_part = u_err_part + ((ux_ori(i)-vortex_particles(i)%results%u(1))**2 + (uy_ori(i)-vortex_particles(i)%results%u(2))**2)
         end if
      end if


      if (ispecial == 11) then
         r1 = sqrt(dot_product(vortex_particles(i)%x,vortex_particles(i)%x))
         if (r1 > 0) then
            ux_ori(i) = -(0.1D01)/((0.24D02)*(r1**2))*(0.1D01-&
                 exp(-0.12D02*(r1**2)))*vortex_particles(i)%x(2)
            uy_ori(i) = (0.1D01)/((0.24D02)*(r1**2))*(0.1D01-&
                 exp(-0.12D02*(r1**2)))*vortex_particles(i)%x(1)
            uz_ori(i) = 0.
            u_err_part_max = max(u_err_part_max,sqrt((ux_ori(i)-vortex_particles(i)%results%u(1))**2 + (uy_ori(i)-vortex_particles(i)%results%u(2))**2))
            u_max_part = max(u_max_part,sqrt(ux_ori(i)**2 + uy_ori(i)**2))
            u_err_rel_part = u_err_rel_part + ux_ori(i)**2 + uy_ori(i)**2
            u_err_part = u_err_part + ((ux_ori(i)-vortex_particles(i)%results%u(1))**2 + (uy_ori(i)-vortex_particles(i)%results%u(2))**2)
         end if
      end if

      if (ispecial == 15) then
        ! TODO: use alpha notation here
         r1 = sqrt(dot_product(vortex_particles(i)%x,vortex_particles(i)%x))
         wx_ori(i) = 0.
         wy_ori(i) = 0.
         wz_ori(i) = 0.1D01/(0.4D01*pi*nu*itime)*exp(-r1**2/(0.4D01*nu*itime))
         !vortex_particles(i)%az = wz_ori(i)
         !wz_ori(i) = 0.1D01/(0.4D01*pi*nu*(itime-(vortex_particles(i)%sigma**2)/(2*nu)))*exp(-r1**2/(0.4D01*nu*(itime-(vortex_particles(i)%sigma**2)/(2*nu))))
         w_err_part_max = max(w_err_part_max,sqrt((wx_ori(i)-vortex_particles(i)%data%alpha(1))**2 + (wy_ori(i)-vortex_particles(i)%data%alpha(2))**2 + (wz_ori(i)-vortex_particles(i)%data%alpha(3))**2))
         w_max_part = max(w_max_part,sqrt(wx_ori(i)**2 + wy_ori(i)**2 + wz_ori(i)**2))
         w_err_rel_part = w_err_rel_part + wx_ori(i)**2 + wy_ori(i)**2 + wz_ori(i)**2
         !w_err_rel_part = w_err_rel_part + vortex_particles(i)%ax**2 + vortex_particles(i)%ay**2 + vortex_particles(i)%az**2
         w_err_part = w_err_part + ((wx_ori(i)-vortex_particles(i)%data%alpha(1))**2 + (wy_ori(i)-vortex_particles(i)%data%alpha(2))**2 + (wz_ori(i)-vortex_particles(i)%data%alpha(3))**2)
      end if
   end do

   if (ispecial == 11) then
      call MPI_ALLREDUCE(u_max_part,u_max_all,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(u_err_part_max,u_err_all_max,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(u_err_part,u_err_all,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(u_err_rel_part,u_err_rel_all,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

      u_err_all_max = u_err_all_max/u_max_all
      u_err_rel_all = sqrt(u_err_all)/sqrt(u_err_rel_all)
      u_err_all = sqrt(u_err_all)*h

      if (my_rank == 0) then
         open(666,file='error_u_isp11.dat',STATUS='UNKNOWN', POSITION = 'APPEND')
         write(666,*) itime, u_err_all, u_err_all/U_norm, u_err_rel_all, u_err_all_max, u_err_all*n*h
         close(666)
      end if
   end if

   if (ispecial == 10) then
      call MPI_ALLREDUCE(u_max_part,u_max_all,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(u_err_part_max,u_err_all_max,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(u_err_part,u_err_all,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(u_err_rel_part,u_err_rel_all,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

      u_err_all_max = u_err_all_max/u_max_all
      u_err_rel_all = sqrt(u_err_all)/sqrt(u_err_rel_all)
      u_err_all = sqrt(u_err_all/n)

      if (my_rank == 0) then
         open(666,file='error_u_isp10.dat',STATUS='UNKNOWN', POSITION = 'APPEND')
         write(666,*) itime, u_err_all, u_err_all/U_norm, u_err_rel_all, u_err_all_max, u_err_all*n*h
         close(666)
      end if
   end if

   if (ispecial == 15) then
   ! TODO: use alpha notation here
      call MPI_ALLREDUCE(w_max_part,w_max_all,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(w_err_part_max,w_err_all_max,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(w_err_part,w_err_all,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(w_err_rel_part,w_err_rel_all,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

      w_err_all_max = w_err_all_max/w_max_all
      w_err_rel_all = sqrt(w_err_all)/sqrt(w_err_rel_all)
      w_err_all = sqrt(w_err_all)/n

      if (my_rank == 0) then
         open(666,file='error_w_isp15.dat',STATUS='UNKNOWN', POSITION = 'APPEND')
         write(666,*) itime, w_err_all, w_err_all*n, w_err_rel_all, w_err_all_max, w_err_all*n*h
         close(666)
      end if
   end if

end subroutine error_norms


subroutine linear_diagnostics(itime,trun)

   use physvars
   implicit none
   include 'mpif.h'

   integer, intent(in) :: itime
   real, intent(in) :: trun

   real*8 :: omega(3), sendbuf_O(3), linear(3), sendbuf_I(3), angular(3), sendbuf_A(3)
   integer :: ierr, i

   ! local total vorticity
   sendbuf_O = 0
   do i = 1,np
      sendbuf_O(1:3) = sendbuf_O(1:3) + vortex_particles(i)%data%alpha(1:3)
   end do

   ! global total vortiticy
   omega = 0.
   call MPI_ALLREDUCE(sendbuf_O,omega,3,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

   ! local linear impulse
   sendbuf_I = 0
   do i = 1,np
      sendbuf_I(1) = sendbuf_I(1) + 0.5*vortex_particles(i)%x(2)*vortex_particles(i)%data%alpha(3)-vortex_particles(i)%x(3)*vortex_particles(i)%data%alpha(2)
      sendbuf_I(2) = sendbuf_I(2) + 0.5*vortex_particles(i)%x(3)*vortex_particles(i)%data%alpha(1)-vortex_particles(i)%x(1)*vortex_particles(i)%data%alpha(3)
      sendbuf_I(3) = sendbuf_I(3) + 0.5*vortex_particles(i)%x(1)*vortex_particles(i)%data%alpha(2)-vortex_particles(i)%x(2)*vortex_particles(i)%data%alpha(1)
   end do

   ! global linear impulse
   linear = 0.
   call MPI_ALLREDUCE(sendbuf_I,linear,3,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

   ! local angular impulse
   sendbuf_A = 0
   do i = 1,np
      sendbuf_A(1) = sendbuf_A(1) + 0.5*(vortex_particles(i)%x(2)*(vortex_particles(i)%x(1)*vortex_particles(i)%data%alpha(2)-vortex_particles(i)%x(2)*vortex_particles(i)%data%alpha(1))-&
           vortex_particles(i)%x(3)*(vortex_particles(i)%x(3)*vortex_particles(i)%data%alpha(1)-vortex_particles(i)%x(1)*vortex_particles(i)%data%alpha(3))-eps**2*omega(1))
      sendbuf_A(2) = sendbuf_A(2) + 0.5*(vortex_particles(i)%x(3)*(vortex_particles(i)%x(2)*vortex_particles(i)%data%alpha(3)-vortex_particles(i)%x(3)*vortex_particles(i)%data%alpha(2))-&
           vortex_particles(i)%x(1)*(vortex_particles(i)%x(1)*vortex_particles(i)%data%alpha(2)-vortex_particles(i)%x(2)*vortex_particles(i)%data%alpha(1))-eps**2*omega(2))
      sendbuf_A(3) = sendbuf_A(3) + 0.5*(vortex_particles(i)%x(1)*(vortex_particles(i)%x(3)*vortex_particles(i)%data%alpha(1)-vortex_particles(i)%x(1)*vortex_particles(i)%data%alpha(3))-&
           vortex_particles(i)%x(2)*(vortex_particles(i)%x(2)*vortex_particles(i)%data%alpha(3)-vortex_particles(i)%x(3)*vortex_particles(i)%data%alpha(2))-eps**2*omega(3))
   end do

   ! global angular impulse
   angular = 0.
   call MPI_ALLREDUCE(sendbuf_A,angular,3,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

   ! std and file output
   if (my_rank == 0) then
      write(*,*) '=============Linear diagnostics============='
      write(*,*) 'Omega:           ', sqrt(omega(1)**2+omega(2)**2+omega(3)**2), omega(1)+omega(2)+omega(3), omega(1), omega(2), omega(3)
      write(*,*) 'Linear Impulse:  ', sqrt(linear(1)**2+linear(2)**2+linear(3)**2), linear(1)+linear(2)+linear(3), linear(1), linear(2), linear(3)
      write(*,*) 'Angular Impulse: ', sqrt(angular(1)**2+angular(2)**2+angular(3)**2), angular(1)+angular(2)+angular(3), angular(1), angular(2), angular(3)
      write(*,*) '============================================'
      write(66,*) itime, trun, sqrt(omega(1)**2+omega(2)**2+omega(3)**2),' ', omega(1)+omega(2)+omega(3),' ', omega(1),' ', omega(2),' ', omega(3) ,' ', &
                               sqrt(linear(1)**2+linear(2)**2+linear(3)**2),' ', linear(1)+linear(2)+linear(3),' ', linear(1),' ', linear(2),' ', linear(3),' ', &
                               sqrt(angular(1)**2+angular(2)**2+angular(3)**2),' ', angular(1)+angular(2)+angular(3),' ', angular(1),' ', angular(2),' ', angular(3),' '
   end if


end subroutine linear_diagnostics

end module diagnostics
