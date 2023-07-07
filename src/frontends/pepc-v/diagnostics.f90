! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
!
! Copyright (C) 2002-2023 Juelich Supercomputing Centre,
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

module diagnostics
   implicit none
contains

   subroutine error_norms(itime)

      use physvars
      use mpi
      implicit none

      real, intent(in) :: itime

      real(kind_physics), parameter :: U_norm = 0.1504581298
      integer :: ierr
      integer(kind_particle) :: i
      real(kind_physics) :: r1, u_err_part, u_err_all, u_err_rel_all, u_err_rel_part, u_err_part_max, &
                            u_err_all_max, u_max_part, u_max_all, w_err_part, w_err_all, w_err_rel_all, w_err_rel_part, &
                            w_err_part_max, w_err_all_max, w_max_part, w_max_all
      real(kind_physics), dimension(1:np) :: ux_ori, uy_ori, uz_ori, wx_ori, wy_ori, wz_ori

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

      do i = 1, np

         if (ispecial .eq. 10) then
            r1 = sqrt(dot_product(vortex_particles(i)%x, vortex_particles(i)%x))
            if (r1 .gt. 0) then
               if (r1 .lt. 1) then
                  ux_ori(i) = -1.0 / (8.0 * r1) * (1.0 - (1.0 - r1)**4) * vortex_particles(i)%x(2)
                  uy_ori(i) = 1.0 / (8.0 * r1) * (1.0 - (1.0 - r1)**4) * vortex_particles(i)%x(1)
                  uz_ori(i) = 0.
               else
                  ux_ori(i) = -1.0 / (8.0 * r1) * vortex_particles(i)%x(2)
                  uy_ori(i) = 1.0 / (8.0 * r1) * vortex_particles(i)%x(1)
                  uz_ori(i) = 0.
               end if
               u_err_part_max = max(u_err_part_max, sqrt((ux_ori(i) - vortex_particles(i)%results%u(1))**2 + (uy_ori(i) - vortex_particles(i)%results%u(2))**2))
               u_max_part = max(u_max_part, sqrt(ux_ori(i)**2 + uy_ori(i)**2))
               u_err_rel_part = u_err_rel_part + ux_ori(i)**2 + uy_ori(i)**2
               u_err_part = u_err_part + ((ux_ori(i) - vortex_particles(i)%results%u(1))**2 + (uy_ori(i) - vortex_particles(i)%results%u(2))**2)
            end if
         end if

         if (ispecial .eq. 11) then
            r1 = sqrt(dot_product(vortex_particles(i)%x, vortex_particles(i)%x))
            if (r1 .gt. 0) then
               ux_ori(i) = -(0.1D01) / ((0.24D02) * (r1**2)) * (0.1D01 - &
                                                                exp(-0.12D02 * (r1**2))) * vortex_particles(i)%x(2)
               uy_ori(i) = (0.1D01) / ((0.24D02) * (r1**2)) * (0.1D01 - &
                                                               exp(-0.12D02 * (r1**2))) * vortex_particles(i)%x(1)
               uz_ori(i) = 0.
               u_err_part_max = max(u_err_part_max, sqrt((ux_ori(i) - vortex_particles(i)%results%u(1))**2 + (uy_ori(i) - vortex_particles(i)%results%u(2))**2))
               u_max_part = max(u_max_part, sqrt(ux_ori(i)**2 + uy_ori(i)**2))
               u_err_rel_part = u_err_rel_part + ux_ori(i)**2 + uy_ori(i)**2
               u_err_part = u_err_part + ((ux_ori(i) - vortex_particles(i)%results%u(1))**2 + (uy_ori(i) - vortex_particles(i)%results%u(2))**2)
            end if
         end if

         if (ispecial .eq. 15) then
            ! TODO: use alpha notation here
            r1 = sqrt(dot_product(vortex_particles(i)%x, vortex_particles(i)%x))
            wx_ori(i) = 0.
            wy_ori(i) = 0.
            wz_ori(i) = 0.1D01 / (0.4D01 * pi * nu * itime) * exp(-r1**2 / (0.4D01 * nu * itime))
            !vortex_particles(i)%az = wz_ori(i)
            !wz_ori(i) = 0.1D01/(0.4D01*pi*nu*(itime-(vortex_particles(i)%sigma**2)/(2*nu)))*exp(-r1**2/(0.4D01*nu*(itime-(vortex_particles(i)%sigma**2)/(2*nu))))
            w_err_part_max = max(w_err_part_max, sqrt((wx_ori(i) - vortex_particles(i)%data%alpha(1))**2 + (wy_ori(i) - vortex_particles(i)%data%alpha(2))**2 + (wz_ori(i) - vortex_particles(i)%data%alpha(3))**2))
            w_max_part = max(w_max_part, sqrt(wx_ori(i)**2 + wy_ori(i)**2 + wz_ori(i)**2))
            w_err_rel_part = w_err_rel_part + wx_ori(i)**2 + wy_ori(i)**2 + wz_ori(i)**2
            !w_err_rel_part = w_err_rel_part + vortex_particles(i)%ax**2 + vortex_particles(i)%ay**2 + vortex_particles(i)%az**2
            w_err_part = w_err_part + ((wx_ori(i) - vortex_particles(i)%data%alpha(1))**2 + (wy_ori(i) - vortex_particles(i)%data%alpha(2))**2 + (wz_ori(i) - vortex_particles(i)%data%alpha(3))**2)
         end if
      end do

      if (ispecial .eq. 11) then
         call MPI_ALLREDUCE(u_max_part, u_max_all, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
         call MPI_ALLREDUCE(u_err_part_max, u_err_all_max, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
         call MPI_ALLREDUCE(u_err_part, u_err_all, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
         call MPI_ALLREDUCE(u_err_rel_part, u_err_rel_all, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

         u_err_all_max = u_err_all_max / u_max_all
         u_err_rel_all = sqrt(u_err_all) / sqrt(u_err_rel_all)
         u_err_all = sqrt(u_err_all) * h

         if (my_rank .eq. 0) then
            open (666, file='error_u_isp11.dat', STATUS='UNKNOWN', POSITION='APPEND')
            write (666, *) itime, u_err_all, u_err_all / U_norm, u_err_rel_all, u_err_all_max, u_err_all * n * h
            close (666)
         end if
      end if

      if (ispecial .eq. 10) then
         call MPI_ALLREDUCE(u_max_part, u_max_all, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
         call MPI_ALLREDUCE(u_err_part_max, u_err_all_max, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
         call MPI_ALLREDUCE(u_err_part, u_err_all, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
         call MPI_ALLREDUCE(u_err_rel_part, u_err_rel_all, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

         u_err_all_max = u_err_all_max / u_max_all
         u_err_rel_all = sqrt(u_err_all) / sqrt(u_err_rel_all)
         u_err_all = sqrt(u_err_all / n)

         if (my_rank .eq. 0) then
            open (666, file='error_u_isp10.dat', STATUS='UNKNOWN', POSITION='APPEND')
            write (666, *) itime, u_err_all, u_err_all / U_norm, u_err_rel_all, u_err_all_max, u_err_all * n * h
            close (666)
         end if
      end if

      if (ispecial .eq. 15) then
         ! TODO: use alpha notation here
         call MPI_ALLREDUCE(w_max_part, w_max_all, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
         call MPI_ALLREDUCE(w_err_part_max, w_err_all_max, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
         call MPI_ALLREDUCE(w_err_part, w_err_all, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
         call MPI_ALLREDUCE(w_err_rel_part, w_err_rel_all, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

         w_err_all_max = w_err_all_max / w_max_all
         w_err_rel_all = sqrt(w_err_all) / sqrt(w_err_rel_all)
         w_err_all = sqrt(w_err_all) / n

         if (my_rank .eq. 0) then
            open (666, file='error_w_isp15.dat', STATUS='UNKNOWN', POSITION='APPEND')
            write (666, *) itime, w_err_all, w_err_all * n, w_err_rel_all, w_err_all_max, w_err_all * n * h
            close (666)
         end if
      end if

   end subroutine error_norms

   subroutine linear_diagnostics(itime, trun)

      use physvars
      use files, only: diag_unit
      use mpi
      implicit none

      integer, intent(in) :: itime
      real, intent(in) :: trun

      real(kind_physics)     :: omega(3), sendbuf_O(3), linear(3), sendbuf_I(3), angular(3), sendbuf_A(3)
      integer                :: ierr
      integer(kind_particle) :: i

      ! local total vorticity
      sendbuf_O = 0
      do i = 1, np
         sendbuf_O(1:3) = sendbuf_O(1:3) + vortex_particles(i)%data%alpha(1:3)
      end do

      ! global total vortiticy
      omega = 0.
      call MPI_ALLREDUCE(sendbuf_O, omega, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

      ! local linear impulse
      sendbuf_I = 0
      do i = 1, np
         sendbuf_I(1) = sendbuf_I(1) + 0.5 * vortex_particles(i)%x(2) * vortex_particles(i)%data%alpha(3) - vortex_particles(i)%x(3) * vortex_particles(i)%data%alpha(2)
         sendbuf_I(2) = sendbuf_I(2) + 0.5 * vortex_particles(i)%x(3) * vortex_particles(i)%data%alpha(1) - vortex_particles(i)%x(1) * vortex_particles(i)%data%alpha(3)
         sendbuf_I(3) = sendbuf_I(3) + 0.5 * vortex_particles(i)%x(1) * vortex_particles(i)%data%alpha(2) - vortex_particles(i)%x(2) * vortex_particles(i)%data%alpha(1)
      end do

      ! global linear impulse
      linear = 0.
      call MPI_ALLREDUCE(sendbuf_I, linear, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

      ! local angular impulse
!~    sendbuf_A = 0
!~    do i = 1,np
!~       sendbuf_A(1) = sendbuf_A(1) + 0.5*(vortex_particles(i)%x(2)*(vortex_particles(i)%x(1)*vortex_particles(i)%data%alpha(2)-vortex_particles(i)%x(2)*vortex_particles(i)%data%alpha(1))-&
!~            vortex_particles(i)%x(3)*(vortex_particles(i)%x(3)*vortex_particles(i)%data%alpha(1)-vortex_particles(i)%x(1)*vortex_particles(i)%data%alpha(3))-eps**2*omega(1))
!~       sendbuf_A(2) = sendbuf_A(2) + 0.5*(vortex_particles(i)%x(3)*(vortex_particles(i)%x(2)*vortex_particles(i)%data%alpha(3)-vortex_particles(i)%x(3)*vortex_particles(i)%data%alpha(2))-&
!~            vortex_particles(i)%x(1)*(vortex_particles(i)%x(1)*vortex_particles(i)%data%alpha(2)-vortex_particles(i)%x(2)*vortex_particles(i)%data%alpha(1))-eps**2*omega(2))
!~       sendbuf_A(3) = sendbuf_A(3) + 0.5*(vortex_particles(i)%x(1)*(vortex_particles(i)%x(3)*vortex_particles(i)%data%alpha(1)-vortex_particles(i)%x(1)*vortex_particles(i)%data%alpha(3))-&
!~            vortex_particles(i)%x(2)*(vortex_particles(i)%x(2)*vortex_particles(i)%data%alpha(3)-vortex_particles(i)%x(3)*vortex_particles(i)%data%alpha(2))-eps**2*omega(3))
!~    end do

      ! global angular impulse
      angular = 0.
!~      call MPI_ALLREDUCE(sendbuf_A, angular, 3, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

      ! std and file output
      if (my_rank .eq. 0) then
         write (*, *) '=============Linear diagnostics============='
         write (*, *) 'Omega:           ', sqrt(omega(1)**2 + omega(2)**2 + omega(3)**2), omega(1) + omega(2) + omega(3), omega(1), omega(2), omega(3)
         write (*, *) 'Linear Impulse:  ', sqrt(linear(1)**2 + linear(2)**2 + linear(3)**2), linear(1) + linear(2) + linear(3), linear(1), linear(2), linear(3)
         write (*, *) 'Angular Impulse: ', sqrt(angular(1)**2 + angular(2)**2 + angular(3)**2), angular(1) + angular(2) + angular(3), angular(1), angular(2), angular(3)
         write (*, *) '============================================'
         write (diag_unit, *) itime, trun, sqrt(omega(1)**2 + omega(2)**2 + omega(3)**2), ' ', omega(1) + omega(2) + omega(3), ' ', omega(1), ' ', omega(2), ' ', omega(3), ' ', &
            sqrt(linear(1)**2 + linear(2)**2 + linear(3)**2), ' ', linear(1) + linear(2) + linear(3), ' ', linear(1), ' ', linear(2), ' ', linear(3), ' ', &
            sqrt(angular(1)**2 + angular(2)**2 + angular(3)**2), ' ', angular(1) + angular(2) + angular(3), ' ', angular(1), ' ', angular(2), ' ', angular(3), ' '
      end if

   end subroutine linear_diagnostics

   subroutine divergence_diag()

      use physvars
      use mpi
      implicit none

      integer :: ierr
      integer(kind_particle) :: i
      real(kind_physics) :: div_max_local, div_min_local, div_mean_local, div_max, div_mean, div_min

      div_mean_local = 0.

      div_max_local = maxval(vortex_particles(1:np)%results%div)
      div_min_local = minval(vortex_particles(1:np)%results%div)
      do i = 1, np
         div_mean_local = div_mean_local + vortex_particles(i)%results%div**2
      end do

      call MPI_ALLREDUCE(div_max_local, div_max, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(div_min_local, div_min, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(div_mean_local, div_mean, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
      div_mean = sqrt(div_mean / n)

      if (my_rank .eq. 0) write (*, *) 'Divergence (max/min/mean/denorm):', div_max, div_min, div_mean, div_mean * n

   end subroutine divergence_diag

   subroutine verify_direct()

      use physvars
      use module_interaction_specific_types, only: t_particle_results
      use manipulate_particles
      use mpi
      implicit none

      real(kind_physics), dimension(np) :: rel_error_u, rel_error_af

      integer :: ierr
      integer(kind_particle) :: i
      type(t_particle_results) :: direct_results(1:np), pepc_results(1:np)
      real(kind_physics) :: diff_u_mean_local, diff_u_mean, diff_af_mean_local, diff_af_mean, &
                            diff_u_max_local, diff_af_max_local, diff_u_max, diff_af_max, t1, &
                            L2_u_mean_local, L2_u_mean, L2_af_mean_local, L2_af_mean

      rel_error_u = 0.
      rel_error_af = 0.

      if (my_rank .eq. 0) write (*, *) 'Starting direct sum ...'
      t1 = MPI_WTIME()
      call direct_sum(np, vortex_particles, direct_results, my_rank, n_cpu)
      if (my_rank .eq. 0) write (*, *) '                    ... done in', MPI_WTIME() - t1, 'sec.'

      if (my_rank .eq. 0) write (*, *) 'Starting post-processing ...'
      t1 = MPI_WTIME()

      do i = 1, np
         pepc_results(i)%u(1:3) = vortex_particles(i)%results%u(1:3)
         pepc_results(i)%af(1:3) = vortex_particles(i)%results%af(1:3)
         direct_results(i)%u(1:3) = direct_results(i)%u(1:3) * force_const
         direct_results(i)%af(1:3) = direct_results(i)%af(1:3) * force_const
      end do

      diff_u_max_local = 0.
      diff_u_mean_local = 0.
      diff_af_max_local = 0.
      diff_af_mean_local = 0.
      L2_u_mean_local = 0.
      L2_af_mean_local = 0.

      do i = 1, np
         rel_error_u(i) = sqrt((direct_results(i)%u(1) - pepc_results(i)%u(1))**2 + &
                               (direct_results(i)%u(2) - pepc_results(i)%u(2))**2 + &
                               (direct_results(i)%u(3) - pepc_results(i)%u(3))**2) / sqrt(dot_product(direct_results(i)%u, direct_results(i)%u))
         rel_error_af(i) = sqrt((direct_results(i)%af(1) - pepc_results(i)%af(1))**2 + &
                                (direct_results(i)%af(2) - pepc_results(i)%af(2))**2 + &
                                (direct_results(i)%af(3) - pepc_results(i)%af(3))**2) / sqrt(dot_product(direct_results(i)%af, direct_results(i)%af))

         diff_u_max_local = max(diff_u_max_local, rel_error_u(i))
         diff_u_mean_local = diff_u_mean_local + rel_error_u(i)
         L2_u_mean_local = L2_u_mean_local + ((direct_results(i)%u(1) - pepc_results(i)%u(1))**2 + &
                                              (direct_results(i)%u(2) - pepc_results(i)%u(2))**2 + &
                                              (direct_results(i)%u(3) - pepc_results(i)%u(3))**2) / dot_product(direct_results(i)%u, direct_results(i)%u)

         diff_af_max_local = max(diff_af_max_local, rel_error_af(i))
         diff_af_mean_local = diff_af_mean_local + rel_error_af(i)
         L2_af_mean_local = L2_af_mean_local + ((direct_results(i)%af(1) - pepc_results(i)%af(1))**2 + &
                                                (direct_results(i)%af(2) - pepc_results(i)%af(2))**2 + &
                                                (direct_results(i)%af(3) - pepc_results(i)%af(3))**2) / dot_product(direct_results(i)%af, direct_results(i)%af)

      end do

      call MPI_ALLREDUCE(diff_u_max_local, diff_u_max, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(diff_af_max_local, diff_af_max, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)

      call MPI_ALLREDUCE(diff_u_mean_local, diff_u_mean, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(diff_af_mean_local, diff_af_mean, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(L2_u_mean_local, L2_u_mean, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(L2_af_mean_local, L2_af_mean, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

      diff_u_mean = diff_u_mean / n
      diff_af_mean = diff_af_mean / n
      L2_u_mean = sqrt(L2_u_mean) / n
      L2_af_mean = sqrt(L2_af_mean) / n

      if (my_rank .eq. 0) write (*, *) '                        ... done in', MPI_WTIME() - t1, 'sec.'

      if (my_rank .eq. 0) then
         write (*, *) 'Error in u (rel. L2 mean / rel. mean / rel. max) ', L2_u_mean, diff_u_mean, diff_u_max
         write (*, *) 'Error in af (rel. L2 mean / rel. mean / rel. max)', L2_af_mean, diff_af_mean, diff_af_max
      end if

      !call MPI_ABORT(MPI_COMM_WORLD,1,ierr)

   end subroutine verify_direct

end module diagnostics
