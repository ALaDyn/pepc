! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
!
! Copyright (C) 2002-2017 Juelich Supercomputing Centre,
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>  Encapsulates external particle beam
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_particle_beam

   implicit none
   save
   private

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   public beam_control     !< Particle beam initialisation (interactive)
   public beam_dust        !< Sets up spherical dust particle
   public beam             !< Sets up particle beam on CPU 0

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

   ! ==============================================
   !
   !                BEAM_CONTROL
   !
   !  Interactive control of particle source
   !
   ! ==============================================

   subroutine beam_control

      use module_physvars
      use module_particle_props
      use module_utilities
      use mpi
      implicit none

      integer :: i, p, iseed1, iseed2, ierr
      real :: Volb, dpx, yt, zt, vosc_old, sigma_old, tpulse_old, u_old, theta_old, phi_old
      integer :: lvisit_active = 0
      real :: ct, st, cp, sp, vx_beam, vy_beam, vz_beam, th_beam, xb, yb, zb
      !  logical :: beam_on = .true.
      logical :: beam_debug = .true.

      integer, save :: np_beam_dt  ! current # beam particles

      ! First check for VISIT connection

#ifdef VISIT_NBODY
      !  if (my_rank==0)   call flvisit_spk_check_connection(lvisit_active)
      if (my_rank .eq. 0) call flvisit_nbody2_check_connection(lvisit_active)
      call MPI_BCAST(lvisit_active, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#endif

      if (lvisit_active .eq. 0) then
         if (my_rank .eq. 0) write (*, *) ' No Connection to Visualization'
         return
      end if

      if (itime .eq. 0) then
         np_beam_dt = np_beam  ! initial # beam particles to introduce per timestep
      end if

      !  if (.not. beam_on) return

      !  if (npart + np_beam_dt > npartm .or. max_list_length > nintmax-20) then
      !     if (my_rank==0) write(*,*) 'Array limit reached: switching off beam'
      !     beam_on = .false.
      !     return
      !  endif

      if (beam_config .eq. 5) then
         ! Define beam from laser parameters
         tpulse_old = tpulse
         vosc_old = vosc
         sigma_old = sigma
         u_beam = vosc
         rho_beam = sigma
         r_beam = tpulse

      else if (beam_config .ge. 3 .and. beam_config .le. 4) then
         ! Define beam from laser parameters
         tpulse_old = tpulse
         vosc_old = vosc
         sigma_old = sigma
         u_beam = vosc
         rho_beam = sigma
         th_beam = theta_beam  ! incidence angle instead of pulse duration
         theta_old = th_beam

      else if (scheme .eq. 4 .and. .not. (beam_config .le. 3 .and. beam_config .gt. 0)) then
         ! Temperature clamp mode - laser should be off
         u_beam = Te_kev

      else if (scheme .eq. 5) then
         ! ion crystal eqm  mode:
         !  r_beam is mean ion spacing
         !  u_beam is ion temperature (eV)
         !  rho_beam is potential constant
         r_beam = a_ii
         u_beam = Ti_keV
         rho_beam = log10(bond_const)

      else if (beam_config .eq. 8) then ! dust particle
         u_old = u_beam
         theta_old = theta_beam
         phi_old = phi_beam
      end if

#ifdef VISIT_NBODY
      if (itime .eq. 0 .and. my_rank .eq. 0) then
         !     call flvisit_spk_check_connection(lvisit_active)
         call flvisit_nbody2_check_connection(lvisit_active)
         ! Specify default parameters at beginning of run
         call flvisit_spk_beam_paraminit_send(th_beam, phi_beam, r_beam, rho_beam, u_beam)
      end if

      if (my_rank .eq. 0) then
         !     call flvisit_spk_check_connection(lvisit_active)
         call flvisit_nbody2_check_connection(lvisit_active)

         ! Fetch real-time, user-specified control parameters
         if (lvisit_active .ne. 0) then
            !  TODO:  Need XNBODY equivalent here
            !        call flvisit_spk_beam_param_recv( th_beam,phi_beam,r_beam,rho_beam,u_beam)
         else
            write (*, *) ' No Connection to Visualization'
            return
         end if

      end if

#endif

      ! Broadcast beam parameters to all other PEs
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)   ! Synchronize first
      call MPI_BCAST(lvisit_active, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      if (lvisit_active .ne. 0) then
         call MPI_BCAST(th_beam, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
         call MPI_BCAST(phi_beam, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
         call MPI_BCAST(r_beam, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
         call MPI_BCAST(rho_beam, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
         call MPI_BCAST(u_beam, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      else
         if (my_rank .eq. 0) write (*, *) ' No Connection to Visualization'
         return
      end if

      !  if (rho_beam==0 )then
      !     if (my_rank==0) write(*,*) ' Switching off beam'
      !
      !     return
      !  endif

      if (beam_config .ge. 3 .and. beam_config .le. 5) then
         ! laser standing wave or pond bullet

         !     u_beam = max(abs(u_beam),0.1)
         !     rho_beam = max(abs(rho_beam),0.5)
         !     r_beam = max(abs(r_beam),0.1)

         if (u_beam / vosc .lt. 10.0 .and. u_beam / vosc .gt. 0.1) then
            vosc = u_beam ! limit amplitude change
            if (my_rank .eq. 0 .and. vosc .ne. vosc_old) write (*, *) 'Laser amplitude changed'
         else
            if (my_rank .eq. 0) write (*, *) 'Amplitude change too big - readjust!'
         end if

         if (sigma / rho_beam .gt. 0.1 .and. sigma / rho_beam .lt. 10.0) then
            sigma = rho_beam
            if (my_rank .eq. 0 .and. sigma .ne. sigma_old) write (*, *) 'Laser spot size changed'
         else
            if (my_rank .eq. 0) write (*, *) 'Spot size change too big - readjust!'
         end if
         if (tpulse / r_beam .gt. 0.1 .and. tpulse / r_beam .lt. 10.0) then
            tpulse = r_beam
         else
            if (my_rank .eq. 0) write (*, *) 'Pulse length change too big - readjust!'
            if (my_rank .eq. 0 .and. tpulse .ne. tpulse_old) write (*, *) 'Laser pulse length changed'
         end if

         if (th_beam .ne. theta_beam) then
            theta_beam = th_beam
            if (my_rank .eq. 0 .and. theta_beam .ne. theta_old) write (*, *) 'Incidence angle changed'
         end if

      else if (beam_config .eq. 2) then
         nb_pe = np_beam_dt / n_cpu  ! # beam particles to load per PE

         if (my_rank .eq. 0 .and. beam_debug) then
            write (*, *) 'beam parts ', nb_pe, ' theta ', theta_beam, ' phi ', phi_beam
            write (*, *) 'r ', r_beam, ' rho ', rho_beam, ' u ', u_beam, ' dt ', np_beam_dt
         end if

         Volb = pi * r_beam**2 * x_beam       !  beam cylinder volumy_rank:  r_beam is radius
         qeb = Volb * rho_beam / np_beam    ! charge
         dpx = x_beam / nb_pe           ! x-axis spacing
         iseed2 = -131 - np_beam - 3 * my_rank ! seed
         iseed1 = -333 - np_beam - 3 * my_rank

         if (qeb .lt. 0) then
            mass_beam = 1.  ! electrons
         else
            mass_beam = 1836.  ! protons
         end if

         ! particle beam initialised along x-axis and rotated by theta, phi
         ct = cos(theta_beam)
         st = sin(theta_beam)
         cp = cos(phi_beam)
         sp = sin(phi_beam)
         vz_beam = -u_beam * st
         vx_beam = u_beam * ct * cp
         vy_beam = u_beam * ct * sp

         i = 0

         do while (i .lt. nb_pe)
            yt = r_beam * (2 * rano(iseed2) - 1.)
            zt = r_beam * (2 * rano(iseed1) - 1.)
            if (yt**2 + zt**2 .le. r_beam**2) then
               i = i + 1
               p = np_local + i  ! put them after last particle on this PE
               xb = dpx * i + dpx / n_cpu * my_rank
               yb = yt
               zb = zt

               ! Now rotate disc about Euler angles

               x(p) = xb * ct * cp - yb * ct * sp + zb * st
               y(p) = xb * sp + yb * cp
               z(p) = -xb * st * cp + yb * st * sp + zb * ct

               ! Add starting point
               x(p) = x(p) + start_beam
               y(p) = y(p) + yl / 2.
               z(p) = z(p) + zl / 2.

               q(p) = qeb
               m(p) = abs(qeb) * mass_beam
               ux(p) = vx_beam
               uy(p) = vy_beam
               uz(p) = vz_beam
               pepid(p) = my_rank                ! processor ID
               pelabel(p) = npart_total + my_rank * nb_pe + i  ! labels
               Ex(p) = 0.
               Ey(p) = 0.
               Ez(p) = 0.
               Bx(p) = 0.
               By(p) = 0.
               Bz(p) = 0.
               Ax(p) = 0.
               Ay(p) = 0.
               Az(p) = 0.
               pot(p) = 0.
               work(p) = 0.
            end if
         end do

         ! Augmy_ranknt total # particles - have to limit increase to max array size
         np_local = np_local + nb_pe
         np_beam = np_beam + np_beam_dt
         npart_total = npart_total + np_beam_dt

      else if (beam_config .eq. 8) then
         ! Dust particle - # beam particles constant; infinite mass

         if (my_rank .eq. 0 .and. u_beam .ne. u_old) write (*, *) 'Beam velocity changed'

         if (my_rank .eq. 0 .and. beam_debug) then
            write (*, *) ' theta ', theta_beam, ' phi ', phi_beam, ' u ', u_beam
         end if
         ! dust particle velocity rotated by theta, phi
         ! beam particles could be sitting anywhere now
         Volb = 4 * pi / 3.*r_beam**3
         qeb = Volb * rho_beam / np_beam    ! new charge
         ct = cos(theta_beam)
         st = sin(theta_beam)
         cp = cos(phi_beam)
         sp = sin(phi_beam)
         vy_beam = u_beam * st * cp * vte * 10   ! Scale by thermal velocity
         vx_beam = u_beam * st * sp * vte * 10
         vz_beam = u_beam * ct * vte * 10
         do i = 1, np_local
            if (pelabel(i) .gt. ne + ni) then
               ux(i) = vx_beam
               uy(i) = vy_beam
               uz(i) = vz_beam
               q(i) = qeb
            end if

         end do

      else if (scheme .eq. 5) then
         ! ion crystal eqm  mode:
         !  r_beam is mean ion spacing
         !  u_beam is ion temperature (eV)
         !  rho_beam is potential constant
         !     a_ii = r_beam
         !     Ti_kev = u_beam
         !     bond_const = 10**(rho_beam)
         if (my_rank .eq. 0) write (*, *) 'Steering pars: a_i=', a_ii, ' Ti=', Ti_kev, ' Pot strength=', bond_const
      else if (scheme .eq. 4) then
         ! Electron temp clamp
         Te_kev = u_beam
      end if

   end subroutine beam_control

   ! ==============================================
   !
   !                BEAM
   !
   !  Sets up particle beam on CPU 0
   !
   ! ==============================================

   subroutine beam

      use module_physvars
      use module_particle_props
      use module_utilities
      implicit none
      integer :: i, p, iseed1, iseed2
      real :: Volb, dpx, yt, zt
      real :: vx_beam, vy_beam, vz_beam, ct, cp, st, sp

      !  evaluate system constants  from inputs:
      !  beam cylinder volume:  r_beam is radius

      if (my_rank .eq. 0) then

         write (*, *) 'Setting up particle beam ', r_beam, ' x ', x_beam
         Volb = pi * r_beam**2 * x_beam
         qeb = Volb * rho_beam / np_beam    ! charge
         dpx = x_beam / np_beam           ! x-axis spacing
         iseed2 = -131 - np_beam - 3 * my_rank ! seed
         iseed1 = -333 - np_beam - 3 * my_rank

         ! proton beam: initialised along x-axis

         ! beam initialised along x-axis and rotated by theta, phi

         ct = cos(theta_beam)
         st = sin(theta_beam)
         cp = cos(phi_beam)
         sp = sin(phi_beam)
         vz_beam = u_beam * st
         vx_beam = u_beam * ct * cp
         vy_beam = u_beam * ct * sp

         i = 0

         do while (i .lt. np_beam)
            yt = r_beam * (2 * rano(iseed2) - 1.)
            zt = r_beam * (2 * rano(iseed1) - 1.)
            if (yt**2 + zt**2 .le. r_beam**2) then
               i = i + 1
               p = np_local + i  ! put them after plasma ions
               x(p) = start_beam + dpx * i
               y(p) = yt + yl / 1.9
               z(p) = zt + zl / 1.8
               q(p) = qeb
               m(p) = -qeb * mass_beam
               ux(p) = vx_beam
               uy(p) = vy_beam
               uz(p) = vz_beam
               pepid(p) = my_rank                ! processor ID
               pelabel(p) = npart_total + np_beam + i  ! labels
            end if
         end do
      end if

      np_local = np_local + np_beam

      npart_total = npart_total + np_beam  ! Augmy_ranknt particle numbers for all CPUs
      ! zero fields
      Ex(1:np_local) = 0
      Ey(1:np_local) = 0
      Ez(1:np_local) = 0
      Bx(1:np_local) = 0
      By(1:np_local) = 0
      Bz(1:np_local) = 0
      Ax(1:np_local) = 0
      Ay(1:np_local) = 0
      Az(1:np_local) = 0
      pot(1:np_local) = 0
      work(1:np_local) = 1.   ! set work load balanced initially

   end subroutine beam

   ! ==============================================
   !
   !                BEAM_DUST
   !
   !  Sets up spherical dust particle
   !
   ! ==============================================

   subroutine beam_dust

      use module_physvars
      use module_particle_props
      use module_utilities
      implicit none
      integer :: i, p, iseed1
      real :: Volb, xt, yt, zt
      real :: vx_beam, vy_beam, vz_beam, ct, cp, st, sp

      !  evaluate system constants  from inputs:
      !  beam cylinder volume:  r_beam is radius

      Volb = 4 * pi / 3.*r_beam**3
      qeb = Volb * rho_beam / np_beam    ! charge

      if (my_rank .eq. 0) then
         iseed1 = -333 - np_beam - 3 * my_rank

         ! proton beam: initialised along x-axis

         ! beam initialised along x-axis and rotated by theta, phi
         write (*, *) 'Setting up dust particle'
         write (*, '(a20,f12.5)') 'Beam charge ', qeb
         write (*, '(a20,f12.5)') 'Ion charge ', qi
         ct = cos(theta_beam)
         st = sin(theta_beam)
         cp = cos(phi_beam)
         sp = sin(phi_beam)
         vz_beam = u_beam * st
         vy_beam = u_beam * ct * cp
         vx_beam = u_beam * ct * sp

         i = 0

         ! Put dust particle on root PE
         do while (i .lt. np_beam)
            yt = r_beam * (2 * rano(iseed1) - 1.)
            zt = r_beam * (2 * rano(iseed1) - 1.)
            xt = r_beam * (2 * rano(iseed1) - 1.)
            if (xt**2 + yt**2 + zt**2 .le. r_beam**2) then
               i = i + 1
               p = np_local + i  ! put them after plasma ions
               x(p) = xl / 2.+xt        ! Assume disc or slab geometry and inject from side
               y(p) = yt + start_beam
               z(p) = zt + zl / 2.
               q(p) = qeb
               m(p) = -qeb * mass_beam
               ux(p) = vx_beam
               uy(p) = vy_beam
               uz(p) = vz_beam
               pepid(p) = my_rank                ! processor ID
               pelabel(p) = ne + ni + i  ! labels
            end if
         end do
         np_local = np_local + np_beam

      end if
      npart_total = npart_total + np_beam  ! Augment particle numbers

   end subroutine beam_dust

end module module_particle_beam
