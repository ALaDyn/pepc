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

!>
!> helper module
!>
module helper
   use module_pepc_kinds
   use module_pepc_types
   use module_timings
   use iso_fortran_env
   use particles_resize
   implicit none

   ! timing variables
   integer, parameter :: t_user_total = t_userdefined_first
   integer, parameter :: t_user_init = t_userdefined_first + 1
   integer, parameter :: t_user_step = t_userdefined_first + 2
   integer, parameter :: t_user_directsum = t_userdefined_first + 3
   integer, parameter :: t_user_particleio = t_userdefined_first + 4
   integer, parameter :: t_boris = t_userdefined_first + 5

   ! MPI variables
   integer(kind_pe) :: my_rank, n_ranks, ierr
   logical :: root

   ! time variables
   real*8 :: dt
   integer :: step

   ! control variables
   integer :: nt ! number of timesteps
   integer(kind_particle) :: tnp ! total number of particles
   integer(kind_particle) :: np ! local number of particles
   logical :: particle_output ! turn vtk output on/off
   logical :: domain_output ! turn vtk output on/off
   logical :: particle_test ! check tree code results against direct summation
   logical :: reflecting_walls ! reflect particles at walls
   integer :: diag_interval ! number of timesteps between all diagnostics and IO
   real(kind_physics) :: plasma_dimensions(3) ! size of the simulation box
   real(kind_physics) :: external_e(3)

   integer, parameter :: particle_direct = -1 ! number of particle for direct summation

   ! particle data (position, velocity, mass, charge)
   type(t_particle), allocatable   :: particles(:)
   real(kind_physics), allocatable :: direct_L2(:)

   ! buffer to record newly generated particles & related counters
   type(linked_list_elem), pointer :: buffer, particle_guide
   integer :: electron_num, i, new_particle_cnt, local_electron_num, swapped_num

   ! variables for random number generation
   integer :: dummy
   integer(kind=int32) :: ctr_s(4), key_s(4)
   real(kind_physics):: rand_num(8)

   ! variables related to cross sections and probabilistic collisions
   real(kind_physics), dimension(:), allocatable :: cross_sections_vector
   real(kind_physics) :: abs_max_CS, neutral_density, init_temperature, pressure, &
                         charge_count(2), total_charge_count(2)

   ! lookup tables for cross section data
   character(255) :: file_path
   type(linked_list_CS), pointer :: CS_tables, CS_guide

   ! variables describing external fields
   real(kind_physics) :: major_radius, minor_radius, B0, B_p, V_loop

   ! constants & scaling factors
   real(kind_physics) :: c = 299792458.0 ! m/s
   real(kind_physics) :: e_mass = 510998.9461 ! eV/c^2
   real(kind_physics) :: e = 1.6021766208e-19 ! Coulomb
   real(kind_physics) :: eps_0 = 8.85418781762e-12 ! Coulomb/Vm
   real(kind_physics) :: pi = 3.141592653589793238462643383279502884197
   real(kind_physics) :: kboltzmann = 1.38064852e-23 ! J K^(-1)
   real(kind_physics) :: E_q_dt_m

   interface random
      module procedure random8, random16
   end interface

contains

   subroutine set_parameter()

      use module_pepc
      use module_interaction_specific, only: theta2, eps2, force_law, include_far_field_if_periodic
      implicit none

      integer, parameter :: fid = 12
      character(255)     :: para_file
      logical            :: read_para_file

      namelist /pepcbreakup/ tnp, nt, dt, particle_output, domain_output, reflecting_walls, &
         particle_test, diag_interval, plasma_dimensions, init_temperature, pressure, external_e, &
         major_radius, minor_radius, B0, B_p, V_loop

      ! set default parameter values
      tnp = 10000
      nt = 25
      dt = 1e-2
      particle_test = .false.
      particle_output = .false.
      domain_output = .false.
      reflecting_walls = .false.
      diag_interval = 1
      plasma_dimensions = (/1.0_8, 1.0_8, 1.0_8/)
      external_e = (/0.0_8, 0.0_8, 0.0_8/)
      electron_num = 0

      ! read in namelist file
      call pepc_read_parameters_from_first_argument(read_para_file, para_file)

      if (read_para_file) then
         if (root) write (*, '(a)') " == reading parameter file, section pepcbreakup: ", para_file
         open (fid, file=para_file)
         read (fid, NML=pepcbreakup)
         close (fid)
      else
         if (root) write (*, *) " == no param file, using default parameter "
      end if

      if (root) then
         write (*, '(a,i12)') " == total number of particles           : ", tnp
         write (*, '(a,i12)') " == number of time steps                : ", nt
         write (*, '(a,es12.4)') " == time step                           : ", dt
         write (*, '(a,i12)') " == diag & IO interval                  : ", diag_interval
         write (*, '(a,l12)') " == particle test                       : ", particle_test
         write (*, '(a,l12)') " == particle output                     : ", particle_output
         write (*, '(a,l12)') " == domain output                       : ", domain_output
         write (*, '(a,l12)') " == reflecting walls                    : ", reflecting_walls
         write (*, '(a,3(es12.4))') " == plasma dimensions                   : ", plasma_dimensions
         write (*, '(a,es12.4)') " == initial electron temperature(K)     : ", init_temperature
         write (*, '(a,es12.4)') " == pressure(Pa)                        : ", pressure
         write (*, '(a,3(es12.4))') " == external electric field(V/m)        : ", external_e
         write (*, '(a,es12.4)') " == major radius(m)                     : ", major_radius
         write (*, '(a,es12.4)') " == minor radius(m)                     : ", minor_radius
         write (*, '(a,es12.4)') " == toroidal magnetic field strength (T): ", B0
         write (*, '(a,es12.4)') " == poloidal magnetic field strength (T): ", B_p
         write (*, '(a,es12.4)') " == toroidal loop voltage(V)            : ", V_loop
      end if

      ! NOTE: Scale the appropriate read-in variables!
      !       1. multiply 4*pi*eps_0*(1e-12 sec*light speed)^2/(elementary charge) to electric field (V/m)
      !       2. multiply 4*pi*eps_0*(1e-12 sec*light speed)/(elementary charge) to voltage
      !       3. multiply ((1.e-12*c)**2)/e_mass to B field (Tesla)
      !       4. divide length variables with (1.e-12*c)
      external_e = external_e*4.0*pi*eps_0*((c*1e-12)**2)/e
      V_loop = V_loop*4.0*pi*eps_0*((c*1e-12))/e
      B0 = B0 * ((1.e-12*c)**2)/e_mass
      B_p = B_p * ((1.e-12*c)**2)/e_mass
      plasma_dimensions = plasma_dimensions/(c*1e-12)
      major_radius = major_radius/(c*1e-12)
      minor_radius = minor_radius/(c*1e-12)

      call pepc_prepare(3_kind_dim)
   end subroutine set_parameter

   function torus_geometry() result(pos)
     implicit none
     real(kind_physics) :: pos(3), ran(3)
     real*8 :: theta, phi, l, r

     call random_number(ran)
     r = ran(1) * minor_radius
     phi = ran(2) * 2. * pi
     theta  = ran(3) * 2. * pi
     l = major_radius + r*sin(phi)
     pos(1) = l * sin(theta)
     pos(2) = l * cos(theta)
     pos(3) = r * cos(phi)
   end function torus_geometry

   function thermal_velocity_mag(mass, temp) result(velocity)
     implicit none
     real(kind_physics), intent(in) :: mass, temp
     real(kind_physics) :: velocity, kb_ev

     kb_ev = 8.61733035e-5 !boltzmann constant in eV/K
     velocity = sqrt(temp*kb_ev/(mass*e_mass)) !dimensionless velocity
    !  print *, velocity
   end function thermal_velocity_mag

   subroutine init_particles(p)
      implicit none

      type(t_particle), allocatable, intent(inout) :: p(:)
      integer(kind_particle) :: ip
      integer :: rc
      real*8 :: dummy
      real(kind_physics) :: magnitude

      if (root) write (*, '(a)') " == [init] init particles "

      ! set initially number of local particles
      np = tnp/n_ranks
      if (my_rank < MOD(tnp, 1_kind_particle*n_ranks)) np = np + 1
      !1_kind_particle means a value of 1, enforced into specified precision defined by kind_particle

      allocate (p(np), stat=rc)
      if (rc .ne. 0) write (*, *) " === particle allocation error!"

      allocate (direct_L2(np), stat=rc)
      if (rc .ne. 0) write (*, *) " === direct_L2 allocation error!"
      direct_L2 = -1.0_8

      ! set random seed
      dummy = par_rand(1*my_rank)

      ! setup random qubic particle cloud
      do ip = 1, np
         p(ip)%label = my_rank*(tnp/n_ranks) + ip - 1
         p(ip)%data%q = -1.0_8
         p(ip)%data%m = 1.0_8
         p(ip)%data%species = 0

         p(ip)%data%age = 0.0_8

         call random(p(ip)%x)
         p(ip)%x = p(ip)%x*plasma_dimensions - plasma_dimensions*0.5
         p(ip)%x(3) = 0.0
        !  p(ip)%x = torus_geometry()

        !  call random_gauss(p(ip)%data%v)
        !  p(ip)%data%v = p(ip)%data%v/c * 1e6

         magnitude = thermal_velocity_mag(p(ip)%data%m, 873.15_kind_physics)! 0.0108359158316
         p(ip)%data%v = 0.0
         p(ip)%data%v(3) = -1.0 * magnitude

         p(ip)%work = 1.0_8
      end do
   end subroutine init_particles

   subroutine push_particles(p)
      use module_mirror_boxes
      implicit none

      type(t_particle), allocatable, intent(inout) :: p(:)
      integer(kind_particle) :: ip
      real*8  :: fact

      if (root) write (*, '(a)') " == [pusher] push particles "

      fact = dt

      do ip = 1, np
         p(ip)%data%v = p(ip)%data%v + fact*p(ip)%data%q/p(ip)%data%m*p(ip)%results%e
         p(ip)%x = p(ip)%x + dt*p(ip)%data%v
      end do
   end subroutine push_particles

   subroutine filter_particles(p)
      implicit none
      include 'mpif.h'

      type(t_particle), allocatable, intent(inout) :: p(:)
      integer(kind_particle) :: ip
      integer :: id, ncoll, ncoll_total, ierr

      ncoll = 0

      do ip = 1, np
         do id = 1, 3
            if (p(ip)%x(id) < 0.0_8) then
               p(ip)%x(id) = -p(ip)%x(id)
               p(ip)%data%v(id) = -p(ip)%data%v(id)
               ncoll = ncoll + 1
            else if (p(ip)%x(id) > plasma_dimensions(id)) then
               p(ip)%x(id) = 2*plasma_dimensions(id) - p(ip)%x(id)
               p(ip)%data%v(id) = -p(ip)%data%v(id)
               ncoll = ncoll + 1
            end if
         end do
      end do

      call mpi_reduce(ncoll, ncoll_total, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      if (root) write (*, '(a,i12)') " == [filter] total number of wall collisions      : ", ncoll_total
   end subroutine filter_particles

   subroutine test_particles()
      use module_pepc_types
      use module_directsum
      implicit none
      include 'mpif.h'

      integer(kind_particle), allocatable   :: tindx(:)
      real*8, allocatable                   :: trnd(:)
      type(t_particle_results), allocatable :: trslt(:)
      integer(kind_particle)                :: tn, tn_global, ti
      integer                               :: rc
      real(kind_physics)                    :: L2sum_local, L2sum_global, L2

      call timer_start(t_user_directsum)

      if (allocated(direct_L2)) then
         deallocate (direct_L2)
      end if
      allocate (direct_L2(np))
      direct_L2 = -1.0_8

      if (particle_direct .eq. -1) then
         tn = np
      else
         tn = particle_direct/n_ranks
         if (my_rank .eq. (n_ranks - 1)) tn = tn + MOD(particle_direct, n_ranks)
      endif

      allocate (tindx(tn), trnd(tn), trslt(tn))

      if (particle_direct .eq. -1) then
         do ti = 1, tn
            tindx(ti) = ti
         enddo
      else
         call random(trnd(1:tn))

         tindx(1:tn) = int(trnd(1:tn)*(np - 1)) + 1
      endif

      call directforce(particles, tindx, tn, trslt, MPI_COMM_WORLD)

      L2sum_local = 0.0
      L2sum_global = 0.0
      do ti = 1, tn
         L2 = &
            (particles(tindx(ti))%results%e(1) - trslt(ti)%e(1))**2 + &
            (particles(tindx(ti))%results%e(2) - trslt(ti)%e(2))**2 + &
            (particles(tindx(ti))%results%e(3) - trslt(ti)%e(3))**2
         L2sum_local = L2sum_local + L2
         direct_L2(tindx(ti)) = L2
      end do

      call MPI_ALLREDUCE(tn, tn_global, 1, MPI_KIND_PARTICLE, MPI_SUM, MPI_COMM_WORLD, rc)
      call MPI_ALLREDUCE(L2sum_local, L2sum_global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)

      L2sum_global = sqrt(L2sum_global)/tn_global

      call timer_stop(t_user_directsum)
      if (root) then
         write (*, '(a,i12)') " == [direct test] number tested particles         : ", tn_global
         write (*, '(a,es12.4)') " == [direct test] L2 error in probed particles    : ", L2sum_global
         write (*, '(a,es12.4)') " == [direct test] time in test [s]                : ", timer_read(t_user_directsum)
      end if

      deallocate (tindx)
      deallocate (trnd)
      deallocate (trslt)
   end subroutine test_particles

   integer function vtk_step_of_step(step) result(vtk_step)
      use module_vtk
      implicit none

      integer, intent(in) :: step

      if (step .eq. 0) then
         vtk_step = VTK_STEP_FIRST
      else if (nt - 1 - step < diag_interval) then
         vtk_step = VTK_STEP_LAST
      else
         vtk_step = VTK_STEP_NORMAL
      endif
   end function vtk_step_of_step

   subroutine write_text_output(anode_charge_count, cathode_charge_count, step)
     implicit none
     real(kind_physics), intent(in) :: anode_charge_count, cathode_charge_count
     integer, intent(in) :: step

    !  open(12, file = fname, action='WRITE')
     if (step == 0) then
       write(12, *) 'step ', 'anode_count ', 'cathode_count'
     end if
     write(12, *) step, anode_charge_count, cathode_charge_count
    !  close(12)
   end subroutine write_text_output

   subroutine write_particles(p)
      use module_vtk_helpers
      implicit none

      include 'mpif.h'

      type(t_particle), intent(in) :: p(:)

      integer :: vtk_step

      call timer_start(t_user_particleio)
      vtk_step = vtk_step_of_step(step)
      call vtk_write_particles("particles", MPI_COMM_WORLD, step, dt*step, vtk_step, p, coulomb_and_l2)
      call timer_stop(t_user_particleio)
      if (root) write (*, '(a,es12.4)') " == [write particles] time in vtk output [s]      : ", timer_read(t_user_particleio)

   contains

      subroutine coulomb_and_l2(d, r, vtkf)
         use module_vtk
         use module_interaction_specific_types
         implicit none

         type(t_particle_data), intent(in) :: d(:)
         type(t_particle_results), intent(in) :: r(:)
         type(vtkfile_unstructured_grid), intent(inout) :: vtkf

         call vtk_write_particle_data_results(d, r, vtkf)
         if (particle_test) call vtkf%write_data_array("L2 error", direct_L2(:))
      end subroutine
   end subroutine write_particles

   subroutine write_domain(p)
      use module_vtk
      use module_vtk_helpers
      use module_pepc, only: global_tree
      implicit none

      type(t_particle), allocatable, intent(in) :: p(:)

      integer :: vtk_step

      ! output of tree diagnostics
      vtk_step = vtk_step_of_step(step)
      call vtk_write_branches(step, dt*step, vtk_step, global_tree)
      call vtk_write_leaves(step, dt*step, vtk_step, global_tree)
      call vtk_write_spacecurve(step, dt*step, vtk_step, p)
   end subroutine write_domain

   subroutine random_gauss(list)
      implicit none

      real(kind_physics), intent(inout) :: list(:)

      real(kind_physics) :: v(2), pi, r, p
      integer :: n, i

      pi = 2.0_kind_physics*acos(0.0_kind_physics)
      n = size(list)

      do i = 1, n, 2

         call random(v)

         r = sqrt(-2.0_8*log(v(1)))
         p = 2.0_8*pi*v(2)

         list(i) = r*sin(p)
         if ((i + 1) <= n) list(i + 1) = r*cos(p)

      end do
   end subroutine

   subroutine random8(array)
      implicit none
      real*8 :: array(:)
      integer :: i

      do i = 1, size(array)
         array(i) = par_rand()
      end do
   end subroutine random8

   subroutine random16(array)
      implicit none
      real*16 :: array(:)
      integer :: i

      do i = 1, size(array)
         array(i) = par_rand()
      end do
   end subroutine random16

   !>
   !> portable random number generator, see numerical recipes
   !> check for the random numbers:
   !> the first numbers should be 0.2853809, 0.2533582 and 0.0934685
   !> the parameter iseed is optional
   !>
   function par_rand(iseed)
      implicit none
      real :: par_rand
      integer, intent(in), optional :: iseed

      integer, parameter :: IM1 = 2147483563
      integer, parameter :: IM2 = 2147483399
      real, parameter :: AM = 1.0/IM1
      integer, parameter :: IMM1 = IM1 - 1
      integer, parameter :: IA1 = 40014
      integer, parameter :: IA2 = 40692
      integer, parameter :: IQ1 = 53668
      integer, parameter :: IQ2 = 52774
      integer, parameter :: IR1 = 12211
      integer, parameter :: IR2 = 3791
      integer, parameter :: NTAB = 32
      integer, parameter :: NDIV = 1 + IMM1/NTAB
      real, parameter :: eps_ = 1.2e-7 ! epsilon(eps_)
      real, parameter :: RNMX = 1.0 - eps_

      integer :: j, k
      integer, volatile, save :: idum = -1
      integer, volatile, save :: idum2 = 123456789
      integer, volatile, save :: iy = 0
      integer, volatile, save :: iv(NTAB)

      if (idum <= 0 .or. present(iseed)) then
         if (present(iseed)) then
            idum = iseed
         else
            if (-idum < 1) then
               idum = 1
            else
               idum = -idum
            endif
         endif

         idum2 = idum

         do j = NTAB + 7, 0, -1
            k = idum/IQ1
            idum = IA1*(idum - k*IQ1) - k*IR1
            if (idum < 0) idum = idum + IM1

            if (j < NTAB) iv(j + 1) = idum

         end do
         iy = iv(1)
      end if

      k = idum/IQ1
      idum = IA1*(idum - k*IQ1) - k*IR1
      if (idum < 0) idum = idum + IM1

      k = idum2/IQ2
      idum2 = IA2*(idum2 - k*IQ2) - k*IR2
      if (idum2 < 0) idum2 = idum2 + IM2

      j = iy/NDIV + 1
      iy = iv(j) - idum2
      iv(j) = idum

      if (iy < 1) iy = iy + IMM1
      par_rand = AM*iy
      if (par_rand > RNMX) par_rand = RNMX
   end function par_rand
end module
