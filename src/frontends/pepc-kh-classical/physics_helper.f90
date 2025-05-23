module physics_helper
   use module_pepc_kinds

   implicit none

   type physics_nml_t
      real(kind_physics) :: m_ratio, T_ratio, B0, shear_halfwidth, shear_strength
      real(kind_physics), dimension(3) :: l_plasma
      integer(kind=kind_particle) :: ni
   end type physics_nml_t

   ! 1 / (2 * pi)
   real(kind_physics), parameter :: force_const = 0.1591549430918953357688837633725143620345_kind_physics

   integer, parameter :: file_energy = 70

   integer, parameter :: LABEL_ELECTRON_LEFT  = -1                              !&
   integer, parameter :: LABEL_ELECTRON_RIGHT = +1                              !&
   integer, parameter :: LABEL_ION_LEFT       = -2                              !&
   integer, parameter :: LABEL_ION_RIGHT      = +2                              !&

contains

   subroutine setup_physics(physics_pars, time_pars, p, pepc_pars)
      use module_pepc, only: pepc_prepare
      use module_pepc_types, only: t_particle
      use module_checkpoint
      use module_debug
      use module_mirror_boxes, only: mirror_box_layers, t_lattice_1, t_lattice_2, t_lattice_3, &
                                     periodicity, spatial_interaction_cutoff
      use module_interaction_specific, only: include_far_field_if_periodic
      use encap
      use pepc_helper, only: para_file_available, para_file_name
      implicit none

      type(physics_pars_t), intent(out) :: physics_pars
      type(time_pars_t), intent(in) :: time_pars
      type(t_particle), allocatable, dimension(:), intent(inout) :: p
      type(pepc_pars_t), intent(inout) :: pepc_pars

      type(physics_nml_t) :: physics_nml
      integer :: dummy_nresume
      character(len=255) :: dummy_file_name

      call read_in_physics_params(physics_nml, para_file_available, para_file_name)

      physics_pars%B0 = physics_nml%B0
      physics_pars%l_plasma = physics_nml%l_plasma
      physics_pars%vte = 1.
      physics_pars%vti = physics_pars%vte * sqrt(physics_nml%T_ratio / physics_nml%m_ratio)
      physics_pars%qe = -physics_pars%l_plasma(1) * physics_pars%l_plasma(2) / physics_nml%ni
      physics_pars%qi = abs(physics_pars%qe)
      physics_pars%me = physics_pars%qi
      physics_pars%mi = physics_nml%m_ratio * physics_pars%qi
      physics_pars%shear_halfwidth = physics_nml%shear_halfwidth
      physics_pars%shear_velocity = physics_nml%shear_strength * physics_pars%qi * &
                                    physics_pars%B0 / physics_pars%mi * physics_pars%shear_halfwidth
      physics_pars%ni = physics_nml%ni

      t_lattice_1 = [physics_pars%l_plasma(1), 0._kind_physics, 0._kind_physics]
      t_lattice_2 = [0._kind_physics, physics_pars%l_plasma(2), 0._kind_physics]
      t_lattice_3 = [0._kind_physics, 0._kind_physics, 1._kind_physics]
      if (.not. include_far_field_if_periodic) then
         spatial_interaction_cutoff = huge(1._kind_physics)*[1., 1., 1.]
         spatial_interaction_cutoff(2) = max(1, mirror_box_layers) * physics_pars%l_plasma(2)
      end if
      periodicity = [.false., .true., .false.]

      call pepc_prepare(2_kind_dim)

      if (time_pars%nresume .gt. 0) then
         call read_particles_mpiio(time_pars%nresume, pepc_pars%pepc_comm%mpi_comm, &
                                   dummy_nresume, pepc_pars%np, p, dummy_file_name)

         if (dummy_nresume .ne. time_pars%nresume) then
            DEBUG_ERROR(*, "Resume timestep mismatch, parameter file says: ", time_pars%nresume, " checkpoint file says: ", dummy_nresume)
         end if

         if (pepc_pars%pepc_comm%mpi_rank .eq. 0) then
            print *, "== [resume]"
            print *, "   resuming with ", pepc_pars%np, " particles."
            print *, " "
         end if
      else
         call special_start(pepc_pars, physics_pars, p)
      end if

   end subroutine setup_physics

   subroutine read_in_physics_params(physics_namelist, file_available, file_name)
      use mpi
      implicit none

      type(physics_nml_t), intent(out) :: physics_namelist
      logical, intent(in) :: file_available
      character(len=255), intent(in) :: file_name

      real(kind_physics) :: m_ratio = 100.
      real(kind_physics) :: T_ratio = 1.
      real(kind_physics) :: B0 = 1.
      real(kind_physics), dimension(3) :: l_plasma = [50., 125., 0.]
      real(kind_physics) :: shear_halfwidth = 2.
      real(kind_physics) :: shear_strength = 1.
      integer(kind=kind_particle) :: ni = 0

      namelist /physics_nml/ m_ratio, T_ratio, B0, l_plasma, shear_halfwidth, shear_strength, ni

      integer, parameter :: para_file_id = 10

      if (file_available) then
         open (para_file_id, file=trim(file_name), action='read')
         rewind (para_file_id)
         read (para_file_id, NML=physics_nml)
         close (para_file_id)
      end if

      physics_namelist%m_ratio = m_ratio
      physics_namelist%T_ratio = T_ratio
      physics_namelist%B0 = B0
      physics_namelist%l_plasma = l_plasma
      physics_namelist%shear_halfwidth = shear_halfwidth
      physics_namelist%shear_strength = shear_strength
      physics_namelist%ni = ni

   end subroutine read_in_physics_params

   subroutine write_physics_params(physics_pars, file_name)
      use encap
      implicit none

      type(physics_pars_t), intent(in) :: physics_pars
      character(len=255), intent(in) :: file_name

      real(kind_physics) :: m_ratio
      real(kind_physics) :: T_ratio
      real(kind_physics) :: shear_halfwidth
      real(kind_physics) :: shear_strength
      real(kind_physics) :: B0
      real(kind_physics), dimension(3) :: l_plasma

      namelist /physics_nml/ m_ratio, T_ratio, B0, l_plasma, shear_halfwidth, shear_strength

      integer, parameter :: para_file_id = 10

      m_ratio = physics_pars%mi / physics_pars%me
      T_ratio = m_ratio * (physics_pars%vti / physics_pars%vte)**2
      B0 = physics_pars%B0
      l_plasma = physics_pars%l_plasma
      shear_halfwidth = physics_pars%shear_halfwidth
      shear_strength = physics_pars%shear_velocity * physics_pars%mi / &
                       (physics_pars%qi * physics_pars%B0 * physics_pars%shear_halfwidth)

      open (para_file_id, file=trim(file_name), status='old', position= &
            'append', action='write')
      write (para_file_id, nml=physics_nml)
      close (para_file_id)

   end subroutine write_physics_params

   subroutine special_start(pepc_pars, physics_pars, p)
      use module_pepc_types, only: t_particle
      use module_mirror_boxes, only: constrain_periodic
      use module_debug
      use module_rng
      use encap
      use time_helper
      use mpi
      implicit none

      type(pepc_pars_t), intent(inout) :: pepc_pars
      type(physics_pars_t), intent(in) :: physics_pars
      type(t_particle), allocatable, intent(inout) :: p(:)

      type(t_particle), allocatable :: pt(:)
      integer(kind_particle) :: ni, ne, ip, i, nil, nel, nebefore
      integer(kind_default) :: mpi_rank, mpi_size, mpi_err
      integer, allocatable :: period(:)
      integer :: it, ic
      real(kind_physics) :: lambda, dx, lx, ly, qi, qe, mi, me, vti, vte, B0, rwci, rwce, v0
      real(kind_physics), allocatable :: xi(:, :), vi(:, :), xp(:)
      type(time_pars_t) :: time_pars
      integer(kind_default), parameter :: nhist = 128
      integer(kind_particle), allocatable :: nihist(:), nehist(:)

      real(kind_physics), parameter :: pi = 3.141592653589793238462643383279502884197_kind_physics

      lx = physics_pars%l_plasma(1)
      ly = physics_pars%l_plasma(2)
      qi = physics_pars%qi
      qe = physics_pars%qe
      mi = physics_pars%mi
      me = physics_pars%me
      vti = physics_pars%vti
      vte = physics_pars%vte
      B0 = physics_pars%B0
      lambda = physics_pars%shear_halfwidth
      v0 = physics_pars%shear_velocity
      ni = physics_pars%ni

      rwce = abs(me / (qe * B0))
      rwci = abs(mi / (qi * B0))

      mpi_rank = pepc_pars%pepc_comm%mpi_rank
      mpi_size = pepc_pars%pepc_comm%mpi_size
      nil = ni / mpi_size
      if (mpi_rank .lt. mod(ni, int(mpi_size, kind=kind_particle))) nil = nil + 1

      if (mpi_rank .eq. 0) then
         print *, "== [special_start]"
         print *, "   Lx x Ly = ", lx, " x ", ly
         print *, "   mi / me = ", mi / me
         print *, "   Ti / Te = ", (mi * vti**2) / (me * vte**2)
         print *, "   B0      = ", B0
         print *, " "
      end if

      allocate (pt(nil), xi(nil, 3), vi(nil, 3), xp(nil))

      ! (1) assume a smooth profile for the electric polarization field, hence for the E x B drift velocity

      do ip = 1, nil
         ! (2) load the ions into the simulation domain with a uniform distribution along the x axis
         pt(ip)%work = 1.
         pt(ip)%data%q = qi
         pt(ip)%data%m = mi

         pt(ip)%x(1) = lx * rng_next_real()
         pt(ip)%x(2) = ly * rng_next_real()
         pt(ip)%x(3) = 0.

         xi(ip, :) = pt(ip)%x(:)
         xp(ip) = pt(ip)%x(1)

         ! (3) give the ions random initial velocities in the positive x direction with a Rayleigh distribution at the desired
         ! temperature
         pt(ip)%data%v(1) = max(velocity_1d_rayleigh(vti), real(1.0E-6, kind=kind_physics))
         ! and in the y direction give them the local E x B drift velocity
         pt(ip)%data%v(2) = vdrift(pt(ip)%x(1))
         pt(ip)%data%v(3) = 0.

         vi(ip, :) = pt(ip)%data%v(:)
      end do

      ! (4) advance each ion along its trajectory, in the presence of the constant but nonuniform electric field, for a random
      ! fraction of its orbital period
      allocate (period(nil))
      period = 0

      time_pars%dt = 2.*pi * rwci / 160.

      it = 0
      do
         do ip = 1, nil
            pt(ip)%results%e(1) = estatic(pt(ip)%x(1))
            pt(ip)%results%e(2) = 0.
         end do

         call push_particles(time_pars, physics_pars, pt(:))
         it = it + 1

         do ip = 1, nil
            if (period(ip) .eq. 0 .and. xp(ip) .lt. xi(ip, 1) .and. pt(ip)%x(1) .ge. xi(ip, 1)) period(ip) = it
            xp(ip) = pt(ip)%x(1)
         end do

         if (all(period .ne. 0)) then
            exit
         end if
      end do

      do ip = 1, nil
         pt(ip)%x(:) = xi(ip, :)
         pt(ip)%data%v(:) = vi(ip, :)

         do it = 1, int(rng_next_real() * period(ip))
            pt(ip)%results%e(1) = estatic(pt(ip)%x(1))
            pt(ip)%results%e(2) = 0.

            call push_particles(time_pars, physics_pars, pt(ip:ip))
         end do
      end do

      deallocate (period, xi, vi, xp)

      call constrain_periodic(pt)
      do ip = 1, nil
         if ((pt(ip)%x(1) .gt. lx) .or. (pt(ip)%x(1) .lt. 0.)) then
            pt(ip)%data%v(1) = -pt(ip)%data%v(1)
            pt(ip)%data%v(2) = 2.*vdrift(pt(ip)%x(1)) - pt(ip)%data%v(2)

            pt(ip)%x(1) = lx - modulo(pt(ip)%x(1), lx)
         end if
      end do

      ! (5) calculate the (ensemble average) ion density profile
      allocate (nihist(nhist), nehist(nhist))
      nihist = 0
      dx = lx / nhist

      do ip = 1, nil
         ic = min(ceiling(pt(ip)%x(1) / dx), nhist)
         nihist(ic) = nihist(ic) + 1
      end do

      call mpi_allreduce(MPI_IN_PLACE, nihist, nhist, MPI_KIND_PARTICLE, MPI_SUM, pepc_pars%pepc_comm%mpi_comm, mpi_err)

      !print *, "nihist: ", nihist
      !print *, "sum nihist: ", sum(nihist)

      ! (6), (7) calculate the net charge density by the Poisson equation from the electric field
      ! and calculate the electron density by subtracting the charge density from the ion density
      nehist = 0
      do ic = 1, nhist
         nehist(ic) = nihist(ic) + floor(ly * dx * estaticp(dx * (ic - 0.5)) / qe)
      end do

      ne = sum(nehist)

      !print *, "nehist: ", nehist
      !print *, "sum nehist: ", ne

      pepc_pars%np = ne + ni

      nel = ne / mpi_size
      if (mpi_rank .lt. mod(ne, int(mpi_size, kind=kind_particle))) nel = nel + 1

      allocate (p(nil + nel))
      p(1:nil) = pt(:)
      deallocate (pt)

      nebefore = (ne / mpi_size + 1) * min(int(mpi_rank, kind=kind_particle), mod(ne, int(mpi_size, kind=kind_particle))) + &
                 (ne / mpi_size) * max(0_kind_particle, mpi_rank - mod(ne, int(mpi_size, kind=kind_particle)))

      do ic = 1, nhist
         if (nebefore .gt. nehist(ic)) then
            nebefore = nebefore - nehist(ic)
            nehist(ic) = 0
         else
            nehist(ic) = nehist(ic) - nebefore
            exit
         end if
      end do

      ip = nil
      do ic = 1, nhist
         do i = 1, nehist(ic)
            ip = ip + 1

            p(ip)%work = 1.
            p(ip)%data%q = qe
            p(ip)%data%m = me

            ! (8) load the electrons into the simulation domain in such a way as to yield the required electron density
            p(ip)%x(1) = dx * (ic - 1.+rng_next_real())
            p(ip)%x(2) = ly * rng_next_real()
            p(ip)%x(3) = 0.

            ! (9) give the electrons random initial velocities with a displaced Maxwellian distribution at the modified temperature
            p(ip)%data%v(1:2) = velocity_2d_maxwell(vte / sqrt(sqrt(1 + vdriftp(p(ip)%x(1)) * rwce)))
            p(ip)%data%v(2) = p(ip)%data%v(2) + vdrift(p(ip)%x(1))
            p(ip)%data%v(3) = 0.

            if (ip .eq. nil + nel) exit
         end do

         if (ip .eq. nil + nel) exit
      end do

      deallocate (nihist, nehist)

      !print *, "num particles: ", ip

      !vte = 0.
      !vti = 0.
      !do ip = 1, nil
      !  vti = vti + p(ip)%data%v(1)**2 + (p(ip)%data%v(2) - vdrift(p(ip)%x(1)))**2
      !end do
      !do ip = nil + 1, nil + nel
      !  vte = vte + p(ip)%data%v(1)**2 + (p(ip)%data%v(2) - vdrift(p(ip)%x(1)))**2
      !end do

      !call mpi_allreduce(MPI_IN_PLACE, vti, 1, MPI_KIND_PHYSICS, MPI_SUM, pepc_pars%pepc_comm%mpi_comm, mpi_err)
      !call mpi_allreduce(MPI_IN_PLACE, vte, 1, MPI_KIND_PHYSICS, MPI_SUM, pepc_pars%pepc_comm%mpi_comm, mpi_err)

      !if (mpi_rank == 0) &
      !  print *, "Ti / Te = ", (mi * vti) / (me * vte)

      ! assign negative labels to all left particles, positive labels to right particles
      do ip = 1, nil
         if (p(ip)%x(1) .le. lx / 2.) then
            p(ip)%label = LABEL_ION_LEFT
         else
            p(ip)%label = LABEL_ION_RIGHT
         end if
      end do

      do ip = nil + 1, size(p, kind=kind(ip))
         if (p(ip)%x(1) .le. lx / 2.) then
            p(ip)%label = LABEL_ELECTRON_LEFT
         else
            p(ip)%label = LABEL_ELECTRON_RIGHT
         end if
      end do

   contains

      function vdrift(x)
         implicit none

         real(kind_physics) :: vdrift
         real(kind_physics), intent(in) :: x

         vdrift = v0 * tanh((x - lx / 2.) / lambda)
      end function

      function vdriftp(x)
         implicit none

         real(kind_physics) :: vdriftp
         real(kind_physics), intent(in) :: x

         vdriftp = v0 * (1.-tanh((x - lx / 2.) / lambda)**2.) / lambda
      end function

      function estatic(x)
         implicit none

         real(kind_physics) :: estatic
         real(kind_physics), intent(in) :: x

         estatic = -B0 * vdrift(x)
      end function

      function estaticp(x)
         implicit none

         real(kind_physics) :: estaticp
         real(kind_physics), intent(in) :: x

         estaticp = -B0 * vdriftp(x)
      end function

   end subroutine special_start

   function velocity_2d_maxwell(vt) result(v)
      use module_rng
      implicit none

      real(kind_physics), intent(in) :: vt
      real(kind_physics), dimension(2) :: v

      real(kind_physics) :: xi, v0
      real(kind_physics), parameter :: pi = 3.141592653589793238462643383279502884197_kind_physics

      xi = rng_next_real()
      v0 = vt * sqrt(-2.*log(max(xi, 10.0_kind_physics**(-15))))
      xi = rng_next_real()
      v(1) = v0 * cos(2 * pi * xi)
      v(2) = v0 * sin(2 * pi * xi)

   end function velocity_2d_maxwell

   function velocity_1d_rayleigh(vt) result(v)
      use module_rng
      implicit none

      real(kind_physics), intent(in) :: vt
      real(kind_physics) :: v

      real(kind_physics) :: xi

      xi = rng_next_real()
      v = vt * sqrt(-2.*log(xi))
   end function velocity_1d_rayleigh

   subroutine physics_dump(pepc_pars, physics_pars, time_pars, step, p)
      use mpi
      use module_pepc_types
      use encap
      implicit none

      type(pepc_pars_t), intent(in) :: pepc_pars
      type(physics_pars_t), intent(in) :: physics_pars
      type(time_pars_t), intent(in) :: time_pars
      integer, intent(in) :: step
      type(t_particle), dimension(:), intent(in) :: p

      integer(kind_particle) :: ip
      integer(kind_default) :: mpi_err
      real(kind_physics) :: e_kin, e_pot, e_kin_g, e_pot_g

      e_kin = 0.
      e_pot = 0.

      do ip = 1, size(p)
         e_kin = e_kin + e_kin_of_particle(p(ip))
         e_pot = e_pot + e_pot_of_particle(p(ip))
      end do

      call mpi_reduce(e_kin, e_kin_g, 1, MPI_KIND_PHYSICS, MPI_SUM, 0, pepc_pars%pepc_comm%mpi_comm, mpi_err)
      call mpi_reduce(e_pot, e_pot_g, 1, MPI_KIND_PHYSICS, MPI_SUM, 0, pepc_pars%pepc_comm%mpi_comm, mpi_err)

      if (pepc_pars%pepc_comm%mpi_rank .eq. 0) then
         call write_to_stdout()
         call write_to_file()
      end if

   contains

      subroutine write_to_stdout()
         implicit none

         print *, "== [physics_dump]"
         print *, "   T : ", e_kin_g
         print *, "   V : ", e_pot_g
         print *, "   H : ", e_kin_g + e_pot_g
         print *, " "

      end subroutine write_to_stdout

      subroutine write_to_file()
         implicit none

         if (step .eq. 0) then
            open (file_energy, file='energy.csv', status='replace', action='write')
            write (file_energy, *) "t, T, Tc, V, H"
         else
            open (file_energy, file='energy.csv', status='old', position='append', action='write')
         end if

         if (.not. (time_pars%nresume .gt. 0 .and. time_pars%nresume .eq. step)) then
            write (file_energy, *) time_pars%dt * step, ", ", e_kin_g, ", ", &
               e_pot_g, ", ", (e_kin_g + e_pot_g)
         end if
         close (file_energy)

      end subroutine write_to_file

   end subroutine physics_dump

   pure function e_kin_of_particle(p)
      use module_pepc_types, only: t_particle
      implicit none

      type(t_particle), intent(in) :: p
      real(kind_physics) :: e_kin_of_particle

      e_kin_of_particle = 0.5 * p%data%m * dot_product(p%data%v, p%data%v)
   end function e_kin_of_particle

   pure function e_pot_of_particle(p)
      use module_pepc_types, only: t_particle
      implicit none

      type(t_particle), intent(in) :: p
      real(kind_physics) :: e_pot_of_particle

      e_pot_of_particle = 0.5 * p%results%pot * p%data%q
   end function e_pot_of_particle

end module physics_helper
