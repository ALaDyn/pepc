module physics_helper

   implicit none

   type physics_nml_t
      real(kind=8) :: m_ratio, T_ratio, B0, shear_halfwidth, shear_strength
      real(kind=8), dimension(3) :: l_plasma
   end type physics_nml_t

   real(kind = 8), parameter :: force_const = 0.159154943D0

   integer, parameter :: file_energy = 70

contains


   subroutine setup_physics(physics_pars, time_pars, p, pepc_pars)
      use module_pepc, only: pepc_prepare
      use module_pepc_types, only: t_particle, kind_dim
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
      character(len = 255) :: dummy_file_name

      call read_in_physics_params(physics_nml, para_file_available, para_file_name)

      physics_pars%B0 = physics_nml%B0
      physics_pars%l_plasma = physics_nml%l_plasma
      physics_pars%vte = 1.0D0
      physics_pars%vti = physics_pars%vte * sqrt(physics_nml%T_ratio / physics_nml%m_ratio)
      physics_pars%qe  = -2.0D0 * physics_pars%l_plasma(1) * physics_pars%l_plasma(2) / pepc_pars%np
      physics_pars%qi  = abs(physics_pars%qe)
      physics_pars%me  = physics_pars%qi
      physics_pars%mi  = physics_nml%m_ratio * physics_pars%qi
      physics_pars%shear_halfwidth = physics_nml%shear_halfwidth
      physics_pars%shear_velocity = physics_nml%shear_strength * physics_pars%qi * &
        physics_pars%B0 / physics_pars%mi * physics_pars%shear_halfwidth

      t_lattice_1 = [ physics_pars%l_plasma(1), 0.0D0, 0.0D0 ]
      t_lattice_2 = [ 0.0D0, physics_pars%l_plasma(2), 0.0D0 ]
      t_lattice_3 = [ 0.0D0, 0.0D0, 1.0D0 ]
      if (.not. include_far_field_if_periodic) then
        spatial_interaction_cutoff = huge(0.0D0) * [ 1.0D0, 1.0D0, 1.0D0 ]
      end if
      periodicity = [ .false., .true., .false. ]

      call pepc_prepare(2_kind_dim)

      if (time_pars%nresume > 0) then
        call read_particles_mpiio(time_pars%nresume, pepc_pars%pepc_comm%mpi_comm, &
          dummy_nresume, pepc_pars%np, p, dummy_file_name)

        if (dummy_nresume .ne. time_pars%nresume) then
          DEBUG_ERROR(*, "Resume timestep mismatch, parameter file says: ", time_pars%nresume, " checkpoint file says: ", dummy_nresume)
        end if

        if (pepc_pars%pepc_comm%mpi_rank == 0) then
          print *, "== [resume]"
          print *, "   resuming with ", pepc_pars%np, " particles."
          print *, ""
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
      character(len = 255), intent(in) :: file_name

      real(kind=8) :: m_ratio = 100.0D0
      real(kind=8) :: T_ratio = 1.0D0
      real(kind=8) :: B0 = 1.0D0
      real(kind=8), dimension(3) :: l_plasma = [ 50.0D0, 125.0D0, 0.0D0 ]
      real(kind=8) :: shear_halfwidth = 2.0D0
      real(kind=8) :: shear_strength = 1.0D0

      namelist /physics_nml/ m_ratio, T_ratio, B0, l_plasma, shear_halfwidth, shear_strength

      integer, parameter :: para_file_id = 10

      if (file_available) then
         open(para_file_id,file=trim(file_name),action='read')
         rewind(para_file_id)
         read(para_file_id, NML=physics_nml)
         close(para_file_id)
      end if

      physics_namelist%m_ratio = m_ratio
      physics_namelist%T_ratio = T_ratio
      physics_namelist%B0 = B0
      physics_namelist%l_plasma = l_plasma
      physics_namelist%shear_halfwidth = shear_halfwidth
      physics_namelist%shear_strength = shear_strength

   end subroutine read_in_physics_params


  subroutine write_physics_params(physics_pars, file_name)
    use encap
    implicit none

    type(physics_pars_t), intent(in) :: physics_pars
    character(len = 255), intent(in) :: file_name

    real(kind=8) :: m_ratio
    real(kind=8) :: T_ratio
    real(kind=8) :: shear_halfwidth
    real(kind=8) :: shear_strength
    real(kind=8) :: B0
    real(kind=8), dimension(3) :: l_plasma

    namelist /physics_nml/ m_ratio, T_ratio, B0, l_plasma, shear_halfwidth, shear_strength

    integer, parameter :: para_file_id = 10

    m_ratio  = physics_pars%mi / physics_pars%me
    T_ratio  = m_ratio * (physics_pars%vti / physics_pars%vte)**2
    B0       = physics_pars%B0
    l_plasma = physics_pars%l_plasma
    shear_halfwidth = physics_pars%shear_halfwidth
    shear_strength = physics_pars%shear_velocity * physics_pars%mi / &
      (physics_pars%qi * physics_pars%B0 * physics_pars%shear_halfwidth)

    open(para_file_id, file = trim(file_name), status = 'old', position = &
      'append', action = 'write')
    write(para_file_id, nml=physics_nml)
    close(para_file_id)

  end subroutine write_physics_params


   subroutine special_start(pepc_pars, physics_pars, p)
      use module_pepc_types, only: t_particle
      use module_mirror_boxes, only: constrain_periodic
      use module_debug
      use module_rng
      use encap
      use time_helper
      implicit none

      type(pepc_pars_t), intent(in) :: pepc_pars
      type(physics_pars_t), intent(in) :: physics_pars
      type(t_particle), allocatable, intent(inout) :: p(:)

      integer(kind_particle) :: np, npp, ni, ne, ip
      integer(kind_default) :: mpi_rank, mpi_size
      integer, allocatable :: period(:)
      integer :: it
      real(kind = 8) :: dx, lx, ly, qi, qe, mi, me, vti, vte, B0, rwce, rwci, v0
      real(kind = 8), allocatable :: xi(:,:), vi(:,:), xp(:)
      type(time_pars_t) :: time_pars

      real(kind=8), parameter :: pi = 3.1415926535897932384626434D0

      lx = physics_pars%l_plasma(1)
      ly = physics_pars%l_plasma(2)
      qi = physics_pars%qi
      qe = physics_pars%qe
      mi = physics_pars%mi
      me = physics_pars%me
      vti = physics_pars%vti
      vte = physics_pars%vte
      B0 = physics_pars%B0
      dx = physics_pars%shear_halfwidth
      v0 = physics_pars%shear_velocity

      rwce = abs(me / (qe * B0))
      rwci = abs(mi / (qi * B0))

      mpi_rank = pepc_pars%pepc_comm%mpi_rank
      mpi_size = pepc_pars%pepc_comm%mpi_size
      np = pepc_pars%np
      ni = np / (2 * mpi_size)
      if (mpi_rank < mod(np, int(mpi_size, kind=kind_particle))) then
        ni = ni + 1
      end if
      ne = ni
      npp = ni + ne

      if (mpi_rank == 0) then
        print *, "== [special_start]"
        print *, "   arranging ", np, " particles."
        print *, "   Lx x Ly = ", lx, " x ", ly
        print *, "   mi / me = ", mi / me
        print *, "   Ti / Te = ", (mi * vti**2) / (me * vte**2)
        print *, "   B0      = ", B0
        print *, ""
      end if

      allocate(p(1:npp), xi(ni,3), vi(ni,3), xp(ni))

      ! (1) assume a smooth profile for the electric polarization field, hence for the E x B drift velocity

      do ip = 1, ni
        ! (2) load the ions into the simulation domain with a uniform distribution along the x axis
        p(ip)%work = 1.0D0
        p(ip)%data%q = qi
        p(ip)%data%m = mi

        p(ip)%x(1) = (lx / ni) * ip - lx / (2.0D0 * ni)
        p(ip)%x(2) = ly * rng_next_real()
        p(ip)%x(3) = 0.0D0

        xi(ip,:) = p(ip)%x(:)
        xp(ip) = p(ip)%x(1)

        ! (3) give the ions random initial velocities in the positive x direction with a Rayleigh distribution at the desired 
        ! temperature
        p(ip)%data%v(1) = velocity_1d_rayleigh(vti)
        ! and in the y direction give them the local E x B drift velocity
        p(ip)%data%v(2) = vdrift(p(ip)%x(1))
        p(ip)%data%v(3) = 0.0D0

        vi(ip,:) = p(ip)%data%v(:)
      end do

      ! (4) advance each ion along its trajectory, in the presence of the constant but nonuniform electric field, for a random 
      ! fraction of its orbital period
      allocate(period(ni))
      period = 0

      time_pars%dt = 2.0D0 * pi * rwci / 160.0D0

      it = 0
      do
        do ip = 1, ni
          p(ip)%results%e(1) = estatic(p(ip)%x(1))
          p(ip)%results%e(2) = 0.0D0
        end do

        call push_particles(time_pars, physics_pars, p(1:ni))
        it = it + 1

        do ip = 1, ni
          if (period(ip) == 0 .and. xp(ip) < xi(ip, 1) .and. p(ip)%x(1) >= xi(ip, 1)) period(ip) = it
          xp(ip) = p(ip)%x(1)
        end do

        if (all(period /= 0)) then
          exit
        end if
      end do

      do ip = 1, ni
        p(ip)%x(:) = xi(ip,:)
        p(ip)%data%v(:) = vi(ip,:)

        do it = 1, int(rng_next_real() * period(ip))
          p(ip)%results%e(1) = estatic(p(ip)%x(1))
          p(ip)%results%e(2) = 0.0D0

          call push_particles(time_pars, physics_pars, p(ip:ip))
        end do
      end do

      deallocate(period, xi, vi, xp)

      call constrain_periodic(p(1:ni))
      do ip = 1, ni
        if ((p(ip)%x(1) .gt. lx) .or. (p(ip)%x(1) .lt. 0.0D0)) then
          p(ip)%data%v(1) = -p(ip)%data%v(1)
          p(ip)%data%v(2) = 2.0D0 * vdrift(p(ip)%x(1)) - p(ip)%data%v(2)

          p(ip)%x(1) = lx - modulo(p(ip)%x(1), lx)
        end if
      end do

      do ip = ni + 1, ni + ne
        p(ip)%work = 1.0D0
        p(ip)%data%q = qe
        p(ip)%data%m = me

        p(ip)%x(1) = lx * rng_next_real()
        p(ip)%x(2) = ly * rng_next_real()
        p(ip)%x(3) = 0.0D0

        p(ip)%data%v(1:2) = velocity_2d_maxwell(vte / sqrt(sqrt(1 + vdriftp(p(ip)%x(1)) * rwce)))
        p(ip)%data%v(2) = p(ip)%data%v(2) + vdrift(p(ip)%x(1))
        p(ip)%data%v(3) = 0.0D0
      end do

      contains

      function vdrift(x)
        implicit none

        real*8 :: vdrift
        real*8, intent(in) :: x

        vdrift = v0 * tanh((x - lx / 2.0D0) / dx)
      end function

      function vdriftp(x)
        implicit none

        real*8 :: vdriftp
        real*8, intent(in) :: x

        vdriftp = v0 * (1.0D0 - tanh((x - lx / 2.0D0) / dx)**2.0D0) / dx
      end function

      function estatic(x)
        implicit none

        real*8 :: estatic
        real*8, intent(in) :: x

        estatic = -v0 * B0 * tanh((x - lx / 2.0D0) / dx)
      end function

   end subroutine special_start


  function velocity_2d_maxwell(vt) result(v)
    use module_rng
    implicit none

    real*8, intent(in) :: vt
    real*8, dimension(2) :: v

    real*8 :: xi, v0
    real(kind=8), parameter :: pi = 3.1415926535897932384626434D0

    xi = rng_next_real()
    v0 = vt * sqrt(-2.0D0 * log(max(xi, 10.0_8**(-15))))
    xi = rng_next_real()
    v(1) = v0 * cos(2 * pi * xi)
    v(2) = v0 * sin(2 * pi * xi)
 
  end function velocity_2d_maxwell


  function velocity_1d_rayleigh(vt) result(v)
    use module_rng
    implicit none

    real*8, intent(in) :: vt
    real*8 :: v

    real*8 :: xi
    
    xi = rng_next_real()
    v = vt * sqrt(-2.0D0 * log(xi))
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
    real(kind = 8) :: e_kin, e_pot, e_kin_g, e_pot_g

    e_kin = 0.0D0
    e_pot = 0.0D0

    do ip = 1, size(p)
      e_kin = e_kin + e_kin_of_particle(p(ip)) 
      e_pot = e_pot + e_pot_of_particle(p(ip))
    end do

    call mpi_reduce(e_kin, e_kin_g, 1, MPI_REAL8, MPI_SUM, 0, &
      pepc_pars%pepc_comm%mpi_comm, mpi_err)
    call mpi_reduce(e_pot, e_pot_g, 1, MPI_REAL8, MPI_SUM, 0, &
      pepc_pars%pepc_comm%mpi_comm, mpi_err)

    if (pepc_pars%pepc_comm%mpi_rank == 0) then
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
      print *, ""

    end subroutine write_to_stdout

    subroutine write_to_file()
      implicit none

      if (step == 0) then
        open(file_energy, file = 'energy.csv', status = 'replace', action = 'write')
        write(file_energy, *) "t, T, Tc, V, H"
      else
        open(file_energy, file = 'energy.csv', status = 'old', position = 'append', action = 'write')
      end if

      if (.not. (time_pars%nresume > 0 .and. time_pars%nresume == step)) then
        write(file_energy, *) time_pars%dt * step, ", ", e_kin_g, ", ", &
          e_pot_g, ", ", (e_kin_g + e_pot_g)
      end if
      close(file_energy)

    end subroutine write_to_file

  end subroutine physics_dump


  pure function e_kin_of_particle(p)
    use module_pepc_types, only: t_particle
    implicit none

    type(t_particle), intent(in) :: p
    real(kind = 8) :: e_kin_of_particle

    e_kin_of_particle = 0.5D0 * p%data%m * dot_product(p%data%v, p%data%v)
  end function e_kin_of_particle


  pure function e_pot_of_particle(p)
    use module_pepc_types, only: t_particle
    implicit none

    type(t_particle), intent(in) :: p
    real(kind = 8) :: e_pot_of_particle

    e_pot_of_particle = 0.5 * p%results%pot * p%data%q
  end function e_pot_of_particle


end module physics_helper
