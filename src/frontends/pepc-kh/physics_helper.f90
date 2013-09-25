module physics_helper

   implicit none

   type physics_nml_t
      real(kind=8) :: m_ratio, T_ratio, B0
      real(kind=8), dimension(3) :: l_plasma
   end type physics_nml_t

   real(kind = 8), parameter :: force_const = 0.159154943D0
   real(kind = 8) :: e_constraint

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

      e_constraint = 0.0D0

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

      namelist /physics_nml/ m_ratio, T_ratio, B0, l_plasma

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

   end subroutine read_in_physics_params


  subroutine write_physics_params(physics_pars, file_name)
    use encap
    implicit none

    type(physics_pars_t), intent(in) :: physics_pars
    character(len = 255), intent(in) :: file_name

    real(kind=8) :: m_ratio
    real(kind=8) :: T_ratio
    real(kind=8) :: B0
    real(kind=8), dimension(3) :: l_plasma

    namelist /physics_nml/ m_ratio, T_ratio, B0, l_plasma

    integer, parameter :: para_file_id = 10

    m_ratio  = physics_pars%mi / physics_pars%me
    T_ratio  = m_ratio * (physics_pars%vti / physics_pars%vte)**2
    B0       = physics_pars%B0
    l_plasma = physics_pars%l_plasma

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
      implicit none

      type(pepc_pars_t), intent(in) :: pepc_pars
      type(physics_pars_t), intent(in) :: physics_pars
      type(t_particle), allocatable, intent(inout) :: p(:)

      integer(kind_particle) :: nx, ny, ipl, ipg, ix, iy, np, npp
      integer(kind_default) :: mpi_rank, mpi_size
      real(kind = 8) :: dx, dy, lx, ly, qi, qe, mi, me, vti, vte, xgc, ygc, B0, rwce, rwci

      lx = physics_pars%l_plasma(1)
      ly = physics_pars%l_plasma(2)
      qi = physics_pars%qi
      qe = physics_pars%qe
      me = physics_pars%me
      mi = physics_pars%mi
      vte = physics_pars%vte
      vti = physics_pars%vti
      B0 = physics_pars%B0

      rwce = abs(me / (qe * B0))
      rwci = abs(mi / (qi * B0))

      ! place ions on a regular grid, np = 2 x ni = 2 x nx x ny
      ! nx / ny = Lx / Ly
      mpi_rank = pepc_pars%pepc_comm%mpi_rank
      mpi_size = pepc_pars%pepc_comm%mpi_size
      np = pepc_pars%np
      npp = np / mpi_size
      if (mpi_rank < mod(np, int(mpi_size, kind=kind_particle))) then
        npp = npp + 1
      end if
      nx = nint(sqrt(lx * np / (2 * ly)))
      ny = np / (2 * nx)

      if (2 * nx * ny /= np) then
        DEBUG_ERROR(*, "No. of particles np = ", np, " cannot be partitioned into nx x ny with nx / ny = Lx / Ly. nx = ", nx, " ny = ", ny, " Lx = ", lx, " Ly = ", ly)
      end if

      if (mpi_rank == 0) then
        print *, "== [special_start]"
        print *, "   arranging ", np, " particles on a grid."
        print *, "   nx x ny = ", nx, " x ", ny, " ion/electron pairs."
        print *, "   Lx x Ly = ", lx, " x ", ly
        print *, "   mi / me = ", mi / me
        print *, "   Ti / Te = ", (mi * vti**2) / (me * vte**2)
        print *, "   B0      = ", physics_pars%B0
        print *, ""
      end if

      dx = lx / nx
      dy = ly / ny

      allocate(p(1:npp))

      do ipl = 1, npp
        ipg = ipl + min(int(mpi_rank, kind=kind_particle), mod(np, int(mpi_size, kind=kind_particle))) * (np / mpi_size + 1) + &
          max(0_kind_particle, mpi_rank - mod(np, int(mpi_size, kind=kind_particle))) * (np / mpi_size)

        ! ipg = 1, 3, ... are ions, ipg = 2, 4, ... are electrons, distribute on the same lattice
        ix = mod((ipg - 1) / 2, nx) + 1
        iy = (ipg - 1) / (2 * nx) + 1

        xgc = (ix - 0.5D0) * dx
        ygc = (iy - 0.5D0) * dy

        p(ipl)%work = 1.0D0
        p(ipl)%label = ipg

        if (mod(ipg, 2_kind_particle) == 0) then ! this is an electron
          p(ipl)%data%q = qe
          p(ipl)%data%m = me

          p(ipl)%data%v = velocity_2d_maxwell(vte)

          p(ipl)%x(1) = xgc + p(ipl)%data%v(2) * rwce
          p(ipl)%x(2) = ygc - p(ipl)%data%v(1) * rwce
          p(ipl)%x(3) = 0.0D0

        else ! this is an ion
          p(ipl)%data%q = qi
          p(ipl)%data%m = mi

          p(ipl)%data%v = velocity_2d_maxwell(vti)

          p(ipl)%x(1) = xgc - p(ipl)%data%v(2) * rwci
          p(ipl)%x(2) = ygc + p(ipl)%data%v(1) * rwci
          p(ipl)%x(3) = 0.0D0

          if (p(ipl)%x(1) < 0) then
            p(ipl)%x(1) = modulo(-p(ipl)%x(1), lx)
            p(ipl)%data%v(2) = -p(ipl)%data%v(2)
          end if

        end if

        if (p(ipl)%x(1) > lx) then
          p(ipl)%x(1) = lx - modulo(p(ipl)%x(1), lx)
          p(ipl)%data%v(2) = -p(ipl)%data%v(2)
        end if

      end do

      call constrain_periodic(p)

   end subroutine special_start


  function velocity_2d_maxwell(vt) result(v)
    use module_rng
    implicit none

    real*8, intent(in) :: vt
    real*8, dimension(3) :: v

    real*8 :: xi, v0
    real(kind=8), parameter :: pi = 3.1415926535897932384626434D0

    xi = rng_next_real()
    v0 = vt * sqrt(-2.0D0 * log(max(xi, 10.0_8**(-15))))
    xi = rng_next_real()
    v(1) = v0 * cos(2 * pi * xi)
    v(2) = v0 * sin(2 * pi * xi)
    v(3) = 0.0D0
 
  end function velocity_2d_maxwell


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
    real(kind = 8) :: e_kin, e_pot, e_kin_g, e_pot_g, e_constraint_g

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
    call mpi_reduce(e_constraint, e_constraint_g, 1, MPI_REAL8, MPI_SUM, 0, &
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
      print *, "   Tc: ", e_constraint_g
      print *, "   V : ", e_pot_g
      print *, "   H : ", e_kin_g + e_constraint + e_pot_g
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
          e_constraint_g, ", ", e_pot_g, ", ", (e_kin_g + e_constraint_g + e_pot_g)
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
