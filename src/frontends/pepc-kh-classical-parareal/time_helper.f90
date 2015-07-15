module time_helper
  use module_pepc_kinds
  implicit none

contains

  function time_of_step(step, time_p)
    use encap
    implicit none
    integer, intent(in) :: step
    type(time_pars_t), intent(in) :: time_p

    real(kind = 8) :: time_of_step

    time_of_step = time_p%tresume + time_p%dt * (step - time_p%nresume)
  end function


  subroutine setup_time(file_name, pepc_comm, time_pars)
    use encap

    character(*), intent(in) :: file_name
    type(pepc_comm_t), intent(in) :: pepc_comm
    type(time_pars_t), intent(out) :: time_pars

    call read_in_time_params(file_name, time_pars)

    if (pepc_comm%mpi_rank == 0) then
        print *, "== [setup_time]"
        print *, "   dt      = ", time_pars%dt
        print *, "   tresume = ", time_pars%tresume
        print *, "   nsteps  = ", time_pars%nsteps
        print *, "   nresume = ", time_pars%nresume
        print *, ""
    end if
  end subroutine setup_time


  subroutine read_in_time_params(file_name, time_pars)
    use encap
    implicit none
    include 'mpif.h'

    character(*), intent(in) :: file_name
    type(time_pars_t), intent(out) :: time_pars

    real(kind=8) :: dt = 0.
    real(kind=8) :: tresume = 0.
    integer :: nsteps = 0
    integer :: nresume = 0

    namelist /time_nml/ dt, tresume, nsteps, nresume

    integer, parameter :: param_file_id = 10

    open(param_file_id,file=trim(file_name),action='read')
    rewind(param_file_id)
    read(param_file_id, NML=time_nml)
    close(param_file_id)

    time_pars%dt = dt
    time_pars%tresume = tresume
    time_pars%nsteps = nsteps
    time_pars%nresume = nresume
  end subroutine read_in_time_params


  subroutine write_time_params(file_name, time_pars, step)
    use encap
    implicit none

    character(*), intent(in) :: file_name
    type(time_pars_t), intent(in) :: time_pars
    integer, intent(in) :: step

    real(kind=8) :: dt = 0.
    real(kind=8) :: tresume = 0.
    integer :: nsteps = 0
    integer :: nresume = 0

    namelist /time_nml/ dt, tresume, nsteps, nresume

    integer, parameter :: param_file_id = 10

    dt = time_pars%dt
    tresume = time_of_step(step, time_pars)
    nsteps = time_pars%nsteps
    nresume = step

    open(param_file_id, file = trim(file_name), status = 'old', position =&
      'append', action = 'write')
    write(param_file_id, NML=time_nml)
    close(param_file_id)
  end subroutine write_time_params


  subroutine push_particles(time_pars, physics_pars, p)
    use module_pepc_types
    use encap
    implicit none

    type(time_pars_t), intent(in) :: time_pars
    type(physics_pars_t), intent(in) :: physics_pars
    type(t_particle), intent(inout) :: p(:)

    integer(kind_particle) :: ip
    real*8 :: beta, gam
    real*8, dimension(3) :: uminus, uprime, uplus, t, s

    real*8, dimension(3) :: B0

    B0(1:2) = 0.0D0
    B0(3) = physics_pars%B0

    do ip = 1, size(p)
      ! charge/mass*time-constant
      beta   = p(ip)%data%q / (2. * p(ip)%data%m) * time_pars%dt
      ! first half step with electric field
      uminus(1:2) = p(ip)%data%v(1:2) + beta * p(ip)%results%e(1:2)
      uminus(3)   = 0
      ! gamma factor
      !gam    = sqrt( 1.0 + ( dot_product(uminus, uminus) ) / unit_c2 )
      gam    = 1.0D0
      ! rotation with magnetic field
      t      = beta/gam * B0
      uprime = uminus + cross_product(uminus, t)
      s      = 2. * t / (1 + dot_product(t, t))
      uplus  = uminus + cross_product(uprime, s)
      ! second half step with electric field
      p(ip)%data%v(1:2) = uplus(1:2) + beta * p(ip)%results%e(1:2)
      p(ip)%data%v(3)   = 0

      ! gam = sqrt(1.0D0 + dot_product(p(ip)%data%v * p(ip)%data%v) / unit_c2)
      gam = 1.0D0
      p(ip)%x = p(ip)%x + p(ip)%data%v / gam * time_pars%dt
    end do

    contains

    pure function cross_product(a, b)
      implicit none

      real*8, dimension(3), intent(in) :: a, b
      real*8, dimension(3) :: cross_product

      cross_product(1) = a(2) * b(3) - a(3) * b(2)
      cross_product(2) = a(3) * b(1) - a(1) * b(3)
      cross_product(3) = a(1) * b(2) - b(2) * a(1)
    end function cross_product
  end subroutine push_particles


  subroutine constrain_particles(physics_pars, p)
    use module_pepc_types
    use module_mirror_boxes

    use encap
    implicit none

    type(physics_pars_t), intent(in) :: physics_pars
    type(t_particle), intent(inout) :: p(:)

    integer(kind_particle) :: ip
    real(kind=8) :: vte, vti, lx

    vte = physics_pars%vte
    vti = physics_pars%vti
    lx  = physics_pars%l_plasma(1)

    do ip = 1, size(p)
      if ((p(ip)%x(1) .gt. lx) .or. (p(ip)%x(1) .lt. 0.0D0)) then

        p(ip)%x(1) = lx - modulo(p(ip)%x(1), lx)
        p(ip)%data%v(1) = -p(ip)%data%v(1)

      end if
    end do

    call constrain_periodic(p)
  end subroutine constrain_particles
end module time_helper
