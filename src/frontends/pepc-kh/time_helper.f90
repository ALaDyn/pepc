module time_helper
   use module_pepc_kinds
   implicit none

   type time_nml_t
      real(kind_physics) :: te
      integer :: nsteps, nresume
   end type time_nml_t

contains

   subroutine setup_time(time_pars, pepc_comm)
      use pepc_helper, only: para_file_available, para_file_name
      use encap

      type(time_pars_t), intent(out) :: time_pars
      type(pepc_comm_t), intent(in) :: pepc_comm

      type(time_nml_t) :: time_nml

      call read_in_time_params(time_nml, para_file_available, para_file_name)

      time_pars%te = time_nml%te
      time_pars%nsteps = time_nml%nsteps
      time_pars%nresume = time_nml%nresume
      time_pars%dt = time_pars%te / time_pars%nsteps

      if (pepc_comm%mpi_rank .eq. 0) then
         print *, "== [setup_time]"
         print *, "   te      = ", time_pars%te
         print *, "   nsteps  = ", time_pars%nsteps
         print *, "   nresume = ", time_pars%nresume
         print *, "   dt      = ", time_pars%dt
         print *, ""
      end if

   end subroutine setup_time

   subroutine read_in_time_params(time_namelist, file_available, file_name)
      use mpi
      implicit none

      type(time_nml_t), intent(out) :: time_namelist
      logical, intent(in) :: file_available
      character(len=255), intent(in) :: file_name

      real(kind_physics) :: te = 0.
      integer :: nsteps = 0
      integer :: nresume = 0

      namelist /time_nml/ te, nsteps, nresume

      integer, parameter :: para_file_id = 10

      if (file_available) then
         open (para_file_id, file=trim(file_name), action='read')
         rewind (para_file_id)
         read (para_file_id, NML=time_nml)
         close (para_file_id)
      end if

      time_namelist%te = te
      time_namelist%nsteps = nsteps
      time_namelist%nresume = nresume

   end subroutine read_in_time_params

   subroutine write_time_params(time_pars, step, file_name)
      use encap
      implicit none

      type(time_pars_t), intent(in) :: time_pars
      integer, intent(in) :: step
      character(len=255), intent(in) :: file_name

      real(kind_physics) :: te = 0.
      integer :: nsteps = 0
      integer :: nresume = 0

      namelist /time_nml/ te, nsteps, nresume

      integer, parameter :: para_file_id = 10

      te = time_pars%te
      nsteps = time_pars%nsteps
      nresume = step

      open (para_file_id, file=trim(file_name), status='old', position= &
            'append', action='write')
      write (para_file_id, NML=time_nml)
      close (para_file_id)

   end subroutine write_time_params

   subroutine push_particles(time_pars, physics_pars, p)
      use module_pepc_types
      use encap
      implicit none

      type(time_pars_t), intent(in) :: time_pars
      type(physics_pars_t), intent(in) :: physics_pars
      type(t_particle), intent(inout) :: p(:)

      integer(kind_particle) :: ip
      real(kind_physics) :: beta, gam
      real(kind_physics), dimension(3) :: uminus, uprime, uplus, t, s

      real(kind_physics), dimension(3) :: B0

      B0(1:2) = 0.
      B0(3) = physics_pars%B0

      do ip = 1, size(p)
         ! charge/mass*time-constant
         beta = p(ip)%data%q / (2.*p(ip)%data%m) * time_pars%dt
         ! first half step with electric field
         uminus(1:2) = p(ip)%data%v(1:2) + beta * p(ip)%results%e(1:2)
         uminus(3) = 0
         ! gamma factor
         !gam    = sqrt( 1.0 + ( dot_product(uminus, uminus) ) / unit_c2 )
         gam = 1.
         ! rotation with magnetic field
         t = beta / gam * B0
         uprime = uminus + cross_product(uminus, t)
         s = 2.*t / (1 + dot_product(t, t))
         uplus = uminus + cross_product(uprime, s)
         ! second half step with electric field
         p(ip)%data%v(1:2) = uplus(1:2) + beta * p(ip)%results%e(1:2)
         p(ip)%data%v(3) = 0

         ! gam = sqrt(1. + dot_product(p(ip)%data%v * p(ip)%data%v) / unit_c2)
         gam = 1.
         p(ip)%x = p(ip)%x + p(ip)%data%v / gam * time_pars%dt
      end do

   contains

      pure function cross_product(a, b)
         implicit none

         real(kind_physics), dimension(3), intent(in) :: a, b
         real(kind_physics), dimension(3) :: cross_product

         cross_product(1) = a(2) * b(3) - a(3) * b(2)
         cross_product(2) = a(3) * b(1) - a(1) * b(3)
         cross_product(3) = a(1) * b(2) - b(2) * a(1)
      end function cross_product

   end subroutine push_particles

   subroutine constrain_particles(physics_pars, p)
      use module_pepc_types
      use module_mirror_boxes

      use module_rng
      use encap
      use pepc_helper
      implicit none

      type(physics_pars_t), intent(in) :: physics_pars
      type(t_particle), intent(inout) :: p(:)

      integer(kind_particle) :: ip
      real(kind_physics) :: vte, vti, lx

      vte = physics_pars%vte
      vti = physics_pars%vti
      lx = physics_pars%l_plasma(1)

      do ip = 1, size(p)
         if (p(ip)%x(1) .gt. lx) then

            p(ip)%x(1) = lx - modulo(p(ip)%x(1), lx)
            p(ip)%data%v(1) = -p(ip)%data%v(1)

         end if
      end do

      call constrain_periodic(p)

   end subroutine constrain_particles

end module time_helper
