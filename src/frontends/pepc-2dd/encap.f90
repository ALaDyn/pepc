module mpi
  implicit none
  include 'mpif.h'
end module



module encap
   use module_pepc_kinds
   use module_pepc_types
   use iso_c_binding
   implicit none

    ! variables for MPI within pepc
   type, bind(c) :: pepc_comm_t
      integer(c_int)        :: mpi_size, mpi_rank
      integer(kind_default) :: mpi_comm
   end type pepc_comm_t

   type, bind(c) :: pepc_pars_t
      integer(c_int)        :: pdump, fdump, cdump
      integer(c_int64_t)    :: np
      type(pepc_comm_t)     :: pepc_comm
   end type pepc_pars_t

  type :: physics_pars_t
    real(kind_physics)               :: vte, vti, qe, qi, me, mi, shear_halfwidth, shear_velocity
    real(kind_physics), dimension(3) :: l_plasma
    integer(kind = kind_particle)    :: ni
  end type physics_pars_t

   type :: time_pars_t
      real(kind_physics) :: te, dt
      integer            :: nsteps, nresume
   end type time_pars_t

  type field_grid_t
    integer(kind = kind_particle), dimension(3)       :: n
    integer(kind = kind_particle)                     :: ntot, nl
    real(kind_physics), dimension(3)                  :: offset, extent, dx
    type(t_particle), dimension(:), allocatable       :: p
    real(kind_physics), dimension(:,:,:), allocatable :: ne, ni,nb, qe, qi,qb, vex, vey, vez, vix, viy, viz, vbx, vby, vbz!,&
!                                                         Jex, Jey, Jez, Jix, Jiy, Jiz, Jbx, Jby, Jbz                      ,&
!                                                         Aex, Aey, Aez, Aix, Aiy, Aiz, Abx, Aby, Abz                      ,&
!                                                         Bex, Bey, Bez, Bix, Biy, Biz, Bbx, Bby, Bbz                      ,&
!                                                         Eex, Eey, Eez, Eix, Eiy, Eiz, Ebx, Eby, Ebz      
  end type field_grid_t

end module encap
