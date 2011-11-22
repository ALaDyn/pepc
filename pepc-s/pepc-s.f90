
! ==============================================================
!
!
!                  PEPC-S
!
!    Parallel Efficient Parallel Coulomb-solver: Single Call Version
!
!  ==============================================================

subroutine pepc(nparts, npart_tot, pos_x, pos_y, pos_z, charge, mass, Ex, Ey, Ez, pot, lat_x, lat_y, lat_z, lat_period, lat_corr)

  use treetypes
  use physvars
  use module_fmm_framework
  use module_pepcfields
  implicit none
  include 'mpif.h'

  integer, intent(in) :: nparts, npart_tot
  real*8, intent(in), dimension(nparts) :: pos_x, pos_y, pos_z
  real*8, intent(in), dimension(nparts) :: charge, mass
  real*8, intent(out), dimension(nparts) :: Ex, Ey, Ez, pot

  real*8, intent(in), dimension(3) :: lat_x, lat_y, lat_z
  logical, intent(in), dimension(3) :: lat_period
  logical, intent(in) :: lat_corr
  type(t_calc_force_params) ::cf_par

  integer :: ip
  integer :: ierr

  t_lattice_1 = lat_x
  t_lattice_2 = lat_y
  t_lattice_3 = lat_z

  periodicity = lat_period 
  do_extrinsic_correction = lat_corr

  np_local = nparts
  npart_total = npart_tot

  ! Get the id number of the current task
  call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)

  ! Get the number of MPI tasks
  call MPI_COMM_size(MPI_COMM_WORLD, n_cpu, ierr)

  ! Each CPU gets copy of initial data
  call pepc_setup()

  ! Allocate array space for tree
  call libpepc_setup(my_rank,n_cpu,db_level)

  do ip=1, nparts
     pelabel(ip) = ip
  end do

  itime = 1

  ! initialize calc force params
  cf_par%theta       = theta
  cf_par%mac         = mac
  cf_par%eps         = eps
  cf_par%force_const = force_const
  cf_par%force_law   = 3


  call pepc_fields_coulomb_wrapper(np_local, npart_total, &
       pos_x, pos_y, pos_z, &
       charge, work, pelabel, &
       ex, ey, ez, pot, &
       np_mult, cf_par, &
       itime, weighted, curve_type, &
       num_neighbour_boxes, neighbour_boxes, .false.)
  
  ! cleanup of lpepc static data
  call libpepc_finalize()

  ! deallocate array space for particles
  call pepc_cleanup(my_rank,n_cpu)

end subroutine pepc



