! =============================================
!
!                CONFIGURE
!
!  Sets up physical system: particle positions, velocities
!
! ==============================================

subroutine configure

  use physvars
  !  use utils
  implicit none
  include 'mpif.h'

  integer :: i, ipe, idummy=0, ierr, ifile,mac_init
  real :: t_walk, t_walkc, t_force
  real :: t0, t1, t_key, t_domain=0., t_build=0., t_branches=0., t_fill=0., t_properties=0., t_restore=0.
  real :: t_begin=0.
  real :: t_stuff_1=0., t_stuff_2=0., t_en, t_mpi=0., t_end=0.
  real :: t_push, t_diag, t_start_push, t_prefetch=0., Tpon, ttot, t_laser, t_all = 0.
  integer :: npnew

  integer, allocatable :: indxl(:),irnkl(:), label_tmp(:)
  integer, allocatable :: islen(:),irlen(:)
  integer, allocatable :: fposts(:),gposts(:)

  npnew = np_local

  config: select case(system_config)

  case(1)              ! Set up particles according to geometry
     call randion         


     if (vte > 0) then
        call maxwell1(ux,nppm,1,nep,vte)
        call maxwell1(uy,nppm,1,nep,vte)
        call maxwell1(uz,nppm,1,nep,vte)
        call scramble_v(1,nep)   ! remove x,y,z correlations
     else
        call cold_start(1,nep)
     endif

     if (vti > 0) then
        call maxwell1(ux,nppm,nep+1,nip,vti)
        call maxwell1(uy,nppm,nep+1,nip,vti)
        call maxwell1(uz,nppm,nep+1,nip,vti)
        call scramble_v(nep+1,nip) ! remove x,y,z correlations
     else
        call cold_start(nep+1,nip)
     endif

  case(2)
     call special_start(ispecial)

  case default     ! Default = 0 - no plasma target
     if (my_rank==0) write (6,*) 'Warning: no particles set up'
     npart_total=0
!     npp = 0
  end select config

 ! Compute initial field values - need these to get vec. pots consistent with velocities
  mac_init=0  ! Use standard s/d mac for first step
  if (my_rank==0) write(*,*) 'Computing initial fields'
  allocate(indxl(nppm),irnkl(nppm),islen(n_cpu),irlen(n_cpu),fposts(n_cpu+1),gposts(n_cpu+1),label_tmp(nppm))
  
  call pepc_fields_p1(np_local,x(1:np_local),y(1:np_local),z(1:np_local), &
       q(1:np_local),m(1:np_local),work(1:np_local),pelabel(1:np_local), &
       ex,ey,ez,pot,npnew,indxl,irnkl,islen,irlen,fposts,gposts,0, 0.3, eps, force_const, err_f, xl, yl, zl, itime, &
       t_begin,t_domain,t_build,t_branches,t_fill,t_properties,t_prefetch, &
       t_stuff_1,t_stuff_2,t_walk,t_walkc,t_force,t_restore,t_mpi,t_end,t_all)
  
  

end subroutine configure

