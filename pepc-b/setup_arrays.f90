subroutine setup_arrays(init_mb)
  use module_physvars
  implicit none
  include 'mpif.h'
  integer :: mem_fields, npmax, ntm
  integer, intent(out) :: init_mb

  npmax = npart_total/n_cpu*3

  ! array allocation  - TODO
  init_mb = 15*npmax*8

  allocate ( xslice(npmax), yslice(npmax), zslice(npmax), uxslice(npmax), uyslice(npmax), uzslice(npmax), & 
       qslice(npmax), mslice(npmax) )    ! Reserve slice particle array space N/NPE


  mem_fields = npmax*(8*8)
  ntm = nt

  allocate (rhoi(0:ngx+1,0:ngy+1,0:ngz+1),rhoe(0:ngx+1,0:ngy+1,0:ngz+1))  ! global field arrays

! local field arrays for cycle-averages
  allocate (rhoe_loc(0:ngx+1,0:ngy+1,0:ngz+1), rhoi_loc(0:ngx+1,0:ngy+1,0:ngz+1), &
      ex_loc(0:ngx+1,0:ngy+1,0:ngz+1), ey_loc(0:ngx+1,0:ngy+1,0:ngz+1),  ez_loc(0:ngx+1,0:ngy+1,0:ngz+1), &
      bx_loc(0:ngx+1,0:ngy+1,0:ngz+1), by_loc(0:ngx+1,0:ngy+1,0:ngz+1),  bz_loc(0:ngx+1,0:ngy+1,0:ngz+1), &
      Te_loc(0:ngx+1,0:ngy+1,0:ngz+1), Ti_loc(0:ngx+1,0:ngy+1,0:ngz+1),  &
      g_ele(0:ngx+1,0:ngy+1,0:ngz+1), g_ion(0:ngx+1,0:ngy+1,0:ngz+1),  &
      jxe_loc(0:ngx+1,0:ngy+1,0:ngz+1), jye_loc(0:ngx+1,0:ngy+1,0:ngz+1), jze_loc(0:ngx+1,0:ngy+1,0:ngz+1) )
  allocate (ex_ave(0:ngav+1),ey_ave(0:ngav+1,3),ez_ave(0:ngav+1) )

 allocate (rhoe2d_loc(0:ngx+1,0:ngy+1), rhoi2d_loc(0:ngx+1,0:ngy+1) )

! Fields for Helmholtz solver
  allocate (rho_helm(0:nxh+1),az_helm(0:nxh+1) )


  mem_fields = mem_fields + ngx*ngy*ngz * (8*13)


! particle buffers for visualisation

  allocate (vbuffer(0:attrib_max-1,nbuf_max), vbuf_local(0:attrib_max-1,nbuf_max))

  if (my_rank==0) then
     write(*,'(//a/)') 'Initial memory allocation:'
     write(*,'(1(a15,f12.3,a3/)/)') 'Fields:',mem_fields/1.e6,' MB'
  endif

  ex_ave=0.
  ey_ave=0.
  ez_ave=0.

end subroutine setup_arrays






