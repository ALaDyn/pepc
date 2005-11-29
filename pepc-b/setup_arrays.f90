subroutine setup_arrays
  use physvars
  implicit none
  include 'mpif.h'
  integer :: mem_fields, npmax

  npmax = npart_total/n_cpu*3
  allocate ( xslice(npmax), yslice(npmax), zslice(npmax), uxslice(npmax), uyslice(npmax), uzslice(npmax), & 
       qslice(npmax), mslice(npmax) )    ! Reserve slice particle array space N/NPE

  mem_fields = npmax*(8*8)

  allocate (rhoi(0:ngx+1,0:ngy+1,0:ngz+1),rhoe(0:ngx+1,0:ngy+1,0:ngz+1))  ! global field arrays

! local field arrays for cycle-averages
  allocate (rhoe_loc(0:ngx+1,0:ngy+1,0:ngz+1), rhoi_loc(0:ngx+1,0:ngy+1,0:ngz+1), &
      ex_loc(0:ngx+1,0:ngy+1,0:ngz+1), ey_loc(0:ngx+1,0:ngy+1,0:ngz+1),  ez_loc(0:ngx+1,0:ngy+1,0:ngz+1), &
      bx_loc(0:ngx+1,0:ngy+1,0:ngz+1), by_loc(0:ngx+1,0:ngy+1,0:ngz+1),  bz_loc(0:ngx+1,0:ngy+1,0:ngz+1), &
      Te_loc(0:ngx+1,0:ngy+1,0:ngz+1), Ti_loc(0:ngx+1,0:ngy+1,0:ngz+1),  &
      g_ele(0:ngx+1,0:ngy+1,0:ngz+1), g_ion(0:ngx+1,0:ngy+1,0:ngz+1),  &
      jxe_loc(0:ngx+1,0:ngy+1,0:ngz+1), jye_loc(0:ngx+1,0:ngy+1,0:ngz+1), jze_loc(0:ngx+1,0:ngy+1,0:ngz+1) )
  mem_fields = mem_fields + ngx*ngy*ngz * (8*13)


! particle buffers for visualisation

  allocate (vbuffer(0:attrib_max-1,nbuf_max), vbuf_local(0:attrib_max-1,nbuf_max))

  if (my_rank==0) then
     write(*,'(//a/)') 'Initial memory allocation:'
     write(*,'(1(a15,f12.3,a3/)/)') 'Fields:',mem_fields/1.e6,' MB'
  endif

end subroutine setup_arrays






