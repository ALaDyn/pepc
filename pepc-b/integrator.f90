
subroutine integrator

  use treevars
  use physvars
  implicit none
  integer :: i
  ! Integrator

 if (me==0) then
	write(6,*) 'Push with scheme',scheme
	write(24,*) 'Push with scheme',scheme
 endif

  pusher: select case(scheme)

  case(1,2,3,4,5)
     call velocities(1,npp,dt)  ! pure ES, NVT ensembles
     call push_x(1,npp,dt)  ! update positions

  case(6) 
     call push_full3v(1,npp,dt)  ! full EM pusher (all E, B components)
     call push_x(1,npp,dt)  ! update positions

  case(7)
     call velocities(1,npp,dt)  ! nonrelativistic push
     call push_nonrel(1,npp,dt)
  case default
     ! do nothing!

  end select pusher



  boundaries: select case(particle_bcs)

  case(2)
     call constrain   ! relective particle bcs for temperature-clamp mode
  case(3)
     call earth_plate  ! special bcs for grounded target end 
  case default
     ! do nothing

  end select boundaries


end subroutine integrator
