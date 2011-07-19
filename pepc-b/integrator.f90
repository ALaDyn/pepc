
subroutine integrator

  use module_physvars
  use module_integration_scheme
  use module_particle_boundaries
  implicit none

  ! Integrator

 if (my_rank==0) then
	write(6,*) 'Push with scheme',scheme
	write(24,*) 'Push with scheme',scheme
 endif

  pusher: select case(scheme)

  case(1,2,3,4,5)
     call velocities(1,np_local,dt)  ! pure ES, NVT ensembles
     call push_x(1,np_local,dt)  ! update positions

  case(6) 
     call push_full3v(1,np_local,dt)  ! full EM pusher (all E, B components)
     call push_x(1,np_local,dt)  ! update positions

  case(7)
     call velocities(1,np_local,dt)  ! nonrelativistic push
     call push_nonrel(1,np_local,dt)
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
