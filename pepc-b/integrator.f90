
subroutine integrator

  use treevars
  use physvars
  implicit none

  ! Integrator
  if (scheme == 6 ) then
     call push_full3v(1,npp,dt)  ! full EM pusher (all E, B components)
  else
     call velocities(1,npp,dt)  ! pure ES, NVT ensembles
  endif

  call push_x(1,npp,dt)  ! update positions


  boundaries: select case(particle_bcs)

  case(2)
     call constrain   ! relective particle bcs for temperature-clamp mode
  case(3)
     call earth_plate  ! special bcs for grounded target end 
  case default
     ! do nothing

  end select boundaries

end subroutine integrator
