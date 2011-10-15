
subroutine integrator

  use module_physvars
  use module_integration_scheme
  use module_particle_boundaries
  implicit none

  ! Integration schemes (selected via scheme)
!      1 = NVE - total energy conserved                                               
!      2 = NVT - global Te, Ti cons                                                      
!      3 = global NVT electron Te conserved; ions frozen                                     
!      4 = local NVT: each PE keeps Te clamped; ions frozen                                 
!      5 = local NVT, ions only; electrons added at end of run                             
!      6 = Full EM pusher (all E, B components), total energy conserved     
!      7 = ES, Nonrelativistic push, total energy conserved  
!      8 = 2v EM push, TE fields (Ex, Ey, Bz)

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

  case(8) 
      call push_TE(1,np_local,dt)  ! full EM pusher (all E, B components)
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
