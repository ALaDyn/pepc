subroutine closefiles
  use physvars
  implicit none

  if (my_rank == 0) then
     close(15)
     close(81)  ! particle dump 
     close(70)
     close(75)
  endif
  if (debug_level>1) close(20)
  close(80)  ! initial particle data


end subroutine closefiles
