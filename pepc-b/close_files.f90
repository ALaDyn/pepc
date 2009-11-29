subroutine closefiles
  use physvars
  implicit none

  if (my_rank == 0) then
     close(15)
     close(24)
     close(70)
     close(71)
     close(75)
     close(90)
  endif
  if (debug_level>2) close(20)
  close(80)  ! initial particle data


end subroutine closefiles
