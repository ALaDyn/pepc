subroutine closefiles
  use physvars

  if (my_rank == 0) then
     close(15)
     close(81)  ! particle dump 
     close(70)
  endif
  close(20)
  close(80)  ! initial particle data


end subroutine closefiles
