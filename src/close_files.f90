subroutine closefiles
  use treevars

  if (me == 0) then
     close(15)
     close(81)  ! particle dump 
     close(70)
     close(75)
  endif
  close(20)
  close(80)  ! initial particle data


end subroutine closefiles
