! ===========================================
!
!           CHECK_TABLE
!
!  Do some quick checks on the tree structure 

!   
!
! ===========================================

subroutine check_table(callpoint)

  use treevars
  use utils

  implicit none
  character*20 :: callpoint
  integer :: nleaf_check, ntwig_check, nleaf_me_check, ntwig_me_check


  ntwig_check = count(mask =  htable%node <0 )
  nleaf_check = count(mask =  htable%node >0 )
  nleaf_me_check = count(mask = htable%owner==me .and. htable%node >0 )
  ntwig_me_check = count(mask = htable%owner==me .and. htable%node <0 )


  if (nleaf /= nleaf_check) then
     write(*,'(3a,i4)') 'Table check called ',callpoint,' by PE',me
     write(*,*) '# leaves in table = ',nleaf_check,'vs ',nleaf,'accumulated'
     write(*,*) 'Fixing and continuing for now..'
!     nleaf = nleaf_check
  endif
  if (ntwig /= ntwig_check) then
     write(*,'(3a,i4)') 'Table check called ',callpoint,' by PE',me
     write(*,*) ' # twigs in table = ',ntwig_check,'vs ',ntwig,'accumulated'
     write(*,*) 'Fixing and continuing for now..'
!     ntwig = ntwig_check
  endif

  if (nleaf_me /= nleaf_me_check) then
     write(*,'(3a,i4)') 'Table check called ',callpoint,' by PE',me
     write(*,*) ' # own leaves in table = ',nleaf_me_check,'vs ',nleaf_me,'accumulated'
     write(*,*) 'Fixing and continuing for now..'
     nleaf_me = nleaf_me_check
  endif
  if (ntwig_me /= ntwig_me_check) then
     write(*,'(3a,i4)') 'Table check called ',callpoint,' by PE',me
     write(*,*) ' # own twigs in table = ',ntwig_me_check,'vs ',ntwig_me,'accumulated'
     write(*,*) 'Fixing and continuing for now..'
     ntwig_me = ntwig_me_check
  endif

end subroutine check_table
