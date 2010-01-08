
module treetypes

  type particle
     real*8 :: x    ! coords
     real*8 :: y
     real*8 :: z
     real*8 :: ux    ! momenta
     real*8 :: uy
     real*8 :: uz 
     real*8 :: q     ! charge
     real*8 :: m     ! mass
     real*8 :: work  ! work load from force sum
     real*8 :: ax   ! 'forces' 
     real*8 :: ay
     real*8 :: az
     integer*8 :: key           ! Key
     integer :: label    ! label
     integer :: pid      ! owner
  end type particle

end module treetypes
