module bindings

  contains
    ! these routines behave like functions in python due to
    ! the final parameter being intent(out)
    ! WARNING: Do not use functions here, since
    ! their return value is sometimes not passed correctly
    ! to python, e.g. a complex*16 function always (!) returns 1.9j
  
  
    subroutine my_bfunc(n, l, m, rv, maxR, out)
      use math
      implicit none
      integer, intent(in) :: n, l, m
      real*8, intent(in)  :: rv(3), maxR
      complex*16, intent(out) :: out
      
      out = Bfunc(n, l, m, rv, maxR)
    end subroutine

    subroutine my_pfunc(l,m,costheta, out)
      use math
      implicit none
      integer, intent(in) :: l,m
      real*8, intent(in) :: costheta
      real*8, intent(out) :: out

      out = Pfunc(l,m,costheta)
    end subroutine      
    
    subroutine my_mmfunc(l,m, out)
      use math
      implicit none
      integer, intent(in) :: l,m
      real*8, intent(out) :: out
      
      out = MMFunc(l,m)      
    end subroutine
    
    subroutine my_efunc(m,Phi, out)
      use math
      implicit none
      integer, intent(in) :: m
      real*8, intent(in) :: Phi
      complex*16, intent(out) :: out

      out = Efunc(m, Phi)
    end subroutine    
    
    subroutine my_rfunc(n,l,a,r, out)
      use math
      implicit none
      integer, intent(in) :: n, l
      real*8, intent(in) :: r, a
      real*8, intent(out) :: out
      
      out = Rfunc(n,l,a,r)
    end subroutine    
    
    
end module
