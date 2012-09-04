module bindings

  contains
  
    subroutine my_bfunc(n, l, m, rv, maxR, re, im)
      use math
      implicit none
      integer, intent(in) :: n, l, m
      real*8, intent(in)  :: rv(3), maxR
      real*8, intent(out) :: re, im
      complex*16 :: tmp
      
      tmp = Bfunc(n, l, m, rv, maxR)
      re  = real(tmp)
      im  = aimag(tmp)
    end subroutine

    real*8 function my_pfunc(l,m,costheta)
      use math
      implicit none
      integer, intent(in) :: l,m
      real*8, intent(in) :: costheta

      my_pfunc = Pfunc(l,m,costheta)
    end function      
    
    real*8 function my_mmfunc(l,m)
      use math
      implicit none
      integer, intent(in) :: l,m
      
      my_mmfunc = MMFunc(l,m)      
    end function
    
    subroutine my_efunc(m,Phi, re, im)
      use math
      implicit none
      integer, intent(in) :: m
      real*8, intent(in) :: Phi
      real*8, intent(out) :: re, im
      complex*16 :: tmp

      tmp = Efunc(m, Phi)
      re  = real(tmp)
      im  = aimag(tmp)
    end subroutine    
    
    real*8 function my_rfunc(n,l,a,r)
      use math
      implicit none
      integer, intent(in) :: n, l
      real*8, intent(in) :: r, a
      
      my_rfunc = Rfunc(n,l,a,r)
    end function    
    
    
end module
