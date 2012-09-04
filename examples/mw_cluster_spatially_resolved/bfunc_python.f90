module bindings

  contains
  
    complex*16 function my_bfunc(n, l, m, rv, maxR)
      use math
      implicit none
      integer, intent(in) :: n, l, m
      real*8, intent(in)  :: rv(3), maxR
      
      my_bfunc = Bfunc(n, l, m, rv, maxR)
    end function

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
    
    complex*16 function my_efunc(m,Phi)
      use math
      implicit none
      integer, intent(in) :: m
      real*8, intent(in) :: Phi

      my_efunc = Efunc(m, Phi)
    end function    
    
    real*8 function my_rfunc(n,l,a,r)
      use math
      implicit none
      integer, intent(in) :: n, l
      real*8, intent(in) :: r, a
      
      my_rfunc = Rfunc(n,l,a,r)
    end function    
    
    
end module
