module bindings

  contains
  
    complex*16 function my_bfunc(n, l, m, rv, use_Raitza_definition, maxR)
      use math
      implicit none
      integer, intent(in) :: n, l, m
      real*8, intent(in)  :: rv(3), maxR
      logical, intent(in) :: use_Raitza_definition
      
      my_bfunc = Bfunc(n, l, m, rv, use_Raitza_definition, maxR)
      
    end function
    
end module
