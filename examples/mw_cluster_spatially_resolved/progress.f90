module progress_bar
  implicit none
  
  logical, private, parameter :: showprogress            = .true.
  integer, private, parameter :: barlen                  = 100
  character, private, parameter :: progress_char         = '='
  character, private, parameter :: progress_char_special = '+'

  contains
  
  subroutine progress(j,k,special)  
    implicit none  
    integer, intent(in):: j,k  
    logical, intent(in), optional :: special
    character(len=18) :: barfmt = "(a1,I3,a,x,???a,$)"
    
    integer :: value, mybarlen
    integer, save :: lastvalue = -1
    
    character, save :: bar(0:barlen+1)
    character :: pchar
    
    if (.not. showprogress) return
    
    mybarlen = min(barlen, k)
    value    = (mybarlen*j)/k
    
    if (value .eq. lastvalue) return
    lastvalue = value
    
    pchar = progress_char
    
    if (present(special)) then
      if (special) then
        pchar = progress_char_special
      endif
    endif
    
    if (value==0) then
      bar(0)          = "|"
      bar(1:mybarlen) = " "
      bar(mybarlen+1) = "|"
    else
      bar(value)        = pchar
    endif
    

    write(barfmt(12:14), '(I3.3)') mybarlen+2
    
    ! print the progress bar.  
    write(6, barfmt) char(13), (j*100)/k, "%", bar(0:mybarlen+1)
    
  end subroutine progress  

end module
