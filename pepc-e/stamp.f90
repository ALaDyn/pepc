
!     =========================
!     
!     Time stamp 
!     
!     =========================
      
subroutine stamp(istream,ibegin)
  character :: cdate*8, ctime*10, czone*5
  integer :: ibegin

     !      call DATE_AND_TIME(cdate,ctime,czone,vals)
     call DATE_AND_TIME(cdate,ctime,czone)

     if (ibegin.eq.1) then

        write(istream,'(//a20,a12/a20,a12/a20,a12//)') 'PEPC run on ' & 
             ,cdate(7:8)//'/'//cdate(5:6)//'/'//cdate(1:4) &
             ,'Time: ',ctime(1:2)//':'//ctime(3:4)//':'//ctime(5:6),' GMT+',czone

     else 
        write(istream,'(a,a9)') 'Finished run at time: ',ctime(1:2)//':'//ctime(3:4)//':'//ctime(5:6)
     endif
end subroutine stamp


