module files
  implicit none
  private

    public openfiles
    public closefiles
    public stamp


  contains



    subroutine openfiles
      use physvars

      character(30) :: cfile

      if (my_rank == 0) then
         !  master diagnostics output
         open(15,file='run.out')
         open(70,file='domains.dat')
     endif

     if (db_level > 0) then
       write(cfile,'(a,i6.6,a)') "diag_", my_rank, ".dat"
       open(20, file=cfile,STATUS='UNKNOWN', POSITION = 'APPEND')
     endif

    end subroutine openfiles




    subroutine closefiles
      use physvars

      if (my_rank == 0) then
         close(15)
         close(70)
      endif

      close(20)
      close(80)  ! initial particle data

    end subroutine closefiles


	subroutine stamp(istream,ibegin)
        character :: cdate*8, ctime*10, czone*5
        integer :: ibegin
        integer :: istream

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



end module files
