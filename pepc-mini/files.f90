module files
  implicit none
  private

    public openfiles
    public closefiles



  contains



    subroutine openfiles
      use physvars

      character(30) :: cfile

      if (my_rank == 0) then
         !  master diagnostics output
         open(15,file='run.out')
         open(70,file='domains.dat')

         open(59,file='memory.dat')
         write(59,*) "# ID               bytes                  KB                  MB           location"
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
         close(59)
      endif

      close(20)
      close(80)  ! initial particle data

    end subroutine closefiles


end module files
