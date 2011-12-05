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
     endif

     !if (db_level > 0) then
     !  call system("mkdir -p " // "diag")
     !  write(cfile,'("diag/diag_",i6.6,".dat")') my_rank
     !  open(20, file=trim(cfile),STATUS='UNKNOWN', POSITION = 'APPEND')
     !endif

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


end module files
