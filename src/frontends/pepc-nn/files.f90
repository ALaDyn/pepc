module files
  implicit none
  private

    public openfiles
    public closefiles



  contains



    subroutine openfiles
      use physvars

      if (my_rank == 0) then
         !  master diagnostics output
         open(15,file='run.out')
         open(70,file='domains.dat')
     endif

    end subroutine openfiles




    subroutine closefiles
      use physvars

      if (my_rank == 0) then
         close(15)
         close(70)
      endif

      close(80)  ! initial particle data

    end subroutine closefiles


end module files