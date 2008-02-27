! ==============================================
!
!                SPECIAL_START
!
!  Initialise set of particles with zero velocity
!
! ==============================================

subroutine special_start(iconf)

  use physvars
  implicit none

  integer, intent(in) :: iconf  ! Configuration switch
  integer :: max_num
  real*8 :: yt,zt,xt,qt,mt,c_status
  character(30) :: cfile 
  character(66) :: dfile
  character(25) :: cme
  integer*8 :: p,i,j,plt

  max_num = n_cpu


  config: select case(iconf)

  case(1)  ! electron disc at x=0 with 2pi phase spread

  case(2)
     
     if (my_rank == 0) write(*,*) "Special start 2: Reading initial particle positions from FCS..."
     p = 0
     write(cme,*) my_rank+1
     cme = ADJUSTL(cme)
     dfile = "sort_data/data_25x4/parts_info_"//TRIM(cme)//".dat"
     open(90,file=dfile,STATUS='OLD')
     
     do while (p<np_local)
        read(90,*) i,xt,yt,zt
        p=p+1
        x(p) = xt
        y(p) = yt
        z(p) = zt
        ux(p) = 0.
        uy(p) = 0.
        uz(p) = 0.
        q(p) = qt
        pelabel(p) = i
        if (mod(p,2_8) == 0) then
           m(p) = mass_e
           q(p) = qe
        else 
           m(p) = mass_i
           q(p) = qi
        end if   
     end do
!     write(*,*) my_rank,p,i

     close(90)
     goto 112
  case(3)

     if (my_rank == 0) write(*,*) "Using special start... homogeneous slices in a unity cube\n"

     cme = "sort_data/parts_info_" &   
       // achar(my_rank/1000_8+48) &
       // achar(mod(my_rank/100_8,10_8)+48) &
       // achar(mod(my_rank/10_8,10_8)+48) &
       // achar(mod(my_rank/1_8,10_8)+48)  ! Convert 4-digit PE number into character string
     cfile=cme//".dat"
   
!     open(90,file=cfile,STATUS='REPLACE')

     p=0
     call srand(2.0*my_rank+1)
     do while (p < np_local)
        xt = 0.
        yt = 0.
        zt = 0.
        do while ((xt .le. 1.0*my_rank/(1.0*n_cpu)) .or. (xt .ge. 1.0*(my_rank+1)/(1.0*n_cpu)) .or. (xt == 0))
           xt = rand()
        end do
        
        do while (yt == 0)
           yt = rand()
        end do

        do while (zt == 0)
           zt = rand()
        end do

        p = p + 1

        ux(p) = 0.
        uy(p) = 0.
        uz(p) = 0.

        z(p) = zt 
        y(p) = yt
        x(p) = xt
        if (p.le.nep) then 
!             write(90,*) p,xt,yt,zt,qe,mass_e,my_rank*nep+p      ! Electrons
        else
!             write(90,*) p,xt,yt,zt,qi,mass_i,ne+my_rank*nip+p-nep      ! Ions
        end if
     end do

!    close(90)

  case(4)

     if (my_rank == 0) write(*,*) "Special start 4: Reading initial particle positions..."
     j = 0
     do i = my_rank*max_num/n_cpu , (my_rank+1)*max_num/n_cpu-1
        cme = "sort_data/parts_info_" &   
             // achar(i/1000+48) &
             // achar(mod(i/100,10_8)+48) &
             // achar(mod(i/10,10_8)+48) &
             // achar(mod(i/1,10_8)+48)  ! Convert 4-digit PE number into character string
        cfile=cme//".dat"
        open(90,file=cfile)
        p = 0
        do while (p < npart_total/max_num)
           read(90,*) p,xt,yt,zt,qt,mt,plt
           ux(p+j*npart_total/max_num) = 0.
           uy(p+j*npart_total/max_num) = 0.
           uz(p+j*npart_total/max_num) = 0.
           z(p+j*npart_total/max_num) = zt 
           y(p+j*npart_total/max_num) = yt
           x(p+j*npart_total/max_num) = xt
           q(p+j*npart_total/max_num) = qt
           m(p+j*npart_total/max_num) = mt
           pelabel(p+j*npart_total/max_num) = plt
        end do    

        j = j+1
        close(90)
     end do
     goto 112
  end select config


  q(1:nep)             = qe        ! plasma electrons
  q(nep + 1:np_local)       = qi        ! plasma ions (need Z* here)
  m(1:nep)             = mass_e    ! electron mass
  m(nep + 1:np_local)       = mass_i    ! ion mass
  pelabel(1:nep)       = my_rank * nep + (/(i, i = 1, nep)/)      ! Electron labels
  pelabel(nep + 1:np_local) = ne + my_rank * nip + (/(i, i = 1, nip)/) ! Ion labels

112  ex(1:np_local) = 0.
  ey(1:np_local) = 0.
  ez(1:np_local) = 0.
  pot(1:np_local) = 0.
  work(1:np_local) = 1.

end subroutine special_start



