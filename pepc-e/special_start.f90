! ==============================================
!
!                SPECIAL_START
!
!  Initialise set of particles with zero velocity
!
! ==============================================

subroutine special_start(iconf)

  use physvars
  use utils
  implicit none

  integer, intent(in) :: iconf  ! Configuration switch
  integer :: p,i,j,k,l,n1,iseed1,iseed2,iseed3,face_nr,iseed=-17,decomp,max_num,plt
  real*8 :: gamma0,yt,zt,xt,qt,mt,c_status
  character(35) :: cinfile, cdump, cfile
  character(50) :: dfile
  character(25) :: cme

  iseed1 = -1 - my_rank
  iseed2 = -1011 - my_rank
  iseed3 = -200111 - my_rank
  r_beam=sigma/2.
  gamma0=sqrt(1+vosc**2/2.)

  max_num = 4


  config: select case(iconf)
  case(1)  ! electron disc at x=0 with 2pi phase spread

     do while (i < np_local)
        yt = r_beam*(2*rano(iseed)-1.)

        zt = r_beam*(2*rano(iseed)-1.)
        if (yt**2 + zt**2 <= r_beam**2 ) then
           i = i+1
           x(i) = 2*pi*rano(iseed)
           y(i) = yt + focus(2)
           z(i) = zt + focus(3)
           ux(i) = 0.
           uy(i) = 0.
           uz(i) = 0.

        endif
     end do

  case(2)
!     x_crit = focus(1)
! electron disc at laser focus with 2pi phase spread
!     do while (i < np_local)
!        yt = r_beam*(2*rano(iseed)-1.)
!
!        zt = r_beam*(2*rano(iseed)-1.)
!        if (yt**2 + zt**2 <= r_beam**2 ) then
!           i = i+1
!           x(i) = focus(1)
!           y(i) = yt + focus(2)
!           z(i) = zt + focus(3)
!           ux(i) = 0.
!           uy(i) = 0.
!           uz(i) = 0.
!
!        endif
!     end do
     
     if (my_rank == 0) write(*,*) "Special start 2: Reading initial particle positions from FCS..."
     p = 0
     write(cme,*) my_rank+1
     cme = ADJUSTL(cme)
     dfile = "data/parts_info_"//TRIM(cme)//".dat"
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
        pelabel(p) = i
        if (mod(p,2) == 0) then
           m(p) = mass_e
           q(p) = qe
        else 
           m(p) = mass_i
           q(p) = qi
        end if   
     end do
!     write(*,*) my_rank,np_local,p,i

     close(90)
     goto 112

  case(3)

     if (my_rank == 0) write(*,*) "Using special start... case 3"

     cme = "pepc_data/parts_info_" &   
       // achar(my_rank/1000+48) &
       // achar(mod(my_rank/100,10)+48) &
       // achar(mod(my_rank/10,10)+48) &
       // achar(mod(my_rank/1,10)+48)  ! Convert 4-digit PE number into character string
     cfile=cme//".dat"
   
!     open(90,file=cfile,STATUS='REPLACE')

     p=0
     call srand(2.0*my_rank+1)
!     decomp = 2
     do while (p < np_local)
        xt = 0.
        yt = 0.
        zt = 0.
        do while ((xt .le. 1.0*my_rank/(1.0*n_cpu)) .or. (xt .ge. 1.0*(my_rank+1)/(1.0*n_cpu)) .or. (xt == 0))
           xt = rand()
!*(1.0*(mod(my_rank,decomp)+1)/(1.0*decomp)-1.0*mod(my_rank,decomp)/(1.0*decomp))+1.0*mod(my_rank,decomp)/(1.0*decomp)
!           xt = .5 * x_plasma * (2 * rano(iseed1) - 1.) + plasma_centre(1)
        end do
        
        do while (yt == 0)
           yt = rand()
        end do
!*(1.0*(ishft(iand(my_rank,decomp),-1)+1)/(1.0*decomp)-1.0*ishft(iand(my_rank,decomp),-1)/(1.0*decomp)) &
!             +1.0*ishft(iand(my_rank,decomp),-1)/(1.0*decomp)
!        yt = .5 * y_plasma * (2 * rano(iseed1) - 1.) + plasma_centre(2)
        do while (zt == 0)
           zt = rand()
        end do
!        zt = .5 * z_plasma * (2 * rano(iseed1) - 1.) + plasma_centre(3)

        p = p + 1

        ux(p) = 0.
        uy(p) = 0.
        uz(p) = 0.

!        do face_nr = 1, number_faces
!           call face((/xt, yt, zt/), c_status, face_nr)
!           if (c_status .ge. 0.) exit
!        end do
!       if (c_status .lt. 0) then

           z(p) = zt 
           y(p) = yt
           x(p) = xt
           if (p.le.nep) then 
!              write(90,*) p,xt,yt,zt,qe,mass_e,my_rank*nep+p      ! Electrons
           else
!              write(90,*) p,xt,yt,zt,qi,mass_i,ne+my_rank*nip+p-nep      ! Ions
           end if
!        end if
     end do

    ! scramble to remove correlations
!    iseed2 = -17 - 4 * my_rank
!    n1 = np_local    

    !  exclude odd one out
!    if (mod(np_local,2) .ne. 0) then
!        n1 = n1 - 1
!    end if
    
!    do i = 1, n1
!        k = n1 * rano(iseed2) + 1
!        xt = x(i)
!        yt = y(i)
!        zt = z(i)
!        x(i) = x(k)
!        y(i) = y(k)
!        z(i) = z(k)
!        x(k) = xt
!        y(k) = yt
!        z(k) = zt
!    end do

!    close(90)

  case(4)

     if (my_rank == 0) write(*,*) "Special start 4: Reading initial particle positions..."
     cme = "pepc_data/parts_info_" &   
          // achar(my_rank/1000+48) &
          // achar(mod(my_rank/100,10)+48) &
          // achar(mod(my_rank/10,10)+48) &
          // achar(mod(my_rank/1,10)+48)  ! Convert 4-digit PE number into character string
     cfile=cme//".dat"
     open(90,file=cfile)
     p = 0
     do while (p < np_local)
        read(90,*) k,xt,yt,zt,qt,mt,plt
        p = p+1
        ux(p) = 0.
        uy(p) = 0.
        uz(p) = 0.
        z(p) = zt 
        y(p) = yt
        x(p) = xt
        q(p) = qt
        m(p) = mt
        pelabel(p) = plt
     end do
     
     close(90)

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



