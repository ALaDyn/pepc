! ==============================================
!
!                REZONE
!
!   Reset particles
!
! ==============================================

subroutine rezone

  use treevars
  use utils
  implicit none
  integer :: i, p, iseed1, iseed2, iseed3, n1, j, k
  real :: Vdisc, qe_disc, qi_disc, me_disc, mi_disc, dpx, yt, zt, xb, yb, zb, xt
  integer :: lvisit_active
  integer :: ne_reset, ni_reset ! total # electrons, ions to rezone
  integer :: ne_rez, ni_rez ! local rezone totals
  integer, dimension(npp) :: prezi, preze  ! local lists for reset
  integer, dimension(num_pe) :: nelec_pe, nion_pe
  logical, dimension(npp) :: remove

  remove = .false.
! Find 'spent' particles
  ne_rez=0
  ni_rez=0
  do i=1,npp
     if (x(i) < window_min+dt) then
        remove(i) = .true.
        if (q(i)<0) then
           ne_rez = ne_rez+1
           preze(ne_rez) = i
        else
           ni_rez = ni_rez+1
           prezi(ni_rez) = i
        endif
     endif
  end do

!  Remove spent particles from local arrays
  nremove = ne_rez+ni+rez
  npp_new = npp-nremove
  x(1:npp_new) = pack(x(1:npp),mask=remove(1:npp))

  call MPI_ALLGATHER( ne_rez, 1, MPI_INTEGER, nelec_pe,1, MPI_INTEGER, MPI_COMM_WORLD, ierr )
  recv_strides(1:num_pe) = (/ 0,( SUM( nelec_pe(1:i-1) ),i=2,num_pe ) /)
  ne_reset = SUM(nelec_pe)  ! total # electrons to reset

  call MPI_ALLGATHER( ni_rez, 1, MPI_INTEGER, nion_pe,1, MPI_INTEGER, MPI_COMM_WORLD, ierr )
  recv_strides(1:num_pe) = (/ 0,( SUM( nion_pe(1:i-1) ),i=2,num_pe ) /)
  ni_reset = SUM(nion_pe)  ! total # ions to reset


     Vdisc = pi*r_sphere**2*dt       !  disc volume:  r_sphere is radius
     qe_disc = -Vdisc/ne_reset    ! charge
     qi_disc = Vdisc/ni_reset    ! charge
     me_disc = -qe_disc
     mi_disc = me_disc*mass_ratio

  iseed1 = -11 -itime - me      ! Select seed depending on PE
  iseed2 = -10011 -itime - me
  iseed3 = -30013 -itime - me

     do i=1,ni_rez
        yt=r_sphere
        zt=r_sphere
!         
        do while (yt**2 + zt**2 > r_sphere**2 )
           yt = r_sphere*(2*rano(iseed1)-1.)          
           zt = r_sphere*(2*rano(iseed1)-1.)  
        end do

        p = prezi(i)  ! replace ion on this PE
        xb= rano(iseed1)*dt
       
           ! Add starting point
        x(p) = window_min+x_plasma + xb
        y(p) = yt + plasma_centre(2)
        z(p) = zt + plasma_centre(3)
        q(p) = qi
        m(p) = mass_i
        ux(p)= 0.
        uy(p)= 0.
        uz(p)= 0.
        pepid(p) = me                ! processor ID
!        pelabel(p) =  ! keep label
        ax(p) = 0.
        ay(p) = 0.
        az(p) = 0.
!       write (ipefile,'(a20,5f10.3)') 'New ion positions:',q(p),m(p),x(p),y(p),z(p)
     end do

  ! scramble to remove correlations
  iseed3=-17-4*me
  n1=ni_rez
  !  exclude odd one out
  if (mod(ni_rez,2).ne.0) then
     n1=n1-1
  endif

  do  i=1,n1
     p = prezi(i)
     k=prezi(n1*rano(iseed3)+1)
     xt=x(p)
     yt=y(p)
     zt=z(p)
     x(p)=x(k)
     y(p)=y(k)
     z(p)=z(k)
     x(k)=xt
     y(k)=yt
     z(k)=zt
  end do


! same for electrons
     do i=1,ne_rez
        yt=r_sphere
        zt=r_sphere
!         
        do while (yt**2 + zt**2 > r_sphere**2 )
           yt = r_sphere*(2*rano(iseed2)-1.)          
           zt = r_sphere*(2*rano(iseed2)-1.)  
        end do

        p = preze(i)  ! replace electron on this PE
        xb= rano(iseed2)*dt
     
           ! Add starting point
        x(p) = window_min+x_plasma + xb
        y(p) = yt + plasma_centre(2)
        z(p) = zt + plasma_centre(3)
        q(p) = qe
        m(p) = mass_e
        ux(p)= 0.
        uy(p)= 0.
        uz(p)= 0.
        pepid(p) = me                ! processor ID
!        pelabel(p) =  ! keep label
        ax(p) = 0.
        ay(p) = 0.
        az(p) = 0.
!        write (ipefile,'(a20,5f10.3)') 'New electron positions:',q(p),m(p),x(p),y(p),z(p)
     end do

  ! scramble to remove correlations
  iseed3=-17-4*me
  n1=ne_rez
  !  exclude odd one out
  if (mod(ne_rez,2).ne.0) then
     n1=n1-1
  endif

  do  i=1,n1
     p = preze(i)
     k=preze(n1*rano(iseed3)+1)
     xt=x(p)
     yt=y(p)
     zt=z(p)
     x(p)=x(k)
     y(p)=y(k)
     z(p)=z(k)
     x(k)=xt
     y(k)=yt
     z(k)=zt
  end do

!  Now displace all particles back by cdt - stay in laser pulse window
    x(1:npp) = x(1:npp) - dt

end subroutine rezone

