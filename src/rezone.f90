! ==============================================
!
!                REZONE
!
!   Reset particles
!
! ==============================================

subroutine rezone

  use treevars
  use physvars
  use utils
  implicit none

  integer :: i, p, iseed2, iseed3, n1, j, k, ierr
  real :: Vdisc, qe_disc, qi_disc, me_disc, mi_disc, dpx, yt, zt, xb, yb, zb, xt
  integer recv_strides(num_pe)
  integer :: lvisit_active
  integer :: ne_reset, ni_reset ! total # electrons, ions to rezone
  integer :: ne_rez, ni_rez, npp_new, nremove, nremove_tot ! local rezone totals
  integer, dimension(npp) :: spent_label, spent_loc  ! local lists for reset
  integer, dimension(num_pe) :: nelec_pe, nremove_pe
  logical, dimension(npp) :: remove
  integer, save :: iseed1

  remove = .false.
  ! Find 'spent' particles
  nremove=0
  ne_rez=0
  ni_rez=0
  do i=1,npp
     if (x(i) < window_min+dt) then
        remove(i) = .true.
        nremove = nremove+1
        spent_loc(nremove) = pelabel(i)
        if (q(i)>0) then
           ni_rez=ni_rez+1
        else
           ne_rez=ne_rez+1
        endif

     endif
  end do

  call MPI_ALLGATHER( nremove, 1, MPI_INTEGER, nremove_pe,1, MPI_INTEGER, MPI_COMM_WORLD, ierr )
  recv_strides(1:num_pe) = (/ 0,( SUM( nremove_pe(1:i-1) ),i=2,num_pe ) /)
  nremove_tot = (SUM(nremove_pe))  ! total # particles to send

  if (me==0) write(*,*) 'total # particles removed ',nremove_tot

  ! Gather all spent labels on root=last_pe
  call MPI_GATHERV( spent_loc, nremove, MPI_INTEGER, spent_label, nremove_pe, recv_strides, MPI_INTEGER, num_pe-1, MPI_COMM_WORLD, ierr )



  if (me==num_pe-1) then
     ! Replenish label supply if fewer spent particles than new ones
     if (nremove_tot < 2*nslice) then
        spent_label(nremove_tot+1:2*nslice) = (/ (new_label+i,i=1,2*nslice-nremove_tot) /)
        new_label=new_label+2*nslice-nremove_tot
     endif
     !     write(ipefile,'(a/(i8))') 'Recyled labels for rezone slice: ',(spent_label(i),i=1,nremove_tot)
     !     if (nremove_tot < 2*nslice) write(ipefile,'(a/(i8))') 'New labels for rezone slice: ',(spent_label(i),i=nremove_tot+1,2*nslice)

  endif

  !  Remove spent particles from local arrays

  npp_new = npp-nremove
  nep = nep-ne_rez
  nip = nip-ni_rez

  x(1:npp_new) = pack(x(1:npp),mask=.not.remove(1:npp))
  y(1:npp_new) = pack(y(1:npp),mask=.not.remove(1:npp))
  z(1:npp_new) = pack(z(1:npp),mask=.not.remove(1:npp))
  ux(1:npp_new) = pack(ux(1:npp),mask=.not.remove(1:npp))
  uy(1:npp_new) = pack(uy(1:npp),mask=.not.remove(1:npp))
  uz(1:npp_new) = pack(uz(1:npp),mask=.not.remove(1:npp))
  q(1:npp_new) = pack(q(1:npp),mask=.not.remove(1:npp))
  m(1:npp_new) = pack(m(1:npp),mask=.not.remove(1:npp))
  pepid(1:npp_new) = pack(pepid(1:npp),mask=.not.remove(1:npp))
  pelabel(1:npp_new) = pack(pelabel(1:npp),mask=.not.remove(1:npp))
  npp = npp_new  ! Reset # live particles

  ! Add new slice of particles to LAST CPU - domain closest to new particles

  if (idim==3) then
     nslice = pi*r_sphere**2*dt/qi  ! # ions to load per disc
  else
     nslice = y_plasma*dt/qi
     zt = 0.
  endif

  if (me==0) write(*,*) 'total # particles added ',2*nslice

  if (me==num_pe-1) then
     iseed1=-17999-itime

     ! initialize new ion disc

     do i=1,nslice
        xt = dt*rano(iseed1)
        yt=r_sphere
        zt=r_sphere
        do while (yt**2 + zt**2 > r_sphere**2) 
           yt = r_sphere*(2*rano(iseed1)-1.)          
           zt = r_sphere*(2*rano(iseed1)-1.)  
        end do
        !        write(*,*) 'seed, delta-r:',iseed1,xt,yt,zt
        x(npp+i) = window_min + x_plasma + xt 
        y(npp+i) = plasma_centre(2) + yt 
        if (idim==3) then
           z(npp+i) = plasma_centre(3) + zt 
        else
           z(npp+i)=0.
        endif
     end do

     ux(npp+1:npp+nslice) = 0. 
     uy(npp+1:npp+nslice) = 0. 
     uz(npp+1:npp+nslice) = 0. 
     q(npp+1:npp+nslice) = qi 
     m(npp+1:npp+nslice) = mass_i 
     pepid(npp+1:npp+nslice) = me 
     pelabel(npp+1:npp+nslice) = spent_label(1:nslice) 
     work(npp+1:npp+nslice) = 1.   ! set work load balanced initially
     ax(npp+1:npp+nslice) = 0.
     ay(npp+1:npp+nslice) = 0.
     az(npp+1:npp+nslice) = 0.
     pot(npp+1:npp+nslice) = 0.

     npp = npp + nslice  ! particles will get resorted later - for now keep them here

     ! add same # electrons
     q(npp+1:npp+nslice) = qe        ! plasma electrons
     m(npp+1:npp+nslice) = mass_e      ! electron mass
     pepid(npp+1:npp+nslice) = me                ! processor ID


     pelabel(npp+1:npp+nslice) =  spent_label(nslice+1:2*nslice)  ! Electron labels: 1->ne copied from ions

     ! zero accelerations - should really compute these for electrons
     ax(npp+1:npp+nslice) = 0.
     ay(npp+1:npp+nslice) = 0.
     az(npp+1:npp+nslice) = 0.
     pot(npp+1:npp+nslice) = 0.


!  Initialize new electrons

     do i=1,nslice
        xt = .2*a_ii*(2*rano(iseed1)-1.)

        x(npp+i) = x(npp-nslice+i)+xt  ! Place close to ions 

        yb=r_sphere
        zb=r_sphere
!        do while (yb**2 + zb**2 > r_sphere**2) 
           yt = .05*a_ii*(2*rano(iseed1)-1.)          
           zt = .05*a_ii*(2*rano(iseed1)-1.)
           yb = y(npp-nslice+i)+yt 
           zb = z(npp-nslice+i)+zt 
!        end do
        y(npp+i) = yb 
        if (idim==3) then
           z(npp+i) = zb
        else
           z(npp+i) = z(npp-nslice+i)  ! 2D
        endif
     end do


     !  Give electrons a thermal distribution
     if (vte > 0) then
        call maxwell1(ux,nppm,npp+1,nslice,vte)
        call maxwell1(uy,nppm,npp+1,nslice,vte)
        call maxwell1(uz,nppm,npp+1,nslice,vte)
        call scramble_v(npp+1,nslice)   ! remove x,y,z correlations
     else
        call cold_start(npp+1,nslice)
     endif

     work(npp+1:npp+nslice) = 1.   ! set work load balanced initially

     ! update particle #s
     npp = npp + nslice
     nep = nep + nslice
     nip = nip + nslice

  endif

  !  Now displace all particles back by cdt - stay in laser pulse window
  x(1:npp) = x(1:npp) - dt

  ! Update total particle #
  call MPI_ALLREDUCE( npp, npart, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr )
  call MPI_ALLREDUCE( nep, ne, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr )
  call MPI_ALLREDUCE( nip, ni, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr )

end subroutine rezone

