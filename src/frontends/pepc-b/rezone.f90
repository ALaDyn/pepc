! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2016 Juelich Supercomputing Centre, 
!                         Forschungszentrum Juelich GmbH,
!                         Germany
! 
! PEPC is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! PEPC is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
! 
! You should have received a copy of the GNU Lesser General Public License
! along with PEPC.  If not, see <http://www.gnu.org/licenses/>.
!

! ==============================================
!
!                REZONE
!
!   Reset particles
!
! ==============================================

subroutine rezone

  use module_physvars
  use module_particle_props
  use module_velocity_setup
  use module_utilities
  implicit none
  include 'mpif.h'

  integer :: i, ierr
  real :: yt, zt, yb, zb, xt
  integer recv_strides(n_cpu)
  integer :: ne_rez, ni_rez, npp_new, nremove, nremove_tot ! local rezone totals
  integer, dimension(np_local) :: spent_label, spent_loc  ! local lists for reset
  integer, dimension(n_cpu) :: nremove_pe
  logical, dimension(np_local) :: remove
  integer, save :: iseed1

  remove = .false.
  ! Find 'spent' particles
  nremove=0
  ne_rez=0
  ni_rez=0
  do i=1,np_local
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
  recv_strides(1:n_cpu) = (/ 0,( SUM( nremove_pe(1:i-1) ),i=2,n_cpu ) /)
  nremove_tot = (SUM(nremove_pe))  ! total # particles to send

  if (my_rank==0) write(*,*) 'total # particles removed ',nremove_tot

  ! Gather all spent labels on root=last_pe
  call MPI_GATHERV( spent_loc, nremove, MPI_INTEGER, spent_label, nremove_pe, recv_strides, &
			MPI_INTEGER, n_cpu-1, MPI_COMM_WORLD, ierr )



  if (my_rank==n_cpu-1) then
     ! Replenish label supply if fewer spent particles than new ones
     if (nremove_tot < 2*nslice) then
        spent_label(nremove_tot+1:2*nslice) = (/ (new_label+i,i=1,2*nslice-nremove_tot) /)
        new_label=new_label+2*nslice-nremove_tot
     endif
     !     write(ipefile,'(a/(i8))') 'Recyled labels for rezone slice: ',(spent_label(i),i=1,nremove_tot)
     !     if (nremove_tot < 2*nslice) write(ipefile,'(a/(i8))') 'New labels for rezone slice: ',(spent_label(i),i=nremove_tot+1,2*nslice)

  endif

  !  Remove spent particles from local arrays

  npp_new = np_local-nremove
  nep = nep-ne_rez
  nip = nip-ni_rez

  x(1:npp_new) = pack(x(1:np_local),mask=.not.remove(1:np_local))
  y(1:npp_new) = pack(y(1:np_local),mask=.not.remove(1:np_local))
  z(1:npp_new) = pack(z(1:np_local),mask=.not.remove(1:np_local))
  ux(1:npp_new) = pack(ux(1:np_local),mask=.not.remove(1:np_local))
  uy(1:npp_new) = pack(uy(1:np_local),mask=.not.remove(1:np_local))
  uz(1:npp_new) = pack(uz(1:np_local),mask=.not.remove(1:np_local))
  q(1:npp_new) = pack(q(1:np_local),mask=.not.remove(1:np_local))
  m(1:npp_new) = pack(m(1:np_local),mask=.not.remove(1:np_local))
  pepid(1:npp_new) = pack(pepid(1:np_local),mask=.not.remove(1:np_local))
  pelabel(1:npp_new) = pack(pelabel(1:np_local),mask=.not.remove(1:np_local))
  np_local = npp_new  ! Reset # live particles

  ! Add new slice of particles to LAST CPU - domain closest to new particles

  if (idim==3) then
     nslice = pi*r_sphere**2*dt/qi  ! # ions to load per disc
  else
     nslice = y_plasma*dt/qi
     zt = 0.
  endif

  if (my_rank==0) write(*,*) 'total # particles added ',2*nslice

  if (my_rank==n_cpu-1) then
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
        x(np_local+i) = window_min + x_plasma + xt 
        y(np_local+i) = plasma_centre(2) + yt 
        if (idim==3) then
           z(np_local+i) = plasma_centre(3) + zt 
        else
           z(np_local+i)=0.
        endif
     end do

     ux(np_local+1:np_local+nslice) = 0. 
     uy(np_local+1:np_local+nslice) = 0. 
     uz(np_local+1:np_local+nslice) = 0. 
     q(np_local+1:np_local+nslice) = qi 
     m(np_local+1:np_local+nslice) = mass_i 
     pepid(np_local+1:np_local+nslice) = my_rank 
     pelabel(np_local+1:np_local+nslice) = spent_label(1:nslice) 
     work(np_local+1:np_local+nslice) = 1.   ! set work load balanced initially

     pot(np_local+1:np_local+nslice) = 0.

     np_local = np_local + nslice  ! particles will get resorted later - for now keep them here

     ! add same # electrons
     q(np_local+1:np_local+nslice) = qe        ! plasma electrons
     m(np_local+1:np_local+nslice) = mass_e      ! electron mass
     pepid(np_local+1:np_local+nslice) = my_rank                ! processor ID


     pelabel(np_local+1:np_local+nslice) =  spent_label(nslice+1:2*nslice)  ! Electron labels: 1->ne copied from ions

     ! zero potentials - should really compute these for electrons

     pot(np_local+1:np_local+nslice) = 0.


!  Initialize new electrons

     do i=1,nslice
        xt = .2*a_ii*(2*rano(iseed1)-1.)

        x(np_local+i) = x(np_local-nslice+i)+xt  ! Place close to ions 

        yb=r_sphere
        zb=r_sphere
!        do while (yb**2 + zb**2 > r_sphere**2) 
           yt = .05*a_ii*(2*rano(iseed1)-1.)          
           zt = .05*a_ii*(2*rano(iseed1)-1.)
           yb = y(np_local-nslice+i)+yt 
           zb = z(np_local-nslice+i)+zt 
!        end do
        y(np_local+i) = yb 
        if (idim==3) then
           z(np_local+i) = zb
        else
           z(np_local+i) = z(np_local-nslice+i)  ! 2D
        endif
     end do


     !  Give electrons a thermal distribution
     if (vte > 0) then
        call maxwell1(ux,np_local+1,nslice,vte)
        call maxwell1(uy,np_local+1,nslice,vte)
        call maxwell1(uz,np_local+1,nslice,vte)
        call scramble_v(ux,uy,uz,np_local+1,nslice)   ! remove x,y,z correlations
     else
        call cold_start(np_local+1,nslice)
     endif

     work(np_local+1:np_local+nslice) = 1.   ! set work load balanced initially

     ! update particle #s
     np_local = np_local + nslice
     nep = nep + nslice
     nip = nip + nslice

  endif

  !  Now displace all particles back by cdt - stay in laser pulse window
  x(1:np_local) = x(1:np_local) - dt

  ! Update total particle #
  call MPI_ALLREDUCE( np_local, npart_total, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr )
  call MPI_ALLREDUCE( nep, ne, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr )
  call MPI_ALLREDUCE( nip, ni, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr )

end subroutine rezone

