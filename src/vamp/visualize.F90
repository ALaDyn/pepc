! ======================
!
!   VISUALIZE
!
!   Send particle data to VISIT for real-time visualisation
!
!
! ======================

subroutine visualize


  use treevars
!VAMPINST include
      INCLUDE 'VTcommon.h'
      INTEGER VTIERR
!
  implicit none   

  integer, parameter :: npart_visit_max = 120000  ! Max 25k data points for VIS

  real, dimension(npart_visit_max) :: xvis,yvis,zvis,vx,vy,vz,qvis,mvis
  integer, dimension(npart_visit_max) :: ppid, plabel

  integer, dimension(num_pe) :: nparts_pe  ! array of npp on each PE
  integer :: icolour(npartm)
  integer :: lvisit_active
  integer :: i, j, k, ioffset,ixd, iyd, izd, ilev, lcount

  real :: s, simtime, dummy, xd,yd,zd, dx, dz, dy, epond_max, box_max, epondx, epondy, epondz,phipond

!VAMPINST subroutine_start
       CALL VTENTER(IF_visualize,VTNOSCL,VTIERR)
!      write(*,*) 'VT: visualize S>',VTIERR,
!     *    IF_visualize,ICLASSH
!
  simtime = dt*(itime+itime_start)

  call MPI_ALLGATHER( npp, one, MPI_INTEGER8, nparts_pe, one, MPI_INTEGER8, MPI_COMM_WORLD, ierr )

  ! work out stride lengths so that partial arrays placed sequentially in global array

  recv_strides(1:num_pe) = (/ 0,( SUM( nparts_pe(1:i-1) ),i=2,num_pe ) /)
  nbuf = npp

  ! Gather particle data onto root PE
  call MPI_GATHERV( pepid, nbuf, MPI_INTEGER8, ppid, nparts_pe, recv_strides, MPI_INTEGER8, root, MPI_COMM_WORLD, ierr )
  call MPI_GATHERV( pelabel, nbuf, MPI_INTEGER8, plabel, nparts_pe, recv_strides, MPI_INTEGER8, root, MPI_COMM_WORLD, ierr) 
  call MPI_GATHERV( x, nbuf, MPI_REAL8, xvis, nparts_pe, recv_strides, MPI_REAL8, root, MPI_COMM_WORLD, ierr )
  call MPI_GATHERV( y, nbuf, MPI_REAL8, yvis, nparts_pe, recv_strides, MPI_REAL8, root, MPI_COMM_WORLD, ierr )
  call MPI_GATHERV( z, nbuf, MPI_REAL8, zvis, nparts_pe, recv_strides, MPI_REAL8, root, MPI_COMM_WORLD, ierr )

  call MPI_GATHERV( ux, nbuf, MPI_REAL8, vx, nparts_pe, recv_strides, MPI_REAL8, root, MPI_COMM_WORLD, ierr )
  call MPI_GATHERV( uy, nbuf, MPI_REAL8, vy, nparts_pe, recv_strides, MPI_REAL8, root, MPI_COMM_WORLD, ierr )
  call MPI_GATHERV( uz, nbuf, MPI_REAL8, vz, nparts_pe, recv_strides, MPI_REAL8, root, MPI_COMM_WORLD, ierr )
  call MPI_GATHERV( q, nbuf, MPI_REAL8, qvis, nparts_pe, recv_strides, MPI_REAL8, root, MPI_COMM_WORLD, ierr )

  if (me.eq.0) then
!     call flvisit_spk_check_connection(lvisit_active)

!     call flvisit_spk_info_send(npart,xl,yl,zl,zl,ne,ni,np_beam,itime+itime_start)

 !    call flvisit_spk_particles_send(simtime,xvis,yvis,zvis,vx,vy,vz,qvis,ppid,plabel,npart)


! ship pond force as 'domain boxes' on grid
     if (beam_config ==4) then

        box_max = xl/25.
        epond_max = sqrt(1.+vosc**2/2.)/2.
        lcount=0
        dx = xl/15.
        dz = zl/9.
        dy = yl/9.
        do i=1,20
           do j=1,9
              do k=1,9
                 lcount=lcount+1
                 xvis(lcount) = -xl+(i-0.5)*dx
                 yvis(lcount) = (j-0.5)*dy
                 zvis(lcount) = (k-0.5)*dz
                 xd = xvis(lcount) - focus(1) ! position relative to laser focus
                 yd = yvis(lcount) - focus(2)
                 zd = zvis(lcount) - focus(3)
                 call fpond(simtime,tpulse,sigma,vosc,omega,-xd,yd,zd,epondx,epondy,epondz,phipond)  ! evaluate pond force
                 mvis(lcount) = max(0.02, abs(epondx)/epond_max*box_max)  ! scale box size for plot
                 icolour(lcount) = phipond/epond_max*box_max*200

              end do
           end do
        end do
  !      call flvisit_spk_domains_send( simtime, xvis,yvis,zvis,mvis,icolour,xvis,lcount)
                 
     else

! ship branch nodes to show domains
        do j=1,nbranch_sum
           ilev = log( 1.*branch_key(j) )/log(2.**idim)
           ixd = SUM( (/ (2**i*ibits( branch_key(j),idim*i,1 ), i=0,ilev-1) /) )
           iyd = SUM( (/ (2**i*ibits( branch_key(j),idim*i+1,1 ), i=0,ilev-1) /) )
           izd = SUM( (/ (2**i*ibits( branch_key(j),idim*i+2,1 ), i=0,ilev-1) /) )
           mvis(j) = boxsize/2**(ilev)          !  box length
           xvis(j)=ixd*mvis(j) + xmin
           yvis(j)=iyd*mvis(j) + ymin
           zvis(j)=izd*mvis(j) + zmin

        end do
 !       call flvisit_spk_domains_send( simtime, xvis,yvis,zvis,mvis,branch_owner,xvis,nbranch_sum)

     endif

  endif
  call MPI_BARRIER( MPI_COMM_WORLD, ierr)  ! Wait for everyone to catch up


!VAMPINST subroutine_end
       CALL VTLEAVE(ICLASSH,VTIERR)
!      write(*,*) 'VT: visualize S<',VTIERR,ICLASSH
!
end subroutine visualize




