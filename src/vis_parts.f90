! ======================
!
!   VIS_PARTS
!
!   Send particle data to VISIT for real-time visualisation
!
!
! ======================

subroutine vis_parts


  use treevars
  implicit none   

  integer, parameter :: npart_visit_max = 12000  ! Max 25k data points for VIS

  real, dimension(npart_visit_max) :: xvis,yvis,zvis,vx,vy,vz,qvis,mvis
  integer, dimension(npart_visit_max) :: ppid, plabel

  integer, dimension(num_pe) :: nparts_pe  ! array of npp on each PE
  integer :: icolour(npart_visit_max)
  integer :: lvisit_active, nskip
  integer :: i, j, k, ioffset,ixd, iyd, izd, ilev, lcount

  real :: s, simtime, dummy, xd,yd,zd, dx, dz, dy, epond_max, box_max, epondx, epondy, epondz,phipond

  simtime = dt*(itime+itime_start)
  nskip = npart/npart_visit_max + 1
  nbuf = npp/nskip

  call MPI_ALLGATHER( nbuf, 1, MPI_INTEGER, nparts_pe,1, MPI_INTEGER, MPI_COMM_WORLD, ierr )

  ! work out stride lengths so that partial arrays placed sequentially in global array

  recv_strides(1:num_pe) = (/ 0,( SUM( nparts_pe(1:i-1) ),i=2,num_pe ) /)

  ! Gather particle data onto root PE
  call MPI_GATHERV( pepid, nbuf, MPI_INTEGER, ppid, nparts_pe, recv_strides, MPI_INTEGER, root, MPI_COMM_WORLD, ierr )
  call MPI_GATHERV( pelabel, nbuf, MPI_INTEGER, plabel, nparts_pe, recv_strides, MPI_INTEGER, root, MPI_COMM_WORLD, ierr) 
  call MPI_GATHERV( x, nbuf, MPI_REAL8, xvis, nparts_pe, recv_strides, MPI_REAL8, root, MPI_COMM_WORLD, ierr )
  call MPI_GATHERV( y, nbuf, MPI_REAL8, yvis, nparts_pe, recv_strides, MPI_REAL8, root, MPI_COMM_WORLD, ierr )
  call MPI_GATHERV( z, nbuf, MPI_REAL8, zvis, nparts_pe, recv_strides, MPI_REAL8, root, MPI_COMM_WORLD, ierr )

  call MPI_GATHERV( ux, nbuf, MPI_REAL8, vx, nparts_pe, recv_strides, MPI_REAL8, root, MPI_COMM_WORLD, ierr )
  call MPI_GATHERV( uy, nbuf, MPI_REAL8, vy, nparts_pe, recv_strides, MPI_REAL8, root, MPI_COMM_WORLD, ierr )
  call MPI_GATHERV( uz, nbuf, MPI_REAL8, vz, nparts_pe, recv_strides, MPI_REAL8, root, MPI_COMM_WORLD, ierr )
  call MPI_GATHERV( q, nbuf, MPI_REAL8, qvis, nparts_pe, recv_strides, MPI_REAL8, root, MPI_COMM_WORLD, ierr )

  if (me.eq.0) then
     call flvisit_spk_check_connection(lvisit_active)
     call flvisit_spk_info_send(npart/nskip,xl,yl,zl,zl,ne,ni,np_beam,itime+itime_start)
     call flvisit_spk_particles_send(simtime,xvis,yvis,zvis,vx,vy,vz,qvis,ppid,plabel,npart/nskip)
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
     call flvisit_spk_domains_send( simtime, xvis,yvis,zvis,mvis,branch_owner,xvis,nbranch_sum)

  endif




  call MPI_BARRIER( MPI_COMM_WORLD, ierr)  ! Wait for everyone to catch up


end subroutine vis_parts






