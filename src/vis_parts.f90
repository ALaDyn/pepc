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

  integer, parameter :: npart_visit_max = 350000  ! Max 250k data points for VIS

  real, dimension(npart_visit_max) :: xvis,yvis,zvis,vx,vy,vz,qvis,mvis
  integer, dimension(npart_visit_max) :: ppid, plabel
  real, dimension(nppm) :: xv, yv, zv, uxv,uyv,uzv,qv,mv
  integer, dimension(nppm) :: pepidv, pelabelv

  integer, dimension(num_pe) :: nparts_pe  ! array of npp on each PE
  integer :: icolour(npart_visit_max)
  integer :: lvisit_active, nskip, nproot
  integer :: i, j, k, ioffset,ixd, iyd, izd, ilev, lcount, wfdatai

  real :: s, simtime, dummy, xd,yd,zd, dx, dz, dy, epond_max, box_max, epondx, epondy, epondz,phipond
  real :: xl_mu, yl_mu, zl_mu, t_fs, wfdatar, convert_mu

  convert_mu=1.
  simtime = dt*(itime+itime_start)
  nskip = npart/npart_visit_max + 1
  !  nbuf = npp/nskip

  ! Filter out particles for shipping
  nbufe = 0
  nbufi = 0
  nbuf=0
  do i=1,npp
     if (q(i)<0 ) then  !  pick out even labels
!      if (mod(pelabel(i),5).eq.0) then    ! select 1:5
        nbuf=nbuf+1
        if (q(i)<0) then
           nbufe=nbufe+1
        else
           nbufi=nbufi+1
        endif
        ! scaled coordinates
        xv(nbuf) = convert_mu*x(i)
        yv(nbuf) = convert_mu*y(i)
        zv(nbuf) = convert_mu*z(i)
        uxv(nbuf) = ux(i)
        uyv(nbuf) = uy(i)
        uzv(nbuf) = uz(i)
        qv(nbuf) = q(i)
        !        mv(nbuf) = m(i)
        pepidv(nbuf) = pepid(i)
        pelabelv(nbuf) = pelabel(i)
     endif
  end do

  call MPI_ALLGATHER( nbuf, 1, MPI_INTEGER, nparts_pe,1, MPI_INTEGER, MPI_COMM_WORLD, ierr )
  recv_strides(1:num_pe) = (/ 0,( SUM( nparts_pe(1:i-1) ),i=2,num_pe ) /)
  npart_buf = (SUM(nparts_pe))  ! total # particles to send
  call MPI_ALLREDUCE( nbufe, ne_buf, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr )
  call MPI_ALLREDUCE( nbufi, ni_buf, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr )

  if (me==0) write(*,*) '# particles shipped ',npart_buf

  if (npart_buf<npart_visit_max) then
     ! send  particle data from all PEs
     if (me==0) then
       call flvisit_spk_check_connection(lvisit_active)
!       call flvisit_spk_info_send(npart,xl,yl,zl,zl,ne,ni,np_beam,itime+itime_start)
        wfdatai=int(tlaser)
        wfdatar=tlaser
        call flvisit_spk_info_send(npart,xl,yl,zl, zl, &
             x_crit, vosc, sigma, tpulse, &
             ne,ni,npart,wfdatai)

     endif
!
! visit info block format in spk4
! int flvisit_spk_info_send(int *npart,double *xl,double *yl,double *zl,double *sbox,double *boxoffsetx,double *vosc,double *sigma,double *tphase,int *ne,int *ni,int *np,int *itime);


     ! work out stride lengths so that partial arrays placed sequentially in global array


     ! Gather particle data onto root PE
     call MPI_GATHERV( pepidv, nbuf, MPI_INTEGER, ppid, nparts_pe, recv_strides, MPI_INTEGER, root, MPI_COMM_WORLD, ierr )
!     call MPI_GATHERV( pelabelv, nbuf, MPI_INTEGER, plabel, nparts_pe, recv_strides, MPI_INTEGER, root, MPI_COMM_WORLD, ierr) 
     plabel(1:npart_buf) = (/ (i,i=1,npart_buf) /) ! Relabel 
     call MPI_GATHERV( xv, nbuf, MPI_REAL8, xvis, nparts_pe, recv_strides, MPI_REAL8, root, MPI_COMM_WORLD, ierr )
     call MPI_GATHERV( yv, nbuf, MPI_REAL8, yvis, nparts_pe, recv_strides, MPI_REAL8, root, MPI_COMM_WORLD, ierr )
     call MPI_GATHERV( zv, nbuf, MPI_REAL8, zvis, nparts_pe, recv_strides, MPI_REAL8, root, MPI_COMM_WORLD, ierr )

     call MPI_GATHERV( uxv, nbuf, MPI_REAL8, vx, nparts_pe, recv_strides, MPI_REAL8, root, MPI_COMM_WORLD, ierr )
     call MPI_GATHERV( uyv, nbuf, MPI_REAL8, vy, nparts_pe, recv_strides, MPI_REAL8, root, MPI_COMM_WORLD, ierr )
     call MPI_GATHERV( uzv, nbuf, MPI_REAL8, vz, nparts_pe, recv_strides, MPI_REAL8, root, MPI_COMM_WORLD, ierr )
     call MPI_GATHERV( qv, nbuf, MPI_REAL8, qvis, nparts_pe, recv_strides, MPI_REAL8, root, MPI_COMM_WORLD, ierr )


     if (me.eq.0) then

        call flvisit_spk_particles_send(tlaser,xvis,yvis,zvis,vx,vy,vz,qvis,ppid,plabel,npart_buf)
     endif

  else
     ! Just send particles on root
     if (me==0) then
        nproot = 0.8*npart/num_pe ! fixed # parts close to npp

        call flvisit_spk_check_connection(lvisit_active)
!        call flvisit_spk_info_send(npp,xl,yl,zl,zl,nep,nip,np_beam,itime+itime_start)
        call flvisit_spk_particles_send(simtime,x,y,z,ux,uy,uz,q,pepid,pelabel,npp)

     endif
  endif


  if (me==0) then

     ! ship branch nodes to show domains
     do j=1,nbranch_sum
        ilev = log( 1.*branch_key(j) )/log(2.**idim)
        ixd = SUM( (/ (2**i*ibits( branch_key(j),idim*i,1 ), i=0,ilev-1) /) )
        iyd = SUM( (/ (2**i*ibits( branch_key(j),idim*i+1,1 ), i=0,ilev-1) /) )
        izd = SUM( (/ (2**i*ibits( branch_key(j),idim*i+2,1 ), i=0,ilev-1) /) )
        mvis(j) = boxsize/2**(ilev)          !  box length
        xvis(j)=(xmin + ixd*mvis(j))*convert_mu
        yvis(j)=(ymin + iyd*mvis(j))*convert_mu
        zvis(j)=(zmin + izd*mvis(j))*convert_mu 

     end do
     call flvisit_spk_domains_send( tlaser, xvis,yvis,zvis,mvis,branch_owner,xvis,nbranch_sum)

  endif


  call MPI_BARRIER( MPI_COMM_WORLD, ierr)  ! Wait for everyone to catch up


end subroutine vis_parts

