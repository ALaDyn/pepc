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

  integer, parameter :: npart_visit_max = 50000  ! Max 250k data points for VIS

  real, dimension(npart_visit_max) :: xvis,yvis,zvis,vx,vy,vz,qvis,mvis
  integer, dimension(npart_visit_max) :: ppid, plabel
  real, dimension(nppm) :: xv, yv, zv, uxv,uyv,uzv,qv,mv
  integer, dimension(nppm) :: pepidv, pelabelv

  integer, dimension(num_pe) :: nparts_pe  ! array of npp on each PE
  integer :: icolour(npart_visit_max)
  integer :: lvisit_active, nskip, nproot, ne_buf, ni_buf, npart_buf, nbufe, nbufi
  integer :: i, j, k, ioffset,ixd, iyd, izd, ilev, lcount, wfdatai

  real :: s, simtime, dummy, xd,yd,zd, dx, dz, dy, epond_max, box_max, epondx, epondy, epondz,phipond
  real :: xl_mu, yl_mu, zl_mu, t_fs, wfdatar, u2, amp_las
  integer :: nship

  convert_mu=1.
  simtime = dt*(itime+itime_start)
  nskip = npart/npart_visit_max + 1
  amp_las = vosc*min(1.,simtime/tpulse)

  if (nskip>1 .and. uthresh==-1) uthresh=vte**2

  ! Filter out particles for shipping
  nbufe = 0
  nbufi = 0
  nship=0


  do i=1,npp
     u2=0.5*0.511*mass_ratio*(ux(i)**2+uy(i)**2+uz(i)**2) ! in MeV
    if (u2>uthresh .and. q(i)>0) then
        nship=nship+1
        nbufi=nbufi+1

        !        else if (q(i)>0 .and. mod(pelabel(i),nskip).eq.0) then
        !           nbufi=nbufi+1
        !          nship=nship+1
     ! scaled coordinates
     xv(nship) = convert_mu*x(i)
     yv(nship) = convert_mu*y(i)
     zv(nship) = convert_mu*z(i)
     uxv(nship) = ux(i)
     uyv(nship) = uy(i)
     uzv(nship) = uz(i)
     qv(nship) = q(i)
     !        mv(nship) = m(i)
     pepidv(nship) = pepid(i)
     pelabelv(nship) = pelabel(i)
   endif
  end do

  call MPI_ALLGATHER( nship, 1, MPI_INTEGER, nparts_pe,1, MPI_INTEGER, MPI_COMM_WORLD, ierr )
  recv_strides(1:num_pe) = (/ 0,( SUM( nparts_pe(1:i-1) ),i=2,num_pe ) /)
  npart_buf = (SUM(nparts_pe))  ! total # particles to send
  call MPI_ALLREDUCE( nbufe, ne_buf, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr )
  call MPI_ALLREDUCE( nbufi, ni_buf, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr )

  if (me==0) write(*,*) '# particles shipped ',npart_buf

  if (npart_buf>0 .and. npart_buf < npart_visit_max) then
     ! send  particle data from all PEs
     ! increase momentum threshold if close to max ship # 
     if (1.0*npart_buf/npart_visit_max > 0.9) uthresh=uthresh*1.1
     if (me==0) then
        call flvisit_spk_check_connection(lvisit_active)
        !       call flvisit_spk_info_send(npart_buf,xl,yl,zl,zl,ne,ni,np_beam,itime+itime_start)
        wfdatai=int(tlaser)
        wfdatar=tlaser
        call flvisit_spk_info_send(npart_buf,xl,yl,zl, zl, &
             x_crit, amp_las, sigma, tpulse, &
             ne,ni,npart,wfdatai)

     endif
     !
     ! visit info block format in spk4
     ! int flvisit_spk_info_send(int *npart,double *xl,double *yl,double *zl,double *sbox,double *boxoffsetx,double *vosc,double *sigma,double *tphase,int *ne,int *ni,int *np,int *itime);


     ! work out stride lengths so that partial arrays placed sequentially in global array


     ! Gather particle data onto root PE
     call MPI_GATHERV( pepidv, nship, MPI_INTEGER, ppid, nparts_pe, recv_strides, MPI_INTEGER, root, MPI_COMM_WORLD, ierr )
     !     call MPI_GATHERV( pelabelv, nship, MPI_INTEGER, plabel, nparts_pe, recv_strides, MPI_INTEGER, root, MPI_COMM_WORLD, ierr) 
     plabel(1:npart_buf) = (/ (i,i=1,npart_buf) /) ! Relabel 
     call MPI_GATHERV( xv, nship, MPI_REAL8, xvis, nparts_pe, recv_strides, MPI_REAL8, root, MPI_COMM_WORLD, ierr )
     call MPI_GATHERV( yv, nship, MPI_REAL8, yvis, nparts_pe, recv_strides, MPI_REAL8, root, MPI_COMM_WORLD, ierr )
     call MPI_GATHERV( zv, nship, MPI_REAL8, zvis, nparts_pe, recv_strides, MPI_REAL8, root, MPI_COMM_WORLD, ierr )

     call MPI_GATHERV( uxv, nship, MPI_REAL8, vx, nparts_pe, recv_strides, MPI_REAL8, root, MPI_COMM_WORLD, ierr )
     call MPI_GATHERV( uyv, nship, MPI_REAL8, vy, nparts_pe, recv_strides, MPI_REAL8, root, MPI_COMM_WORLD, ierr )
     call MPI_GATHERV( uzv, nship, MPI_REAL8, vz, nparts_pe, recv_strides, MPI_REAL8, root, MPI_COMM_WORLD, ierr )
     call MPI_GATHERV( qv, nship, MPI_REAL8, qvis, nparts_pe, recv_strides, MPI_REAL8, root, MPI_COMM_WORLD, ierr )

 
     if (me.eq.0) then

        call flvisit_spk_particles_send(tlaser,xvis,yvis,zvis,vx,vy,vz,qvis,ppid,plabel,npart_buf)
     endif

  else
     ! Just send particles on root
     if (me==0) then
        if (npart_buf>npart_visit_max) then
           write(*,*) 'Too many particles to ship - sending ones on root'
           uthresh = uthresh*2
           write(*,*) 'Increasing momentum threshold to: ',sqrt(abs(uthresh))
        else if (npart_buf==0 .and. ne > 0) then
           uthresh = vte*2
           write(*,*) 'Reducing momentum threshold to: ',sqrt(abs(uthresh))
        endif
           nproot = 0.8*npart/num_pe ! fixed # parts close to npp

           call flvisit_spk_check_connection(lvisit_active)
           !        call flvisit_spk_info_send(npp,xl,yl,zl,zl,nep,nip,np_beam,itime+itime_start)
           wfdatai=int(tlaser)
           wfdatar=tlaser
           call flvisit_spk_info_send(nproot,xl,yl,zl, zl, &
                x_crit, amp_las, sigma, tpulse, &
                ne,ni,npart,wfdatai)

           call flvisit_spk_particles_send(tlaser,x,y,z,ux,uy,uz,q,pepid,pelabel,npp)

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
           xvis(j)=(xmin + ixd*mvis(j))*convert_mu ! centres?
           yvis(j)=(ymin + iyd*mvis(j))*convert_mu
           zvis(j)=(zmin + izd*mvis(j))*convert_mu 
           mvis(j) = mvis(j)*0.8
        end do
        call flvisit_spk_domains_send( tlaser, xvis,yvis,zvis,mvis,branch_owner,xvis,nbranch_sum)

     endif

! Make sure everyone else knows about new momentum threshold
     call MPI_BCAST( uthresh, 1, MPI_REAL, root, MPI_COMM_WORLD,ierr)

     call MPI_BARRIER( MPI_COMM_WORLD, ierr)  ! Wait for everyone to catch up


   end subroutine vis_parts


