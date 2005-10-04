! ======================
!
!   VIS_PARTS
!
!   Send particle data to VISIT for real-time visualisation
!
!
! ======================

subroutine vis_parts_nbody


  use physvars
  use treevars
  implicit none
  include 'mpif.h'

  integer, parameter :: npart_visit_max = 100000  ! Max 250k data points for VIS
  integer, parameter :: ship_max = 100000, attrib_max=22
  real*4, dimension(0:attrib_max-1,npart_visit_max) :: vbuffer
  !  real*4, dimension(0:attrib_max-1,npart_total) :: vbuffer
  real*4, dimension(0:attrib_max-1,npart_visit_max) :: vbuf_local

  integer, dimension(num_pe) :: nparts_pe, recv_strides, nbuf_pe  ! array of npp on each PE
  integer :: icolour(npart_visit_max)
  integer :: lvisit_active, nskip, nproot, ne_buf, ni_buf, npart_buf, nbufe, nbufi
  integer :: i, j, k, ioffset,ixd, iyd, izd, ilev, lcount, wfdatai

  real :: s, simtime, dummy, xd,yd,zd, dx, dz, dy, epond_max, box_max, epondx, epondy, epondz,phipond
  real :: plasma1, plasma2, plasma3, t_display, wfdatar, u2, amp_las, box_x, box_y, box_z
  integer :: nship, ierr, cpuid
  integer :: type
  integer :: vbufcols = 22, incdf, ndom_vis, ivisdom
  real :: lbox, work_ave

  convert_mu=1.
  simtime = dt*(itime+itime_start)

  nskip = npart/ship_max + 1
  if (beam_config==4) then
     amp_las = vosc*min(1.,simtime/tpulse)
  else 
     amp_las = vosc
  endif

  !  if (beam_config>=3 .and. beam_config<=5) then
  !     t_display = tlaser*convert_fs
  !  else
  t_display = simtime
  !  endif

  if (beam_config<=3) then
     plasma1 = Ukine*convert_keV
     plasma2 = Ukini*convert_keV
  else
     plasma1=sigma
     plasma2=tpulse
  endif

  if (nskip>1 .and. uthresh==-1) uthresh=vte**2

  nbufe = 0
  nbufi = 0

  if (me.eq.0) then
     nship=1  ! leave gap on root for info block (filled in later)
  else
     nship = 0
  endif

  ! Filter out particles for shipping
  type = 1
     ! TODO:  separate interesting ions and electrons
     do i=1,npp
        u2=0.5*0.511*mass_ratio*(ux(i)**2+uy(i)**2+uz(i)**2) ! in MeV

        !    if ((npart>=100000 .and. (u2>uthresh .and. q(i)<0)) .or. (npart<100000 .and. mod(pelabel(i),nskip).eq.0)) then
        nship=nship+1
        nbufe=nbufe+1
        !        nbufi=nbufi+1

        !        else if (q(i)>0 .and. mod(pelabel(i),nskip).eq.0) then
        !           nbufi=nbufi+1
        !          nship=nship+1
        ! Store attributes for visualizing
        vbuf_local(0,nship) = t_display
        vbuf_local(1,nship) = dt
        ! scaled coordinates
        vbuf_local(2,nship) = convert_mu*x(i)
        vbuf_local(3,nship) = convert_mu*y(i)
        vbuf_local(4,nship) = convert_mu*z(i)
        ! velocities
        vbuf_local(5,nship) = ux(i)
        vbuf_local(6,nship) = uy(i)
        vbuf_local(7,nship) = uz(i)
        ! E-field
        vbuf_local(8,nship) = Ex(i)
        vbuf_local(9,nship) = Ey(i)
        vbuf_local(10,nship) = Ez(i)
        vbuf_local(11,nship) = 0.
        vbuf_local(12,nship) = 0.
        vbuf_local(13,nship) = 0.
        vbuf_local(14,nship) = pelabel(i)
        vbuf_local(15,nship) = 0.
        vbuf_local(16,nship) = type ! particle type
        vbuf_local(17,nship) = q(i)
        vbuf_local(18,nship) = m(i)
        vbuf_local(19,nship) = pepid(i)
        vbuf_local(20,nship) = 0.
        vbuf_local(21,nship) = 0.
        !   endif
     end do

     ! Find # particles to collect from each cpu
     call MPI_ALLGATHER( nship, 1, MPI_INTEGER, nparts_pe,1, MPI_INTEGER, MPI_COMM_WORLD, ierr )

     ! Total buffer lengths/pe
     nbuf_pe=nparts_pe*attrib_max

     recv_strides(1:num_pe) = (/ 0,( SUM( nbuf_pe(1:i-1) ),i=2,num_pe ) /)

     npart_buf = (SUM(nparts_pe))  ! global total # particles to visualize

     call MPI_ALLREDUCE( nbufe, ne_buf, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr )
     call MPI_ALLREDUCE( nbufi, ni_buf, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr )


     if (npart_buf>0 .and. npart_buf < ship_max) then
        ! send  particle data from all PEs
        ! increase momentum threshold if close to max ship # 
        if (1.0*npart_buf/ship_max > 0.9) uthresh=uthresh*1.1

        call MPI_GATHERV( vbuf_local, nship*attrib_max, MPI_REAL, vbuffer, nbuf_pe, recv_strides, MPI_REAL, 0, MPI_COMM_WORLD, ierr )


        ! Ship to visualization
        if (me.eq.0) then

           ! Stick info block in front of particle data on root

           type=8

              vbuffer(0,1) = t_display
              vbuffer(1,1) = npart_buf-1
              ! scaled coordinates
              vbuffer(2,1) = t_display
              vbuffer(3,1) = 0
              vbuffer(4,1) = 1
              ! velocities
              vbuffer(5,1) = npart_buf-1
              vbuffer(6,1) = 0.
              vbuffer(7,1) = 0.
              ! E-field
              vbuffer(8,1) = 0.
              vbuffer(9,1) = 0.
              vbuffer(10,1) = 0.
              vbuffer(11,1) = 0.
              vbuffer(12,1) = 0.
              vbuffer(13,1) = 0.
              vbuffer(14,1) = 0.
              vbuffer(15,1) = 0.
              vbuffer(16,1) = type ! particle type
              vbuffer(17,1) = 0.
              vbuffer(18,1) = 0.
              vbuffer(19,1) = 0.
              vbuffer(20,1) = 0.
              vbuffer(21,1) = 0.

              work_ave=SUM(work_loads)/num_pe


             ! Add branch nodes to show domains

              ndom_vis=0

! First count how many to be shipped

              do k=1,nbranch_sum
             
                 ilev = log( 1.*branch_key(k) )/log(8.)
                 ixd = SUM( (/ (2**i*ibits( branch_key(k),3*i,1 ), i=0,ilev-1) /) )
                 iyd = SUM( (/ (2**i*ibits( branch_key(k),3*i+1,1 ), i=0,ilev-1) /) )
                 izd = SUM( (/ (2**i*ibits( branch_key(k),3*i+2,1 ), i=0,ilev-1) /) )
                 lbox = boxsize/2**(ilev)          !  box length
                 box_x =  xmin + lbox*(ixd+.5) ! box centres
                 box_y =  ymin + lbox*(iyd+.5) ! box centres
                 box_z =  zmin + lbox*(izd+.5) ! box centres

                 if (box_z < domain_cut) ndom_vis=ndom_vis+1
              end do

              ivisdom = 0
              do k=1,nbranch_sum
             
                 ilev = log( 1.*branch_key(k) )/log(8.)
                 ixd = SUM( (/ (2**i*ibits( branch_key(k),3*i,1 ), i=0,ilev-1) /) )
                 iyd = SUM( (/ (2**i*ibits( branch_key(k),3*i+1,1 ), i=0,ilev-1) /) )
                 izd = SUM( (/ (2**i*ibits( branch_key(k),3*i+2,1 ), i=0,ilev-1) /) )
                 lbox = boxsize/2**(ilev)          !  box length
                 box_x =  xmin + lbox*(ixd+.5) ! box centres
                 box_y =  ymin + lbox*(iyd+.5) ! box centres
                 box_z =  zmin + lbox*(izd+.5) ! box centres

                 if (box_z < domain_cut) then   ! only show domains below z=domain_cut
                    ivisdom=ivisdom+1
                    j=npart_buf+1+ivisdom
                    ! Store attributes for visualizing
                    vbuffer(0,j) = t_display
                    vbuffer(1,j)= box_x
                    vbuffer(2,j)= box_y
                    vbuffer(3,j)= box_z 
                    vbuffer(4,j) = lbox
                    vbuffer(5,j) = branch_owner(k)  ! cpu id of branch node
                    cpuid = branch_owner(k)+1
                    vbuffer(6,j) = work_loads(cpuid)/work_ave   ! total work load of branch node/cpu
                    vbuffer(7,j) = npps(cpuid) ! total # particles on cpu
                    vbuffer(8,j) = 0.
                    vbuffer(9,j) = 0.
                    vbuffer(10,j) = 0.
                    vbuffer(11,j) = 0.
                    vbuffer(12,j) = 0.
                    vbuffer(13,j) = 0.
                    vbuffer(14,j) = 0.
                    vbuffer(15,j) = 0.
                    vbuffer(16,j) = 16    ! Domain type
                    vbuffer(17,j) = ndom_vis  ! Total # branch nodes
                    vbuffer(18,j) = 0.
                    vbuffer(19,j) = 0.
                    vbuffer(20,j) = 0.
                    vbuffer(21,j) = 0.
                 endif
                 !        write (*,'(7f13.4)') vbuffer(1,j),vbuffer(2,j), vbuffer(4,j), vbuffer(5,j), & 
                 !             vbuffer(6,j), vbuffer(7,j), vbuffer(17,j)
              end do

! Fill out dummy values for netcdf
              do j=npart_buf+1+ndom_vis,nbuf_max
                 vbuffer(0:attrib_max-1,j) = 0.
              end do


! ---- Preprocess VISIT setup -----------
#ifdef VISIT_NBODY

     if (me==0) then 
        write(*,*) '# particles shipped ',npart_buf,nship
        write(*,*) '# branches shipped ',ndom_vis, '/', nbranch_sum
        write(*,*) 'Total # objects shipped :',nbranch_sum+1+npart_buf,' /',nbuf_max
        write(*,*) 'u_thresh: (MeV)     ',uthresh
     endif

              call flvisit_nbody2_check_connection(lvisit_active)
! send particles and branch boxes together
              call flvisit_nbody2_partstep_send(vbuffer,npart_buf+ndom_vis+1,attrib_max)
! netcdf needs fixed buffer size, so take max used for initialisation

             call ncnbody_put(ncid,vbuffer,nbuf_max,attrib_max,incdf)

              !
              !        write (90,*) 'local',i
              !        write (90,'((22(f12.5/)//))') vbuf_local(0:attrib_max-1,i)
              !        write (90,*) 'global',i
              !        write (90,'((22(f15.8/)//))') vbuffer(0:attrib_max-1,i)
              !        write (90,*) vbuffer(0:attrib_max-1,i)
              !              end do
#else
!  --- No VISIT installed ----
#endif
!  --- end preprocess --------

           endif

        else
           ! Just send particles on root

! TODO: this needs removing or fixing

           if (me==0) then
              if (npart_buf>ship_max) then
                 write(*,*) 'Too many particles to ship - reduce numbers'
                 uthresh = uthresh*2
                 write(*,*) 'Increasing momentum threshold to: ',sqrt(abs(uthresh))
              else if (npart_buf==0 .and. ne > 0) then
                 uthresh = vte*2
                 write(*,*) 'Reducing momentum threshold to: ',sqrt(abs(uthresh))
              endif
              nproot = 0.8*npart/num_pe ! fixed # parts close to npp

              !           call flvisit_spk_check_connection(lvisit_active)
!              call flvisit_nbody2_check_connection(lvisit_active)
              !        call flvisit_spk_info_send(npp,xl,yl,zl,zl,nep,nip,np_beam,itime+itime_start)
              wfdatai=int(tlaser)
              wfdatar=tlaser
              !           call flvisit_spk_info_send(nproot,xl,yl,zl, t_display, &
              !                x_crit, amp_las, sigma, tpulse, &
              !                ne,ni,npart,wfdatai)
              !           call flvisit_nbody2_info_send(t_display,xl,yl,zl, xl, &
              !             t_display, amp_las, plasma1, plasma2, xl, &
              !             xl, xl, xl, xl, xl, & 
              !             xl, xl, xl, xl, xl )
              !           call flvisit_spk_particles_send(t_display,x,y,z,ux,uy,uz,q,pepid,pelabel,npp)
              !        call flvisit_nbody2_particles_send(t_display,xvis,yvis,zvis,vx,vy,vz,qvis,ppid,plabel,npart_buf)
           endif
        endif


        ! Make sure everyone else knows about new momentum threshold
        call MPI_BCAST( uthresh, 1, MPI_REAL8, 0, MPI_COMM_WORLD,ierr)

        call MPI_BARRIER( MPI_COMM_WORLD, ierr)  ! Wait for everyone to catch up


      end subroutine vis_parts_nbody


