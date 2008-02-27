! ======================
!
!   VIS_PARTS
!
!   Send particle data to VISIT for real-time visualisation
!
!
! ======================

subroutine vis_parts_nbody(vcount)


  use physvars
!  use visit
  implicit none
  include 'mpif.h'

! Buffer arrays vbuffer, vbuf_local allocated in setup_arrays 

  integer, dimension(n_cpu) :: nparts_pe, recv_strides, nbuf_pe  ! array of npp on each PE
  integer :: icolour(nbuf_max)
  integer :: lvisit_active, nskip, nproot, ne_buf, ni_buf, npart_buf, nbufe, nbufi
  integer :: i, j, k, ioffset,ixd, iyd, izd, ilev, lcount, wfdatai

  real :: s, simtime, dummy, xd,yd,zd, dx, dz, dy, epond_max, box_max, epondx, epondy, epondz,phipond
  real :: plasma1, plasma2, plasma3, t_display, wfdatar, u2, amp_las, box_x, box_y, box_z
  integer :: nship, ierr, cpuid
  integer :: type
  integer :: vbufcols = 22, incdf, ndom_vis, ivisdom, ndomain_vis=1000
  real :: lbox, work_ave, upmax, uproton_max, uxmax
  logical :: vis_debug=.false.
  logical :: pick=.false.
  integer :: vcount


  convert_mu=1.
  simtime = dt*(itime+itime_start)

!  pick_particles: select case(vis_select)
!  case(1,3) 
!    nproot = npart
!  case(2)
!    nproot = n_layer(1)
!  case default
    nproot=0
!  end select pick_particles
 
!  if ( nproot > nbuf_max-ndomain_vis) then
!        nskip = nproot/nbuf_max + 2
!  else
	nskip = 1
!  endif

!  if (me==0 .and. nproot> nbuf_max-ndomain_vis) then
!	write(*,*) "VISNB | # particles > vis nbuf_max - reducing number shipped"
!	write(*,*) "VISNB | nbuf_max=",nbuf_max," nskip=",nskip
!  endif

!  if (beam_config==4) then
!     amp_las = vosc*min(1.,simtime/tpulse)
!  else 
!     amp_las = vosc
!  endif

!    if (.not. launch) then
       t_display = vcount
!    else if (beam_config>=3 .and. beam_config<=6) then
!       t_display = tlaser*convert_fs
!    else
!      t_display = simtime
!    endif

!  if (beam_config<=3) then
!     plasma1 = Ukine*convert_keV
!     plasma2 = Ukini*convert_keV
!  else
!     plasma1=sigma
!     plasma2=tpulse
!  endif

!  if (nskip>1 .and. uthresh==-1) uthresh=vte**2

  nbufe = 0
  nbufi = 0

  if (my_rank.eq.0) then
     nship = 1  ! leave gap on root for info block (filled in later)
  else
     nship = 0
  endif

  ! Filter out particles for shipping
  type = 1
!  upmax = 0.
!  uxmax = 0.
!  uproton_max=0.
     ! TODO:  separate interesting ions and electrons
     do i=1,np_local
!        u2=0.5*0.511*mratio_layer(1)*(ux(i)**2+uy(i)**2+uz(i)**2) ! in MeV
!       pick_part: select case(vis_select)
!	case(1)  ! select from bulk
!	  pick = ( npart<nbuf_max .or. (npart > nbuf_max .and. mod(pelabel(i),nskip).eq.0)) 
!	case(2)  ! select from 2nd layer
!          pick = (  mod(pelabel(i),nskip).eq.0 .and. pelabel(i)>npart-n_layer(1) )
!	case(3)  ! select above energy threshold 
!  	  pick = (u2>uthresh)
!	case default
!	  pick = .false.
!	end select pick_part

      pick = .true.
      if (pick) then 
        nship=nship+1
        if (q(i)<0) then
	  nbufe=nbufe+1
        else
          nbufi=nbufi+1
        endif
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
        vbuf_local(19,nship) = my_rank
!        vbuf_local(19,nship) = pepid(i)
        vbuf_local(20,nship) = 0.
        vbuf_local(21,nship) = 0.
!        if (mod(i,2**7) == 0) write(*,*) my_rank,"Im Here: End"
       endif
     end do

     ! Find # particles to collect from each cpu
     call MPI_ALLGATHER( nship, 1, MPI_INTEGER, nparts_pe,1, MPI_INTEGER, MPI_COMM_WORLD, ierr )
!     call MPI_ALLREDUCE( upmax, uproton_max, 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr )

     ! Total buffer lengths/pe
     nbuf_pe=nparts_pe*attrib_max

     recv_strides(1:n_cpu) = (/ 0,( SUM( nbuf_pe(1:i-1) ),i=2,n_cpu ) /)

     npart_buf = (SUM(nparts_pe))  ! global total # particles to visualize

     call MPI_ALLREDUCE( nbufe, ne_buf, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr )
     call MPI_ALLREDUCE( nbufi, ni_buf, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr )


     if (npart_buf>0 .and. npart_buf < nbuf_max-ndomain_vis) then
        ! send  particle data from all PEs
        ! increase momentum threshold if close to max ship # 
!       if (1.0*npart_buf/nbuf_max > 0.9) uthresh=uthresh*1.1

        call MPI_GATHERV( vbuf_local, nship*attrib_max, MPI_REAL, vbuffer, nbuf_pe, recv_strides, MPI_REAL, 0, MPI_COMM_WORLD, ierr )


        ! Ship to visualization
        if (my_rank.eq.0) then

           ! Stick info block in front of particle data on root

           type=8

              vbuffer(0,1) = t_display
!              vbuffer(1,1) = npart_buf-1
              vbuffer(1,1) = npart_total
              ! scaled coordinates
              vbuffer(2,1) = t_display
              vbuffer(3,1) = 0
              vbuffer(4,1) = 1
              ! velocities
              vbuffer(5,1) = npart_buf-1
              vbuffer(6,1) = npart_total
              vbuffer(7,1) = ops_per_sec
              vbuffer(8,1) = uproton_max
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

!             work_ave=SUM(work_loads)/n_cpu


             ! Add branch nodes to show domains
              ndom_vis=0

! Fill out dummy values for netcdf
              do j=npart_buf+1+ndom_vis,nbuf_max
                 vbuffer(0:attrib_max-1,j) = 0.
              end do

! ---- Preprocess VISIT setup -----------
#ifdef VISIT_NBODY

              if (my_rank==0) then 
                 write(*,*) 'VISNB | Display time/count ',t_display
                 write(*,*) 'VISNB | # particles shipped ',npart_buf,nship
!                 write(*,*) 'VISNB | Total # objects shipped :',ndom_vis+1+npart_buf,' /',nbuf_max
!                 write(*,*) 'VISNB | ux, up_max: (MeV)     ',uxmax, upmax, uproton_max
                 write(*,*) 'VISNB | vis_select=  ',vis_select
              endif
              
              write(*,*) my_rank,npart_buf+1
              call flvisit_nbody2_check_connection(lvisit_active)
              ! send particles and branch boxes together
              write(*,*) my_rank,npart_buf+1
              call flvisit_nbody2_partstep_send(vbuffer,npart_buf+1,attrib_max)
! netcdf needs fixed buffer size, so take max used for initialisation
!#ifdef NETCDFLIB
!	      if (netcdf) then
!		call ncnbody_put(ncid,vbuffer,nbuf_max,attrib_max,incdf)
!	        write(*,*) "VIS_PARTS | Writing particles to netcdf"
!          	endif
!#endif

#else
!  --- No VISIT installed ----
#endif
!  --- end preprocess --------


           endif

        endif


        ! Make sure everyone else knows about new momentum threshold
!       call MPI_BCAST( uthresh, 1, MPI_REAL, 0, MPI_COMM_WORLD,ierr)

        call MPI_BARRIER( MPI_COMM_WORLD, ierr)  ! Wait for everyone to catch up

!	vcount=vcount+1

      end subroutine vis_parts_nbody


