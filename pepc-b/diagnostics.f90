!  ================================
!
!         DIAGNOSTICS
!
!     Perform diagnostics
!
!  $Revision 1.15$
!
!  ================================


subroutine diagnostics

  use treevars
  use physvars
  include 'mpif.h'

  implicit none
  integer :: i,lvisit_active, ifile, ierr
  integer :: max_fetches, max_reqs, max_sum_fetches, max_local_f, max_local_r, sum_local_f

  ! Interface to VISIT (Online visualisation)


  if (u_beam>0 .and. beam_config==3) scheme=1  ! Switch off Te control if beam on

  if ( vis_on ) then
     !     if ( mod(itime,ivis) ==0 ) call vis_parts       
     if ( mod(itime,ivis) ==0 ) call vis_parts_nbody       
!     if ( mod(itime,ivis_domains) ==0 ) call vis_domains_nbody       
     if ( mod(itime,ivis)==0 .and. steering) call beam_control
     if ( mod(itime,ivis_fields)==0 ) then
        !     call pot_grid
        call densities
        call sum_fields
        call vis_fields
     endif
  endif

  !  if (target_geometry.eq.4) then
  !    do i=1,npp
  !	if (pelabel(i)==60) then
  !	if (itime.eq.0) write(90,'(a)') '! t x ux Ex Ax Axo -dA/dt Bx' 
  !       write (90,'(8(1pe15.4))') itime*dt,x(i),ux(i),Ex(i),Ax(i),Axo(i),(Axo(i)-Ax(i))/dt,Bx(i)
  !	endif
  !    enddo
  !  endif

  if (mod(itime,idump) >= idump-navcycle) call sum_fields    ! Accumulate cycle-averaged fields on grid
  !  - assume for now that idump > navcycle

  if (beam_config == 4 .and. mod(itime,itrack)==0 ) call track_nc          ! Gather densities and track critical surface 



  if (itime_start>0 .and. itime==0) return  ! Avoid over-writing restart data
  call energy_cons(Ukine,Ukini,Umagnetic,Ubeam)       ! Compute energy balance

  if ( dump_tree .and. mod(itime,iprot) ==0 ) then
     call diagnose_tree   ! Printed tree info (htable etc)
     call draw_tree2d(xl,yl)     ! Draw PE-trees
     call draw_lists      ! Draw interaction lists
     call draw_domains(itime+itime_start)   ! Domains
  endif

  if ( idump>0 .and. (mod(itime+itime_start,idump)==0 .or. itime==nt) ) then
     call dump(itime+itime_start)     ! Dump complete set of particle data
     call dump_fields(itime+itime_start)  ! Field data
     if (vis_on)  call vis_fields
  endif

!  if (debug_level.ge.2) then
     max_local_f =  maxval(nfetch_total)
     sum_local_f =  sum(nfetch_total)  
     call MPI_ALLREDUCE(max_local_f, max_fetches, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr )  
     call MPI_ALLREDUCE(sum_local_f, max_sum_fetches, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr )  
! max_reqs, max_fetches should be equal
!  endif

  if (my_rank.eq.debug_rank .and. debug_level.ge.2) then

     do ifile = 6,15,9
        write(ifile,'(/a,i4)') 'Tree stats for CPU ',debug_rank
        write(ifile,'(a50,2i8,a3,i8)') 'new npp, npart, (max): ',npp,npart,'/',nppm
        write(ifile,'(a50,3i8)') 'local # leaves, twigs, keys: ',nleaf_me,ntwig_me,nleaf_me+ntwig_me
        write(ifile,'(a50,3i8)') 'non-local # leaves, twigs, keys: ',nleaf-nleaf_me,ntwig-ntwig_me,nleaf+ntwig-nleaf_me-ntwig_me
        write(ifile,'(a50,3i8,f12.1,a6,i8)') 'final # leaves, twigs, keys, (max): ',nleaf,ntwig,nleaf+ntwig, &
             (nleaf+ntwig)/(.01*maxaddress),' % of ',maxaddress
        write(ifile,'(a50,2i8,a3,i7)') 'local, global # branches, (max): ',nbranch,nbranch_sum,'/',nbranch_max
        write (ifile,'(a50,i8,a3,i7)') 'Max length of all interaction lists: ',max_list_length,' / ',nintmax
        write (ifile,'(a50,i8)') 'Max # traversals ',maxtraverse

        write (ifile,'(a50,i8)') 'Max # multipole ships/cpu/walk ',max_fetches-max_prefetches
        write (ifile,'(a50,i8)') 'Max # multipole ships/cpu/prefetch ',max_prefetches
        write (ifile,'(a50,i8,a3,i7)') 'Max # multipole ships/cpu/iteration / limit',max_fetches,'/',size_fetch
        write (ifile,'(a50,2i8)') 'Total multipole ships/prefetch, global max ',sum_local_f,max_sum_fetches
!        write (ifile,*) ' cumulative # requested keys:  ',nreqs_total
!        write (ifile,*) ' cumulative # fetched keys:    ',nfetch_total

     end do
  endif

! Array bound check

  if (nleaf+ntwig > .95*maxaddress) then
     write (6,'(a,i4)') '*** WARNING:  hash table >95% full on CPU ',my_rank 
     call cleanup
  endif
  if (max_local_f > .95*size_fetch) then
     write (6,'(a,i4,a1,i8,a3,i8)') '*** WARNING:  # fetches >95% max on CPU ',my_rank,':',max_local_f,'/',size_fetch 
     call cleanup
  endif
  if (npp > nppm) then
     write (6,'(a,i4)') '*** WARNING:  particle arrays full on CPU ',my_rank 
     write (6,'(a,i4)') '*** WARNING:  npp, nppm:',npp, nppm 
     call cleanup
  endif
end subroutine diagnostics





