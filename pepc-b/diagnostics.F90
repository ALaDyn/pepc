!  ================================
!
!         DIAGNOSTICS
!
!     Perform diagnostics
!
!
!  ================================


subroutine diagnostics

  use treevars
  use physvars
  implicit none
  include 'mpif.h'

  integer :: i,lvisit_active, ifile, ierr
  integer :: max_fetches, max_reqs, max_sum_fetches, max_sum_ships,max_local_f, max_local_r, sum_local_f
  integer :: max_nbranch, gmax_leaves, gmax_twigs
  integer :: vcount=0

! Tree diagnostics
! If interaction lists needed, must ensure that intlist() is large enough to contain all lists
! - will otherwise just get last pass of tree walk

  if ( dump_tree .and. mod(itime,iprot) ==0 ) then
     call diagnose_tree   ! Printed tree info (htable etc)
     call draw_tree2d(xl,yl)     ! Draw PE-trees
     call draw_lists      ! Draw interaction lists
     call draw_domains(itime+itime_start)   ! Domains
  endif


 if ( mod(itime,ivis_fields)==0 ) then
   call densities
   call sum_fields
 endif 


  ! Routines for particle and field processing
  ! for VISIT (Online visualisation) and/or Netcdf file

  if ( vis_on ) then
     !     if ( mod(itime,ivis) ==0 ) call vis_parts       
     if ( mod(itime,ivis) ==0 ) then
        call vis_parts_nbody(vcount)
        vcount = vcount + 1
     endif
     if ( mod(itime,min(ivis,ivis_fields))==0 .and. steering) call vis_control
     if ( mod(itime,ivis_fields)==0 ) then
        !     call pot_grid
        call vis_fields_nbody(itime+itime_start)
        call vis_vecfields_nbody(itime+itime_start)
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

!  if (mod(itime,idump) >= idump-navcycle .or. nt.lt.idump) then
!	call sum_fields    ! Accumulate cycle-averaged fields on grid
!	call sum_fieldave
!  endif

  !  - assume for now that idump > navcycle
  !  - should really use running average to be compatible with online vis.


  if (beam_config == 4 .and. mod(itime,itrack)==0 ) call track_nc          ! Gather densities and track critical surface 
  if( mod(itime+itime_start,ivis_fields)==0 .and. mod(target_geometry,10)==1) then
!     call sum_radial(itime+itime_start)  ! Radial moments
  endif


  if (itime_start>0 .and. itime==0) return  ! Avoid over-writing restart data

  call energy_cons(Ukine,Ukini,Umagnetic,Ubeam)       ! Compute energy balance
  call laser_hist

  if ( idump>0 .and. (mod(itime+itime_start,idump)==0 .or. itime==nt) ) then
     call dump(itime+itime_start)     ! Dump complete set of particle data
     call dump_fields(itime+itime_start)  ! Field data
  endif


!  if (debug_level.ge.2) then
     sum_local_f =  sum(nfetch_total)  
     call MPI_ALLREDUCE(sum_fetches, max_sum_fetches, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr )  
     call MPI_ALLREDUCE(sum_ships, max_sum_ships, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr )  
     call MPI_ALLREDUCE(sum_local_f, max_fetches, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr )  
     call MPI_ALLREDUCE(nbranch, max_nbranch, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr )  
  nnodes=nleaf+ntwig
  call MPI_ALLREDUCE(nleaf, gmax_leaves, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr )  
  call MPI_ALLREDUCE(ntwig, gmax_twigs, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr )  
! max_reqs, max_fetches should be equal
!  endif

  if ((my_rank.eq.debug_rank .and. debug_tree > 0) .or. (nleaf+ntwig > .98*maxaddress)) then

     do ifile = 6,15,9
        write(ifile,'(/a,i4)') 'Tree stats for CPU ',debug_rank
        write(ifile,'(a50,2i10,a3,i8)') 'new npp, npart, (max): ',npp,npart,'/',nppm
        write(ifile,'(a50,3i8)') 'local # leaves, twigs, keys: ',nleaf_me,ntwig_me,nleaf_me+ntwig_me
        write(ifile,'(a50,3i8)') 'non-local # leaves, twigs, keys: ',nleaf-nleaf_me,ntwig-ntwig_me,nleaf+ntwig-nleaf_me-ntwig_me
        write(ifile,'(a50,3i8,f12.1,a6,i8)') 'final # leaves, twigs, keys, (max): ',nleaf,ntwig,nleaf+ntwig, &
             (nleaf+ntwig)/(.01*maxaddress),' % of ',maxaddress
        write(ifile,'(a50,2i8,a6,2f12.1)') 'Global max # leaves, twigs: ',gmax_leaves, gmax_twigs, &
             ' %s: ',gmax_leaves/(.01*maxleaf),gmax_twigs/(.01*maxtwig)
        write(ifile,'(a50,3i8,a3,2i7)') 'Local, max local, global # branches, (max): ',nbranch,max_nbranch,nbranch_sum,'/',nbranch_local_max,nbranch_max
        write (ifile,'(a50,i8,a3,i7)') 'Max length of all interaction lists: ',max_list_length,' / ',nintmax
        write (ifile,'(a50,i8)') 'Max # traversals ',maxtraverse

        write (ifile,'(a50,i8)') 'Max # multipole fetches/cpu/walk ',max_fetches-max_prefetches
        write (ifile,'(a50,i8)') 'Max # multipole prefetches/cpu/prefetch ',max_prefetches
        write (ifile,'(a50,2i8)') 'Local #  multipole fetches & ships/iteration ',sum_fetches,sum_ships
        if (debug_tree.ge.3) write (ifile,*) ' cumulative # requested keys:  ',nreqs_total
        if (debug_tree.ge.3) write (ifile,*) ' cumulative # fetched keys:    ',nfetch_total
        write (ifile,'(a50,i8,a3,i7)') 'Global max # multipole fetches/cpu/iteration / limit',max_sum_fetches,'/',size_fetch
        write (ifile,'(a50,i8,a3,i7)') 'Global max # multipole ships/cpu/iteration / limit',max_sum_ships,'/',size_fetch
	write (ifile,'(a50,3f12.3)') 'Load imbalance percent,min,max: ',work_imbal,work_imbal_min,work_imbal_max
	write (ifile,'(a50,f15.3,2i15)') 'Particle imbalance ave,min,max: ',part_imbal,part_imbal_min,part_imbal_max
     end do
  endif

  if (part_imbal>5 .or. work_imbal_min < 0.1) then
   if (my_rank.eq.debug_rank) write(ifile,'(a)') 'WARNING: Load imbalance too large - switching to debug mode'
   debug_tree = 2
  endif

  if (debug_tree.ge.1 .and. mod(itime,iprot)==0 ) then
    call tree_stats(itime+itime_start)
  endif

! Array bound check
  if (nleaf+ntwig > .98*maxaddress) then
     write (6,'(a,i4)') '*** WARNING:  hash table >98% full on CPU ',my_rank 
  endif


!  if (max_local_f > .95*size_fetch) then
!     write (6,'(a,i4,a1,i8,a3,i8)') '*** WARNING:  # fetches >95% max on CPU ',my_rank,':',max_local_f,'/',size_fetch 
!     call cleanup
!  endif
  if (npp > nppm) then
     write (6,'(a,i4)') '*** WARNING:  particle arrays full on CPU ',my_rank 
     write (6,'(a,i4)') '*** WARNING:  npp, nppm:',npp, nppm 
     call MPI_ABORT(MPI_COMM_WORLD,ierr)
     stop 
  endif
end subroutine diagnostics





