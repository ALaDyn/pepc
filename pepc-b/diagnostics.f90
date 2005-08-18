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
  integer :: max_fetches, max_reqs, max_local_f, max_local_r

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

  if ((mod(itime+itime_start,idump)==0 .or. itime==nt) ) then
     call dump(itime+itime_start)     ! Dump complete set of particle data
     call dump_fields(itime+itime_start)  ! Field data
     if (vis_on)  call vis_fields
  endif

  if (debug_level.ge.2) then
     max_local_f =  maxval(nfetch_total)
     max_local_r =  maxval(nreqs_total)  
     call MPI_ALLREDUCE(max_local_f, max_fetches, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr )  
!     call MPI_ALLREDUCE(max_local_r, max_reqs, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr )  
! max_reqs, max_fetches should be equal
  endif

  if (my_rank.eq.0 .and. debug_level.ge.2) then

     do ifile = 6,15,9
        write(ifile,'(/a)') 'Tree stats:'
        write(ifile,'(a50,2i8,a3,i8,a1)') 'new npp, npart, (max): ',npp,npart,'(',nppm,')'
        write(ifile,'(a50,2i8)') 'local # leaves, twigs: ',nleaf_me,ntwig_me
        write(ifile,'(a50,3i8,a3,i8,a1)') 'final # leaves, twigs, keys, (size_tree): ',nleaf,ntwig,nleaf+ntwig,'(',size_tree,')'
        write(ifile,'(a50,2i8,a3,i8,a1)') 'local, global # branches, (max): ',nbranch,nbranch_sum,'(',nbranch_max,')'
        write (ifile,'(a50,i8,a3,i8)') 'Max length of all interaction lists: ',max_list_length,' / ',nintmax
        write (ifile,'(a50,i8)') 'Max # traversals ',maxtraverse

        write (ifile,'(a50,i8)') 'Max # multipole ships/iteration ',max_fetches
        write (ifile,'(a50,i8)') 'Max # multipole ships/prefetch ',max_prefetches
        write (ifile,'(a50,i8)') 'Array limit ',size_fetch
        write (ifile,*) ' cumulative # requested keys:  ',nreqs_total
        write (ifile,*) ' cumulative # fetched keys:    ',nfetch_total

     end do
  endif


end subroutine diagnostics





