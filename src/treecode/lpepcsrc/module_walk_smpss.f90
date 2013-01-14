! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2013 Juelich Supercomputing Centre, 
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

module module_walk

  use module_walk_smpss_utils, only: chunk_size_default
  implicit none

  ! IDs for internal timing measurement
  !integer, public, parameter :: TIMING_COMMLOOP = 1
  !integer, public, parameter :: TIMING_RECEIVE  = 2
  !integer, public, parameter :: TIMING_SENDREQS = 3
  
  namelist /walk_para_smpss/ chunk_size_default

  public tree_walk
  public tree_walk_finalize
  public tree_walk_prepare
  public tree_walk_statistics
  public tree_walk_read_parameters
  public tree_walk_write_parameters

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>
  !> writes walk-specific data to file steam ifile
  !>
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine tree_walk_statistics(ifile, perform_output)
    implicit none
    integer, intent(in) :: ifile !< file stream to write to
    logical, intent(in) :: perform_output !< if set to .false., output in this routine is prevented (e.g. for MPI ranks that shall not write any statistics)
    ! TODO: fill with life here
  end subroutine


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>
  !> reads walk specific parameters from file
  !>
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine tree_walk_read_parameters(filehandle)
    use module_debug, only: pepc_status
    implicit none
    integer, intent(in) :: filehandle

    call pepc_status("READ PARAMETERS, section walk_para_smpss")
    read(filehandle, NML=walk_para_smpss)

  end subroutine


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>
  !> writes walk specific parameters to file
  !>
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine tree_walk_write_parameters(filehandle)
    use module_debug, only: pepc_status
    implicit none
    integer, intent(in) :: filehandle

    write(filehandle, NML=walk_para_smpss)

  end subroutine


  subroutine tree_walk_finalize()
    implicit none
  end subroutine tree_walk_finalize

  subroutine tree_walk_prepare()
    use treevars, only : me
    implicit none

    if (me == 0) then
      write(*,'("MPI-MPSs walk")')
    endif
    
  end subroutine
  
  subroutine tree_walk(nparticles_, particles_, &
       twalk, twalk_loc_, vbox_, tcomm)
    use module_pepc_types
    use module_timings
    
    use module_walk_smpss_utils
    
    implicit none
    include 'mpif.h'
    
    integer, intent(in) :: nparticles_
    type(t_particle), target, intent(in) :: particles_(:)
    real*8, intent(in) :: vbox_(3) !< real space shift vector of box to be processed                           
    real*8, target, intent(inout) :: twalk, twalk_loc_
    real*8, target, intent(out), dimension(3) :: tcomm

    integer :: ierr
    integer :: dummy

    !!! chunk counter
    integer :: ccnt, pcnt, cpart, icnt

    !!! prefetch request list
    integer, parameter :: prefetch_req_size = 1000
    integer*8, dimension(:,:), allocatable :: prefetch_reqs

    interface

       !$CSS TASK
       subroutine tree_walk_smpss_walk_and_interact(particles, &
            status, requests, rlvl1, size) 
         
         implicit none
         integer,   intent(in)                     :: size
         integer,   intent(inout), dimension(size) :: particles
         integer*8, intent(inout), dimension(size) :: status, requests, rlvl1
         
       end subroutine tree_walk_smpss_walk_and_interact
       
       !$CSS TASK TARGET(COMM_THREAD) 
       !HIGHPRIORITY
       subroutine tree_walk_smpss_communicate(requests, rsize) 
         
         implicit none
         integer,   intent(in)                      :: rsize
         integer*8, intent(inout), dimension(rsize) :: requests
         
       end subroutine tree_walk_smpss_communicate

       !$CSS TASK
       subroutine tree_walk_smpss_local_finish(status, cs1, cs2)

         implicit none
         integer*8, intent(inout), dimension(cs1,cs2) :: status
         integer,   intent(in)                        :: cs1, cs2

       end subroutine tree_walk_smpss_local_finish

       !$CSS TASK 
       !TARGET(COMM_THREAD) 
       subroutine tree_walk_smpss_comm_serve(dummy)
         implicit none
         integer, intent(in) :: dummy
       end subroutine tree_walk_smpss_comm_serve

       !$CSS TASK HIGHPRIORITY
        subroutine tree_walk_smpss_walk_prefetch(x, requests, size)
         
          implicit none
          integer,   intent(in)                     :: size
          real*8,    intent(in)                     :: x(1:3)
          integer*8, intent(inout), dimension(size) :: requests
        end subroutine tree_walk_smpss_walk_prefetch

    end interface

    twalk = MPI_WTIME()

    !!!!! init constants in common module
    call tree_walk_smpss_setconst(nparticles_, particles_, vbox_)

    !!!!! setup chunks 
    call tree_walk_smpss_setup_chunks

    allocate(prefetch_reqs(prefetch_req_size, chunk_number))

    prefetch_reqs(:,:) = 0

    call MPI_BARRIER(MPI_COMM_lpepc, ierr)

    walk_status = WALK_STILL_RUNNING

    call tree_walk_smpss_comm_serve(dummy)

    do pcnt=1, 5

       do ccnt=1, chunk_number !chunk_number
          !cpart = MOD(pcnt, chunk_sizes(ccnt)) + 1
          cpart = pcnt
          call tree_walk_smpss_walk_prefetch(particles(chunk_particles(cpart,ccnt))%x, &
               prefetch_reqs(:,ccnt), prefetch_req_size)
          call tree_walk_smpss_communicate(prefetch_reqs(:,ccnt), prefetch_req_size)
       end do

    end do


    pcnt = 0
    icnt = 0
    do while (walk_status .ne. WALK_ALL_FINISHED)       

       icnt = icnt + 1

       pcnt = pcnt + 1

       !write(*,*) "start walk on rank", me

       do ccnt=1, chunk_number

          if(walk_status.eq.WALK_STILL_RUNNING) then
             call tree_walk_smpss_walk_and_interact(chunk_particles(:,ccnt), &
                  chunk_status(:,ccnt), chunk_requests(:,ccnt), chunk_rlvl1(:,ccnt), chunk_sizes(ccnt))
             
             !cpart = MOD(pcnt, chunk_sizes(ccnt)) + 1
             !call tree_walk_smpss_walk_prefetch(x(chunk_particles(cpart,ccnt)), y(chunk_particles(cpart,ccnt)), z(chunk_particles(cpart,ccnt)), prefetch_reqs(:,ccnt), prefetch_req_size)
             !call tree_walk_smpss_communicate(prefetch_reqs(:,ccnt), prefetch_req_size)
             
          end if
          
          call tree_walk_smpss_communicate(chunk_requests(:,ccnt), chunk_sizes(ccnt))

       end do


       if(MOD(icnt,10).eq.0) then
       do ccnt=1, chunk_number
          !$CSS WAIT ON(chunk_status(:,ccnt))
       end do
       call tree_walk_smpss_local_finish(chunk_status(:,:), chunk_size_default, chunk_number)
       end if

    end do

    !$CSS BARRIER

    !!!!! free chunk memory
    call tree_walk_smpss_free_chunks

    deallocate(prefetch_reqs)

    twalk = MPI_WTIME() - twalk

  end subroutine tree_walk

end module module_walk


!$CSS TASK
subroutine tree_walk_smpss_local_finish(status, cs1, cs2)
  
  use module_walk_smpss_utils

  implicit none
  integer*8, intent(inout), dimension(cs1,cs2) :: status
  integer,   intent(in)                        :: cs1, cs2

  integer :: c1, c2
  logical :: finished

  !if(wstatus(1).eq.WALK_STILL_RUNNING .and. all(status(:,:) .eq. -1)) wstatus(1) = WALK_IAM_FINISHED
  !if(wstatus(1).eq.WALK_STILL_RUNNING .and. all(status .eq. -1)) wstatus(1) = WALK_IAM_FINISHED

  !write(*,*) "evaluate local finish", cs1, cs2

  if(walk_status.eq.WALK_STILL_RUNNING) then

     finished = .true.
     do c1=1, cs1
        do c2=1, cs2
           if(status(c1, c2) .ne. -1_8) finished = .false.
        end do
     end do

     if(finished) walk_status = WALK_IAM_FINISHED

  end if
  
end subroutine tree_walk_smpss_local_finish


!$CSS TASK TARGET(COMM_THREAD) 
!HIGHPRIORITY
subroutine tree_walk_smpss_communicate(full_requests, full_size)
  
  use treevars, only: me
  use module_walk_smpss_utils
  
  implicit none
  
  integer,   intent(in)                          :: full_size
  integer*8, intent(inout), dimension(full_size) :: full_requests
  
!!! compressed request list
  integer*8, dimension(full_size) :: short_requests
  integer,   dimension(full_size) :: short_req_owner
  integer                         :: short_size
  
!!! global node exchange 
  integer, dimension(0:max_rank) :: send_req_number, send_node_number
  integer, dimension(0:max_rank) :: recv_req_number, recv_node_number
  
!!! counter
  integer :: rcnt, scnt
  
!!! tmp variables
  integer*8 :: tkey
  
!!! flags
  logical :: flog = .false.
  
  !if(flog) call tree_walk_smpss_print_keylist(full_requests, full_size, "full req list")
  
!!!!! compress the full request list
  call tree_walk_smpss_unique_keylist(full_requests, full_size, short_requests, short_req_owner, short_size)
  
!!!!! call tree_walk_smpss_gather_request_numbers(short_requests, short_size, recv_req_number, send_req_number)
  
  !if(flog) call tree_walk_smpss_print_keylist(short_requests, short_size, "short req list")
  
  
  if(walk_status .ne. WALK_ALL_FINISHED) call tree_walk_smpss_comm_loop_inner(short_requests, short_req_owner, short_size)
  
  full_requests = 0_8

  end subroutine tree_walk_smpss_communicate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!$CSS TASK
subroutine tree_walk_smpss_walk_and_interact(particle_list, particle_status, particle_requests, rlvl1, size)

  use treevars
  use module_walk_smpss_utils
  use module_htable
  use module_interaction_specific
  use module_spacefilling

  implicit none

  integer,   intent(in)                     :: size
  integer,   intent(inout), dimension(size) :: particle_list
  integer*8, intent(inout), dimension(size) :: particle_status, particle_requests, rlvl1

  integer :: pcnt, wcnt

!!! current (inside inner loop) values/copies of freq. used variables
  integer*8 :: creqs, cstat, cnext
  integer   :: cnode, caddr, cpart, clevl

!!! flags
  logical :: flog = .false.
  logical :: fmac
  logical :: fsamenode

!!! tmp variables
  integer*8 :: tkey
  real*8    :: tdist(3), tdist2

!!! tmp child variables
  integer*8 :: ch_addr(8)
  integer   :: ch_num

!!!!! loop over all particles in chunk
  do pcnt=1, size

     if(flog) write(*,*) "*** traverse particle ", pcnt

     creqs = 0
     cstat = particle_status(pcnt)
     cpart = particle_list(pcnt)

!!!!! traverse as long as possible (no request), or finished (status.eq.0)
     do while ((creqs .eq. 0) .and. (cstat .ge. 0))

        if(flog) write(*,'(a,2o21)') "** traverse status & request ", cstat, creqs

        caddr = key2addr(cstat, 'SMPSS-WALK')
        cnode = htable( caddr ) % node
        cnext = get_next_node_key( cstat )
        clevl = level_from_key( cstat )

!!!!! not myself in central box
        fsamenode = (fcentral .and. (cnode .eq. cpart))

        tdist  = particles(cpart)%x - (tree_nodes(cnode)%coc + vbox)
        tdist2 = tdist(1)**2 + tdist(2)**2 + tdist(3)**2 

        fmac = mac(particles(cpart), cnode, tdist2, boxlength2(clevl))

        if(flog) write(*,'(a,l2,i5,o21)') "** mac & cnode & cnext", fmac, cnode, cnext

!!!!! mac fits, or node is a leaf (node > 0), but not myself
        if ( (fmac .or. cnode.gt.0) .and. .not.fsamenode ) then

           if(flog) write(*,'(a)') "** compute force"

!!!!! compute force
           call calc_force_per_interaction(particles(cpart), tree_nodes(cnode), cstat, tdist, tdist2, vbox, cnode > 0)

!!!!! need to traverse deeper, cnode is a twig (node < 0)
        else if ( .not.fmac .and. cnode .lt. 0 ) then

!!!!! are the nodes children available? yes -> set first child as next node
           if ( children_available( caddr ) ) then
              if(flog) write(*,*) "** traverse, child available"
              !if(cnode .lt. -1) then
              !if(.true.) then
              call get_childkeys(caddr, ch_num, ch_addr)
              cstat = ch_addr(1)
              cnext = ch_addr(1)
              !else

!               do while(rlvl1(pcnt).lt.8 .and. (.not. btest(htable(caddr)%childcode, MOD(cpart+rlvl1(pcnt),8_8))))
!                  rlvl1(pcnt) = rlvl1(pcnt) + 1
!               end do
              
!               rlvl1(pcnt) = rlvl1(pcnt) + 1
              
!               cstat = 8 + MOD(cpart+rlvl1(pcnt),8_8)
!               cnext = cstat
              
!               if(rlvl1(pcnt).eq.8) cnext = -1

              !end if

!!!!! no -> stop traversal and put it on request list
           else
              if(flog) write(*,'(a,o21)') "** traverse, child NOT available -> request list entry: ", cstat
              cnext = cstat
              creqs = cstat
           end if

        end if

        cstat = cnext

!!!!! if next key is 1 (root) -> traversal is finished!
        if (cnext.eq.1) cstat = -1

     end do

!!!!! store back current status and request in chunk array
     particle_requests(pcnt) = creqs
     particle_status(pcnt)   = cstat

  end do

  !if(all(status(:).eq.-1))   write(*,*) "finised walk for 1. particle", particles(1), " on rank", me


end subroutine tree_walk_smpss_walk_and_interact

!$CSS TASK
subroutine tree_walk_smpss_walk_prefetch(x, particle_requests, size)

  use treevars
  use module_walk_smpss_utils
  use module_htable
  use module_interaction_specific
  use module_spacefilling

  implicit none

  integer,   intent(in)                     :: size
  real*8,    intent(in)                     :: x(1:3)
  integer*8, intent(inout), dimension(size) :: particle_requests

  integer :: pcnt, wcnt

!!! current (inside inner loop) values/copies of freq. used variables
  integer*8 :: creqs, cstat, cnext
  integer   :: cnode, caddr, cpart, clevl

!!! flags
  logical :: flog = .false.
  logical :: fmac
  logical :: fsamenode

!!! tmp variables
  integer*8 :: tkey
  real*8    :: tdist(3), tdist2

!!! tmp child variables
  integer*8 :: ch_addr(8)
  integer   :: ch_num

  cstat = 8 + int(rand()*8)

!!!!! loop over all particles in chunk
  do pcnt=1, size-1

     !if(flog) write(*,'(a,i3,2o21)') "*** traverse prefetch particle ", pcnt, cstat, cnext

     creqs = 0

!!!!! traverse as long as possible (no request), or finished (status.eq.0)
     do while ((creqs .eq. 0) .and. (cstat .ge. 0))

        caddr = key2addr(cstat, 'SMPSS-WALK-PREFETCH')
        cnode = htable( caddr ) % node
        cnext = get_next_node_key( cstat )
        clevl = level_from_key( cstat )

!!!!! not myself in central box
        fsamenode = (fcentral .and. (cnode .eq. cpart))

        tdist = x - (tree_nodes(cnode)%coc + vbox(1))
        tdist2 = tdist(1)**2 + tdist(2)**2 + tdist(3)**2 

        fmac = (theta2 * tdist2 > boxlength2(clevl))

!!!!! mac fits, or node is a leaf (node > 0), but not myself
        if ( (fmac .or. cnode.gt.0) .and. .not.fsamenode ) then
           
           !!!!! do nothing
           fmac = .true.

!!!!! need to traverse deeper, cnode is a twig (node < 0)
        else if ( .not.fmac .and. cnode .lt. 0 ) then

!!!!! are the nodes children available? yes -> set first child as next node
           if ( children_available( caddr ) ) then
              
              call get_childkeys(caddr, ch_num, ch_addr)
              cstat = ch_addr(1)
              cnext = ch_addr(1)
              if(flog) write(*,'(a,2o21,i5)') "** traverse, child available: ", cstat, ch_addr(1), ch_num

!!!!! no -> stop traversal and put it on request list
           else
              if(flog) write(*,'(a,2o21)') "** traverse, child NOT available -> request list entry: ", cstat, cnext
              creqs = cstat
           end if

        end if

        cstat = cnext
        if(flog) write(*,'(a,o21)') "** next node: ", cnext

!!!!! if next key is 1 (root) -> traversal is finished!
        if (cnext.eq.1) cstat = -1

     end do

!!!!! store back current status and request in chunk array
     particle_requests(pcnt) = creqs

  end do

  !write(*,*) "prefetch list: ", requests

end subroutine tree_walk_smpss_walk_prefetch

!$CSS TASK 
!TARGET(COMM_THREAD) 
subroutine tree_walk_smpss_comm_serve(dummy)

  use module_walk_smpss_utils

  implicit none

  integer, intent(in) :: dummy

  integer :: nc = 0
  integer, parameter :: nc_max = 10000
  save nc

  nc = nc + 1

  call tree_walk_smpss_comm_serve_inner()

  if((walk_status .eq. WALK_STILL_RUNNING) .and. &
       (nc .lt. nc_max) ) then
     !call css_fortran_restart_
  else 
     nc = 0
  end if
  
end subroutine tree_walk_smpss_comm_serve

