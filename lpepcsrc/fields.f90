!  ===================================================================
!
!                              FIELDS
!
!
!   Calculate fields and potential from coordinates x,y,z:
!
!
!   ** Returns fields Ex, Ey, Ez and potential pot excluding external terms **
!
!
!  ===================================================================


subroutine pepc_fields(np_local,npart_total,nppm_ori,p_x, p_y, p_z, p_q, p_m, p_w, p_label, &
     p_Ex, p_Ey, p_Ez, p_pot, np_mult_,&
     mac, theta, eps, force_const, itime, choose_sort,weighted, &
     num_neighbours, neighbours)

  use treevars
  use timings
  use module_calc_force
  use module_fmm_framework
  use tree_walk_utils
  use tree_walk_communicator
  implicit none
  include 'mpif.h'

  integer, intent(in) :: np_local, npart_total, nppm_ori  ! # particles on this CPU
  real, intent(in) :: theta, np_mult_       ! multipole opening angle
  real, intent(in) :: force_const       ! scaling factor for fields & potential
  real, intent(in) :: eps         ! potential softening distance
  integer, intent(in) :: itime  ! timestep
  integer, intent(in) :: mac, choose_sort, weighted
  real*8, intent(in), dimension(np_local) :: p_x, p_y, p_z  ! coords and velocities: x1,x2,x3, y1,y2,y3, etc 
  real*8, intent(in), dimension(np_local) :: p_q, p_m ! charges, masses
  integer, intent(in), dimension(np_local) :: p_label  ! particle label 
  real*8, intent(out), dimension(np_local) :: p_ex, p_ey, p_ez, p_pot  ! fields and potential to return
  integer, intent(in) :: num_neighbours !< number of shift vectors in neighbours list (must be at least 1 since [0, 0, 0] has to be inside the list)
  integer, intent(in) :: neighbours(3, num_neighbours) !< list with shift vectors to neighbour boxes that shall be included in interaction calculation, at least [0, 0, 0] should be inside this list
  real*8, dimension(np_local) :: p_w ! work loads
  
  integer :: npnew, npold

  integer :: indxl(nppm_ori),irnkl(nppm_ori)
  integer :: islen(num_pe),irlen(num_pe)
  integer :: fposts(num_pe+1),gposts(num_pe+1)

  integer :: i
  real*8 :: ttrav, tfetch, tcomm(3) ! timing integrals

  integer :: timestamp

  integer :: ibox
  real*8 :: vbox(3)
  character(30) :: cfile

  allocate(ex_tmp(nppm_ori), ey_tmp(nppm_ori), ez_tmp(nppm_ori), pot_tmp(nppm_ori), w_tmp(nppm_ori))

  call timer_start(t_all)
  call timer_start(t_fields_begin)

  np_mult = np_mult_
  npart   = npart_total

  npp = np_local

  if (force_debug) then
     write (*,'(a7,a50/2i5,4f15.2)') 'PEPC | ','Params: itime, mac, theta, eps, force_const:', &
		itime, mac, theta, eps, force_const
     write (*,'(a7,a20/(i16,4f15.3,i8))') 'PEPC | ','Initial buffers: ',(p_label(i), p_x(i), p_y(i), p_z(i), p_q(i), &
		p_label(i),i=1,npp) 
  endif

 ! Copy particle buffers to tree arrays
  do i=1,npp
     x(i) = p_x(i)
     y(i) = p_y(i)
     z(i) = p_z(i)

     ux(i) = 0.  ! No B-fields for now
     uy(i) = 0.
     uz(i) = 0.
     q(i) = p_q(i)
     m(i) = p_m(i)
!     if (p_label(i) <=0) then
!        ! Trap bad particle labels
!        write (*,*) '*** Error: particle labels must be positive integers (1,2,3,...)! '
!        write (*,*) p_label(1:20)
!        call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
!        stop
!     else
        pelabel(i) = p_label(i)     
!     endif
     pepid(i) = me
     if (num_pe==1 .or. p_w(i)<1.) then
        work(i) = 1.
     else
       work(i) = p_w(i)
     endif
     ax(i) = 0.
     ay(i) = 0.
     az(i) = 0.
  end do

  call timer_stop(t_fields_begin)
  call timer_start(t_fields_tree)

  ! Domain decomposition: allocate particle keys to PEs
  call tree_domains(indxl,irnkl,islen,irlen,fposts,gposts,npnew,npold, choose_sort, weighted)
  call tree_allocate(theta)

  ! calculate spherical multipole expansion of central box
  call fmm_framework_timestep

  ! build local part of tree
  call timer_stamp(t_stamp_before_local)
  call tree_local
  ! exchange branch nodes
  call timer_stamp(t_stamp_before_exchange)
  call tree_exchange
  ! build global part of tree
  call timer_stamp(t_stamp_before_global)
  call tree_global

  call timer_stop(t_fields_tree)
  call timer_reset(t_walk)
  call timer_reset(t_comm_total)
  call timer_reset(t_comm_recv)
  call timer_reset(t_comm_sendreqs)

  max_req_list_length = 0
  cum_req_list_length = 0
  comm_loop_iterations = 0
  sum_fetches=0      ! total # multipole fetches/iteration
  sum_ships=0      ! total # multipole shipments/iteration

  call timer_start(t_fields_passes)

  pot_tmp = 0.
  ex_tmp  = 0.
  ey_tmp  = 0.
  ez_tmp  = 0.
  work    = 1.

  call timer_stamp(t_stamp_before_walkloop)

  do ibox = 1,num_neighbours ! sum over all boxes within ws=1

    vbox = lattice_vect(neighbours(:,ibox))

    ! tree walk finds interaction partners and calls interaction routine for particles on short list
    call tree_walk(npp,theta,eps,force_const,itime,mac,ttrav,tfetch, vbox, work, tcomm)

    call timer_add(t_walk, ttrav)    ! traversal time (serial)
    call timer_add(t_comm_total,    tcomm(TIMING_COMMLOOP))
    call timer_add(t_comm_recv,     tcomm(TIMING_RECEIVE))
    call timer_add(t_comm_sendreqs, tcomm(TIMING_SENDREQS))

  end do ! ibox = 1,num_neighbours

  call timer_stamp(t_stamp_after_walkloop)

  w_tmp(1:npp) = work(1:npp)  ! send back work load for next iteration

  ! add lattice contribution
  call timer_start(t_lattice)
  call calc_force_per_particle(eps, force_const)
  call timer_stop(t_lattice)

  ! restore initial particle order specified by calling routine to reassign computed forces
  ! notice the swapped order of the index-fields -> less changes in restore.f90 compared to tree_domains.f90

  call timer_stop(t_fields_passes)
  call timer_start(t_restore)

  call restore(npnew,npold,nppm_ori,irnkl,indxl,irlen,islen,gposts,fposts, &
       pot_tmp(1:npnew),ex_tmp(1:npnew),ey_tmp(1:npnew),ez_tmp(1:npnew),w_tmp(1:npnew),p_pot,p_ex,p_ey,p_ez,p_w)    

  call timer_stop(t_restore)
  call timer_start(t_fields_stats)

  timestamp = itime

  nkeys_total = nleaf+ntwig  

  if (force_debug) then
     write (ipefile,101)
     write (*,101)
     write (ipefile,101) force_const

     do i=1,npp
        write (ipefile,102) pelabel(i), & 
             q(i), m(i), ux(i), p_pot(i), p_ex(i)
        write (*,102) pelabel(i), x(i), & 
             q(i), m(i), ux(i), p_pot(i)
     end do

101  format('Tree forces:'/'   p    q   m   ux   pot  ',f8.2)
102  format(1x,i7,5(1pe14.5))

  endif

  call tree_stats(itime)

  call timer_stop(t_fields_stats)

  call timer_start(t_deallocate)
  call tree_deallocate(nppm_ori)
  call timer_stop(t_deallocate)

  call timer_stop(t_all)

  write(cfile,'(a,i6.6,a)') "load_", me, ".dat"  
  open(60, file=cfile,STATUS='UNKNOWN', POSITION = 'APPEND')
  write(60,'(i5,2f20.10, i12)') itime,interactions_local, mac_evaluations_local,npp
  close(60)   

  deallocate(ex_tmp, ey_tmp, ez_tmp, pot_tmp, w_tmp)

end subroutine pepc_fields









