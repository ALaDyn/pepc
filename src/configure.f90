! =============================================
!
!                CONFIGURE
!
!  Sets up physical system: particle positions, velocities
!
! ==============================================

subroutine configure

  use physvars
  use treevars
  use utils
  implicit none
  integer :: i, ipe, idummy=0, ierr
  real :: t_walk, t_force

  if (restart) then

     call predef_parts    ! Predefined particle properties read from peXX/parts_dump.NNNNNN
     !  Find critical/shelf densities if laser switched on

     if (beam_config ==4 ) call track_nc 


  else
     config: select case(plasma_config)

     case(1)              ! Set up particles according to geometry
        call randion         

        if (scheme /= 5 .and. ramp) then
           call add_ramp     ! add exponential ramp to target
        endif

        if (vte > 0) then
           call maxwell1(ux,nppm,1,nep,vte)
           call maxwell1(uy,nppm,1,nep,vte)
           call maxwell1(uz,nppm,1,nep,vte)
           call scramble_v(1,nep)   ! remove x,y,z correlations
        else
           call cold_start(1,nep)
        endif

        if (vti > 0) then
           call maxwell1(ux,nppm,nep+1,nip,vti)
           call maxwell1(uy,nppm,nep+1,nip,vti)
           call maxwell1(uz,nppm,nep+1,nip,vti)
           call scramble_v(nep+1,nip) ! remove x,y,z correlations
        else
           call cold_start(nep+1,nip)
        endif
     
     case(2)
        call special_start(ispecial)

     case default     ! Default = 0 - no plasma target
        if (me==0) write (6,*) 'Warning: no plasma set up'
        npart=0
        npp = 0
     end select config

  endif

  if (target_dup) then
     if (me==0) write(*,*) 'DOUBLE-TARGET CONFIGURATION' 
     call double_target   ! Special target-duplication mode
  endif

  if (mc_init) call mc_config  ! Do MC min-PE initialisation depending on config

  beamconf: select case(beam_config)

  case(1)
     call beam           ! Fixed beam
     call beam_control   ! Display default parameters

  case(2)
     call beam_control   ! Constant particle source

  case(3)
     call beam_dust   ! Dust particle

  case(4:) ! laser on
     if (me==0) write(*,*) 'Laser switched on' 
 
  end select beamconf

  ! Do tree-build for initial P.E. value

  if (me==0) then
     do ifile = 6,15,9
        write(ifile,'(//a,i8,a,f10.5)') 'Timestep ',itime,' t=',itime*dt
     end do
  endif

  ! Initial tree construction and force computation

  call make_domains(xl,yl,zl)    ! Domain decomposition: allocate particle keys to PEs
  call MPI_BARRIER( MPI_COMM_WORLD, ierr)  ! Wait for everyone to catch up
  call tree_build      ! Build trees from local particle lists
  call make_branches   ! Determine and concatenate branch nodes
  call tree_fill       ! Fill in remainder of local tree
  call tree_properties ! Compute multipole moments for local tree

!  call diagnose_tree
!    call MPI_FINALIZE(ierr)
!    call closefiles
!   stop
  if (coulomb .or. lenjones) then
     call forces(1,npp,dt,t_walk,t_force)          ! Calculate initial potentials and forces
  endif
  if (.not. perf_anal) call diagnostics
end subroutine configure


