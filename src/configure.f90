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
  integer :: i, ipe, idummy=0
  real :: t_walk, t_force

  if (restart) then
     call predef_parts    ! Predefined particle properties read from peXX/parts_dump.NNNNNN
!  Find critical/shelf densities

     if (beam_config ==4 )call track_nc 

  else
     config: select case(initial_config)

     case(:8)
        call randion         ! Set up particles according to geometry

        if (scheme<>5 .and. ramp) then
           call add_ramp     ! add exponential ramp
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
     
     case(10)
        call special_start(ispecial)

     case default
        call randion         ! Random by default
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
  end select beamconf

  ! Do tree-build for initial P.E. value

  if (me==0) then
     do ifile = 6,15,9
        write(ifile,'(//a,i8,a,f10.5)') 'Timestep ',itime,' t=',itime*dt
     end do
  endif

  ! Initial tree construction and force computation

  call make_domains(xl,yl,zl)    ! Domain decomposition: allocate particle keys to PEs
  call tree_build      ! Build trees from local particle lists
  call make_branches   ! Determine and concatenate branch nodes
!  if (me>=250) call  diagnose_tree
!  close(ipefile)
  call tree_fill       ! Fill in remainder of local tree

  call tree_properties ! Compute multipole moments for local tree
!    call MPI_FINALIZE(ierr)
!    call closefiles
!   stop
  if (coulomb .or. lenjones) then
     call forces(1,npp,dt,t_walk,t_force)          ! Calculate initial potentials and forces
  endif
  if (.not. perf_anal) call diagnostics
end subroutine configure


