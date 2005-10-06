! =============================================
!
!                CONFIGURE
!
!  Sets up physical system: particle positions, velocities
!
! ==============================================

subroutine configure

  use physvars   ! Use internal particle arrays from lpepc
  use treevars
  implicit none
  include 'mpif.h'

  integer :: i, ipe, idummy=0, ierr, ifile
  real :: t_walk, t_walkc, t_force, t_domain,t_build,t_prefetch

  npp = nep+nip  ! initial # particles on this CPU

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
        if (me==0) then
           call special_start(ispecial)
        else
           npp=0
        endif

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



  beamconf: select case(mod(beam_config,10))  ! Configure laser or particle beam

  case(1)
     call beam           ! Fixed beam
     if (steering) call beam_control   ! Display default parameters

  case(2)
     if (steering) call beam_control   ! Constant particle source
     if (me==0) write(*,'(//a)') '===> Particle beam switched on' 

  case(8)
     call beam_dust   ! Dust particle

  case(3:6) ! laser on

     if (me==0) write(*,'(//a)') '===> Laser switched on' 
     if (steering) call beam_control 

  end select beamconf

! local and global # particles for  main routine
  np_local = npp
  npart_total = npart

! Compute initial field values - need these to get vec. pots consistent with velocities

  if (me==0) write(*,*) 'Computing initial fields'
  call pepc_fields_p(np_local, mac, theta, ifreeze, eps, force_tolerance, balance, force_const, bond_const, &
          dt, xl, yl, zl, 0, &
          coulomb, bfields, bonds, lenjones, &
          t_domain,t_build,t_prefetch,t_walk,t_walkc,t_force, iprot)   

!  Initialise vec. pots. to avoid jump in induced E-field
  Axo(1:npp) = Ax(1:npp)
  Ayo(1:npp) = Ay(1:npp)
  Azo(1:npp) = Az(1:npp)

end subroutine configure

