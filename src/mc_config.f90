! ==============================================
!
!                MC_CONFIG
!
!  Perform Monte-Carlo minimisation of PE by adjusting particle positions
!  according to constraints given by system geometry
!
! ==============================================

subroutine mc_config

  use treevars
  use utils
  implicit none

  integer :: nmove,i, ipar, j, pe_move
  integer :: iseed0, i_count, i_rate, i_max, mc_dump
  real :: epot, etrial, delta_E, prob1, epold, r 
  real :: xold, yold, zold, xt, yt, zt, xs, ys, zs
  real :: Qplas, Tplas
  real :: r_limit, x_limit, y_limit
  logical :: constrained


  !  idump = max(1,mc_steps/100)
  mc_dump = 10
  !  epot is p.e. of current configuration
  !  etrial is p.e. of new trial config

  call make_domains    ! Domain decomposition: allocate particle keys to PEs
  call tree_build      ! Build trees from local particle lists
  call make_branches   ! Determine and concatenate branch nodes
  call tree_fill       ! Fill in remainder of local tree
  call tree_properties ! Compute multipole moments for local tree
  call potenergy(epot) ! Compute potential energy
  call vis_parts


  if (me==0) then
     do i=6,15,9
        write (i,*) 'Configuring particle positions with MC routine ...'
        write (i,*) 'Initial p.e.=',epot
     end do
     open(76,file='energy_mc.dat')
     write(76,'(i6,1pe14.5)') 0,epot   ! Initial PE
  endif

     Qplas = abs(qe)*ne  ! Norm factor for plasma kinetic energy
     Tplas = Te_keV/511.*Qplas  ! convert temperature to code units

  call system_clock(i_count,i_rate,i_max)
  iseed0= -i_count/10 + 100*me  ! random number seed

  !  maximum step size as fraction of inter-particle spacing
  !      delta_l=0.5*min(aii,aee)
  write(ipefile,*) 'MC step length ',delta_mc


  do j=1, mc_steps
!     call MPI_BARRIER( MPI_COMM_WORLD, ierr)  ! Wait for active PE to catch up

     epold=epot
     if (me==0) then 
	! root picks PE for move
	pe_move = num_pe*rano(iseed0)
	pe_move = max(0,min(num_pe,pe_move))
     endif

     call MPI_BCAST(pe_move,one,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)  ! Pass result to other PEs

     if (pe_move == me) then
	i = npp*rano(iseed0)+1  ! select particle
	i = max(1,min(npp,i))
	!  keep old position
	xold=x(i)
	yold=y(i)
	zold=z(i)
	!  shift particle to new position
	constrained = .false.
	do while (.not. constrained) 
	   xs = (rano(iseed0)*2.-1.)*delta_mc
	   ys = (rano(iseed0)*2.-1.)*delta_mc
	   zs = (rano(iseed0)*2.-1.)*delta_mc

	   x(i) = xold + xs
	   y(i) = yold + ys
	   z(i) = zold + zs
	   xt = x(i) - plasma_centre(1)
	   yt = y(i) - plasma_centre(2)
	   zt = z(i) - plasma_centre(3)
	   if (q(i) <0) then 
	      r_limit = r_sphere +10*vte   ! allow for Debye sheath
	      x_limit = x_plasma/2. + 10*vte
	      y_limit = y_plasma/2. + 10*vte
	   else
	      r_limit = r_sphere
	      x_limit = x_plasma/2.
	      y_limit = y_plasma/2.
	   endif

	   if (initial_config==1) then
	      ! sphere
	      constrained = (xt**2 + yt**2 + zt**2 <= r_limit**2) 
	   else if (initial_config==2) then
	      ! disc
	      constrained = (xt >= -x_limit .and. xt<= x_limit .and. yt**2 + zt**2 <= r_limit**2) 
	   else if (initial_config==3) then
	      ! wire
	      constrained = (zt >= -x_limit .and. zt<= x_limit .and. yt**2 + xt**2 <= r_limit**2)   
	   else
	      ! slab
	      constrained = (     xt >= -x_limit .and. xt <= x_limit &
		   .and. yt >= -y_limit .and. yt <= y_limit &
		   .and. zt >= -y_limit .and. zt <= y_limit )
	   endif
	end do


	!  check for wrapping
	!   call boundry(i,i,x,y,z,ux,uy,uz)
     else
	! do nothing
     endif
     call MPI_BARRIER( MPI_COMM_WORLD, ierr)  ! Wait for active PE to catch up


     call tree_properties   ! rebuild multipole moments with same tree
     call potenergy(etrial) ! Estimate new potential energy


     delta_E = etrial - epold

     prob1=exp(-amin1(delta_E/Tplas, 30.))
     r = rano(iseed0)
     call MPI_BCAST(r,one,MPI_REAL8,root,MPI_COMM_WORLD,ierr)  ! Pass random no. to other PEs

     if (Te_keV == 0.) prob1=0.

     if ( delta_E < 0 .or.  r < prob1 ) then
        ! accept move either i) if it decreases pot. energy 
        !               or   ii) if it increases with random prob exp(-delta_E/kT)
        call make_domains    ! rebuild tree from scratch
        call tree_build     
        call make_branches   
        call tree_fill  
!	call tree_properties   ! rebuild multipole moments 
!	call potenergy(epot) ! Compute new (accurate) potential energy
        epot = etrial

     else
        !             epot unchanged
        ! restore particle to old position - tree unchanged
	if (pe_move == me) then
	   x(i)=xold
	   y(i)=yold
	   z(i)=zold
	endif

        epot = epold   ! Keep old tree and P.E. value
     endif

     !  distributions
     if (mod(j,mc_dump).eq.0 .and. me==pe_move) then
	do ifile=6,15,9
           write(ifile,'(a5,i5,a,i4,3(a8,1pe12.3))') &
                'Step ',j,' particle ',i, &
		' P.E.=',epot/Qplas*511.,' dE ', delta_E/Qplas*511.,' prob ',prob1	
	enddo
     endif
     if (vis_on .and. mod(j,mc_steps/10) == 0 )  call vis_parts

     if (mod(j,mc_dump)==0) then
	call system_clock(i_count,i_rate,i_max)
	iseed0= -i_count/10 + 100*me  ! generate new sequence
     endif

     write(76,'(i6,1pe14.5)') j,epot
  end do

  close(76) 

  call dump(mc_steps)
end subroutine mc_config
