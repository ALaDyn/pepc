!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>
!>  Encapsulates diagnostic routines
!>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module module_diagnostics
 
      implicit none
      save
      private

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  public subroutine declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	public param_dump	!< Write out inputs to run protocol
	public energy_cons	!< Energy conservation I/O	
	public kinenergy	!< Calculate kinetic energies
	public potenergy	!< Calculate electrostatic and magnetic energies
	public error_test	!< Estimate force errors
	public force_direct	!< PP force-sum


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!  private variable declarations  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      contains


subroutine param_dump

  use module_physvars
  use module_utilities
!  use tree_walk_pthreads
  implicit none
  integer :: ifile


  if (my_rank==0) then
     do ifile = 6,24,18

        write (ifile,'(a20,i4,a11)') ' Plasma config: ',plasma_config,plasma_configs(plasma_config)
        write (ifile,'(a20,i4,a11)') ' Target geometry: ',target_geometry,geometries(mod(target_geometry,10))
        write (ifile,'(a20,i4,a8)') ' Beam config: ',beam_config,beam_configs(beam_config)
        write (ifile,'(a20,i4)') ' type: ',beam_config_in
        write (ifile,'(a20,i4,a8)') ' Scheme: ',scheme,schemes(scheme)
        write (ifile,'(a20,1pe12.3)') ' Plasma volume: ',Vplas
        write (ifile,'(a20,f12.3)') ' Sphere radius: ',r_sphere
        write (ifile,'(a20,f12.3)') ' Plasma length: ',x_plasma
        write (ifile,'(a20,f12.3)') ' Plasma width: ',y_plasma
        write (ifile,'(a20,f12.3)') ' Plasma height: ',z_plasma
        write (ifile,'(a20,1pe12.3)') ' Electron charge: ',qe
        write (ifile,'(a20,1pe12.3)') ' Electron mass: ',mass_e
        write (ifile,'(a20,1pe12.3)') ' Ion Z:',Zion
        write (ifile,'(a20,1pe12.3)') ' Ion m_i/m_p:',mass_ratio
        write (ifile,'(a20,1pe12.3)') ' Ion charge: ',qi
        write (ifile,'(a20,1pe12.3)') ' Ion mass: ',mass_i
        write (ifile,'(a20,f12.3)') ' Te: ',Te_keV
        write (ifile,'(a20,f12.3)') ' Ti: ',Ti_keV
        write (ifile,'(a20,1pe12.3)') ' Debye length: ',vte
        write (ifile,'(a20,f12.3)') ' n_e/n_c: ',1./omega**2
        if (ramp) then
           write (ifile,'(a20,f12.3)') ' n_min: ', rho_min
           write (ifile,'(a20,f12.3)') ' k_p L: ', lolam
           write (ifile,'(a20,f12.3)') ' x_sol: ', plasma_centre(1)-x_plasma/2.+lolam*(1-rho_min)

        endif

        write (ifile,'(a20,f12.4)') ' Ion spacing a: ',a_ii
        write (ifile,'(a20,f12.3)') ' Cloud radius R: ',eps
        
        write (ifile,'(a20,f12.3)') ' Collision freq.: ',nu_ei
        write (ifile,'(a20,f12.3)') ' Electron mfp: ',vte/nu_ei
        write (ifile,'(a20,f12.3)') ' Conductivity: ',sigma_e
        write (ifile,'(a20,f12.3)') ' Laser amplitude: ',vosc
        write (ifile,'(a20,f12.3)') ' Pulse width: ',sigma
        write (ifile,'(a20,f12.3)') ' Pulse duration: ',tpulse
        write (ifile,'(a20,f12.3)') ' Wavelength: ',lambda
        write (ifile,'(a20,f12.3)') ' Incidence angle: ',theta_beam
        write (ifile,'(a20,f12.3)') ' Laser intensity: ',intensity


        write (ifile,'(a20,i12)') ' # dt per cycle: ',navcycle

        write (ifile,'(a20,1pe12.3)') ' N_D: ',4*pi/3.*(vte/a_ii)**3
        write (ifile,'(a20,f12.3)') ' N_c: ',4*pi/3.*(eps/a_ii)**3
        write (ifile,'(a20,f12.3)') ' R/lambda_De: ',eps/(max(vte,1.e-8))
        write (ifile,'(a20,f12.3)') ' Timestep: ',dt
        write (ifile,'(a20,f12.3)') ' Max timestep: ',0.45*sqrt(3.)*eps**2/abs(qe)*vte

        write (ifile,'(a20,1pe12.3)') ' Neighbour search radius: ',r_neighbour
        write (ifile,'(a20,f12.3)') ' MAC theta: ',theta
        write (ifile,'(a20,1pe12.3)') ' Particle # ratio: ',4.e6*lambda*omega*abs(qe)

        write (ifile,'(a,f9.2,a3,f9.2,a3,f9.2/)') ' Graphics box: ',xl,' x ',yl,' x ',zl
        write (ifile,'(a,3f9.2/)') ' Laser focus: ',focus(1:3)


        write (ifile,*) ' Electrons: ', ne
        write (ifile,*) ' Ions: ', ni
        write (ifile,*) ' Beam particles: ', np_beam
 	write (ifile,*) ' Beam angles ',theta_beam, phi_beam
        write (ifile,*) ' Particles per PE: ', np_local


        write (ifile,'(/a/a)') ' Switches:','--------'
        write (ifile,'(a20,l3)') ' Coulomb forces: ',coulomb
        write (ifile,'(a20,l3)') ' Lennard-Jones forces: ',lenjones

        write (ifile,'(a20,l3)') ' load balance: ',balance
        write (ifile,'(a20,l3)') ' restart: ',restart
        write (ifile,'(a20,l3)') ' ramp: ',ramp

! Tree stuff
!        write (ifile,'(a20,l3)') ' domain debug: ',domain_debug
!        write (ifile,'(a20,l3)') ' walk debug: ',walk_debug
!        write (ifile,'(a20,l3)') ' walk summary: ',walk_summary
!        write (ifile,'(a20,l3)') ' dump tree: ',dump_tree
!        write (ifile,'(a20,l3)') ' performance: ',perf_anal
!        write (ifile,'(a20,l3)') ' visit: ',vis_on
!        write (ifile,'(a20,l3/)') ' steering: ',steering

        if (debug_level>=2) then
           write (ifile,'(/a)') 'Other inputs:'
           write(ifile,NML=pepcdata)
        else 
           write (24,'(/a)') 'Other inputs:'
           write(24,NML=pepcdata)
        endif


        !  ibig = 2**63 - 1
        !  write (*,'(i20,b64/o24,z20)') ibig,ibig,ibig,ibig
     end do
  endif

end subroutine param_dump

!  ================================
!
!         ENERGY_CONS
!
!     Find potential, kinetic energies
!
!
!  ================================


subroutine energy_cons(ekine,ekini,emag,ebeam)

  use module_physvars

  implicit none

  integer :: ifile
  real*8 :: epot, ekine, ekini, ebeam,  etot
  real*8 :: emag, eproton_max

  call potenergy(epot,emag)
  call kinenergy(ekine, ekini, ebeam, eproton_max)

  etot = epot + emag + ekine + ekini + ebeam




  laser_energy: select case(beam_config)
  case(4)
     elaser = 3./8.*omega**2*sigma**2*vosc**2*tlaser
  case(5)
     elaser = 3./8.*omega**2*sigma**2*vosc**2*tpulse
  case default
     elaser = 0
  end select laser_energy

  if ( my_rank == 0 .and. debug_level.ge.1 ) then
     do ifile = 6,24,18
        write (ifile,'(20x,3a20/6(a20,1pe18.8,8x,1pe12.5,8x,1pe12.5/),a20,1pe18.8)') &
	     'norm  ','keV ','erg ', &
	     ' P.E. = ',epot, epot*convert_keV, epot*convert_erg, &
	     ' Magnetic E. = ',emag, emag*convert_keV, emag*convert_erg, &
	     ' Electron K.E. norm/erg: ',ekine, ekine*convert_keV, ekine*convert_erg, &
             ' Ion K.E. norm/erg: ',ekini, ekini*convert_keV, ekini*convert_erg, &
	     ' Beam K.E.  = ',ebeam, ebeam*convert_keV, ebeam*convert_erg, &
	     ' Total: ',etot, etot*convert_keV, etot*convert_erg, &
             ' Laser energy = ',elaser

        write (ifile,'(3(a20,f12.5/))') 'Plasma Te (keV):',convert_kev*ekine/max(1,ne), &
					'Ti (keV):',convert_kev*ekini/max(1,ni), &
					'Max Eproton (MeV):',convert_kev*1.e-3*eproton_max
     end do
     ! Write out to energy.dat file
     if (itime.eq.1)  write(75,'(a)') '! time  Upot  Umag  Ukin_e Ukin_i Ukin_beam Utot Tpon xc'
     write (75,'(f12.5,5(1pe12.3),1pe13.5)') &
          trun, convert_kev*epot, convert_kev*emag, convert_kev*ekine, convert_kev*ekini,&
          convert_kev*ebeam, convert_kev*etot
  endif
end subroutine energy_cons



!  ===================================================================
!
!                              POTENERGY
!
!   Calculate electrostatic and magnetic potential energies:
!
!  ===================================================================

subroutine potenergy(epot_total,emag_total)
  use module_physvars
  use module_particle_props
  implicit none
  include 'mpif.h'


  integer :: p, ierr 

  real*8 :: upartial, umag, gamma
  real*8, intent(out) :: epot_total, emag_total
  logical :: pot_debug=.false.

  epot_total = 0.  ! Global potential energy
  emag_total = 0.  ! Global magnetic energy

  !  Sum Potential Energy

  upartial = 0. ! Single PE partial potential energy sum
  umag = 0. ! Single PE partial potential energy sum

  do p=1, np_local
     upartial = upartial + 0.5*q(p)*pot(p)
     if (bfields) then
        gamma = sqrt(1.0+ux(p)**2+uy(p)**2+uz(p)**2)
        umag = umag + 0.5*q(p)/gamma*( ux(p)*ax(p) + uy(p)*ay(p) + uz(p)*az(p) )
     endif
     if (pot_debug) then
        write (ipefile,'(a,i5,a,i5,3f10.3,a,f12.4,a,3f12.5)') & 
	'local particle ',p,' label ',pelabel(p),x(p),y(p),q(p),' pot ',pot(p),' mag ',ax(p), ay(p), az(p)
     endif
  end do

  if (pot_debug) write (ipefile,'(a,1pe11.4)') 'partial PE sum',upartial


  call MPI_ALLREDUCE(upartial, epot_total,1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(umag, emag_total, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)


end subroutine potenergy



!  ===================================================================
!
!                              KINENERGY
!
!   Calculate kinetic energy
!
!  ===================================================================

subroutine kinenergy(ekine,ekini,ebeam,ebeam_max)
  use module_physvars
  use module_particle_props

  implicit none
  include 'mpif.h'

  integer :: p, ierr
  real*8 :: ekine, ekini, ebeam, sum_plas_e, sum_plas_i, sum_beam, gamma
  real*8 :: ebeam_max, ebm_local, u2, dUk
  real*8, dimension(np_local) :: uhx, uhy, uhz

  sum_plas_e = 0.
  sum_plas_i = 0.

  sum_beam = 0.
  ebeam_max=0.
  ebm_local=0.

  do p=1, np_local
  ! Velocities at previous 1/2-step to synch with P.E.
    uhx(p) = ux(p)-dt*q(p)*Ex(p)/m(p)/2. 
    uhy(p) = uy(p)-dt*q(p)*Ey(p)/m(p)/2.
    uhz(p) = uz(p)-dt*q(p)*Ez(p)/m(p)/2.
    u2 = uhx(p)**2 + uhy(p)**2 + uhz(p)**2
    gamma = sqrt(1.0 + uhx(p)**2 + uhy(p)**2 + uhz(p)**2)

! KE contrib
    if (scheme.eq.7 .or. scheme.eq.8) then
      dUk = 0.5*m(p)*u2
    else
      dUk = u2*m(p)/(gamma+1)  ! equivalent to mc^2(g-1) but less prone to rounding error
    endif

    if (pelabel(p) <= ne) then
!  Sum local plasma electron kinetic energy
	sum_plas_e = sum_plas_e+dUk

!  Sum beam energy - assumed to be protons
    else if (nproton>0 .and. pelabel(p) >= proton_label .and. pelabel(p) <= proton_label+nproton) then
        sum_beam = sum_beam + dUk
	ebm_local = max(ebm_local,dUk)

!  Plasma ion energies 
    else if (pelabel(p) <=ne+ni) then
        sum_plas_i = sum_plas_i + dUk
    endif
 end do

! Gather partial sums together for global energies
 
  call MPI_ALLREDUCE(sum_plas_e, ekine, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(sum_plas_i, ekini, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
  call MPI_ALLREDUCE(sum_beam, ebeam, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
! Max proton energy
  call MPI_ALLREDUCE(ebm_local, ebeam_max, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)


end subroutine kinenergy

!  =========================
!
!  Estimate force errors
!
!  =========================

      subroutine error_test(ntest)

      use module_physvars
      use module_particle_props
      use module_utilities
      use module_gle
      use module_htable
      use module_pepc_wrappers

      implicit none
!      integer, parameter :: ntest = 3
      integer, intent(in) ::  ntest
      integer :: listerr(ntest)
      real*8, dimension(ntest) ::  potd, exd, eyd, ezd
      real*8 :: dfx2, dfy2, dfz2, fxs, fys, fzs, dpot, spot
      real*8 :: errfx, errfy, errfz, errf_ave, err_pot
      integer :: i, iseed = -317, isamp


! direct force evaluation
! TODO: make random list of particles for sample

	write (*,*) "Doing error analysis .."

      do i=1,ntest
         isamp = np_local*rano(iseed)+1
         listerr(i) = isamp

      end do

! Get direct forces, potential
      if (force_law==3) then
	call force_direct(np_local,ntest,x(1:np_local),y(1:np_local),z(1:np_local),q(1:np_local),listerr(1:ntest),eps,force_const,exd,eyd,ezd,potd)

      else if (force_law==2) then
	call force_direct_2d(np_local,ntest,x(1:np_local),y(1:np_local),q(1:np_local),listerr(1:ntest),eps,force_const,exd,eyd,potd)
      else
! no force
	exd(1:np_local)=0
	eyd(1:np_local)=0
	potd(1:np_local)=0
      endif

! find rms error

      dfx2 = 0.
      dfy2 = 0.
      dfz2 = 0.
      fxs = 0.
      fys = 0.
      fzs = 0.
      dpot=0.
      spot=0.

      do i=1,ntest
        dfx2 = dfx2 + (exd(i) - ex(listerr(i)))**2
        dfy2 = dfy2 + (eyd(i) - ey(listerr(i)))**2
        if (force_law==3) dfz2 = dfz2 + (ezd(i) - ez(listerr(i)))**2
        dpot = dpot + (potd(i) - pot(listerr(i)))**2
	spot = spot+ potd(i)**2
        fxs = fxs + exd(i)**2
        fys = fys + eyd(i)**2
        if (force_law==3) fzs = fzs + ezd(i)**2
      end do


      err_pot = sqrt(dpot/spot)
      errfx = sqrt(dfx2/fxs)
      errfy = sqrt(dfy2/fys)
      if (force_law==3) then
	errfz = sqrt(dfz2/fzs) 
      else
	errfz=0

      endif 
      errf_ave = (errfx+errfy+errfz)/idim

      open(60,file='forces.dat')
      write (*,*) 'Writing forces, potentials to forces.dat'

      if (force_law==3) then
        write (60,*) '    i     list   nlist,  pot_tree       pot_direct       ex_tree       ex_direct', &
	'   ey_tree        ey_direct       ez_tree        ez_direct'
!      write (60,'((2i8,i6,8(1pe14.4)))') (i,listerr(i),nterm(i),pot(listerr(i)),potd(i), &
!	ex(listerr(i)),exd(i), ey(listerr(i)), eyd(i), ez(listerr(i)), ezd(i), i=1,ntest)
! Currently no access to nterm(i)
        write (60,'((2i8,8(1pe14.4)))') (i,listerr(i),pot(listerr(i)),potd(i), &
	ex(listerr(i)),exd(i), ey(listerr(i)), eyd(i), ez(listerr(i)), ezd(i), i=1,ntest)

  
        write (6,'(a/a20,1pe13.6/a20,3(1pe13.6)/a20,1pe13.6)') 'Relative rms errors:','Potential ',err_pot, &
	'Forces (x,y,z) ',errfx,errfy,errfz,'Average force',errf_ave
        write (60,'(a/a20,1pe13.6/a20,3(1pe13.6)/a20,1pe13.6)') 'Relative rms errors:','Potential ',err_pot, &
	'Forces (x,y,z) ',errfx,errfy,errfz,'Average force',errf_ave
!      write (6,'(a,i8,a1,i8)') 'Ave. list length',SUM(nterm(1:np_local))/np_local,'/',np_local
!      write (60,'(a,i8,a1,i8)') 'Ave. list length',SUM(nterm(1:np_local))/np_local,'/',np_local


      else if (force_law==2) then
        write (60,*) '    i     list   nlist,  pot_tree       pot_direct       ex_tree       ex_direct', &
	'   ey_tree        ey_direct'
!      write (60,'((2i8,i6,8(1pe14.4)))') (i,listerr(i),nterm(i),pot(listerr(i)),potd(i), &
!	ex(listerr(i)),exd(i), ey(listerr(i)), eyd(i), ez(listerr(i)), ezd(i), i=1,ntest)
! Currently no access to nterm(i)
        write (60,'((2i8,6(1pe14.4)))') (i,listerr(i),pot(listerr(i)),potd(i), &
	ex(listerr(i)),exd(i), ey(listerr(i)), eyd(i), i=1,ntest)

  
        write (6,'(a/a20,1pe13.6/a20,2(1pe20.12)/a20,1pe20.12)') 'Relative rms errors:','Potential ',err_pot, &
	'Forces (x,y) ',errfx,errfy,'Average force',errf_ave
        write (60,'(a/a20,1pe13.6/a20,2(1pe20.12)/a20,1pe20.12)') 'Relative rms errors:','Potential ',err_pot, &
	'Forces (x,y) ',errfx,errfy,'Average force',errf_ave
!      write (6,'(a,i8,a1,i8)') 'Ave. list length',SUM(nterm(1:np_local))/np_local,'/',np_local
!      write (60,'(a,i8,a1,i8)') 'Ave. list length',SUM(nterm(1:np_local))/np_local,'/',np_local
      else
      endif

      close(60)

! Tree diagnostics
! If interaction lists needed, must ensure that intlist() is large enough to contain all lists
! - will otherwise just get last pass of tree walk

        call diagnose_tree(particles)   ! Printed tree info (htable etc)
        call draw_tree2d(xl)     ! Draw PE-trees
!        call draw_lists      ! Draw interaction lists
        call draw_domains()   ! Domains
       

      end subroutine error_test


!  ====================================================================
!
!                              FORCE_DIRECT
!
!   Direct PP force sum for error estimates
!
!  ====================================================================

subroutine force_direct(n,ntest,x,y,z,q, list, eps, const, ex, ey, ez, pot)
  implicit none
  real, intent(in) :: eps, const
  integer, intent(in) :: n, ntest
  integer, intent(in), dimension(n) :: list 
  real*8, intent(in), dimension(n) :: x, y, z, q  ! coords and charge 
  real*8,  dimension(ntest) :: ex, ey, ez, pot  ! fields and potential to return

  real :: eps2, d, d3, dx, dy, dz
  integer :: i,j,k
!  write (*,*) 'Direct sum params: ',eps,const
!  write (*,'(a10,a20/(i6,4f15.3))') 'DIRECT | ','Initial buffers: ',(i, x(i), y(i), z(i), q(i),i=1,n) 
 eps2=eps**2
  ex(1:ntest) = 0.
  ey(1:ntest) = 0.
  ez(1:ntest) = 0.
  pot(1:ntest) = 0.

  !  PP contribution from simulation volume

  do  k=1,ntest
     i=list(k)
!     write (*,*) i,'x_i=',x(i)
     do  j=1,n
        if (j.ne.i) then
           dx=x(i)-x(j)
           dy=y(i)-y(j)
           dz=z(i)-z(j)
           d=sqrt(dx**2+dy**2+dz**2+eps2)
           d3=d**3
           ex(k) = ex(k) + const*q(j)*dx/d3
           ey(k) = ey(k) + const*q(j)*dy/d3
           ez(k) = ez(k) + const*q(j)*dz/d3
           pot(k) = pot(k) + const*q(j)/d
!          write (*,'(i5,a5,f12.3,a5,f12.3,a5,f12.3)') j,'q_j=',q(j),' x_j=',x(j), ' d=',d
        endif
     end do
  end do


end subroutine force_direct

!  ====================================================================
!
!                              FORCE_DIRECT_2D
!
!   Direct PP force sum for error estimates - 2D force law
!
!  ====================================================================

subroutine force_direct_2d(n,ntest,x,y,q, list, eps, const, ex, ey, pot)
  implicit none
  real, intent(in) :: eps, const
  integer, intent(in) :: n, ntest
  integer, intent(in), dimension(n) :: list 
  real*8, intent(in), dimension(n) :: x, y, q  ! coords and charge 
  real*8,  dimension(ntest) :: ex, ey, pot  ! fields and potential to return

  real :: eps2, d2, dx, dy
  integer :: i,j,k
!  write (*,*) 'Direct sum params: ',eps,const
!  write (*,'(a10,a20/(i6,4f15.3))') 'DIRECT | ','Initial buffers: ',(i, x(i), y(i), z(i), q(i),i=1,n) 
 eps2=eps**2
  ex(1:ntest) = 0.
  ey(1:ntest) = 0.
  pot(1:ntest) = 0.

  !  PP contribution from simulation volume

  do  k=1,ntest
     i=list(k)
!     write (*,*) i,'x_i=',x(i)
     do  j=1,n
        if (j.ne.i) then
           dx=x(i)-x(j)
           dy=y(i)-y(j)
           d2 = dx**2+dy**2+eps2
           ex(k) = ex(k) + const*q(j)*dx/d2
           ey(k) = ey(k) + const*q(j)*dy/d2
           pot(k) = pot(k) - const*0.5*q(j)*log(d2)
!          write (*,'(i5,a5,f12.3,a5,f12.3,a5,f12.3)') j,'q_j=',q(j),' x_j=',x(j), ' d=',d
        endif
     end do
  end do


end subroutine force_direct_2d
end module module_diagnostics
