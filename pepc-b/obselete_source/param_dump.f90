subroutine param_dump

  use physvars
  use treevars
  use utils
  implicit none
  integer :: ifile
  include 'namelist.h'


  if (me==0) then
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
        write (ifile,*) ' Particles per PE: ', npp


        write (ifile,'(/a/a)') ' Switches:','--------'
        write (ifile,'(a20,l3)') ' Coulomb forces: ',coulomb
        write (ifile,'(a20,l3)') ' Lennard-Jones forces: ',lenjones

        write (ifile,'(a20,l3)') ' load balance: ',load_balance
        write (ifile,'(a20,l3)') ' restart: ',restart
        write (ifile,'(a20,l3)') ' ramp: ',ramp

! Tree stuff
        write (ifile,'(a20,l3)') ' domain debug: ',domain_debug
        write (ifile,'(a20,l3)') ' walk debug: ',walk_debug
        write (ifile,'(a20,l3)') ' walk summary: ',walk_summary
        write (ifile,'(a20,l3)') ' dump tree: ',dump_tree
        write (ifile,'(a20,l3)') ' performance: ',perf_anal
        write (ifile,'(a20,l3)') ' visit: ',vis_on
        write (ifile,'(a20,l3/)') ' steering: ',steering

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