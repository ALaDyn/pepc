subroutine param_dump

  use physvars
  implicit none
  integer :: ifile


  if (my_rank==0) then
     do ifile = 6,15,9

        write (ifile,'(a20,i4,a11)') ' System config: ',system_config,plasma_configs(system_config)
        write (ifile,'(a20,i4,a11)') ' Target geometry: ',target_geometry,geometries(target_geometry)
        write (ifile,'(a20,1pe12.3)') ' System volume: ',Vplas
        write (ifile,'(a20,f12.3)') ' Sphere radius: ',r_sphere
        write (ifile,'(a20,f12.3)') ' System length: ',x_plasma
        write (ifile,'(a20,f12.3)') ' System width: ',y_plasma
        write (ifile,'(a20,f12.3)') ' System height: ',z_plasma
        write (ifile,'(a20,1pe12.3)') ' Electron charge: ',qe
        write (ifile,'(a20,1pe12.3)') ' Electron mass: ',mass_e
        write (ifile,'(a20,1pe12.3)') ' Ion mass: ',mass_i
        write (ifile,'(a20,f12.3)') ' Te: ',Te_keV
        write (ifile,'(a20,f12.3)') ' Ti: ',Ti_keV
        write (ifile,'(a20,1pe12.3)') ' Debye length: ',vte

        write (ifile,'(a20,f12.4)') ' Ion spacing a: ',a_ii
        write (ifile,'(a20,f12.3)') ' Cloud radius R: ',eps
       
        write (ifile,'(a20,1pe12.3)') ' N_D: ',4*pi/3.*(vte/a_ii)**3
        write (ifile,'(a20,f12.3)') ' N_c: ',4*pi/3.*(eps/a_ii)**3
        write (ifile,'(a20,f12.3)') ' R/lambda_De: ',eps/vte
        write (ifile,'(a20,f12.3)') ' Timestep: ',dt
        write (ifile,'(a20,f12.3)') ' Max timestep: ',0.45*sqrt(3.)*eps**2/abs(qe)*vte

        write (ifile,'(a20,f12.3)') ' MAC theta: ',theta

        write (ifile,'(a,f9.2,a3,f9.2,a3,f9.2/)') ' Graphics box: ',xl,' x ',yl,' x ',zl


        write (ifile,*) ' Electrons: ', ne
        write (ifile,*) ' Ions: ', ni
        write (ifile,*) ' Particles per PE: ', np_local

        if (debug_level>=1) then
           write (ifile,'(/a)') 'Other inputs:'
!           write(ifile,NML=pepcdata)
        else 
           write (15,'(/a)') 'Other inputs:'
!           write(15,NML=pepcdata)
        endif


        !  ibig = 2**63 - 1
        !  write (*,'(i20,b64/o24,z20)') ibig,ibig,ibig,ibig
     end do
  endif

end subroutine param_dump
