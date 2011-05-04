! =============================================
!
!                CONFIGURE
!
!  Sets up physical system: particle positions, velocities
!
! ==============================================

subroutine configure

  use physvars
  use module_velocity_setup
  implicit none
  real*8 :: vte_, vti_

  vte_ = real(vte)
  vti_ = real(vti)

  config: select case(system_config)

  case(1)              ! Set up particles according to geometry
     call randion         


     if (vte > 0) then
        call maxwell1(ux,nppm,1,nep,vte_)
        call maxwell1(uy,nppm,1,nep,vte_)
        call maxwell1(uz,nppm,1,nep,vte_)
        call scramble_v(ux,uy,uz,nppm,1,nep)   ! remove x,y,z correlations
     else
        call cold_start(ux,uy,uz,nppm,1,nep)
     endif

     if (vti > 0) then
        call maxwell1(ux,nppm,nep+1,nip,vti_)
        call maxwell1(uy,nppm,nep+1,nip,vti_)
        call maxwell1(uz,nppm,nep+1,nip,vti_)
        call scramble_v(ux,uy,uz,nppm,nep+1,nip) ! remove x,y,z correlations
     else
        call cold_start(ux,uy,uz,nppm,nep+1,nip)
     endif

  case(2)
     call special_start(ispecial)

  case default     ! Default = 0 - no plasma target
     if (my_rank==0) write (6,*) 'Warning: no particles set up'
     npart_total=0
  end select config

end subroutine configure

