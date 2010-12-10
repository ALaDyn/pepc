! =============================================
!
!                CONFIGURE
!
!  Sets up physical system: particle positions, velocities
!
! ==============================================

subroutine configure

  use physvars
  !  use utils
  implicit none
  include 'mpif.h'

  config: select case(system_config)

  case(1)              ! Set up particles according to geometry
     call randion         


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
     if (my_rank==0) write (6,*) 'Warning: no particles set up'
     npart_total=0
  end select config

end subroutine configure

