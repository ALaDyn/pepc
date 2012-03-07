! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2012 Juelich Supercomputing Centre, 
!                         Forschungszentrum Juelich GmbH,
!                         Germany
! 
! PEPC is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! PEPC is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
! 
! You should have received a copy of the GNU Lesser General Public License
! along with PEPC.  If not, see <http://www.gnu.org/licenses/>.
!

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

