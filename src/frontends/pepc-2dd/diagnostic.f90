! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
!
! Copyright (C) 2002-2014 Juelich Supercomputing Centre,
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

!!!!!!!!!!!!!!!!!!!!
!! helper module
!!!!!!!!!!!!!!!!!!!!


module module_diagnostic

!  use module_pepc_kinds
!  use module_pepc_types
!  use module_interaction_Specific_types
  use module_globals    !,only: root
  use module_shortcut   ,only: zero,one,half
  use module_tool       ,only: cross_product
  implicit none

  contains

    subroutine hamiltonian(np,p,pold,t)
      use module_globals,only: flag_classic,lorentz_tilde
      implicit none
      include 'mpif.h'
      type(t_particle), allocatable, intent(in) :: p(:),pold(:)
      real(kind_particle)          , intent(in) :: t
      integer(kind_particle)       , intent(in) :: np
      integer(kind_particle)                    :: ip,rc=201,rd=202
      real(kind_particle)                       :: upot,uekin,uikin,udar,gam,upot_loc,uekin_loc,uikin_loc,udar_loc,&
                                                   v2,ploc(3),pglo(3),p0(3),q0(3),pn(3),qn(3),   &
                                                   cross(3),cross_glo(3),crossloc(1:3)

      upot     = zero
      uekin    = zero
      uikin    = zero
      udar     = zero
      upot_loc = zero
      uekin_loc= zero
      uikin_loc= zero
      udar_loc = zero
      ploc     = -one
      crossloc = -one 
      pglo     = -one

      do ip = 1,np

      !!!!   ENERGY !!!!!!!!!!
        v2       = dot_product( p(ip)%data%v/vtilde, p(ip)%data%v/vtilde )
        gam      = v2/( p(ip)%data%g + one ) 
!        if (flag_classic) gam         =  half*v2
        
        upot_loc = upot_loc + half*p(ip)%data%q*p(ip)%results%pot
        udar_loc = udar_loc + half/lorentz_tilde*p(ip)%data%q*dot_product( p(ip)%results%A(1:3) , p(ip)%data%v(1:3) )/p(ip)%data%g
!        udar_loc = udar_loc + half/lorentz_tilde*p(ip)%data%q*( p(ip)%results%A(3)*p(ip)%data%v(3) )
        
        if ( p(ip)%data%q .lt. zero ) then
            uekin_loc = uekin_loc +                    p(ip)%data%m*gam*vtilde**2
        else 
            uikin_loc = uikin_loc +                    p(ip)%data%m*gam*vtilde**2
        endif

      !!!!!! CANONICAL MOMENTUM - Max Norm !!!!!!
        ploc(1:3)= max(abs( p(ip)%data%m*( p(ip)%data%v(1:3) - pold(ip)%data%v(1:3) ) +                         &
                p(ip)%data%q/lorentz_tilde*( p(ip)%results%A(1:3) - pold(ip)%results%A(1:3) ) ) , ploc(1:3) )

      !!!!!!!!!!!!!! Symplecity !!!!!!!!!!!!!
        p0(1:3)  = pold(ip)%x(1:3)
        pn(1:3)  = p(ip)%x(1:3)
        
        q0(1:3)  = p(ip)%data%m*pold(ip)%data%v(1:3) + p(ip)%data%q/lorentz_tilde*pold(ip)%results%A(1:3)
        qn(1:3)  = p(ip)%data%m*p(ip)%data%v(1:3)        + p(ip)%data%q/lorentz_tilde*p(ip)%results%A(1:3)


        cross        = cross_product( p0 , q0 ) - cross_product( pn , qn )
        cross        = abs(cross)
        crossloc(1:3)= max( cross(1:3)  , crossloc(1:3) )

      enddo

      call MPI_ALLREDUCE(upot_loc , upot     , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
      call MPI_ALLREDUCE(uekin_loc, uekin    , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
      call MPI_ALLREDUCE(uikin_loc, uikin    , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
      call MPI_ALLREDUCE(udar_loc , udar     , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
      call MPI_ALLREDUCE(ploc     , pglo     , 3, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, rc)
      call MPI_ALLREDUCE(crossloc , cross_glo, 3, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, rc)

      if (root) then
        open(unit=rc,file=trim(folder)//trim("energy_")//trim(adjustl(ischeme))//".dat",form='formatted',status='unknown',position='append')
        open(unit=rd,file=trim(folder)//trim("momentum_")//trim(adjustl(ischeme))//".dat",form='formatted',status='unknown',position='append')
        write(rc,*) t,uekin,uikin,upot,udar,upot+udar+uekin+uikin
        write(rd,*) t,pglo,cross
        close (rc )
        close (rd )
      endif

      end subroutine hamiltonian
      
      
      
      subroutine hamiltonian_weibel(np,p,pold,t)
      use module_globals,only: flag_classic,lorentz_tilde
      implicit none
      include 'mpif.h'
      type(t_particle), allocatable, intent(in) :: p(:),pold(:)
      real(kind_particle)          , intent(in) :: t
      integer(kind_particle)       , intent(in) :: np
      integer(kind_particle)                    :: ip,rc=201,rd=202
      real(kind_particle)                       :: uepot,ubpot,uipot,uekin,uikin,uedar,ubdar,uidar,gam,uekin_loc,uikin_loc,    &
                                                   uedar_loc,ubdar_loc,uidar_loc,v2,ploc(3),pglo(3),p0(3),q0(3),pn(3),qn(3),   &
                                                   cross(3),cross_glo(3),crossloc(1:3),ubkin,ubkin_loc,                        &
                                                   uepot_loc,uipot_loc,ubpot_loc,utot,dar,keav,kiav,kbav,                      &
                                                   ke2,ki2,kb2,ke2_tot,ki2_tot,kb2_tot,rtnp

      uekin     = zero
      ubkin     = zero
      uikin     = zero
      uedar     = zero
      ubdar     = zero
      uidar     = zero
      uepot     = zero
      ubpot     = zero
      uipot     = zero
      uepot_loc = zero
      ubpot_loc = zero
      uipot_loc = zero
      uekin_loc = zero
      ubkin_loc = zero
      uikin_loc = zero
      uedar_loc = zero
      ubdar_loc = zero
      uidar_loc = zero
      ploc      = -one
      crossloc  = -one 
      pglo      = -one
      keav      =  zero
      kiav      =  zero
      kbav      =  zero
      ke2       =  zero
      ki2       =  zero
      kb2       =  zero
      ke2_tot   =  zero
      ki2_tot   =  zero
      kb2_tot   =  zero
      
      rtnp            = real(tnp, kind=kind_particle)   
      
      do ip = 1,np

      !!!!   ENERGY !!!!!!!!!!
        v2       = dot_product( p(ip)%data%v/vtilde, p(ip)%data%v/vtilde )
        gam      = v2/( p(ip)%data%g + one ) 
        dar      = dot_product( p(ip)%results%A(1:3) , p(ip)%data%v(1:3) )
        
        if ( p(ip)%label .eq. -1 ) then
            uekin_loc = uekin_loc +                    p(ip)%data%m*gam*vtilde**2
            ke2       = ke2       +                  ( p(ip)%data%m*gam*vtilde**2 )**2 
            uepot_loc = uepot_loc + half*p(ip)%data%q*p(ip)%results%pot
            uedar_loc = uedar_loc + half/lorentz_tilde*p(ip)%data%q*dar/p(ip)%data%g
        elseif ( p(ip)%label .eq. 1 ) then 
            uikin_loc = uikin_loc +                    p(ip)%data%m*gam*vtilde**2
            ki2       = ki2       +                  ( p(ip)%data%m*gam*vtilde**2 )**2
            uipot_loc = uipot_loc + half*p(ip)%data%q*p(ip)%results%pot
            uidar_loc = uidar_loc + half/lorentz_tilde*p(ip)%data%q*dar/p(ip)%data%g
        elseif ( p(ip)%label .eq. 0 ) then 
            ubkin_loc = ubkin_loc +                    p(ip)%data%m*gam*vtilde**2
            kb2       = kb2       +                  ( p(ip)%data%m*gam*vtilde**2 )**2
            ubpot_loc = ubpot_loc + half*p(ip)%data%q*p(ip)%results%pot
            ubdar_loc = ubdar_loc + half/lorentz_tilde*p(ip)%data%q*dar/p(ip)%data%g
        endif

      !!!!!! CANONICAL MOMENTUM - Max Norm !!!!!!
        ploc(1:3)= max(abs( p(ip)%data%m*( p(ip)%data%v(1:3) - pold(ip)%data%v(1:3) )                   +&
                p(ip)%data%q/lorentz_tilde*( p(ip)%results%A(1:3) - pold(ip)%results%A(1:3) ) ) , ploc(1:3) )

      !!!!!!!!!!!!!! Symplecity !!!!!!!!!!!!!
        p0(1:3)  = pold(ip)%x(1:3)
        pn(1:3)  = p(ip)%x(1:3)
        
        q0(1:3)  = p(ip)%data%m*pold(ip)%data%v(1:3) + p(ip)%data%q/lorentz_tilde*pold(ip)%results%A(1:3)
        qn(1:3)  = p(ip)%data%m*p(ip)%data%v(1:3)        + p(ip)%data%q/lorentz_tilde*p(ip)%results%A(1:3)


        cross        = cross_product( p0 , q0 ) - cross_product( pn , qn )
        cross        = abs(cross)
        crossloc(1:3)= max( cross(1:3)  , crossloc(1:3) )

      enddo
      
      
      call MPI_ALLREDUCE(uepot_loc, uepot   , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
      call MPI_ALLREDUCE(ubpot_loc, ubpot   , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
      call MPI_ALLREDUCE(uipot_loc, uipot   , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
      call MPI_ALLREDUCE(uekin_loc, uekin   , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
      call MPI_ALLREDUCE(ubkin_loc, ubkin   , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
      call MPI_ALLREDUCE(uikin_loc, uikin   , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
      call MPI_ALLREDUCE(ki2      , ki2_tot , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
      call MPI_ALLREDUCE(ke2      , ke2_tot , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
      call MPI_ALLREDUCE(kb2      , kb2_tot , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
      call MPI_ALLREDUCE(uedar_loc, uedar   , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
      call MPI_ALLREDUCE(ubdar_loc, ubdar   , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
      call MPI_ALLREDUCE(uidar_loc, uidar   , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
      call MPI_ALLREDUCE(ploc     , pglo    , 3, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, rc)
      call MPI_ALLREDUCE(crossloc, cross_glo, 3, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, rc)
      
      keav      =    uekin/(rtnp)
      kiav      =    uikin/(rtnp)
      kbav      =    ubkin/(rtnp)
      ki2_tot   =   rtnp/( rtnp - one )*( ki2_tot/(rtnp) - kiav**2  )
      ke2_tot   =   rtnp/( rtnp - one )*( ke2_tot/(rtnp) - keav**2 )
      kb2_tot   =   rtnp/( rtnp - one )*( kb2_tot/(rtnp) - kbav**2  )

      if (root) then
        open(unit=rc,file=trim(folder)//trim("energy_")//trim(adjustl(ischeme))//".dat",form='formatted',status='unknown',position='append')
        open(unit=rd,file=trim(folder)//trim("momentum_")//trim(adjustl(ischeme))//".dat",form='formatted',status='unknown',position='append')
        utot = uekin+uikin+ubkin+uepot+uipot+ubpot+uedar+uidar+ubdar
        write(rc,*) t,uekin,uikin,ubkin,keav,kiav,kbav,ke2_tot,ki2_tot,kb2_tot,uepot,uipot,ubpot,uedar,uidar,ubdar,utot
        write(rd,*) t,pglo,cross
        close (rc )
        close (rd )
      endif

      end subroutine hamiltonian_weibel
      
      
      
      subroutine densities_weibel(np,p,t)
      implicit none
      include 'mpif.h'
      type(t_particle), allocatable, intent(in) :: p(:)
      real(kind_particle)          , intent(in) :: t
      integer(kind_particle)       , intent(in) :: np
      integer(kind_particle)                    :: ip,j,rc=201
      real(kind_particle)                       :: Eeloc(3),Eeglo(3),Eiloc(3),Eiglo(3),Ebloc(3),Ebglo(3),&
                                                   Beloc(3),Beglo(3),Biloc(3),Biglo(3),Bbloc(3),Bbglo(3)

      Eeloc = zero
      Beloc = zero
      Eeglo = zero
      Beglo = zero
      Eiloc = zero
      Biloc = zero
      Eiglo = zero
      Biglo = zero
      Ebloc = zero
      Bbloc = zero
      Ebglo = zero
      Bbglo = zero

      do ip = 1,np
        if ( p(ip)%label .eq. -1 ) then
            do j =1,3
                Eeloc(j) = Eeloc(j) + p(ip)%results%E(j)**2
                Beloc(j) = Beloc(j) + p(ip)%results%B(j)**2
            enddo
        else if ( p(ip)%label .eq. 1 ) then
            do j =1,3
                Eiloc(j) = Eiloc(j) + p(ip)%results%E(j)**2
                Biloc(j) = Biloc(j) + p(ip)%results%B(j)**2
            enddo
        else if ( p(ip)%label .eq. 0 ) then
            do j =1,3
                Ebloc(j) = Ebloc(j) + p(ip)%results%E(j)**2
                Bbloc(j) = Bbloc(j) + p(ip)%results%B(j)**2
            enddo
        endif
      enddo
      
      call MPI_ALLREDUCE(Eeloc, Eeglo, 3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
      call MPI_ALLREDUCE(Eiloc, Eiglo, 3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
      call MPI_ALLREDUCE(Ebloc, Ebglo, 3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
      call MPI_ALLREDUCE(Beloc, Beglo, 3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
      call MPI_ALLREDUCE(Biloc, Biglo, 3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
      call MPI_ALLREDUCE(Bbloc, Bbglo, 3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)

      if (root) then
        open(unit=rc,file=trim(folder)//trim("density_")//trim(adjustl(ischeme))//".dat",form='formatted',status='unknown',position='append')
        write(rc,*) t,Eeglo(1:3),Eiglo(1:3),Ebglo(1:3),Beglo(1:3),Biglo(1:3),Bbglo(1:3)
        close (rc )

      endif

      end subroutine densities_weibel


      subroutine densities(np,p,t)
      implicit none
      include 'mpif.h'
      type(t_particle), allocatable, intent(in) :: p(:)
      real(kind_particle)          , intent(in) :: t
      integer(kind_particle)       , intent(in) :: np
      integer(kind_particle)                    :: ip,j,jp,rc=201
      real(kind_particle)                       :: Eloc(3),Eglo(3),Bloc(3),Bglo(3)

      Eloc = zero
      Bloc = zero
      Eglo = zero
      Bglo = zero


      do ip = 1,np
!        if ( MOD(p(ipl)%label,2_kind_particle) .eq. one ) then
            jp = ip! 2*ip-1
            do j =1,3
                Eloc(j) = Eloc(j) + p(jp)%results%E(j)**2
                Bloc(j) = Bloc(j) + p(jp)%results%B(j)**2
            enddo
!        endif
      enddo

      call MPI_ALLREDUCE(Eloc, Eglo, 3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
      call MPI_ALLREDUCE(Bloc, Bglo, 3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)

      if (root) then
        open(unit=rc,file=trim(folder)//trim("density_")//trim(adjustl(ischeme))//".dat",form='formatted',status='unknown',position='append')
        write(rc,*) t,Eglo(1:3),Bglo(1:3)
        close (rc )

      endif

      end subroutine


      subroutine interpolated_densities(field,t)
      use encap, only: field_grid_t
      implicit none
      include 'mpif.h'
      type(field_grid_t)           ,   intent(in) :: field
      real(kind_particle)          , intent(in)   :: t
      integer(kind_particle)                      :: np,ip,j,jp,rc=201
      real(kind_particle)                         :: Eloc(3),Eglo(3),Bloc(3),Bglo(3)

      Eloc = zero
      Bloc = zero
      Eglo = zero
      Bglo = zero

      np   = field%nl

      do ip = 1,np
!        if ( MOD(p(ipl)%label,2_kind_particle) .eq. one ) then
            jp = ip! 2*ip-1
            do j =1,3
                Eloc(j) = Eloc(j) + field%p(jp)%results%E(j)**2
                Bloc(j) = Bloc(j) + field%p(jp)%results%B(j)**2
            enddo
!        endif
      enddo

      call MPI_ALLREDUCE(Eloc, Eglo, 3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
      call MPI_ALLREDUCE(Bloc, Bglo, 3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)

      if (root) then
        open(unit=rc,file=trim(folder)//trim("interpolated_density_")//trim(adjustl(ischeme))//".dat",form='formatted',status='unknown',position='append')
        write(rc,*) t,Eglo(1:3),Bglo(1:3)
        close (rc )

      endif

      end subroutine


      subroutine beam_rnv(tnp,p,t)
      implicit none
      include 'mpif.h'
      type(t_particle), allocatable, intent(in) :: p(:)
      real(kind_particle)          , intent(in) :: t
      integer(kind_particle)       , intent(in) :: tnp
      integer(kind_particle)                    :: ip,jp,np,rc=201
      real(kind_particle)                       :: rloc,rglo,vloc(1:3),vglo(1:3),vrloc,vrglo,rtnp

      rloc      = zero
      vloc(1:3) = zero
      vrloc     = zero
      rglo      = zero
      vglo(1:3) = zero      
      vrglo     = zero

      np    = size(p, kind=kind_particle)
      rtnp  = real(tnp, kind=kind_particle)

      do ip = 1,np
            rloc  = rloc   + dot_product(p(ip)%x , p(ip)%x)
            vrloc = vrloc  + dot_product(p(ip)%data%v , p(ip)%x)/p(ip)%data%g
            do jp = 1,3
                vloc(jp)  = vloc(jp)   + p(ip)%data%v(jp)**2/p(ip)%data%g**2
            enddo
      enddo

      call MPI_ALLREDUCE(rloc , rglo , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
      call MPI_ALLREDUCE(vloc , vglo , 3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
      call MPI_ALLREDUCE(vrloc, vrglo, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)

      rglo       = sqrt(rglo/rtnp)
      vglo(1:3)  = sqrt(vglo(1:3)/rtnp)
      vrglo      = vrglo/rtnp/rglo


      if (root) then
        open(unit=rc,file=trim(folder)//trim("beam_rnv_")//trim(adjustl(ischeme))//".dat",form='formatted',status='unknown',position='append')
        write(rc,*) t,rglo,vglo(1:3),vrglo
        close (rc )
      endif

      end subroutine



end module
