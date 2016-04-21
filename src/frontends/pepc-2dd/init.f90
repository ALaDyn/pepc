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
!! module init
!!!!!!!!!!!!!!!!!!!!

module module_init

  use module_pepc_kinds
  use module_pepc_types
  use module_shortcut  ! , only:two,one,zero,pi
  use module_tool       , only: random,random_boltzmann,random_gauss
  use module_globals    , only: me,mi
  implicit none


  contains

  subroutine read_restart_2d(p,filename)
    use module_globals, only: my_rank
    implicit none
    character(*)                 , intent(in)    :: filename
    type(t_particle), allocatable, intent(inout) :: p(:)
    integer(kind_particle)                       :: ip,rc,jp,size_tmp,io_stat,open_status,np
    real(kind_particle), allocatable             :: tmp(:,:)
    character(255)                               :: str_proc


    np       = size(p, kind = kind_particle )
    size_tmp = 6

    if (allocated(tmp)) deallocate(tmp)
    allocate(tmp(np,size_tmp), stat = rc)

    write( str_proc , '(i10)' ) my_rank
    open (unit=my_rank,file=trim(filename)//trim(adjustl(str_proc))//".dat",action="read", position='rewind')

    do ip = 1,np
        read(my_rank,*) tmp(ip,1:size_tmp)

        p(ip)%label       = int(tmp(ip,1), kind= kind_particle)
        p(ip)%data%q      = tmp(ip,2)
        p(ip)%data%m      = 1.0_8

        if(p(ip)%data%q .gt. 0.0) p(ip)%data%m = p(ip)%data%m * mi/me

        p(ip)%x(1)        = tmp(ip,3)
        p(ip)%x(2)        = tmp(ip,4)
        p(ip)%x(3)        = 0.0_8

        p(ip)%data%v(1)   = tmp(ip,5)
        p(ip)%data%v(2)   = tmp(ip,6)
        p(ip)%data%v(3)   = 0.0_8


    end do


    close(unit=my_rank)
    deallocate(tmp)

  end subroutine read_restart_2d


!  subroutine langmuir_waves(p,field_grid)
!    use module_globals, only: nppd,extent,offset,my_rank,n_ranks,nsp,veth,vith,vedrift,vidrift,Volume,tnp,nsp
!    use zufall, only: random_gaussian_flux
!    use encap
!    implicit none
!
!    type(t_particle), allocatable, intent(inout) :: p(:)
!    type(field_grid_t)           , intent(in)    :: field_grid
!    integer(kind_particle)                       :: ipl,ipg,ix,iy,iz,np,ni(2),jp!,tnp,rc
!    real(kind_particle)                          :: dx(3),L(3)
!
!    real(kind_particle),parameter                :: v0 = 1.0e-2
!
!    np = size(p, kind=kind_particle)
!
!    dx(1:2) = zero
!    ni(1:2) = nppd(1:2)/nsp
!    dx(1:2) = extent(1:2)/ni(1:2)
!    
!    L  = extent  + offset
!
!
!    Volume = L(1)
!    if ( dx(2) .gt. zero ) Volume = Volume*L(2)
!    if ( dx(3) .gt. zero ) Volume = Volume*L(3)
!
!    do ipl = 1, np
!
!      p(ipl)%data%q                      = -one
!      if (nsp .eq. 2) p(ipl)%data%q      = (-one + two*MOD(ipl,2_kind_particle))
!      p(ipl)%data%m                      = one
!      if(p(ipl)%data%q .lt. zero)   then
!
!        call random(p(ipl)%x)
!
!        jp          = mod(ipl,field_grid%nl) + 1
!
!        p(ipl)%x(1) = field_grid%dx(1)*p(ipl)%x(1) + field_grid%p(jp)%x(1)
!        p(ipl)%x(2) = field_grid%dx(2)*p(ipl)%x(2) + field_grid%p(jp)%x(2)
!        p(ipl)%x(3) = zero
!        call random_boltzmann(p(ipl)%data%v,veth,vedrift)
!
!
!            !!! perturbations
!
!          p(ipl)%data%v(1)    = p(ipl)%data%v(1) +  v0*cos( two*pi/L(1)*p(ipl)%x(1)  )
!    !      p(ipl)%data%v(2)    = zero!p(ipl)%data%v(2) +  v0*cos( two*pi/L(2)*p(ipl)%x(2)  )
!          p(ipl)%data%v(3)    = zero
!          
!      else
!          
!        ipg = ipl + min(int(my_rank, kind=kind_particle), mod(tnp, int(n_ranks, kind=kind_particle))) * (tnp / n_ranks + 1) + &
!          max(0_kind_particle, my_rank - mod(tnp, int(n_ranks, kind=kind_particle))) *(tnp / n_ranks)
!
!        ix = mod(ipg - 1, ni(1) ) + 1
!        iy = mod( (ipg - 1) / ni(1) , ni(2) )  + 1
!        iz = (ipg - 1) / ( ni(1) * ni(2) )  + 1
!
!
!        p(ipl)%x(1) = (ix - half)*dx(1) + offset(1)
!        p(ipl)%x(2) = (iy - half)*dx(2) + offset(2)
!        p(ipl)%x(3) = (iz - half)*dx(3) + offset(3)
!
!        p(ipl)%data%m    = p(ipl)%data%m * mi/me
!        p(ipl)%data%v    = zero
!        
!      endif
!
!    end do
!
!  end subroutine langmuir_waves

  subroutine langmuir_waves(p,field_grid)
    use module_globals, only: nppd,extent,offset,my_rank,n_ranks,nsp,veth,vith,vedrift,vidrift,Volume,tnp
    use zufall        , only: random_gaussian_flux
    use encap         , only: field_grid_t
    implicit none

    type(t_particle), allocatable, intent(inout) :: p(:)
    type(field_grid_t)           , intent(in)    :: field_grid
    integer(kind_particle)                       :: ipl,ipg,ix,iy,iz,np,jp!,tnp,rc
    real(kind_particle)                          :: L(3)!,dx(3)

    real(kind_particle),parameter                :: v0 = 1.0e-2

    np = size(p, kind=kind_particle)

!    dx = zero
!    dx(1:2) = extent(1:2)/nppd(1:2)
    
    L  = extent  + offset


    Volume = L(1)
    if ( field_grid%dx(2) .gt. zero ) Volume = Volume*L(2)
    if ( field_grid%dx(3) .gt. zero ) Volume = Volume*L(3)

    do ipl = 1, np

      p(ipl)%data%q                      = -one
      if (nsp .eq. 2) p(ipl)%data%q      = (-one + two*MOD(ipl,2_kind_particle))
      p(ipl)%data%m                      = one


!      ipg = ipl + min(int(my_rank, kind=kind_particle), mod(tnp, int(n_ranks, kind=kind_particle))) * (tnp / n_ranks + 1) + &
!        max(0_kind_particle, my_rank - mod(tnp, int(n_ranks, kind=kind_particle))) *(tnp / n_ranks)
!
!      ix = mod(ipg - 1, nppd(1) ) + 1
!      iy = mod( (ipg - 1) / nppd(1) , nppd(2) )  + 1
!      iz = (ipg - 1) / ( nppd(1) * nppd(2) )  + 1
!
!
!      p(ipl)%x(1) = (ix - half)*dx(1) + offset(1)
!      p(ipl)%x(2) = (iy - half)*dx(2) + offset(2)
!      p(ipl)%x(3) = (iz - half)*dx(3) + offset(3)
      
      call random(p(ipl)%x)
      
      jp          = mod(ipl,field_grid%nl) + 1
      
      p(ipl)%x(1) = field_grid%dx(1)*p(ipl)%x(1) + field_grid%p(jp)%x(1)
      p(ipl)%x(2) = field_grid%dx(2)*p(ipl)%x(2) + field_grid%p(jp)%x(2)
!      write(*,*)ipl,jp,field_grid%p(jp)%x(1),field_grid%p(jp)%x(2)
      p(ipl)%x(3) = zero
      call random_boltzmann(p(ipl)%data%v,veth,vedrift)
      

        !!! perturbations

      p(ipl)%data%v(1)    = p(ipl)%data%v(1) +  v0*cos( two*pi/L(1)*p(ipl)%x(1)  )
!      p(ipl)%data%v(2)    = zero!p(ipl)%data%v(2) +  v0*cos( two*pi/L(2)*p(ipl)%x(2)  )
      p(ipl)%data%v(3)    = zero

      if(p(ipl)%data%q .gt. zero)   then
        p(ipl)%data%m    = p(ipl)%data%m * mi/me
        p(ipl)%data%v    = zero
      endif

      
      p(ipl)%label       = my_rank * (tnp / n_ranks) + ipl

      p(ipl)%results%E   = zero
      p(ipl)%results%pot = zero
      p(ipl)%results%A   = zero
      p(ipl)%results%dxA = zero
      p(ipl)%results%dyA = zero
      p(ipl)%results%B   = zero
      p(ipl)%results%J   = zero
      p(ipl)%results%Jirr= zero
      p(ipl)%work        = one

    end do

  end subroutine langmuir_waves

  subroutine weibell_instability(p,field_grid)
    use module_globals, only: nppd,extent,offset,my_rank,n_ranks,nsp,veth,vith,vedrift,vidrift,Volume,tnp,root
    use zufall        , only: random_gaussian_flux
    use encap         , only: field_grid_t
    implicit none
    include 'mpif.h'

    type(t_particle), allocatable, intent(inout) :: p(:)
    type(field_grid_t)           , intent(in)    :: field_grid
    integer(kind_particle)                       :: ipl,jp,np,n,npp,nb,ne,ni,ni_tot,ne_tot,nb_tot,pratio,rc=202
    real(kind_particle)                          :: r,r0,theta,alpha,rtnp,rtnb,sumV,sumv_tot,rtne,x(1:2),L(1:3)
    logical                                      :: tmp

!    r0 = two*tentominusthree   
!    r0  = two*tentominusfive
!    r0 = tentominustwo*one/five
    r0          = three!one/five!ten/four!one/five!
    pratio      = 10

    np     = size(p, kind=kind_particle)
    ne     = 0
    ni     = 0
    nb     = 0
    ne_tot = 0
    ni_tot = 0
    nb_tot = 0
    sumV   = zero 
    sumV_tot= zero
    
    rtnp   = real(tnp, kind=kind_particle)
    rtnb   = rtnp/real(pratio, kind=kind_particle)
     
    if(my_rank.eq.(n_ranks-1)) rtnb = rtnb + MOD(tnp/pratio, int(n_ranks, kind=kind_particle))
    if(my_rank.eq.(n_ranks-1)) np = np + MOD(tnp, int(n_ranks, kind=kind_particle))
    
    rtne   = ( half*rtnp -  rtnb )
    
    
!    Volume = pi*r0**2
    L  = extent  + offset


    Volume = L(1)
    if ( field_grid%dx(2) .gt. zero ) Volume = Volume*L(2)
    if ( field_grid%dx(3) .gt. zero ) Volume = Volume*L(3)


    alpha = pi*(three - sqrt(five) )

    do ipl = 1, np
!    ipl     = 0
!    do while (ipl < np)
!
!        call random(x)
!        x           = r0*( two*x - one )
!        
!        if ( x(1)**2 + x(2)**2 .le. r0**2 ) then
!            
!            ipl                = ipl +1
            p(ipl)%data%q      = (-one + two*MOD(ipl,2_kind_particle))
            p(ipl)%data%m      =  mi/me

            p(ipl)%label       = 1
            
            call random(p(ipl)%x)
!      
!            jp          = mod(ipl,field_grid%nl) + 1
      
            p(ipl)%x(1) = extent(1)*p(ipl)%x(1) + offset(1)
            p(ipl)%x(2) = extent(2)*p(ipl)%x(2) + offset(2)
!            
!            theta = real(ipl, kind=kind_particle)*alpha 
!            r     = r0*sqrt( ( real(my_rank *np + ipl , kind=kind_particle) )/rtnp )
!
!            p(ipl)%x(1)   = r*cos(theta)
!            p(ipl)%x(2)   = r*sin(theta)
            
!            p(ipl)%x(1:2) = x
            p(ipl)%x(3)   = zero
            ni            = ni + 1
            p(ipl)%data%v = zero


            if ( (mod(ipl,pratio) .eq. 0 )  ) then

                p(ipl)%data%m    =  p(ipl)%data%m * me/mi
                p(ipl)%data%q    = -one
                call random_boltzmann(p(ipl)%data%v,veth,vedrift)
                p(ipl)%data%v(3) = veth(3)
!                p(ipl)%data%v(1:3) = veth(1:3)
                p(ipl)%label     = 0
                nb               = nb+1
                ni               = ni-1

            else if ( (p(ipl)%data%q .lt. zero).and.( mod(ipl,pratio) .ne. 0 ) ) then

              p(ipl)%data%m    =  p(ipl)%data%m * me/mi
!              p(ipl)%data%v(1:3) = veth(1:3)
!              p(ipl)%data%v(3) = -veth(3)
              call random_boltzmann(p(ipl)%data%v,veth,vedrift)
              p(ipl)%data%v(3) = -rtnb*veth(3)/rtne
              p(ipl)%label     = -1
              ne               = ne+1
              ni               = ni-1

            endif



            p(ipl)%results%E   = zero
            p(ipl)%results%pot = zero
            p(ipl)%results%A   = zero
            p(ipl)%results%dxA = zero
            p(ipl)%results%dyA = zero
            p(ipl)%results%B   = zero
            p(ipl)%results%J   = zero
            p(ipl)%results%Jirr= zero
            p(ipl)%work        = one

            sumV               = sumV + p(ipl)%data%v(3)
        
!        endif


    end do
    
     call MPI_ALLREDUCE(ne, ne_tot     , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
     call MPI_ALLREDUCE(ni, ni_tot     , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
     call MPI_ALLREDUCE(nb, nb_tot     , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
     call MPI_ALLREDUCE(sumV, sumV_tot , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)

     if (root) then
        write(*,*) ne_tot,ni_tot,nb_tot,ne_tot+ni_tot+nb_tot,rtne,rtnb,sumv_tot 
     endif
    

  end subroutine weibell_instability
  
  
  
  subroutine weibell_instability_equilibrium(p,field_grid)
    use module_globals, only: nppd,extent,offset,my_rank,n_ranks,nsp,veth,vith,vedrift,vidrift,Volume,tnp,root
    use zufall        , only: random_gaussian_flux
    use encap         , only: field_grid_t
    implicit none
    include 'mpif.h'

    type(t_particle), allocatable, intent(inout) :: p(:)
    type(field_grid_t)           , intent(in)    :: field_grid
    integer(kind_particle)                       :: ipl,np,rc=202
    real(kind_particle)                          :: x(1:2),L(1:3),v0 = 2.0e-5
    logical                                      :: tmp


    np     = size(p, kind=kind_particle)    
    if(my_rank.eq.(n_ranks-1)) np = np + MOD(tnp, int(n_ranks, kind=kind_particle))

    L  = extent  + offset
    Volume = L(1)
    if ( field_grid%dx(2) .gt. zero ) Volume = Volume*L(2)
    if ( field_grid%dx(3) .gt. zero ) Volume = Volume*L(3)


    do ipl = 1, np

            p(ipl)%data%q                      = -one
            if (nsp .eq. 2) p(ipl)%data%q      = (-one + two*MOD(ipl,2_kind_particle))
            p(ipl)%data%m      =  mi/me

            p(ipl)%label       = 1
            
            call random(p(ipl)%x)
      
            p(ipl)%x(1) = extent(1)*p(ipl)%x(1) + offset(1)
            p(ipl)%x(2) = extent(2)*p(ipl)%x(2) + offset(2)
            p(ipl)%x(3) = zero

            call random_boltzmann(p(ipl)%data%v,vith,vidrift)


            if ( p(ipl)%data%q .lt. zero  ) then

              p(ipl)%label       = -1  
              p(ipl)%data%m      =  p(ipl)%data%m * me/mi
              call random_boltzmann(p(ipl)%data%v,veth,vedrift) 
              p(ipl)%data%v(1)   =  p(ipl)%data%v(1) +  v0*cos( two*pi/L(1)*p(ipl)%x(1)  )
            
            endif



            p(ipl)%results%E   = zero
            p(ipl)%results%pot = zero
            p(ipl)%results%A   = zero
            p(ipl)%results%dxA = zero
            p(ipl)%results%dyA = zero
            p(ipl)%results%B   = zero
            p(ipl)%results%J   = zero
            p(ipl)%results%Jirr= zero
            p(ipl)%work        = one

    end do
    

  end subroutine weibell_instability_equilibrium


    
  subroutine neutral_plasma(p,field_grid)
    use module_globals, only: nppd,extent,offset,my_rank,n_ranks,nsp,veth,vith,vedrift,vidrift,Volume,tnp,root
    use zufall        , only: random_gaussian_flux
    use encap         , only: field_grid_t
    implicit none
    include 'mpif.h'

    type(t_particle), allocatable, intent(inout) :: p(:)
    type(field_grid_t)           , intent(in)    :: field_grid
    integer(kind_particle)                       :: ipl,np,rc=202
    real(kind_particle)                          :: x(1:2),L(1:3),v0 = 2.0e-5
    logical                                      :: tmp


    np     = size(p, kind=kind_particle)    
    if(my_rank.eq.(n_ranks-1)) np = np + MOD(tnp, int(n_ranks, kind=kind_particle))

    L  = extent  + offset
    Volume = L(1)
    if ( field_grid%dx(2) .gt. zero ) Volume = Volume*L(2)
    if ( field_grid%dx(3) .gt. zero ) Volume = Volume*L(3)


    do ipl = 1, np

            p(ipl)%data%q                      = -one
            if (nsp .eq. 2) p(ipl)%data%q      = (-one + two*MOD(ipl,2_kind_particle))
            p(ipl)%data%m      =  mi/me

            p(ipl)%label       = 1
            
            call random(p(ipl)%x)
      
            p(ipl)%x(1) = extent(1)*p(ipl)%x(1) + offset(1)
            p(ipl)%x(2) = extent(2)*p(ipl)%x(2) + offset(2)
            p(ipl)%x(3) = zero

            call random_boltzmann(p(ipl)%data%v,vith,vidrift)


            if ( p(ipl)%data%q .lt. zero  ) then

              p(ipl)%label       = -1  
              p(ipl)%data%m      =  p(ipl)%data%m * me/mi
              call random_boltzmann(p(ipl)%data%v,veth,vedrift) 
            
            endif



            p(ipl)%results%E   = zero
            p(ipl)%results%pot = zero
            p(ipl)%results%A   = zero
            p(ipl)%results%dxA = zero
            p(ipl)%results%dyA = zero
            p(ipl)%results%B   = zero
            p(ipl)%results%J   = zero
            p(ipl)%results%Jirr= zero
            p(ipl)%work        = one

    end do
    

  end subroutine neutral_plasma
  
    subroutine periodic_test(p,field_grid)
    use module_globals, only: nppd,extent,offset,my_rank,n_ranks,nsp,veth,vith,vedrift,vidrift,Volume,tnp
    use module_tool   , only: random_gauss
    use zufall        , only: random_gaussian_flux
    use encap         , only: field_grid_t 
    implicit none

    type(t_particle), allocatable, intent(inout) :: p(:)
    type(field_grid_t)           , intent(in)    :: field_grid
    integer(kind_particle)                       :: ipl,jp,np,n
    real(kind_particle)                          :: r,r0,rtnp,x(1:2)

!    r0 = two*tentominusthree   
!    r0  = two*tentominusfive
!    r0 = tentominustwo*one/five
    r0     = one/five!six/five!one/five!ten/four!one/five!

    np     = size(p, kind=kind_particle)
    rtnp   = real(tnp, kind=kind_particle)
    Volume = pi*r0**2

    ipl     = 0
    do while (ipl < np)

        call random(x)
        x           = r0*( two*x - one )
        
        if ( x(1)**2 + x(2)**2 .le. r0**2 ) then
            
            ipl                = ipl +1

            p(ipl)%data%q                      = -one
            if (nsp .eq. 2) p(ipl)%data%q      = (-one + two*MOD(ipl,2_kind_particle))
            p(ipl)%data%m                      =  one
            
            p(ipl)%x(1:2) = x
            p(ipl)%x(3)   = zero
            
            do jp = 1,3
                if ( p(ipl)%x(jp) .lt. offset(jp) )              p(ipl)%x(jp) = p(ipl)%x(jp) + extent(jp)
                if ( p(ipl)%x(jp) .gt. offset(jp) + extent(jp) ) p(ipl)%x(jp) = mod( p(ipl)%x(jp) , extent(jp) )
            enddo

            p(ipl)%data%v = veth

            if(p(ipl)%data%q .gt. zero) then

                p(ipl)%data%m    = p(ipl)%data%m * mi/me
                call random_gauss(p(ipl)%data%v)!p(ipl)%data%v    = vith

            endif

            p(ipl)%label       = my_rank * (tnp / n_ranks) + ipl

            p(ipl)%results%E   = zero
            p(ipl)%results%pot = zero
            p(ipl)%results%A   = zero
            p(ipl)%results%dxA = zero
            p(ipl)%results%dyA = zero
            p(ipl)%results%B   = zero
            p(ipl)%results%J   = zero
            p(ipl)%results%Jirr= zero
            p(ipl)%work        = one
            
        endif


    end do
    

  end subroutine periodic_test
  
  
  subroutine landau_damping(p,field_grid)
    use module_globals, only: nppd,extent,offset,my_rank,n_ranks,nsp,veth,vith,vedrift,vidrift,Volume,tnp,ixdim
    use zufall        , only: random_gaussian_flux
    use encap         , only: field_grid_t
    implicit none

    type(t_particle), allocatable, intent(inout) :: p(:)
    type(field_grid_t)           , intent(in)    :: field_grid
    integer(kind_particle)                       :: ipl,ipg,ix,iy,iz,n(3),np,jp!,tnp,rc
    real(kind_particle)                          :: dx(3),L(3)

    real(kind_particle),parameter                :: v0 = five*tentominustwo

    np = size(p, kind=kind_particle)

    dx = zero
    dx(1:2) = extent(1:2)/nppd(1:2)
    
    L  = extent  + offset


    Volume = L(1)
    if ( dx(2) .gt. zero ) Volume = Volume*L(2)
    if ( dx(3) .gt. zero ) Volume = Volume*L(3)

    do ipl = 1, np

      p(ipl)%data%q                      = -one
      if (nsp .eq. 2) p(ipl)%data%q      = (-one + two*MOD(ipl,2_kind_particle))
      p(ipl)%data%m                      = one


!      ipg = ipl + min(int(my_rank, kind=kind_particle), mod(tnp, int(n_ranks, kind=kind_particle))) * (tnp / n_ranks + 1) + &
!        max(0_kind_particle, my_rank - mod(tnp, int(n_ranks, kind=kind_particle))) *(tnp / n_ranks)
!
!      ix = mod(ipg - 1, nppd(1) ) + 1
!      iy = mod( (ipg - 1) / nppd(1) , nppd(2) )  + 1
!      iz = (ipg - 1) / ( nppd(1) * nppd(2) )  + 1
!
!
!      p(ipl)%x(1) = (ix - half)*dx(1) + offset(1)
!      p(ipl)%x(2) = (iy - half)*dx(2) + offset(2)
!      p(ipl)%x(3) = (iz - half)*dx(3) + offset(3)
      
      call random(p(ipl)%x)
      
      jp          = mod(ipl,field_grid%nl) + 1
      
      p(ipl)%x(1) = field_grid%dx(1)*p(ipl)%x(1) + field_grid%p(jp)%x(1)
      p(ipl)%x(2) = field_grid%dx(2)*p(ipl)%x(2) + field_grid%p(jp)%x(2)
      p(ipl)%x(3) = zero
      call random_boltzmann(p(ipl)%data%v,veth,vedrift)
      

        !!! perturbations

      p(ipl)%data%v(1)    = p(ipl)%data%v(1) +  v0*cos( two*pi/L(1)*p(ipl)%x(1)  )
!      p(ipl)%data%v(2)    = p(ipl)%data%v(2) +  v0*cos( two*pi/L(2)*p(ipl)%x(2)  )
      p(ipl)%data%v(3)    = zero

      if(p(ipl)%data%q .gt. zero)   then
        p(ipl)%data%m    = p(ipl)%data%m * mi/me
        p(ipl)%data%v    = zero
      endif

      p(ipl)%label       = my_rank * (tnp / n_ranks) + ipl

      p(ipl)%results%E   = zero
      p(ipl)%results%pot = zero
      p(ipl)%results%A   = zero
      p(ipl)%results%dxA = zero
      p(ipl)%results%dyA = zero
      p(ipl)%results%B   = zero
      p(ipl)%results%J   = zero
      p(ipl)%results%Jirr= zero
      p(ipl)%work        = one

      
    end do

  end subroutine landau_damping


  subroutine beam(p)
    use module_globals, only: nppd,extent,offset,my_rank,n_ranks,nsp,veth,vith,vedrift,vidrift,Volume,tnp
    use zufall        , only: random_gaussian_flux
    implicit none

    type(t_particle), allocatable, intent(inout) :: p(:)
    integer(kind_particle)                       :: ipl,ipg,ix,iy,iz,n(3),np!,tnp,rc
    real(kind_particle)                          :: dx(3),L(3)

    real(kind_particle),parameter :: v0     = 1e-1
    real(kind_particle),parameter :: x0     = 1e-1

    np = size(p, kind=kind_particle)
    n  = nppd
    dx = extent/tnp
    L  = extent  !+ offset

    do ipl = 1, np

      p(ipl)%data%q                      = -one
      if (nsp .eq. 2) p(ipl)%data%q      = (-one + two*MOD(ipl,2_kind_particle))
      p(ipl)%data%m                      = one


      ipg = ipl + min(int(my_rank, kind=kind_particle), mod(tnp, int(n_ranks, kind=kind_particle))) * (tnp / n_ranks + 1) + &
        max(0_kind_particle, my_rank - mod(tnp, int(n_ranks, kind=kind_particle))) *(tnp / n_ranks)

      ix = mod(ipg - 1, nppd(1) ) + 1
      iy = mod( (ipg - 1) / nppd(1) , nppd(2) )  + 1
      iz = (ipg - 1) / ( nppd(1) * nppd(2) )  + 1


      p(ipl)%x(1:3) = (/ ix , iy ,iz /) *dx(1:3) + offset(1:3)
      p(ipl)%data%v = veth

      p(ipl)%x(2:3) = zero
      if(p(ipl)%data%q .gt. zero) then

        p(ipl)%data%m    = p(ipl)%data%m * mi/me
        p(ipl)%data%v    = vith

      endif


      p(ipg)%data%v(1:2) = zero
      
      p(ipl)%label       = my_rank * (tnp / n_ranks) + ipl

      p(ipl)%results%E   = zero
      p(ipl)%results%pot = zero
      p(ipl)%results%A   = zero
      p(ipl)%results%dxA = zero
      p(ipl)%results%dyA = zero
      p(ipl)%results%B   = zero
      p(ipl)%results%J   = zero
      p(ipl)%results%Jirr= zero
      p(ipl)%work        = one


    end do

  end subroutine beam



  subroutine beam_disk(p)
    use module_globals, only: nppd,extent,offset,my_rank,n_ranks,nsp,veth,vith,vedrift,vidrift,Volume,tnp
    use module_tool   , only: random_gauss
    use zufall        , only: random_gaussian_flux
    implicit none

    type(t_particle), allocatable, intent(inout) :: p(:)
    integer(kind_particle)                       :: ipl,np,n
    real(kind_particle)                          :: r,r0,theta,alpha, phase,rtnp

!    r0 = two*tentominusthree   
!    r0  = two*tentominusfive
!    r0 = tentominustwo*one/five
    r0     = ten/four!six/five!one/five!ten/four!one/five!

    np     = size(p, kind=kind_particle)
    rtnp   = real(tnp, kind=kind_particle)
    Volume = pi*r0**2

    call random(p(1)%x)

    alpha = pi*(three - sqrt(five) )
    phase = two*pi*p(1)%x(1)

    do ipl = 1, np

      p(ipl)%data%q                      = -one
      if (nsp .eq. 2) p(ipl)%data%q      = (-one + two*MOD(ipl,2_kind_particle))
      p(ipl)%data%m                      =  one

!      theta = ipl*alpha + phase
!      r     = r0*( real(my_rank * (tnp / n_ranks) + ipl , kind=kind_particle) )/rtnp
      
      theta = real(ipl, kind=kind_particle)*alpha 
      r     = r0*sqrt( ( real(my_rank * (tnp / n_ranks) + ipl , kind=kind_particle) )/rtnp )
      
!      call random(p(ipl)%x)
!      theta = two*pi*p(ipl)%x(1)
!      r     = r0*sqrt( pi*p(ipl)%x(2) )

        
      p(ipl)%x(1)   = r*cos(theta)
      p(ipl)%x(2)   = r*sin(theta)
      p(ipl)%x(3)   = zero

      p(ipl)%data%v = veth

      if(p(ipl)%data%q .gt. zero) then

        p(ipl)%data%m    = p(ipl)%data%m * mi/me
        call random_gauss(p(ipl)%data%v)!p(ipl)%data%v    = vith

      endif

      p(ipl)%label       = my_rank * (tnp / n_ranks) + ipl

      p(ipl)%results%E   = zero
      p(ipl)%results%pot = zero
      p(ipl)%results%A   = zero
      p(ipl)%results%dxA = zero
      p(ipl)%results%dyA = zero
      p(ipl)%results%B   = zero
      p(ipl)%results%J   = zero
      p(ipl)%results%Jirr= zero
      p(ipl)%work        = one


    end do


  end subroutine beam_disk
  
  
  subroutine solenoid_infinite(p)
    use module_globals, only: nppd,extent,offset,my_rank,n_ranks,nsp,veth,vith,vedrift,vidrift,Volume,tnp
    use module_tool   , only: random_gauss
    use zufall        , only: random_gaussian_flux
    implicit none

    type(t_particle), allocatable, intent(inout) :: p(:)
    integer(kind_particle)                       :: ipl,np,n
    real(kind_particle)                          :: r,r0,theta,alpha, phase,rtnp

    r0 = one!0.2
    np  = size(p, kind=kind_particle)
    rtnp= real(tnp, kind=kind_particle)
    Volume = pi*r0**2

    call random(p(1)%x)

    alpha = pi*(three - sqrt(five) )
    phase = two*pi*p(1)%x(1)

    do ipl = 1, np

      p(ipl)%data%q                      = -one
      if (nsp .eq. 2) p(ipl)%data%q      = (-one + two*MOD(ipl,2_kind_particle))
      p(ipl)%data%m                      = one

      theta = ipl*alpha + phase

      p(ipl)%x(1)   = r0*cos(theta)
      p(ipl)%x(2)   = r0*sin(theta)
      p(ipl)%x(3)   = zero

      p(ipl)%data%v = veth

      if(p(ipl)%data%q .gt. zero) then

        p(ipl)%data%m    = p(ipl)%data%m * mi/me
        call random_gauss(p(ipl)%data%v)!p(ipl)%data%v    = vith

      endif

      p(ipl)%label       = my_rank * (tnp / n_ranks) + ipl

      p(ipl)%results%E   = zero
      p(ipl)%results%pot = zero
      p(ipl)%results%A   = zero
      p(ipl)%results%dxA = zero
      p(ipl)%results%dyA = zero
      p(ipl)%results%B   = zero
      p(ipl)%results%J   = zero
      p(ipl)%results%Jirr= zero
      p(ipl)%work        = one

      
    end do


  end subroutine solenoid_infinite
  
  
  subroutine solenoid(p)
    use module_globals, only: nppd,extent,offset,my_rank,n_ranks,nsp,veth,vith,vedrift,vidrift,Volume,tnp
    use module_tool   , only: random_gauss
    use zufall        , only: random_gaussian_flux
    implicit none

    type(t_particle), allocatable, intent(inout) :: p(:)
    integer(kind_particle)                       :: ipl,np,n
    real(kind_particle)                          :: r0,rtnp,L,dx

    r0     = half
    L      = 100.0
    np     = size(p, kind=kind_particle)
    rtnp   = real(tnp, kind=kind_particle)
    Volume = L*r0
    dx     = two*L/rtnp
    
    do ipl = 1, np

      p(ipl)%data%q                      = (-one + two*MOD(ipl,2_kind_particle))

      p(ipl)%data%m                      = one
      if(p(ipl)%data%q .gt. zero) p(ipl)%data%m    = p(ipl)%data%m * mi/me

      p(ipl)%x(1)   = -L + ( real(ipl, kind=kind_particle) + (-one + MOD(ipl,2_kind_particle)) )*dx
      p(ipl)%x(2)   = p(ipl)%data%q*r0
      p(ipl)%x(3)   = zero

      p(ipl)%data%v = veth

      p(ipl)%label       = my_rank * (tnp / n_ranks) + ipl

      p(ipl)%results%E   = zero
      p(ipl)%results%pot = zero
      p(ipl)%results%A   = zero
      p(ipl)%results%dxA = zero
      p(ipl)%results%dyA = zero
      p(ipl)%results%B   = zero
      p(ipl)%results%J   = zero
      p(ipl)%results%Jirr= zero
      p(ipl)%work        = one

      
    end do


  end subroutine solenoid



  subroutine thermal(p)
    use module_globals, only: nppd,extent,offset,my_rank,n_ranks,nsp,veth,vith,vedrift,vidrift,Volume,tnp
    use module_tool   , only: random_gauss
    implicit none

    type(t_particle), allocatable, intent(inout) :: p(:)
    integer(kind_particle)                       :: ip,np!,tnp

    real(kind_particle)          , parameter     :: Tkb=1e-3

    np = size(p, kind=kind_particle)
    Volume = one

    do ip=1, np

      p(ip)%data%q                      = one
      if (nsp .eq. 2) p(ip)%data%q      = (-one + two*MOD(ip,2_kind_particle))
      p(ip)%data%m      = one

      if(p(ip)%data%q .gt. 0.0) p(ip)%data%m = p(ip)%data%m * mi/me

      call random(p(ip)%x)

      call random_gauss(p(ip)%data%v)
      p(ip)%data%v      = p(ip)%data%v*sqrt(Tkb * p(ip)%data%m) / p(ip)%data%m

      p(ip)%label       = my_rank * (tnp / n_ranks) + ip

      p(ip)%results%E   = zero
      p(ip)%results%pot = zero
      p(ip)%results%A   = zero
      p(ip)%results%dxA = zero
      p(ip)%results%dyA = zero
      p(ip)%results%B   = zero
      p(ip)%results%J   = zero
      p(ip)%results%Jirr= zero
      p(ip)%work        = one

    end do

  end subroutine thermal

  subroutine thermal2D3V(p)
    use module_globals, only: nppd,extent,offset,my_rank,n_ranks,nsp,veth,vith,vedrift,vidrift,Volume,tnp
    use module_tool   , only: random_gauss
    implicit none

    type(t_particle), allocatable, intent(inout) :: p(:)
    integer(kind_particle)                       :: ip,np!,tnp

    real(kind_particle)          , parameter     :: Tkb=1e-3

    np = size(p, kind=kind_particle)
    Volume = one

    do ip=1, np

      p(ip)%data%q                      = one
      if (nsp .eq. 2) p(ip)%data%q      = (-one + two*MOD(ip,2_kind_particle))
      p(ip)%data%m      = one

      if(p(ip)%data%q .gt. 0.0) p(ip)%data%m = p(ip)%data%m * mi/me

      call random(p(ip)%x)
      p(ip)%x(3)    = zero
!      call random_gauss(p(ip)%data%v)
      p(ip)%data%v      = half!p(ip)%data%v*sqrt(Tkb * p(ip)%data%m) / p(ip)%data%m

      p(ip)%label       = my_rank * (tnp / n_ranks) + ip

      p(ip)%results%E   = zero
      p(ip)%results%pot = zero
      p(ip)%results%A   = zero
      p(ip)%results%dxA = zero
      p(ip)%results%dyA = zero
      p(ip)%results%B   = zero
      p(ip)%results%J   = zero
      p(ip)%results%Jirr= zero
      p(ip)%work        = one

    end do
    
  end subroutine thermal2D3V



end module module_init
