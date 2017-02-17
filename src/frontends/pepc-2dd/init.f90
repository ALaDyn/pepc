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
  use module_globals    , only: me,mi,n_ranks,my_rank
  implicit none


  contains
  
  
  
  subroutine rect(p)
    use module_globals, only: extent,offset,my_rank,n_ranks,nsp,Volume,tnp,root
!    use encap         , only: field_grid_t
    implicit none

    type(t_particle), allocatable, intent(inout) :: p(:)
!    type(field_grid_t)           , intent(in)    :: field_grid
    integer(kind_particle)                       :: ipl,np!,jp
    
    if (root)     write(*,*)            "== ...Loading Particle's Position  "
    np = size(p, kind=kind_particle)
    
    do ipl = 1, np

      call random(p(ipl)%x)
      
!      jp          = mod(ipl,field_grid%nl) + 1
      
!      p(ipl)%x(1) = field_grid%dx(1)*p(ipl)%x(1) + field_grid%p(jp)%x(1)
!      p(ipl)%x(2) = field_grid%dx(2)*p(ipl)%x(2) + field_grid%p(jp)%x(2)
      p(ipl)%x(1:2) = ( extent(1:2) - offset(1:2) )*p(ipl)%x(1:2) + offset(1:2)
      p(ipl)%x(3) = zero



    end do

  end subroutine rect
  
  
  subroutine load_file(p)
    use module_globals, only: folder,root,my_rank,vtilde
    use mpi
    implicit none

    type(t_particle), allocatable, intent(inout) :: p(:)
    integer(kind_particle)                       :: ip,np,rc=125
    character(255)                               :: filename
    character(*), parameter                      :: part_dir = "particles/"
    integer(kind = MPI_OFFSET_KIND)              :: mpi_disp,my_offset
    integer                                      :: fh, mpi_err
    integer         , dimension(MPI_STATUS_SIZE) :: mpi_stat

        
    if (root)     write(*,*)            "== ...Loading from File  "
    np = size(p, kind=kind_particle)

    
    my_offset = np*my_rank
    
    write(filename,'(a,"restart_label.dat")') trim(folder)//trim(part_dir) 
    call read_mpi_int(filename,p(:)%label)
    
    write(filename,'(a,"restart_m.dat")') trim(folder)//trim(part_dir) 
    call read_mpi_real(filename,p(:)%data%q)
    write(filename,'(a,"restart_m.dat")') trim(folder)//trim(part_dir) 
    call read_mpi_real(filename,p(:)%data%m)
    write(filename,'(a,"restart_x.dat")') trim(folder)//trim(part_dir) 
    call read_mpi_real(filename,p(:)%x(1))
    write(filename,'(a,"restart_y.dat")') trim(folder)//trim(part_dir) 
    call read_mpi_real(filename,p(:)%x(2))
    write(filename,'(a,"restart_z.dat")') trim(folder)//trim(part_dir) 
    call read_mpi_real(filename,p(:)%x(3))
    write(filename,'(a,"restart_vx.dat")') trim(folder)//trim(part_dir) 
    call read_mpi_real(filename,p(:)%data%v(1))
    write(filename,'(a,"restart_vy.dat")') trim(folder)//trim(part_dir) 
    call read_mpi_real(filename,p(:)%data%v(2))
    write(filename,'(a,"restart_vz.dat")') trim(folder)//trim(part_dir) 
    call read_mpi_real(filename,p(:)%data%v(3))
    
    
    do ip = 1,np
        p(ip)%data%g = sqrt( one + dot_product( p(ip)%data%v(:), p(ip)%data%v(:) ) )
        call iclean_fields(p(ip))
    enddo
    
        contains
        
        subroutine read_mpi_real(filename,p)
        implicit none

        character(*)                    , intent(in)  :: filename
        real(kind_physics), dimension(:), intent(out) :: p
     
        real(kind = 8)   , dimension(:), allocatable :: real8_buf
        
        
        if ( allocated(real8_buf) ) deallocate(real8_buf)
        allocate( real8_buf( np ) )
                    
             
        call mpi_file_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, fh, mpi_err)
      
!        if (mpi_err .ne. MPI_SUCCESS) DEBUG_ERROR(*, 'write_restart(): I/O error for ', filename)
        call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)       
        call mpi_file_set_view(    fh, 0_MPI_OFFSET_KIND, MPI_DOUBLE_PRECISION , MPI_DOUBLE_PRECISION,'native'           , MPI_INFO_NULL      , mpi_err)  
        call mpi_file_read_at(fh, my_offset        , real8_buf, int(np, kind = kind_default), MPI_DOUBLE_PRECISION, mpi_stat, mpi_err)
        call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)  
        call mpi_file_close(fh, mpi_err)

        p(:)  = real8_buf(:)
        
        deallocate(real8_buf)
        
            
        end subroutine read_mpi_real
        
        
        subroutine read_mpi_int(filename,p)
        implicit none

        character(*)                       , intent(in)  :: filename
        integer(kind_physics), dimension(:), intent(out) :: p
                
        integer(kind_physics)   , dimension(:), allocatable :: real8_buf
        
        
        if ( allocated(real8_buf) ) deallocate(real8_buf)
        allocate( real8_buf( np ) )
            
             
        call mpi_file_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, fh, mpi_err)
      
!        if (mpi_err .ne. MPI_SUCCESS) DEBUG_ERROR(*, 'write_restart(): I/O error for ', filename)
             
        call mpi_file_set_view(    fh, 0_MPI_OFFSET_KIND, MPI_INTEGER8 , MPI_INTEGER8,'native'           , MPI_INFO_NULL      , mpi_err)  
        call mpi_file_read_at(fh, my_offset        , real8_buf, int(np, kind = kind_default), MPI_INTEGER8, mpi_stat, mpi_err)
        call mpi_file_close(fh, mpi_err)

        
        p(:)  = real8_buf(:) 
        
        deallocate(real8_buf)
        
            
        end subroutine read_mpi_int
    
  end subroutine load_file
  
  
  
  
    subroutine cube(p)
    use module_globals, only: extent,offset,root!,my_rank,n_ranks,Volume,tnp,root
!    use encap         , only: field_grid_t1:
    implicit none

    type(t_particle), allocatable, intent(inout) :: p(:)
!    type(field_grid_t)           , intent(in)    :: field_grid
    integer(kind_particle)                       :: ipl,np,jp
    
    if (root)     write(*,*)            "== ...Loading Particle's Position  "
    np = size(p, kind=kind_particle)
    
    
    do ipl = 1, np

      call random(p(ipl)%x)
      
!      jp          = mod(ipl,field_grid%nl) + 1
!      
!      p(ipl)%x(1) = field_grid%dx(1)*p(ipl)%x(1) + field_grid%p(jp)%x(1)
!      p(ipl)%x(2) = field_grid%dx(2)*p(ipl)%x(2) + field_grid%p(jp)%x(2)
!      p(ipl)%x(3) = field_grid%dx(3)*p(ipl)%x(3) + field_grid%p(jp)%x(3)

      
      p(ipl)%x(1:3) = ( extent(1:3) - offset(1:3) )*( p(ipl)%x(1:3) ) + offset(1:3)

    end do

  end subroutine cube
  
  
  subroutine disk(p)
    use module_globals, only: extent,offset,my_rank,n_ranks,Volume,tnp,radius,root
    implicit none
!
    type(t_particle), allocatable, intent(inout) :: p(:)
    integer(kind_particle)                       :: ipl,np,n
    real(kind_particle)                          :: r,r0,theta,rtnp

    if (root)   write(*,*)            "== ...Loading Particle's Position "

    r0     = radius

    np     = size(p, kind=kind_particle)
    rtnp   = real(tnp, kind=kind_particle)
    
      
    ipl    = 1
    
    do while (ipl .le. np)
        
        call random(p(ipl)%x)
        
        p(ipl)%x(1) = r0*( two*p(ipl)%x(1)-one )
        p(ipl)%x(2) = r0*( two*p(ipl)%x(2)-one )
        p(ipl)%x(3) = zero
        
        if ( p(ipl)%x(1)**2 + p(ipl)%x(2)**2 .le. r0**2 ) ipl = ipl+1
	
    end do

  end subroutine disk
  
  subroutine uniform_disk(p)
    use module_globals, only: extent,offset,my_rank,n_ranks,Volume,tnp,radius,root
    implicit none
!
    type(t_particle), allocatable, intent(inout) :: p(:)
    integer(kind_particle)                       :: ipl,np,n!,rc
    real(kind_particle)                          :: r,r0,theta,alpha,r1,rtnp
!    character(100)                                  :: filename_i
    if (root)   write(*,*)            "== ...Loading Particle's Position "

    r0     = radius

    np     = size(p, kind=kind_particle)
    rtnp   = real(tnp, kind=kind_particle)


    alpha = pi*(three - sqrt(five) )

!    write(filename_i,'("part_",i6.6,".dat")')  my_rank
!    open(unit=rc,file=filename_i,form='formatted',status='unknown',position='append')
    do ipl = 1, np

      theta = real(ipl, kind=kind_particle)*alpha 
      call random( p(ipl)%x )
      r1    = dot_product( p(ipl)%x, p(ipl)%x )
      r     = r0*sqrt( r1*( real( (my_rank+1) * (tnp / n_ranks) , kind=kind_particle) )/rtnp )
!      r     = r0*sqrt( ( real(my_rank * (tnp / n_ranks) + ipl , kind=kind_particle) )/rtnp )
!      write(rc,*) ipl,r,sqrt( ( real(my_rank * (tnp / n_ranks) + ipl , kind=kind_particle) )/rtnp ),sqrt( p(ipl)%x(1)*( real( (my_rank+1) * (tnp / n_ranks) , kind=kind_particle) )/rtnp )        
      p(ipl)%x(1)   = r*cos(theta)
      p(ipl)%x(2)   = r*sin(theta)
      p(ipl)%x(3)   = zero


    end do
!    close (rc )

  end subroutine uniform_disk
   
  
  subroutine cylinder(p)
    use module_globals, only: extent,offset,root,tnp,radius!,my_rank,n_ranks,Volume,root
!    use encap         , only: field_grid_t
    implicit none
!
    type(t_particle), allocatable, intent(inout) :: p(:)
!    type(field_grid_t)           , intent(in)    :: field_grid
    integer(kind_particle)                       :: ipl,np,n!,jp
    real(kind_particle)                          :: r,r0,theta,rtnp

    if (root)   write(*,*)            "== ...Loading Particle's Position "

    r0     = radius

    np     = size(p, kind=kind_particle)
    rtnp   = real(tnp, kind=kind_particle)
    
      
    ipl    = 1
    
    do while (ipl .le. np)
        
        call random(p(ipl)%x)
        
        p(ipl)%x(1) = r0*( two*p(ipl)%x(1)-one )
        p(ipl)%x(2) = r0*( two*p(ipl)%x(2)-one )
        
!        jp          = mod(ipl,field_grid%nl) + 1
        
        p(ipl)%x(3) = ( extent(3) - offset(3) )*( p(ipl)%x(3) ) + offset(3)
        
        if ( p(ipl)%x(1)**2 + p(ipl)%x(2)**2 .le. r0**2 ) ipl = ipl+1
	
    end do

  end subroutine cylinder
  
  
  subroutine clean_fields(p)
    use module_globals, only: root
    implicit none
!
    type(t_particle), allocatable, intent(inout) :: p(:)
    integer(kind_particle)                       :: ipl,np
     
    
    if (root)   write(*,*)            "== ...Cleaning Fields  "

    np     = size(p, kind=kind_particle)

    do ipl = 1, np


      p(ipl)%results%E   = zero
      p(ipl)%results%pot = zero
      p(ipl)%results%A   = zero
      p(ipl)%results%dxA = zero
      p(ipl)%results%dyA = zero
      p(ipl)%results%dzA = zero
!      p(ipl)%results%dxxA= zero
!      p(ipl)%results%dxyA= zero
!      p(ipl)%results%dyyA= zero
      p(ipl)%results%B   = zero
      p(ipl)%results%J   = zero
      p(ipl)%results%Jirr= zero
      p(ipl)%work        = one


    end do


  end subroutine clean_fields
  
  
  
  subroutine iclean_fields(p)
    use module_globals, only: root
    implicit none
!
    type(t_particle), intent(inout) :: p
     

      p%results%E   = zero
      p%results%pot = zero
      p%results%A   = zero
      p%results%dxA = zero
      p%results%dyA = zero
      p%results%dzA = zero
!      p%results%dxxA= zero
!      p%results%dxyA= zero
!      p%results%dyyA= zero
      p%results%B   = zero
      p%results%J   = zero
      p%results%Jirr= zero
      p%work        = one


  end subroutine iclean_fields
  
  
  subroutine velocity_profile(p)
   use module_globals, only: nsp,ns_tot,percentages,my_rank,n_ranks,tnp,vth,vdrift,uth,udrift,wth,wdrift,root,nsp
    implicit none
    include 'mpif.h'
!
    type(t_particle), allocatable, intent(inout) :: p(:)
    integer(kind_particle)                       :: ip,np,n,j1,species,count_species(1:nsp),rc
    real(kind_particle)                          :: v_th(1:3),v_drift(1:3),sumV

    np     = size(p, kind=kind_particle)
    if (root)   write(*,*)            "== ...Loading Particle's Velocity  "
    count_species   = 0
    sumV            = zero
    
    do ip = 1,np
                species                            = p(ip)%label
                count_species(species)             = count_species(species)+1
                v_th                               = (/ uth(species),vth(species),wth(species) /)
                v_drift                            = (/ udrift(species),vdrift(species),wdrift(species) /)                
                call random_boltzmann(p(ip)%data%v,v_th,v_drift)
                sumV                               = sumV +  p(ip)%data%v(3)
    end do
    
    call MPI_ALLREDUCE(MPI_IN_PLACE , count_species(1:nsp)  , nsp, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
    call MPI_ALLREDUCE(MPI_IN_PLACE , sumV                  , 1  , MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)

    if ( root ) then
        write(*,*) 'Number of particles per species: ',count_species(1:nsp)
        write(*,*) 'Sum                            : ',sum(count_species)
        write(*,*) 'Sum of currents                : ',sumV       
    endif

  end subroutine velocity_profile
  
  
  subroutine charge_mass_label(p)
    use module_globals, only: nsp,ns_tot,charge_init,mass_init,percentages,my_rank,n_ranks,tracks,tnp,nsp,root
    implicit none
!
    type(t_particle), allocatable, intent(inout) :: p(:)
    integer(kind_particle)                       :: ip,ipl,np,n,delta,j1,species

    np     = size(p, kind=kind_particle)
    
!    rdelta = ( real( np/percentages(ns_tot+1) , kind=kind_particle ) )
    if (root)   write(*,*)            "== ...Loading Particle's Charge  "
!    if ( mod(rdelta,one) .gt. zero ) then
!        write(*,*) my_rank,np,percentages(ns_tot+1),rdelta
!        call exit(1)
!    endif
    
    delta = ( int( np/percentages(ns_tot+1) , kind=kind_particle ) )
    ipl   = 0
    
    do ip = 1, delta

        do species = 1,nsp
            
            do j1 = 1,percentages(species)
                ipl                                = ipl+1
                p(ipl)%data%q                      = real(charge_init(species), kind=kind_particle)
                p(ipl)%data%m                      = real(  mass_init(species), kind=kind_particle)

!                if ( tracks .eq. 0 ) then 
                    ! track single particle
!                    p(ipl)%label       = my_rank * (tnp / n_ranks) + ipl 
!                else if ( tracks .eq. 1 ) then
!                     track per species
                    p(ipl)%label       = species
!                endif

            enddo
            
        enddo
        
    end do


  end subroutine charge_mass_label
  
  subroutine perturbations(p)
    use module_globals, only: nsp,ns_tot,charge_init,mass_init,percentages,my_rank,n_ranks,tracks,x_pert,y_pert,z_pert,&
                              u_pert,v_pert,w_pert,extent,nsp,root
    implicit none
!
    type(t_particle), allocatable, intent(inout) :: p(:)
    integer(kind_particle)                       :: ip,np,n,j1,species
    real(kind_particle)                          :: dL(1:3)

    np     = size(p, kind=kind_particle)
    dL     = extent
    if (root)   write(*,*)            "== ...Loading Perturbation  "
    
    do ip = 1, np
        
                species                           = p(ip)%label 
                
                p(ip)%x(1)                        = p(ip)%x(1) + x_pert(species)*cos( two*pi/dL(1)*p(ip)%x(1)  )
                p(ip)%x(2)                        = p(ip)%x(2) + y_pert(species)*cos( two*pi/dL(2)*p(ip)%x(2)  )
                if ( dL(3) .gt. zero ) then
                p(ip)%x(3)                        = p(ip)%x(3) + z_pert(species)*cos( two*pi/dL(3)*p(ip)%x(3)  )
                endif

                p(ip)%data%v(1)                   = p(ip)%data%v(1) + u_pert(species)*cos( two*pi/dL(1)*p(ip)%x(1)  )
                p(ip)%data%v(2)                   = p(ip)%data%v(2) + v_pert(species)*cos( two*pi/dL(2)*p(ip)%x(2)  )
                if ( dL(3) .gt. zero ) then
                p(ip)%data%v(3)                   = p(ip)%data%v(3) + v_pert(species)*cos( two*pi/dL(3)*p(ip)%x(3)  )
                endif
       
    end do


  end subroutine perturbations
  
  

end module module_init
