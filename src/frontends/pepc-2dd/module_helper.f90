! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
!
! Copyright (C) 2002-2017 Juelich Supercomputing Centre,
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

!!!!!!!!!!!!!!!!!!!!
!! helper module
!!!!!!!!!!!!!!!!!!!!

module helper

  use module_pepc_kinds
  use module_pepc_types
  use module_interaction_Specific_types
  use module_globals
  use module_shortcut
  implicit none

!  ! MPI variables
!  integer :: my_rank, n_ranks
!  logical :: root
!
!  ! time variables
!  real*8 :: dt
!  integer :: step
!
!  ! control variables
!  integer :: nt               ! number of timesteps
!  integer(kind_particle) :: tnp ! total number of particles
!  integer(kind_particle) :: np  ! local number of particles
!  integer(kind_particle) :: initial_setup ! initial setup
!  logical :: particle_output  ! turn vtk output on/off
!  logical :: domain_output    ! turn vtk output on/off
!  logical :: particle_test    ! turn direct summation on/off
!  integer :: diag_interval    ! number of timesteps between all diagnostics and IO

  integer, parameter :: particle_direct = 10000 ! number of particle for direct summation

  type(t_particle), allocatable :: particles(:)
  real*8, allocatable           :: direct_L2(:)
  real*8, allocatable           :: direct_EL2Pot(:)
  real*8, allocatable           :: direct_EL2rho(:)
  real*8, allocatable           :: direct_EL2A(:)
  real*8, allocatable           :: direct_EL2B(:)

  integer(kind_particle)        :: next_label


  contains

  subroutine set_parameter()

    use module_pepc
    use module_interaction_specific , only : theta2, eps2, include_far_field_if_periodic
    use module_mirror_boxes         , only : periodicity, spatial_interaction_cutoff
    use module_shortcut             , only : half,pi   
    use module_utils                , only: create_directory
    implicit none

    integer, parameter :: fid = 12
    character(255)     :: para_file
    logical            :: read_para_file,bool_exp, bool_nk,bool_pic,firstcall  = .true.
    character(12), parameter                        :: part_dir  = 'particles/'
    character(12), parameter                        :: field_dir = 'fields/'

    namelist /pepc2dd/ tnp, nt, dt, particle_output, domain_output, unique_species, &
       diag_interval, periodicity_particles,  spatial_interaction_cutoff,x_distribution,ischeme,&
       restart_file,restart_step,normal,nsp,vth,vdrift,uth,udrift,wth,wdrift,folder,we,vmax,x_pert,y_pert,z_pert,u_pert, &
       v_pert,w_pert,charge_init,mass_init,percentages,tracks,radius,newmark_x,newmark_v,newmark_Es,newmark_Ei,newmark_B, &
       newmark_g,dA_1,dA__1,dA_0,B0,load


    ! set default parameter values
    load            = .false.   
    x_distribution  = 2
    v_distribution  = 1
    tnp             = 1441
    nt              = 20
    dt              = 1e-3
    particle_output = .true.
    ischeme         = "leapfrog"
    folder          = "data"
    periodicity_particles     = .false.
    restart_file    = ""
    restart_step    = 1

    diag_interval   = 2
    nsp             = 1
    unique_species  = 1
    tracks          = 0
    normal          = 0
    
    uth             = one
    vth             = zero
    wth             = zero
    udrift          = zero
    vdrift          = zero
    wdrift          = zero
    x_pert          = zero
    y_pert          = zero
    z_pert          = zero
    u_pert          = zero
    v_pert          = zero
    w_pert          = zero
    !!! Assuming Trapezoidal-Rule
    newmark_x       = one
    newmark_v       = one
    newmark_Es      = one
    newmark_Ei      = half
    newmark_B       = half
    newmark_g       = half
    !! Total Derivative of A
    dA_1            = one 
    dA__1           = zero
    dA_0            = one
    B0              = zero
       
    percentages     = (/ 1,0,0,0,0,1/)
    charge_init     = (/-one,zero,zero,zero,zero  /)
    mass_init       = (/ one,zero,zero,zero,zero  /)
    
    radius          = one
    
    
    etilde          =  one
    phitilde        =  one
    btilde          =  one
    atilde          =  one
    jtilde          =  one
    rhotilde        =  one
    lorentz_tilde   =  one
    vtilde          =  one
    we              =  one
    vmax            =  one
    

    ! read in namelist file
    call pepc_read_parameters_from_first_argument(read_para_file, para_file)

    if (read_para_file) then
      if(root) write(*,'(a)') " == reading parameter file, section pepc-2dd : ", para_file
      open(fid,file=para_file)
      read(fid,NML=pepc2dd)
      close(fid)
    else
      if(root) write(*,*) " == no param file, using default parameter "
    end if

    periodicity = [ .false., .false., .false.]
    if (periodicity_particles) periodicity = [ .true., .true., .false.]
    tnp = tnp

            
    if (firstcall) then
        call create_directory(trim(folder))
        call create_directory(trim(folder)//trim(part_dir))
        call create_directory(trim(folder)//trim(field_dir))
        firstcall = .false.        
    endif
            
            
            
    if(root) then
      write(*,'(a,i12)')    " == total number of particles : ", tnp
      write(*,'(a,i12)')    " == number of time steps      : ", nt
      write(*,'(a,es12.4)') " == time step                 : ", dt
      write(*,*)            "== time integration scheme   : ",  ischeme
      write(*,'(a,es12.4)') " == theta2                    : ", theta2
      write(*,'(a,es12.4)') " == epsilon2                  : ", eps2
      write(*,'(a,l12)')    " == periodicity               : ", periodicity!
      write(*,'(a,i12)')    " == diag & IO interval        : ", diag_interval
      write(*,'(a,l12)')    " == particle output           : ", particle_output
      write(*,'(a,l12)')    " == domain output             : ", domain_output
      write(*,*) ""

    end if


    bool_exp = ischeme .eq. "leapfrog" .or. ischeme .eq. "explicit" .or. ischeme .eq."euler_method".or. ischeme .eq."euler_method_electrostatic"&
    .or. ischeme .eq. "explicit_verlet".or. ischeme .eq. "euler_method3d"
    bool_pic = ischeme .eq. "trapezoidal_picard" .or.ischeme .eq."trapezoidal_picard_electrostatic" 

    adv = -1

    if (bool_exp) then
        adv = 0
    elseif (bool_pic) then
        adv = 1
    elseif (bool_nk) then
        adv = 2
    endif


  end subroutine

  subroutine write_restart(p)
    use mpi
    implicit none
    type(t_particle), allocatable, intent(in)    :: p(:)
    integer(kind_particle)                       :: rc,np
    character(255)                               :: filename
    character(*), parameter                      :: part_dir = "particles/"
    integer(kind = MPI_OFFSET_KIND)              :: mpi_disp,my_offset
    integer                                      :: fh, mpi_err
    integer         , dimension(MPI_STATUS_SIZE) :: mpi_stat
!    real(kind = 8)   , dimension(:), allocatable :: real8_buf
    
    
    
!    if (root)  then
    
        np    = size(particles, kind=kind_particle)
!        if ( allocated(real8_buf) ) deallocate(real8_buf)
!        allocate( real8_buf( np ) )
        my_offset = 0 
      
        write(filename,'(a,"restart_label","_",i6.6,".dat")') trim(folder)//trim(part_dir), my_rank 
        call write_mpi_int(filename,p(:)%label)
        write(filename,'(a,"restart_q","_",i6.6,".dat")') trim(folder)//trim(part_dir), my_rank 
        call write_mpi_real(filename,p(:)%data%q)
        write(filename,'(a,"restart_m","_",i6.6,".dat")') trim(folder)//trim(part_dir), my_rank 
        call write_mpi_real(filename,p(:)%data%m)
        write(filename,'(a,"restart_x","_",i6.6,".dat")') trim(folder)//trim(part_dir), my_rank 
        call write_mpi_real(filename,p(:)%x(1))
        write(filename,'(a,"restart_y","_",i6.6,".dat")') trim(folder)//trim(part_dir), my_rank 
        call write_mpi_real(filename,p(:)%x(2))
        write(filename,'(a,"restart_z","_",i6.6,".dat")') trim(folder)//trim(part_dir), my_rank 
        call write_mpi_real(filename,p(:)%x(3))
        write(filename,'(a,"restart_vx","_",i6.6,".dat")') trim(folder)//trim(part_dir), my_rank 
        call write_mpi_real(filename,p(:)%data%v(1))
        write(filename,'(a,"restart_vy","_",i6.6,".dat")') trim(folder)//trim(part_dir), my_rank 
        call write_mpi_real(filename,p(:)%data%v(2))
        write(filename,'(a,"restart_vz","_",i6.6,".dat")') trim(folder)//trim(part_dir), my_rank 
        call write_mpi_real(filename,p(:)%data%v(3))
        
        contains
        
        subroutine write_mpi_real(filename,p)
        implicit none

        character(*)                    , intent(in) :: filename
        real(kind_physics), dimension(:), intent(in) :: p
     
        real(kind = 8)   , dimension(:), allocatable :: real8_buf
        
        
        if ( allocated(real8_buf) ) deallocate(real8_buf)
        allocate( real8_buf( np ) )
        
        real8_buf(:) = real(p(:), kind = 8)
            
             
        call mpi_file_open(MPI_COMM_WORLD, filename, ior(MPI_MODE_RDWR, MPI_MODE_CREATE), MPI_INFO_NULL, fh, mpi_err)
      
!        if (mpi_err .ne. MPI_SUCCESS) DEBUG_ERROR(*, 'write_restart(): I/O error for ', filename)
        
        real8_buf(:) = real(p(:), kind = 8)
        
        call mpi_file_set_view(    fh, 0_MPI_OFFSET_KIND, MPI_DOUBLE_PRECISION , MPI_DOUBLE_PRECISION,'native'           , MPI_INFO_NULL      , mpi_err)  
        call mpi_file_write_at_all(fh, my_offset        , real8_buf, int(np, kind = kind_default), MPI_DOUBLE_PRECISION, mpi_stat, mpi_err)
        call mpi_file_sync(fh, mpi_err)
        call mpi_file_close(fh, mpi_err)

        deallocate(real8_buf)
        
            
        end subroutine write_mpi_real
        
        
        subroutine write_mpi_int(filename,p)
        implicit none

        character(*)                       , intent(in) :: filename
        integer(kind_physics), dimension(:), intent(in) :: p
                
        integer(kind = 8)   , dimension(:), allocatable :: real8_buf
        
        
        if ( allocated(real8_buf) ) deallocate(real8_buf)
        allocate( real8_buf( np ) )
                    
             
        call mpi_file_open(MPI_COMM_WORLD, filename, ior(MPI_MODE_RDWR, MPI_MODE_CREATE), MPI_INFO_NULL, fh, mpi_err)
      
!        if (mpi_err .ne. MPI_SUCCESS) DEBUG_ERROR(*, 'write_restart(): I/O error for ', filename)
        
        real8_buf(:) = int(p(:), kind = 8)
        
        call mpi_file_set_view(    fh, 0_MPI_OFFSET_KIND, MPI_INTEGER8 , MPI_INTEGER8,'native'           , MPI_INFO_NULL      , mpi_err)  
        call mpi_file_write_at_all(fh, my_offset        , real8_buf, int(np, kind = kind_default), MPI_INTEGER8, mpi_stat, mpi_err)
        call mpi_file_sync(fh, mpi_err)
        call mpi_file_close(fh, mpi_err)

        deallocate(real8_buf)
        
            
        end subroutine write_mpi_int
        
  end subroutine write_restart



subroutine init_particles(p,field_grid)
    use module_interaction_specific, only : force_law
    use module_pepc
    use module_pepc_kinds
    use module_init
    use module_mirror_boxes, only:t_lattice_1,t_lattice_2,t_lattice_3,LatticeOrigin
    use module_tool        , only:par_rand,copy_particle,scramble_particles,icopy_particle
    use module_shortcut    , only:half,pi,oneoverfourpi,oneoverpi
    use module_utils       , only: create_directory
    use module_globals     , only: extent,pold,poldold
    use encap
    
    implicit none
    include 'mpif.h'
    type(t_particle), allocatable, intent(inout) :: p(:)
    type(field_grid_t)           , intent(in)    :: field_grid
    integer(kind_particle)                       :: ip,rc,jp
    real(kind_particle)                          :: dummy,nd,lambda,rtnp,gl
    
    character(100)                               :: filename_i
!    character(12), parameter                     :: part_dir = 'prova/'

    if(my_rank.eq.0) write(*,'(a)') " == [init] init particles "

    ! set initially number of local particles
    np = tnp / n_ranks
!    write(*,*) "Info== ", my_rank, n_ranks,np, tnp
    if(my_rank.eq.(n_ranks-1)) np = np + MOD(tnp, int(n_ranks, kind=kind_particle))

    if (allocated(particles) ) deallocate(particles)
    allocate(particles(np), stat=rc)
    if(rc.ne.0) write(*,*) " === particle allocation error!"
    
    if (allocated(pold) ) deallocate(pold)
    allocate(pold(np), stat=rc)
    if(rc.ne.0) write(*,*) " === particle allocation error!"
    
    if (allocated(poldold) ) deallocate(poldold)
    allocate(poldold(np), stat=rc)
    if(rc.ne.0) write(*,*) " === particle allocation error!"

    allocate(direct_L2(np), stat=rc)
    if(rc.ne.0) write(*,*) " === direct_L2 allocation error!"
    direct_L2 = -1.0_8



    if ( extent(3) .eq. zero ) then
        ixdim                = int( 2, kind = kind_dim )
        force_law            = 2   
    else if ( extent(3) .gt. zero ) then
        ixdim                = int( 3, kind = kind_dim )
        force_law            = 3   
    endif

    call pepc_prepare(ixdim)
    call create_directory(trim(folder))

    ! set random seed
    dummy = par_rand(my_rank)
    
    
    
    t_lattice_1(1)      = extent(1)
    t_lattice_2(2)      = extent(2)
    LatticeOrigin(1:3)  = offset(1:3)

    if ( extent(3).gt. zero )   t_lattice_3(3)      = extent(3)
    
    
    
    Volume          = extent(1)*extent(2)
    if ( extent(3).gt. zero )  Volume          = extent(3)*Volume
    if ( ( x_distribution .eq.2).or.( x_distribution .eq. 3) )  Volume          = radius**2*pi 
    if ( ( x_distribution .eq.5)                             )  Volume          = radius**2*pi*(extent(3) - offset(3)) 

    
    
    rtnp            = real(tnp, kind=kind_particle)         ! total nomber of particles - doubleprecision
    nd              = rtnp/Volume                           ! number density
    
   
    select case (normal)
          case (0)  !
            
            
            qtilde          =  unique_species*we/nd
            mtilde          =  unique_species*we/nd
            
            vtilde          =  one    
            lambda          =  one
         
            rhotilde        =  lambda*two
            jtilde          =  lambda*two
            etilde          =  lambda*half*oneoverpi
            phitilde        =  lambda*half*oneoverpi
            btilde          =  lambda*half*oneoverpi
            atilde          =  lambda*half*oneoverpi
            lorentz_tilde   =  vtilde
            wpe             =  one

          case (1)  !
            
            qtilde          =  unique_species*we/nd
            mtilde          =  unique_species*we/nd
            
            vtilde          =  vmax    
            lambda          =  one
            
            rhotilde        =  lambda*two
            jtilde          =  lambda*two

            etilde          =  lambda*half*oneoverpi
            phitilde        =  lambda*half*oneoverpi
            btilde          =  lambda*half*oneoverpi
            atilde          =  lambda*half*oneoverpi
            lorentz_tilde   =  vtilde
            wpe             =  one
            
          case (2)  !

            qtilde          =  unique_species*we/nd
            mtilde          =  unique_species*we/nd
            
            lambda          =  one
            vtilde          =  one  
            rhotilde        =  one
            jtilde          =  one

            etilde          =  lambda
            phitilde        =  lambda
            btilde          =  lambda
            atilde          =  lambda
            lorentz_tilde   =  vtilde
            wpe             =  one
    
          case (3)  ! CGS System

            qtilde          =  qe
            mtilde          =  me
            if ( unique_species .eq. 2 ) mtilde          =  mi
            
            vtilde          =  c
            
            lambda          =  two
            
            rhotilde        =  lambda
            jtilde          =  lambda

            etilde          =  lambda
            phitilde        =  lambda
            btilde          =  lambda
            atilde          =  lambda
            lorentz_tilde   =  vtilde
            wpe             =  sqrt( four*pi*qtilde**2*nd/me  )           ! electron plasma frequency
            
         case (4)  ! From CGS TO SI
             
            qtilde          =  qe
            mtilde          =  me 
            if ( unique_species .eq. 2 ) mtilde          =  mi

            qtilde          =  qtilde/sqrt( alpha*beta*oneoverfourpi/epsilon0 )
            mtilde          =  mtilde/( beta/alpha**2 )
            vtilde          =  c / alpha             
            lambda          =  two
            
            rhotilde        =  lambda
            jtilde          =  lambda

            etilde          =  lambda*oneoverfourpi/epsilon0 
            phitilde        =  lambda*oneoverfourpi/epsilon0 
            btilde          =  lambda*vtilde*oneoverfourpi*mu0
            atilde          =  lambda*vtilde*oneoverfourpi*mu0
            lorentz_tilde   =  one
            wpe             =  sqrt( qtilde**2*nd/(me*epsilon0) )           ! electron plasma frequency

          case default
              ! CGS System
              
            qtilde          =  qe
            mtilde          =  me
            if ( unique_species .eq. 2 ) mtilde          =  mi

            vtilde          =  c
            
            lambda          =  two
            
            rhotilde        =  lambda
            jtilde          =  lambda

            etilde          =  lambda
            phitilde        =  lambda
            btilde          =  lambda
            atilde          =  lambda
            lorentz_tilde   =  vtilde
            wpe             =  sqrt( four*pi*qtilde**2*nd/me  )           ! electron plasma frequency

    end select

    dt = dt/wpe
    
    
    if ( .not.load ) then
        
        call charge_mass_label(p)

        select case (x_distribution) 
              case (1)  !  Rectangle/Square
                  call rect(p)
              case (2)  ! Disc
                  call disk(p)
              case (3)  ! Disc
                  call uniform_disk(p)
              case (4) ! Random Cube
                  call cube(p)   
              case (5) ! Random Cylinder with axes in z
                  call cylinder(p)   
              case default
                  call rect(p)

        end select
    
        
        call velocity_profile(p)
        call scramble_particles(p)
!        call perturbations(p)
        
        
!        do ip = 1,np
!
!            p(ip)%data%q      = p(ip)%data%q*qtilde
!            p(ip)%data%m      = p(ip)%data%m*mtilde
!            p(ip)%data%v(1:3) = p(ip)%data%v(1:3)*vtilde
!
!            gl           = dot_product( p(ip)%data%v/vtilde, p(ip)%data%v/vtilde )
!            if (gl .gt. one ) write(*,*) "Warning -- Lorentz Factor Bigger than 1!!" 
!            gl           = one/sqrt( one - gl )
!            p(ip)%data%g = gl
!            p(ip)%data%v = p(ip)%data%v*gl
!
!        enddo
        
    else 
        
        call load_file(p)
        
   endif

    

!    call copy_particle(p,pold,np)
!    call copy_particle(p,poldold,np)
!    
!    call clean_fields(pold)
    call clean_fields(p)
!    call clean_fields(poldold)
!    
    
    
    do ip = 1,np
        
        call icopy_particle(p(ip),pold(ip))
        pold(ip)%x(1:2)         = p(ip)%x(1:2) - dt*p(ip)%data%v(1:2)
        pold(ip)%x(3)           = p(ip)%x(3)
!        poldold(ip)%x(1:2)      = p(ip)%x(1:2) - two*dt*p(ip)%data%v(1:2)
                
    enddo
    
!    call write_particles_vtk(p, 0, 0.0_8)
!    call write_particles_vtk(pold, 10, 10.0_8)
    
    call pepc_particleresults_clear(pold)
    call pepc_grow_tree(pold)
    call pepc_traverse_tree(pold)
    call pepc_restore_particles(pold)
    call pepc_timber_tree()
    call normalize(np, pold)
    
    
!    call exit(1)
!    
!    call pepc_particleresults_clear(poldold)
!    call pepc_grow_tree(poldold)
!    call pepc_traverse_tree(poldold)
!    call pepc_restore_particles(poldold)
!    call pepc_timber_tree()
!    call normalize(np, poldold)
    
    
     next_label = tnp+1

  end subroutine init_particles

  subroutine time_step(p)
  use module_shortcut, only: zero,one
  implicit none
  include 'mpif.h'
    type(t_particle), allocatable, intent(in)    :: p(:)
    integer(kind_particle)                       :: ip,np,rc
    real(kind_particle)                          :: wpe,wgf,alpha
    logical                                      :: flag
    
    wpe = zero
    wgf = zero
    
    np = size(p, kind=kind_particle)
    
    do ip = 1,np
        
        wpe     = wpe + p(ip)%data%q**2/p(ip)%data%m 
        wgf     = wgf + dot_product( p(ip)%results%B , p(ip)%results%B ) 
                
    enddo
    
    call MPI_ALLREDUCE(MPI_IN_PLACE, wpe, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
    call MPI_ALLREDUCE(MPI_IN_PLACE, wgf, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        
!    wpe       = abs( wpe/real(tnp, kind=kind_particle) )!sqrt(wpe)
    wgf       = sqrt(wgf)!/real(tnp, kind=kind_particle)
    
    ip      = 1
    flag    = .true.
    alpha   = abs( p(ip)%data%q/p(ip)%data%m )
    
    
    do  while ( flag .and. ( ip .lt. np ) )
        
        if ( p(ip)%data%q .lt. zero ) then 
            flag = .false.
            alpha   = abs( p(1)%data%q/p(1)%data%m )
        endif
        ip  =  ip + 1
                
    enddo

    wgf =  alpha*wgf

    
    if (root) write(*,*) "electron/magnetic frequency =  ", wpe,wgf
    
    if ( wpe .gt. zero )  then 
        wpe = one/wpe
    else
        wpe = zero
    endif
    
    if ( wgf .gt. zero ) then 
        wgf = one/wgf
    else
        wgf = zero
    endif 
    
        
    if (root) write(*,*) "electron/magnetic period = ", wpe,wgf
    
  end subroutine time_step
  
  subroutine periodic_particles(np,p)
!    use module_mirror_boxes
    implicit none

    type(t_particle), allocatable, intent(inout) :: p(:)
    integer(kind_particle)       , intent(in)    :: np
    integer(kind_particle)                       :: ip,jp

    if(root) write(*,'(a)') " ======  periodicity on particles "
    do ip=1, np
        do jp = 1,3
            if ( p(ip)%x(jp) .lt. offset(jp) )              p(ip)%x(jp) = p(ip)%x(jp) + extent(jp)
            if ( p(ip)%x(jp) .gt. offset(jp) + extent(jp) ) p(ip)%x(jp) = mod( p(ip)%x(jp) , extent(jp) )
        enddo
    end do
  end subroutine periodic_particles
  
  subroutine iperiodic_particles(p)
!    use module_mirror_boxes
    implicit none

    type(t_particle)             , intent(inout) :: p
    integer(kind_particle)                       :: jp

        do jp = 1,3
            if ( p%x(jp) .lt. offset(jp) )              p%x(jp) = p%x(jp) + extent(jp)
            if ( p%x(jp) .gt. offset(jp) + extent(jp) ) p%x(jp) = mod( p%x(jp) , extent(jp) )
        enddo
        
  end subroutine iperiodic_particles


  real*8 function get_time()
      implicit none
      include 'mpif.h'

      get_time = MPI_WTIME()

  end function get_time

  subroutine normalize(np,p)
    implicit none

    type(t_particle), allocatable, intent(inout) :: p(:)
    integer(kind_particle)       , intent(in)    :: np

    integer(kind_particle)                       :: ip
!    real(kind_particle) ::err= zero

    do ip = 1,np

        p(ip)%results%E(1:3)    = p(ip)%results%E(1:3)      *etilde
        p(ip)%results%B(1:3)    = p(ip)%results%B(1:3)      *btilde/vtilde
        p(ip)%results%A(1:3)    = p(ip)%results%A(1:3)      *atilde/vtilde
        p(ip)%results%dxA(1:3)  = p(ip)%results%dxA(1:3)    *atilde/vtilde
        p(ip)%results%dyA(1:3)  = p(ip)%results%dyA(1:3)    *atilde/vtilde
        p(ip)%results%dzA(1:3)  = p(ip)%results%dzA(1:3)    *atilde/vtilde
        p(ip)%results%J(1:3)    = p(ip)%results%J(1:3)      *jtilde
        p(ip)%results%Jirr(1:3) = p(ip)%results%Jirr(1:3)   *jtilde
        p(ip)%results%pot       = p(ip)%results%pot         *phitilde
!        p(ip)%results%rho       = p(ip)%results%rho         *rhotilde
        
!        err = err + ( p(ip)%data%v(1)/p(ip)%data%g*p(ip)%results%pot-p(ip)%results%A(1) )**2

    enddo
!    write(*,*) "errore2: ",err
  end subroutine normalize
  
  subroutine inormalize(p)
    implicit none

    type(t_particle), intent(inout) :: p

        p%results%E(1:3)    = p%results%E(1:3)      *etilde
        p%results%B(1:3)    = p%results%B(1:3)      *btilde/vtilde
        p%results%A(1:3)    = p%results%A(1:3)      *atilde/vtilde
        p%results%dxA(1:3)  = p%results%dxA(1:3)    *atilde/vtilde
        p%results%dyA(1:3)  = p%results%dyA(1:3)    *atilde/vtilde
        p%results%dzA(1:3)  = p%results%dzA(1:3)    *atilde/vtilde
        p%results%J(1:3)    = p%results%J(1:3)      *jtilde
        p%results%Jirr(1:3) = p%results%Jirr(1:3)   *jtilde
        p%results%pot       = p%results%pot         *phitilde
!        p(ip)%results%rho       = p(ip)%results%rho         *rhotilde

  end subroutine inormalize
  
  
  subroutine test_particles()
        use module_pepc_types
        use module_directsum
        use module_globals             , only : folder
        use module_interaction_specific, only : theta2
        use module_tool                , only : random,cross_product
        implicit none
        include 'mpif.h'

        integer(kind_particle), allocatable   :: tindx(:)
        real(kind_particle), allocatable      :: trnd(:)
        type(t_particle_results), allocatable :: trslt(:)
        integer(kind_particle)                :: tn, tn_global, ti
        integer                               :: rc
        
        real(kind_particle)                   :: phipp      ,phitree      ,phi_err ,phi_norm
        real(kind_particle)                   :: Epp(1:3)   ,Etree(1:3)   ,E_err   ,E_norm
        real(kind_particle)                   :: Expp(1:3)  ,Extree(1:3)  ,Ex_err  ,Ex_norm
        real(kind_particle)                   :: Eypp(1:3)  ,Eytree(1:3)  ,Ey_err  ,Ey_norm
        real(kind_particle)                   :: Bpp(1:3)   ,Btree(1:3)   ,B_err   ,B_norm
        real(kind_particle)                   :: Jpp(1:3)   ,Jtree(1:3)   ,J_err   ,J_norm
        real(kind_particle)                   :: Jirpp(1:3) ,Jirtree(1:3) ,Jir_err ,Jir_norm
        real(kind_particle)                   :: App(1:3)   ,Atree(1:3)   ,A_err   ,A_norm
        real(kind_particle)                   :: Axpp(1:3)  ,Axtree(1:3)  ,Ax_err  ,Ax_norm
        real(kind_particle)                   :: Aypp(1:3)  ,Aytree(1:3)  ,Ay_err  ,Ay_norm
        real(kind_particle)                   :: Azpp(1:3)  ,Aztree(1:3)  ,Az_err  ,Az_norm
        real(kind_particle)                   :: Axxpp(1:3) ,Axxtree(1:3) ,Axx_err ,Axx_norm
        real(kind_particle)                   :: Axypp(1:3) ,Axytree(1:3) ,Axy_err ,Axy_norm
        real(kind_particle)                   :: Ayypp(1:3) ,Ayytree(1:3) ,Ayy_err ,Ayy_norm
        real(kind_particle)                   :: Fdarpp(1:3),Fdartree(1:3),Fdar_err,Fdar_norm
        real(kind_particle)                   :: FElpp(1:3) ,FEltree(1:3) ,FEl_err ,FEl_norm
        real(kind_particle)                   :: v(1:3)     ,m,q,ta,tb,gradpp(1:3),gradtree(1:3)
        

        ta = get_time()

        phi_err     = zero
        E_err       = zero
        B_err       = zero
        J_err       = zero
        Jir_err     = zero
        A_err       = zero
        Ax_err      = zero
        Ay_err      = zero
        Az_err      = zero
        Axx_err     = zero
        Ayy_err     = zero
        Axy_err     = zero
        Fdar_err    = zero
        FEl_err     = zero
        
        phi_norm    = zero
        E_norm      = zero
        B_norm      = zero
        J_norm      = zero
        Jir_norm    = zero
        A_norm      = zero
        Ax_norm     = zero
        Ay_norm     = zero
        Az_norm     = zero
        Axx_norm    = zero
        Ayy_norm    = zero
        Axy_norm    = zero
        Fdar_norm   = zero
        FEl_norm    = zero
        

        tn = tnp!particle_direct / n_ranks
        if(my_rank.eq.(n_ranks-1)) tn = tn + MOD(particle_direct, n_ranks)

        allocate(tindx(tn), trnd(tn), trslt(tn))

        call random(trnd)

        tindx(1:tn) = int(trnd(1:tn) * (np-1)) + 1

        call directforce(particles, tindx, tn, trslt, MPI_COMM_WORLD)

        do ti = 1, tn

          v(1)        = particles(tindx(ti))%data%v(1)
          v(2)        = particles(tindx(ti))%data%v(2)
          v(3)        = particles(tindx(ti))%data%v(3)

          m           = particles(tindx(ti))%data%m
          q           = particles(tindx(ti))%data%q

          phipp       = trslt(ti)%pot
          phitree     = particles(tindx(ti))%results%pot

          Epp(1)      = trslt(ti)%e(1)
          Epp(2)      = trslt(ti)%e(2)
          Epp(3)      = trslt(ti)%e(3)

          Etree(1)    = particles(tindx(ti))%results%e(1)
          Etree(2)    = particles(tindx(ti))%results%e(2)
          Etree(3)    = particles(tindx(ti))%results%e(3)
          
!          Expp(1)     = trslt(ti)%dxE(1)
!          Expp(2)     = trslt(ti)%dxE(2)
!          Expp(3)     = trslt(ti)%dxE(3)
!
!          Extree(1)   = particles(tindx(ti))%results%dxE(1)
!          Extree(2)   = particles(tindx(ti))%results%dxE(2)
!          Extree(3)   = particles(tindx(ti))%results%dxE(3)
!          
!          
!          Eypp(1)     = trslt(ti)%dyE(1)
!          Eypp(2)     = trslt(ti)%dyE(2)
!          Eypp(3)     = trslt(ti)%dyE(3)
!
!          Eytree(1)   = particles(tindx(ti))%results%dyE(1)
!          Eytree(2)   = particles(tindx(ti))%results%dyE(2)
!          Eytree(3)   = particles(tindx(ti))%results%dyE(3)

          App(1)      = trslt(ti)%A(1)
          App(2)      = trslt(ti)%A(2)
          App(3)      = trslt(ti)%A(3)

          Atree(1)    = particles(tindx(ti))%results%A(1)
          Atree(2)    = particles(tindx(ti))%results%A(2)
          Atree(3)    = particles(tindx(ti))%results%A(3)
          
          Axpp(1)     = trslt(ti)%dxA(1)
          Axpp(2)     = trslt(ti)%dxA(2)
          Axpp(3)     = trslt(ti)%dxA(3)

          Axtree(1)   = particles(tindx(ti))%results%dxA(1)
          Axtree(2)   = particles(tindx(ti))%results%dxA(2)
          Axtree(3)   = particles(tindx(ti))%results%dxA(3)
          
          
          Aypp(1)     = trslt(ti)%dyA(1)
          Aypp(2)     = trslt(ti)%dyA(2)
          Aypp(3)     = trslt(ti)%dyA(3)

          Aytree(1)   = particles(tindx(ti))%results%dyA(1)
          Aytree(2)   = particles(tindx(ti))%results%dyA(2)
          Aytree(3)   = particles(tindx(ti))%results%dyA(3)
          
          Azpp(1)     = trslt(ti)%dzA(1)
          Azpp(2)     = trslt(ti)%dzA(2)
          Azpp(3)     = trslt(ti)%dzA(3)

          Aztree(1)   = particles(tindx(ti))%results%dzA(1)
          Aztree(2)   = particles(tindx(ti))%results%dzA(2)
          Aztree(3)   = particles(tindx(ti))%results%dzA(3)
          
!          Axxpp(1)     = trslt(ti)%dxxA(1)
!          Axxpp(2)     = trslt(ti)%dxxA(2)
!          Axxpp(3)     = trslt(ti)%dxxA(3)
!
!          Axxtree(1)   = particles(tindx(ti))%results%dxxA(1)
!          Axxtree(2)   = particles(tindx(ti))%results%dxxA(2)
!          Axxtree(3)   = particles(tindx(ti))%results%dxxA(3)
!          
!          
!          Axypp(1)     = trslt(ti)%dxyA(1)
!          Axypp(2)     = trslt(ti)%dxyA(2)
!          Axypp(3)     = trslt(ti)%dxyA(3)
!
!          Axytree(1)   = particles(tindx(ti))%results%dxyA(1)
!          Axytree(2)   = particles(tindx(ti))%results%dxyA(2)
!          Axytree(3)   = particles(tindx(ti))%results%dxyA(3)
!          
!          
!          Ayypp(1)     = trslt(ti)%dyyA(1)
!          Ayypp(2)     = trslt(ti)%dyyA(2)
!          Ayypp(3)     = trslt(ti)%dyyA(3)
!
!          Ayytree(1)   = particles(tindx(ti))%results%dyyA(1)
!          Ayytree(2)   = particles(tindx(ti))%results%dyyA(2)
!          Ayytree(3)   = particles(tindx(ti))%results%dyyA(3)

          Jirpp(1)    = trslt(ti)%Jirr(1)
          Jirpp(2)    = trslt(ti)%Jirr(2)
          Jirpp(3)    = trslt(ti)%Jirr(3)


          Jirtree(1)  = particles(tindx(ti))%results%Jirr(1)
          Jirtree(2)  = particles(tindx(ti))%results%Jirr(2)
          Jirtree(3)  = particles(tindx(ti))%results%Jirr(3)

          Jpp(1)       = trslt(ti)%J(1)
          Jpp(2)       = trslt(ti)%J(2)
          Jpp(3)       = trslt(ti)%J(3)

          Jtree(1)     = particles(tindx(ti))%results%J(1)
          Jtree(2)     = particles(tindx(ti))%results%J(2)
          Jtree(3)     = particles(tindx(ti))%results%J(3)

          Bpp(1)       = trslt(ti)%B(1)
          Bpp(2)       = trslt(ti)%B(2)
          Bpp(3)       = trslt(ti)%B(3)

          Btree(1)     = particles(tindx(ti))%results%B(1)
          Btree(2)     = particles(tindx(ti))%results%B(2)
          Btree(3)     = particles(tindx(ti))%results%B(3)


          E_err     = E_err  + dot_product( Epp - Etree  , Epp - Etree )
          E_norm    = E_norm + dot_product( Epp          , Epp )
          
!          Ex_err    = Ex_err  + dot_product( Expp - Extree , Expp - Extree ) !( Expp(1) - Extree(1) )**2!
!          Ex_norm   = Ex_norm + dot_product( Expp , Expp ) !( Expp(1)             )**2!
!          
!          Ey_err    = Ey_err  + dot_product( Eypp - Eytree , Eypp - Eytree ) !( Eypp(1) - Eytree(1) )**2!
!          Ey_norm   = Ey_norm + dot_product( Eypp , Eypp ) !( Eypp(1)             )**2!

          phi_err   =  phi_err    +  ( phipp - phitree )**2
          phi_norm  =  phi_norm   +    phipp**2


          A_err     = A_err  + dot_product( App - Atree        , App - Atree )
          A_norm    = A_norm + dot_product( App          , App )
          
          Ax_err    = Ax_err  + dot_product( Axpp - Axtree     , Axpp - Axtree )
          Ax_norm   = Ax_norm + dot_product( Axpp , Axpp )
          
          Ay_err    = Ay_err  + dot_product( Aypp - Aytree     , Aypp - Aytree )
          Ay_norm   = Ay_norm + dot_product( Aypp , Aypp )
          
          Az_err    = Az_err  + dot_product( Azpp - Aztree     , Azpp - Aztree )
          Az_norm   = Az_norm + dot_product( Azpp , Azpp )
          
!          Axx_err    = Axx_err  + dot_product( Axxpp - Axxtree , Axxpp - Axxtree ) !( Axxpp(2) - Axxtree(2) )**2!
!          Axx_norm   = Axx_norm + dot_product( Axxpp , Axxpp ) !( Axxpp(2)              )**2!
!          
!          Axy_err    = Axy_err  + dot_product( Axypp - Axytree , Axypp - Axytree ) !( Axypp(2) - Axytree(2) )**2!
!          Axy_norm   = Axy_norm + dot_product( Axypp , Axypp ) !( Axypp(2)              )**2!
!          
!          Ayy_err    = Ayy_err  + dot_product( Ayypp - Ayytree , Ayypp - Ayytree ) !( Ayypp(2) - Ayytree(2) )**2!
!          Ayy_norm   = Ayy_norm + dot_product( Ayypp , Ayypp ) !( Ayypp(2)              )**2!
          
          B_err     = B_err  + dot_product( Bpp - Btree , Bpp - Btree )
          B_norm    = B_norm + dot_product( Bpp , Bpp )


          J_err     = J_err  + dot_product( Jpp - Jtree , Jpp - Jtree )
          J_norm    = J_norm + dot_product( Jpp , Jpp )

          Jir_err   = Jir_err  + dot_product( Jirpp - Jirtree , Jirpp - Jirtree )
          Jir_norm  = Jir_norm + dot_product( Jirpp , Jirpp )

          gradpp    = zero   
          gradpp(1) = v(1)*Axpp(1) + v(2)*Axpp(2)  + v(3)*Axpp(3)
          gradpp(2) = v(1)*Aypp(1) + v(2)*Aypp(2)  + v(3)*Aypp(3)
          gradpp(3) = v(1)*Azpp(1) + v(2)*Azpp(2)  + v(3)*Azpp(3)
          gradpp    = gradpp/lorentz_tilde

          gradtree      = zero   
          gradtree(1)   = v(1)*Axtree(1) + v(2)*Axtree(2) + v(3)*Axtree(3)
          gradtree(2)   = v(1)*Aytree(1) + v(2)*Aytree(2) + v(3)*Aytree(3)
          gradtree(2)   = v(1)*Aztree(1) + v(2)*Aztree(2) + v(3)*Aztree(3)
          gradtree      = gradtree/lorentz_tilde

          Fdar_err  = Fdar_err  + (q)**2*( dot_product( Epp - Etree + gradpp - gradtree , Epp - Etree + gradpp - gradtree ) )

          Fdar_norm = Fdar_norm + (q)**2*( dot_product( Epp         + gradpp            , Epp         + gradpp ) )

          FEl_err  = FEl_err    + (q)**2*( dot_product( Epp - Etree                     , Epp - Etree ) )

          FEl_norm = FEl_norm   + (q)**2*( dot_product( Epp                            , Epp          ) )


        end do

        call MPI_ALLREDUCE(MPI_IN_PLACE , phi_err   , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(MPI_IN_PLACE , phi_norm  , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(MPI_IN_PLACE , E_err     , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(MPI_IN_PLACE , E_norm    , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
!        call MPI_ALLREDUCE(MPI_IN_PLACE , Ex_err    , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
!        call MPI_ALLREDUCE(MPI_IN_PLACE , Ex_norm   , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
!        call MPI_ALLREDUCE(MPI_IN_PLACE , Ey_err    , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
!        call MPI_ALLREDUCE(MPI_IN_PLACE , Ey_norm   , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(MPI_IN_PLACE , B_err     , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(MPI_IN_PLACE , B_norm    , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(MPI_IN_PLACE , J_err     , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(MPI_IN_PLACE , J_norm    , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(MPI_IN_PLACE , Jir_err   , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(MPI_IN_PLACE , Jir_norm  , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(MPI_IN_PLACE , A_err     , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(MPI_IN_PLACE , A_norm    , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(MPI_IN_PLACE , Ax_err    , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(MPI_IN_PLACE , Ax_norm   , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(MPI_IN_PLACE , Ay_err    , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(MPI_IN_PLACE , Ay_norm   , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(MPI_IN_PLACE , Az_err    , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(MPI_IN_PLACE , Az_norm   , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
!        call MPI_ALLREDUCE(MPI_IN_PLACE , Axx_err   , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
!        call MPI_ALLREDUCE(MPI_IN_PLACE , Axx_norm  , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
!        call MPI_ALLREDUCE(MPI_IN_PLACE , Axy_err   , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
!        call MPI_ALLREDUCE(MPI_IN_PLACE , Axy_norm  , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
!        call MPI_ALLREDUCE(MPI_IN_PLACE , Ayy_err   , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
!        call MPI_ALLREDUCE(MPI_IN_PLACE , Ayy_norm  , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(MPI_IN_PLACE , Fdar_err  , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(MPI_IN_PLACE , Fdar_norm , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(MPI_IN_PLACE , FEl_err   , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(MPI_IN_PLACE , FEl_norm  , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        
        if ( phi_norm .eq. zero ) phi_norm  = one
        if ( E_norm   .eq. zero ) E_norm    = one
!        if ( Ex_norm  .eq. zero ) Ex_norm   = one
!        if ( Ey_norm  .eq. zero ) Ey_norm   = one
        if ( A_norm   .eq. zero ) A_norm    = one
        if ( Ax_norm  .eq. zero ) Ax_norm   = one
        if ( Ay_norm  .eq. zero ) Ay_norm   = one
        if ( Az_norm  .eq. zero ) Az_norm   = one
!        if ( Axx_norm .eq. zero ) Axx_norm  = one
!        if ( Axy_norm .eq. zero ) Axy_norm  = one
!        if ( Ayy_norm .eq. zero ) Ayy_norm  = one
        if ( B_norm   .eq. zero ) B_norm    = one
        if ( J_norm   .eq. zero ) J_norm    = one
        if ( Jir_norm .eq. zero ) Jir_norm  = one
        if ( Fdar_norm.eq. zero ) Fdar_norm = one
        if ( FEl_norm .eq. zero ) FEl_norm  = one
        
        
        phi_err   = sqrt(phi_err / phi_norm )
        E_err     = sqrt(E_err   / E_norm   )
!        Ex_err    = sqrt(Ex_err  / Ex_norm  )
!        Ey_err    = sqrt(Ey_err  / Ey_norm  )
        A_err     = sqrt(A_err   / A_norm   )
        Ax_err    = sqrt(Ax_err  / Ax_norm  )
        Ay_err    = sqrt(Ay_err  / Ay_norm  )
        Az_err    = sqrt(Az_err  / Az_norm  )
!        Axx_err   = sqrt(Axx_err / Axx_norm )
!        Axy_err   = sqrt(Axy_err / Axy_norm )
!        Ayy_err   = sqrt(Ayy_err / Ayy_norm )
        B_err     = sqrt(B_err   / B_norm   )
        J_err     = sqrt(J_err   / J_norm   )
        Jir_err   = sqrt(Jir_err / Jir_norm )
        Fdar_err  = sqrt(Fdar_err/ Fdar_norm)
        FEl_err   = sqrt(FEl_err / FEl_norm )

        tb = get_time()
        if(root) then
!          write(*,'(a,i12)')    " == [direct test] number tested particles         : ", tn
          write(*,'(a,es12.4)') " == [direct test] Relative error in Pot           : ", phi_err
          write(*,'(a,es12.4)') " == [direct test] Relative error in Eirr          : ", E_err
!          write(*,'(a,es12.4)') " == [direct test] Relative error in dxEirr        : ", Ex_err
!          write(*,'(a,es12.4)') " == [direct test] Relative error in dyEirr        : ", Ey_err
          write(*,'(a,es12.4)') " == [direct test] Relative error in A             : ", A_err
          write(*,'(a,es12.4)') " == [direct test] Relative error in dxA           : ", Ax_err
          write(*,'(a,es12.4)') " == [direct test] Relative error in dyA           : ", Ay_err
          write(*,'(a,es12.4)') " == [direct test] Relative error in dzA           : ", Az_err
!          write(*,'(a,es12.4)') " == [direct test] Relative error in Axx           : ", Axx_err
!          write(*,'(a,es12.4)') " == [direct test] Relative error in Axy           : ", Axy_err
!          write(*,'(a,es12.4)') " == [direct test] Relative error in Ayy           : ", Ayy_err
          write(*,'(a,es12.4)') " == [direct test] Relative error in B             : ", B_err
          write(*,'(a,es12.4)') " == [direct test] Relative error in J             : ", J_err
          write(*,'(a,es12.4)') " == [direct test] Relative error in Jirr          : ", Jir_err
          write(*,'(a,es12.4)') " == [direct test] Relative error in F dar         : ", Fdar_err
          write(*,'(a,es12.4)') " == [direct test] Relative error in F el          : ", FEl_err
          write(*,'(a,es12.4)') " == [direct test] time in test [s]                : ", tb - ta

!          open(unit=rc,file=trim(folder)//trim("monopole3d.dat"),form='formatted',status='unknown',access='append')
!          open(unit=rc,file=trim(folder)//trim("dipole3d.dat"),form='formatted',status='unknown',position='append')
          open(unit=rc,file=trim(folder)//trim("quadrupole3d.dat"),form='formatted',status='unknown',position='append')
          write(rc,*) sqrt(theta2), phi_err, E_err, A_err, Ax_err, Ay_err, Az_err, B_err, J_err, Jir_err, Fdar_err, FEl_err
!          write(rc,*) sqrt(theta2), phi_err, E_err, Ex_err, Ey_err, A_err, Ax_err, Ay_err, Axx_err, Axy_err, Ayy_err,B_err, J_err, Jir_err, Fdar_err, FEl_err
          close(rc)

        end if

        deallocate(tindx)
        deallocate(trnd)
        deallocate(trslt)

      end subroutine test_particles


      subroutine write_field_on_grid_ascii(itime,field_grid)    
            use module_pepc_types
            use module_utils
            use encap
            use module_globals, only: folder
            implicit none
!            integer(kind_pe), intent(in) :: my_rank
            integer(kind_default), intent(in)               :: itime
            type(field_grid_t), intent(in)                  :: field_grid
            character(100)                                  :: filename_i
            integer(kind_particle)                          :: ip,nl
            
            character(12), parameter                        :: part_dir = 'fields/'
            integer, parameter :: filehandle_i = 40
            integer, parameter :: filehandle_e = 41
            integer, parameter :: filehandle_b = 42

                        
            write(filename_i,'(a,"field_",i6.6,"_",i6.6,".dat")') trim(folder)//trim(part_dir), itime, my_rank
            open(filehandle_i, file=trim(filename_i), STATUS='REPLACE')
            
            nl = field_grid%nl
            
            do ip=1, nl
                write(filehandle_i,'(20(f8.3,x) )') field_grid%p(ip)%x(1:2), field_grid%p(ip)%results%E(1:2),    &
                field_grid%p(ip)%results%A(1:3), field_grid%p(ip)%results%B(1:3), field_grid%p(ip)%results%J(1:3),  &
                field_grid%p(ip)%results%dxA(1:3),field_grid%p(ip)%results%dyA(1:3),field_grid%p(ip)%results%pot

            end do
            close(filehandle_i)

        end subroutine
!

      integer function vtk_step_of_step(step) result(vtk_step)
        use module_vtk
!        use pepca_helper
        implicit none

        integer, intent(in) :: step

        if (step .eq. 0) then
          vtk_step = VTK_STEP_FIRST
        else if (step == nt - 1) then
          vtk_step = VTK_STEP_LAST
        else
          vtk_step = VTK_STEP_NORMAL
        endif
      end function vtk_step_of_step

      
      subroutine write_particles_vtk(p, step, realtime)
        use module_vtk_helpers
        use module_pepc_types
!        use pepca_units
        implicit none

        include 'mpif.h'

        type(t_particle), intent(in) :: p(:)
        real*8, intent(in) :: realtime
        integer, intent(in) :: step

        integer :: vtk_step

        vtk_step = vtk_step_of_step(step)
        call vtk_write_particles("particles", MPI_COMM_WORLD, step, realtime, vtk_step, p, vtk_results)

        contains

            subroutine vtk_results(d, r, vtkf)
              use module_vtk
              use module_interaction_specific_types
              implicit none

              type(t_particle_data), intent(in) :: d(:)
              type(t_particle_results), intent(in) :: r(:)
              type(vtkfile_unstructured_grid), intent(inout) :: vtkf

              call vtk_write_particle_data_results(d, r, vtkf)
            end subroutine
            
        end subroutine write_particles_vtk
        
        
        subroutine write_particles_ascii(itime, p)    
            use module_pepc_types
            use module_utils
            use module_globals, only: folder,lorentz_tilde
            implicit none
!            integer(kind_pe), intent(in) :: my_rank
            integer(kind_default), intent(in)               :: itime
            type(t_particle)     , intent(in), dimension(:) :: p
            character(100)                                  :: filename_i,filename_e,filename_b
            integer(kind_particle)                          :: ip
            
            character(12), parameter                        :: part_dir = 'particles/'
            integer, parameter :: filehandle_i = 40
            integer, parameter :: filehandle_e = 41
            integer, parameter :: filehandle_b = 42

                        
            if ( tracks .eq. 0 ) then
                
                write(filename_i,'(a,"particle_",i6.6,".dat")') trim(folder)//trim(part_dir), itime

                open(filehandle_i, file=trim(filename_i), STATUS='REPLACE')

                do ip=1, size(p,kind=kind(ip))
                  
                    write(filehandle_i,'(27(f8.3,x),i12)') p(ip)%x(1:3), p(ip)%data%v(1:3), p(ip)%results%E(1:3),&
                    p(ip)%results%A(1:3), p(ip)%results%B(1:3), p(ip)%results%J(1:3),p(ip)%results%pot, p(ip)%data%g,&
                    p(ip)%results%dxA(1:3),p(ip)%results%dyA(1:3),real(p(ip)%label, kind=kind_particle)
                  
                end do
                close(filehandle_i)

            else if ( tracks .eq. 1 ) then
                
                write(filename_i,'(a,"particle_ions_",i6.6,"_",i6.6,".dat")') trim(folder)//trim(part_dir), itime, my_rank
                write(filename_e,'(a,"particle_elec_",i6.6,"_",i6.6,".dat")') trim(folder)//trim(part_dir), itime, my_rank
                write(filename_b,'(a,"particle_beam_",i6.6,"_",i6.6,".dat")') trim(folder)//trim(part_dir), itime, my_rank



                open(filehandle_i, file=trim(filename_i), STATUS='REPLACE')
                open(filehandle_e, file=trim(filename_e), STATUS='REPLACE')
                open(filehandle_b, file=trim(filename_b), STATUS='REPLACE')

                do ip=1, size(p,kind=kind(ip))
                  if (p(ip)%label .eq. 1)  then

                    write(filehandle_e,'(26(f8.3,x),i12)') p(ip)%x(1:3), p(ip)%data%v(1:3), p(ip)%results%E(1:3),&
                    p(ip)%results%A(1:3), p(ip)%results%B(1:3), p(ip)%results%J(1:3),p(ip)%results%pot,  p(ip)%data%g,&
                    p(ip)%results%dxA(1:3),p(ip)%results%dyA(1:3)

                  else if (p(ip)%label .eq. 2) then

                    write(filehandle_b,'(26(f8.3,x),i12)') p(ip)%x(1:3), p(ip)%data%v(1:3), p(ip)%results%E(1:3),&
                    p(ip)%results%A(1:3), p(ip)%results%B(1:3), p(ip)%results%J(1:3),p(ip)%results%pot,  p(ip)%data%g,&
                    p(ip)%results%dxA(1:3),p(ip)%results%dyA(1:3)

                  else if (p(ip)%label .eq. 3) then

                    write(filehandle_i,'(26(f8.3,x),i12)') p(ip)%x(1:3), p(ip)%data%v(1:3), p(ip)%results%E(1:3),&
                    p(ip)%results%A(1:3), p(ip)%results%B(1:3), p(ip)%results%J(1:3),p(ip)%results%pot, p(ip)%data%g,&
                    p(ip)%results%dxA(1:3),p(ip)%results%dyA(1:3)
                  endif
                end do
                close(filehandle_i)
                close(filehandle_e)
                close(filehandle_b)
                
            endif

        end subroutine write_particles_ascii
        
        
  
!
      subroutine write_domain(p)

        use module_vtk
        use module_vtk_helpers
        use module_pepc, only: global_tree
        implicit none

        type(t_particle), allocatable, intent(in) :: p(:)

        integer :: vtk_step

        ! output of tree diagnostics
        if (step .eq. 0) then
          vtk_step = VTK_STEP_FIRST
        else if (step .eq. nt) then
          vtk_step = VTK_STEP_LAST
        else
          vtk_step = VTK_STEP_NORMAL
        endif
        call vtk_write_branches(step,  dt * step, vtk_step, global_tree)
        call vtk_write_spacecurve(step, dt * step, vtk_step, p)

      end subroutine write_domain

end module
