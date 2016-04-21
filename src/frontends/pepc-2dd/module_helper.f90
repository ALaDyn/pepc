! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
!
! Copyright (C) 2002-2016 Juelich Supercomputing Centre,
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

  ! current measurement
  real*8 :: current_q, current_I
  real*8, parameter :: current_emission = 1.0e4_8
  integer, parameter :: current_file_id = 14

  ! particle data (position, velocity, mass, charge)
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
    use module_interaction_specific, only : theta2, eps2, force_law, include_far_field_if_periodic
    use module_mirror_boxes, only: periodicity, spatial_interaction_cutoff
    use module_shortcut, only: half,pi
!    use module_globals, only:adv
    !use field_helper, only:
    implicit none

    integer, parameter :: fid = 12
    character(255)     :: para_file
    logical            :: read_para_file,bool_exp, bool_nk,bool_pic

    namelist /pepc2dd/ tnp, nt, dt, particle_output, domain_output,  particle_test, &
       diag_interval, periodicity_particles,  spatial_interaction_cutoff,initial_setup,ischeme,&
       restart_file,restart_step,normal,nsp,veth,vedrift,vith,vidrift,nppd,woutput,ivdim,&
       folder,flag_classic,we,vmax

    ! set default parameter values
    initial_setup   = 2
    tnp             = 1441
    nt              = 20
    dt              = 1e-3
    particle_test   = .true.
    particle_output = .true.
    ischeme         = "leapfrog"
    folder          = "data"
    periodicity_particles     = .false.
    restart_file    = ""
    restart_step    = 1

    ivdim           = 3
    diag_interval   = 2
    nsp             = 1
    normal          = 0
    veth            = one
    vedrift         = zero
    vith            = zero
    vidrift         = zero
    nppd            = (/ 1441, 1 , 1 /)

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
    
    flag_classic    = .true.


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

    if(root) then
      write(*,'(a,i12)')    " == total number of particles : ", tnp
      write(*,'(a,i12)')    " == number of time steps      : ", nt
      write(*,'(a,es12.4)') " == time step                 : ", dt
      write(*,*)            "== time integration scheme   : ",  ischeme
      write(*,'(a,es12.4)') " == theta2                    : ", theta2
      write(*,'(a,es12.4)') " == epsilon2                  : ", eps2
      write(*,'(a,l12)')    " == periodicity               : ", periodicity!
      write(*,'(a,i12)')    " == diag & IO interval        : ", diag_interval
      write(*,'(a,l12)')    " == particle test             : ", particle_test
      write(*,'(a,l12)')    " == particle output           : ", particle_output
      write(*,'(a,l12)')    " == domain output             : ", domain_output
      write(*,*) ""

    end if



    ! reste current measurement
    current_q = 0.0_8
    current_I = 0.0_8

    bool_exp = ischeme .eq. "leapfrog"
    bool_nk  = ischeme .eq. "midpoint3D" .or. ischeme .eq. "trapezoidal3D" 
    bool_pic = ischeme .eq. "midpoint_picard" .or. ischeme .eq. "trapezoidal_picard" .or. &
    ischeme .eq."midpoint_picard_electrostatic" .or. ischeme .eq."hamiltonian_boris"!.or. ischeme .eq. "trapezoidal_broyden"
    bool_exp = ischeme .eq. "leapfrog"
    bool_nk  = ischeme .eq. "midpoint3D" .or. ischeme .eq. "trapezoidal3D" 
    bool_pic = ischeme .eq. "midpoint_picard" .or. ischeme .eq. "trapezoidal_picard" .or. &
    ischeme .eq."midpoint_picard_electrostatic" .or. ischeme .eq."hamiltonian_boris".or. ischeme .eq. "trapezoidal_broyden"

    adv = -1

    if (bool_exp) then
        adv = 0
    elseif (bool_pic) then
        adv = 1
    elseif (bool_nk) then
        adv = 2
    endif


  end subroutine

  subroutine write_restart_2d(p,istep)
    implicit none
    type(t_particle), allocatable, intent(in)    :: p(:)
    integer(kind_particle)       , intent(in)    :: istep
    integer(kind_particle)                       :: ip,rc,size_tmp
    real(kind_particle), allocatable             :: tmp(:)
    character(255)                               :: str_istep,str_proc

!    if (root)  then

        size_tmp = 6
        if (allocated(tmp)) deallocate(tmp)
        allocate(tmp(size_tmp), stat = rc)
!        open(unit=rc,file="restart.dat",form='formatted',status='unknown',access='append')
        write( str_istep, '(i10)' ) istep
        write( str_proc , '(i10)' ) my_rank
        open (unit=my_rank,file=trim(folder)//trim("restart/restart_")//trim(adjustl(str_proc))//".dat",action="write",status="replace")
!        open (unit=my_rank,file=trim("restart/restart_")//trim(adjustl(str_istep))//trim("_")//trim(adjustl(str_proc))//".dat",action="write",status="replace")
!        write (rc,*) " "

        do ip = 1, np

              tmp(1)            =   p(ip)%label
              tmp(2)            =   p(ip)%data%q

              tmp(3)            =   p(ip)%x(1)
              tmp(4)            =   p(ip)%x(2)

              tmp(5)            =   p(ip)%data%v(1)
              tmp(6)            =   p(ip)%data%v(2)

              write (my_rank,*) tmp


        end do

        close(unit = my_rank)
        deallocate(tmp)
!    endif



  end subroutine write_restart_2d



subroutine init_particles(p,field_grid)
    use module_pepc
    use module_pepc_kinds
    use module_init
    use module_mirror_boxes, only:t_lattice_1,t_lattice_2,t_lattice_3,LatticeOrigin
    use module_tool        , only:par_rand
    use module_shortcut    , only:half,pi,oneoverfourpi,oneoverpi
    use module_utils       , only: create_directory
    use encap
    
    implicit none
    include 'mpif.h'
    type(t_particle), allocatable, intent(inout) :: p(:)
    type(field_grid_t)           , intent(in)    :: field_grid
!    integer(kind_particle)       , intent(in)    :: initial_setup
    integer(kind_particle)                       :: ip,rc
    real(kind_particle)                          :: dummy,nd,lambda,rtnp,gl

    if(my_rank.eq.0) write(*,'(a)') " == [init] init particles "

    ! set initially number of local particles
    np = tnp / n_ranks
!    write(*,*) "Info== ", my_rank, n_ranks,np, tnp
    if(my_rank.eq.(n_ranks-1)) np = np + MOD(tnp, int(n_ranks, kind=kind_particle))

    allocate(particles(np), stat=rc)
    if(rc.ne.0) write(*,*) " === particle allocation error!"

    allocate(direct_L2(np), stat=rc)
    if(rc.ne.0) write(*,*) " === direct_L2 allocation error!"
    direct_L2 = -1.0_8



    ixdim                = int( 3, kind = kind_dim )
    if (n(3) == 1) ixdim = int( 2, kind = kind_dim )
    if (n(2) == 1) ixdim = int( 1, kind = kind_dim )

    if (nsp .gt. 2 .or. nsp .lt. 1) then
        write(*,*) "Number of species larger than 2 is not yet implemented !!"
        call exit(1)
    endif

    call pepc_prepare(ixdim)
    call create_directory(trim(folder))

    ! set random seed
    dummy = par_rand(my_rank)

    select case (initial_setup)
          case (1)  ! Load State
              call read_restart_2d(p,"restart_0_")
          case (2)  !  Default initial setup - 2D random particles with thermal velocity
              call thermal(p)
          case (3)  ! initial setup - 2D random particles with 3D thermal velocity
              call thermal2D3V(p)
          case (4)  ! 2D Landau Damping
              call landau_damping(p,field_grid)
          case (5)  ! Langmuir waves
              call langmuir_waves(p,field_grid)
          case (6)  ! 1D beam
              call beam(p)
          case (7)  ! 1D beam
              call beam_disk(p)
          case (8)  ! 1D beam
              call weibell_instability(p,field_grid)
          case (9)  ! 1D beam
              call solenoid(p)
          case (10)  ! 1D beam
              call solenoid_infinite(p) 
          case (12)  ! 1D beam 
              call neutral_plasma(p,field_grid)
          case (13)  ! 1D beam 
              call periodic_test(p,field_grid)
          case default
              call thermal(p)
    select case (initial_setup)
          case (1)  ! Load State
              call read_restart_2d(p,"restart_0_")
          case (2)  !  Default initial setup - 2D random particles with thermal velocity
              call thermal(p)
          case (3)  ! initial setup - 2D random particles with 3D thermal velocity
              call thermal2D3V(p)
          case (4)  ! 2D Landau Damping
              call landau_damping(p,field_grid)
          case (5)  ! Langmuir waves
              call langmuir_waves(p,field_grid)
          case (6)  ! 1D beam
              call beam(p)
          case (7)  ! 1D beam
              call beam_disk(p)
          case (8)  ! 1D beam
              call weibell_instability(p,field_grid)
          case (9)  ! 1D beam
              call solenoid(p)
          case (10)  ! 1D beam
              call solenoid_infinite(p) 
          case (12)  ! 1D beam
              call neutral_plasma(p,field_grid)
          case default
              call thermal(p)

    end select

    
    rtnp            = real(tnp, kind=kind_particle)         ! total nomber of particles - doubleprecision
    nd              = rtnp/Volume                           ! number density
    

    select case (normal)
          case (0)  !
            
            
            qtilde          =  nsp*we/nd
            mtilde          =  nsp*we/nd
            
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
            
            qtilde          =  nsp*we/nd
            mtilde          =  nsp*we/nd
            
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

            qtilde          =  nsp*we/nd
            mtilde          =  nsp*we/nd
            
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
            if ( nsp .eq. 2 ) mtilde          =  mi
            
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
            if ( nsp .eq. 2 ) mtilde          =  mi

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
            if ( nsp .eq. 2 ) mtilde          =  mi

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

    
    do ip = 1,np
        
        p(ip)%data%q      = p(ip)%data%q*qtilde
        p(ip)%data%m      = p(ip)%data%m*mtilde
        p(ip)%data%v(1:3) = p(ip)%data%v(1:3)*vtilde
        
        gl           = dot_product( p(ip)%data%v/vtilde, p(ip)%data%v/vtilde )
        gl           = one/sqrt( one - gl )
        p(ip)%data%g = gl
        p(ip)%data%v = p(ip)%data%v*gl
        
    enddo
    
    t_lattice_1(1)      = extent(1)
    t_lattice_2(2)      = extent(2)
    LatticeOrigin(1:3)  = offset(1:3)

    if ( extent(3).gt. zero )   t_lattice_3(3)      = extent(3)

    next_label = tnp+1

  end subroutine init_particles

  subroutine time_step(p)
  use module_shortcut, only: zero,one
  implicit none
  include 'mpif.h'
    type(t_particle), allocatable, intent(in)    :: p(:)
    integer(kind_particle)                       :: ip,jp,np,rc
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
  
  subroutine iperiodic_particles(ip,np,p)
!    use module_mirror_boxes
    implicit none

    type(t_particle), allocatable, intent(inout) :: p(:)
    integer(kind_particle)       , intent(in)    :: ip,np
    integer(kind_particle)                       :: jp

        do jp = 1,3
            if ( p(ip)%x(jp) .lt. offset(jp) )              p(ip)%x(jp) = p(ip)%x(jp) + extent(jp)
            if ( p(ip)%x(jp) .gt. offset(jp) + extent(jp) ) p(ip)%x(jp) = mod( p(ip)%x(jp) , extent(jp) )
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
    integer(kind_particle)                       :: ip

    do ip = 1,np

        p(ip)%results%E(1:3)    = p(ip)%results%E(1:3)      *etilde
        p(ip)%results%B(1:3)    = p(ip)%results%B(1:3)      *btilde/vtilde
        p(ip)%results%A(1:3)    = p(ip)%results%A(1:3)      *atilde/vtilde
        p(ip)%results%dxA(1:3)  = p(ip)%results%dxA(1:3)    *atilde/vtilde
        p(ip)%results%dyA(1:3)  = p(ip)%results%dyA(1:3)    *atilde/vtilde
        p(ip)%results%J(1:3)    = p(ip)%results%J(1:3)      *jtilde
        p(ip)%results%Jirr(1:3) = p(ip)%results%Jirr(1:3)   *jtilde
        p(ip)%results%pot       = p(ip)%results%pot         *phitilde
!        p(ip)%results%rho       = p(ip)%results%rho         *rhotilde
        
!        err = err + ( p(ip)%data%v(1)/p(ip)%data%g*p(ip)%results%pot-p(ip)%results%A(1) )**2
        p(ip)%results%E(1:3)    = p(ip)%results%E(1:3)      *etilde
        p(ip)%results%B(1:3)    = p(ip)%results%B(1:3)      *btilde
        p(ip)%results%A(1:3)    = p(ip)%results%A(1:3)      *atilde
        p(ip)%results%dxA(1:3)  = p(ip)%results%dxA(1:3)    *atilde
        p(ip)%results%dyA(1:3)  = p(ip)%results%dyA(1:3)    *atilde
        p(ip)%results%J(1:3)    = p(ip)%results%J(1:3)      *jtilde
        p(ip)%results%Jirr(1:3) = p(ip)%results%Jirr(1:3)   *jtilde
        p(ip)%results%pot       = p(ip)%results%pot         *phitilde
!        p(ip)%results%rho       = p(ip)%results%rho         *rhotilde

    enddo
!    write(*,*) "errore2: ",err
  end subroutine normalize
  
  subroutine inormalize(ip,np,p)
    implicit none
    enddo

    type(t_particle), allocatable, intent(inout) :: p(:)
    integer(kind_particle)       , intent(in)    :: np,ip
  end subroutine normalize
  
  subroutine inormalize(ip,np,p)
    implicit none

        p(ip)%results%E(1:3)    = p(ip)%results%E(1:3)      *etilde
        p(ip)%results%B(1:3)    = p(ip)%results%B(1:3)      *btilde/vtilde
        p(ip)%results%A(1:3)    = p(ip)%results%A(1:3)      *atilde/vtilde
        p(ip)%results%dxA(1:3)  = p(ip)%results%dxA(1:3)    *atilde/vtilde
        p(ip)%results%dyA(1:3)  = p(ip)%results%dyA(1:3)    *atilde/vtilde
        p(ip)%results%J(1:3)    = p(ip)%results%J(1:3)      *jtilde
        p(ip)%results%Jirr(1:3) = p(ip)%results%Jirr(1:3)   *jtilde
        p(ip)%results%pot       = p(ip)%results%pot         *phitilde
!        p(ip)%results%rho       = p(ip)%results%rho         *rhotilde
    type(t_particle), allocatable, intent(inout) :: p(:)
    integer(kind_particle)       , intent(in)    :: np,ip

  end subroutine inormalize
  
  
        p(ip)%results%E(1:3)    = p(ip)%results%E(1:3)      *etilde
        p(ip)%results%B(1:3)    = p(ip)%results%B(1:3)      *btilde
        p(ip)%results%A(1:3)    = p(ip)%results%A(1:3)      *atilde
        p(ip)%results%dxA(1:3)  = p(ip)%results%dxA(1:3)    *atilde
        p(ip)%results%dyA(1:3)  = p(ip)%results%dyA(1:3)    *atilde
        p(ip)%results%J(1:3)    = p(ip)%results%J(1:3)      *jtilde
        p(ip)%results%Jirr(1:3) = p(ip)%results%Jirr(1:3)   *jtilde
        p(ip)%results%pot       = p(ip)%results%pot         *phitilde
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
        real(kind_particle)                   :: E(1:3),ETilde(1:3),E_norm_loc,E_norm_global,E_local,E_global
!        real(kind_particle)                   :: Et(1:3),EtTilde(1:3),Et_norm_loc,Et_norm_global,Et_local,Et_global
        real(kind_particle)                   :: B(1:3),BTilde(1:3),B_norm_loc,B_norm_global,B_local,B_global
        real(kind_particle)                   :: J(1:3),JTilde(1:3),J_norm_loc,J_norm_global,J_local,J_global
        real(kind_particle)                   :: Jirr(1:3),JirrTilde(1:3),Jirr_norm_loc,Jirr_norm_global,Jirr_local,Jirr_global
        real(kind_particle)                   :: A(1:3),ATilde(1:3),A_norm_loc,A_norm_global,A_local,A_global
        real(kind_particle)                   :: phi,phiTilde,phi_norm_loc,phi_local,phi_global,phi_norm_global,v(3),m,q,rot(3),rot_tilde(3),&
                                                 F_dar_loc,F_dar_den_loc,F_el_loc,F_el_den_loc,devE,devPhi,ta,tb,&
                                                 F_dar_glo,F_dar_den_glo,F_el_glo,F_el_den_glo,dxA(1:3),dyA(1:3),&
                                                 dxA_norm_loc,dxA_norm_global,dxA_local,dxA_global,              &
                                                 dyA_norm_loc,dyA_norm_global,dyA_local,dyA_global,dxATilde(1:3),dyATilde(1:3)

        ta = get_time()

        E_local       = zero
        E_global      = zero
        E_norm_loc    = zero
        E_norm_global = zero

!        Et_norm_loc    = zero
!        Et_norm_global = zero
!        Et_local       = zero
!        Et_global      = zero

        phi_norm_loc    = zero
        phi_norm_global = zero
        phi_local       = zero
        phi_global      = zero

        B_norm_loc    = zero
        B_norm_global = zero
        B_local       = zero
        B_global      = zero

        A_norm_loc    = zero
        A_norm_global = zero
        A_local       = zero
        A_global      = zero

        J_norm_loc    = zero
        J_norm_global = zero
        J_local       = zero
        J_global      = zero

        Jirr_norm_loc    = zero
        Jirr_norm_global = zero
        Jirr_local       = zero
        Jirr_global      = zero

        F_dar_loc     = zero
        F_el_loc      = zero
        F_dar_den_loc = zero
        F_el_den_loc  = zero

        F_dar_glo     = zero
        F_el_glo      = zero
        F_dar_den_glo = zero
        F_el_den_glo  = zero

        devE          = zero
        devPhi        = zero

        tn = tnp!particle_direct / n_ranks
        if(my_rank.eq.(n_ranks-1)) tn = tn + MOD(particle_direct, n_ranks)

        allocate(tindx(tn), trnd(tn), trslt(tn))

        call random(trnd)

        tindx(1:tn) = int(trnd(1:tn) * (np-1)) + 1

        call directforce(particles, tindx, tn, trslt, MPI_COMM_WORLD)

        do ti = 1, tn

          v(1)          = particles(tindx(ti))%data%v(1)
          v(2)          = particles(tindx(ti))%data%v(2)
          v(3)          = particles(tindx(ti))%data%v(3)

          m           = particles(tindx(ti))%data%m
          q           = particles(tindx(ti))%data%q

          phi         = trslt(ti)%pot
          phiTilde    = particles(tindx(ti))%results%pot

          E(1)          = trslt(ti)%e(1)
          E(2)          = trslt(ti)%e(2)
          E(3)          = trslt(ti)%e(3)

          ETilde(1)    = particles(tindx(ti))%results%e(1)
          ETilde(2)    = particles(tindx(ti))%results%e(2)
          ETilde(3)    = particles(tindx(ti))%results%e(3)

!          Et(1)          = trslt(ti)%Et(1)
!          Et(2)          = trslt(ti)%Et(2)
!          Et(3)          = trslt(ti)%Et(3)
!
!          EtTilde(1)    = particles(tindx(ti))%results%Et(1)
!          EtTilde(2)    = particles(tindx(ti))%results%Et(2)
!          EtTilde(3)    = particles(tindx(ti))%results%Et(3)

          A(1)          = trslt(ti)%A(1)
          A(2)          = trslt(ti)%A(2)
          A(3)          = trslt(ti)%A(3)

          ATilde(1)     = particles(tindx(ti))%results%A(1)
          ATilde(2)     = particles(tindx(ti))%results%A(2)
          ATilde(3)     = particles(tindx(ti))%results%A(3)
          
          dxA(1)        = trslt(ti)%dxA(1)
          dxA(2)        = trslt(ti)%dxA(2)
          dxA(3)        = trslt(ti)%dxA(3)

          dxATilde(1)   = particles(tindx(ti))%results%dxA(1)
          dxATilde(2)   = particles(tindx(ti))%results%dxA(2)
          dxATilde(3)   = particles(tindx(ti))%results%dxA(3)
          
          
          dyA(1)        = trslt(ti)%dyA(1)
          dyA(2)        = trslt(ti)%dyA(2)
          dyA(3)        = trslt(ti)%dyA(3)

          dyATilde(1)   = particles(tindx(ti))%results%dyA(1)
          dyATilde(2)   = particles(tindx(ti))%results%dyA(2)
          dyATilde(3)   = particles(tindx(ti))%results%dyA(3)


          JirrTilde(1)  = particles(tindx(ti))%results%Jirr(1)
          JirrTilde(2)  = particles(tindx(ti))%results%Jirr(2)
          JirrTilde(3)  = particles(tindx(ti))%results%Jirr(3)

          Jirr(1)       = trslt(ti)%Jirr(1)
          Jirr(2)       = trslt(ti)%Jirr(2)
          Jirr(3)       = trslt(ti)%Jirr(3)

          JTilde(1)     = particles(tindx(ti))%results%J(1)
          JTilde(2)     = particles(tindx(ti))%results%J(2)
          JTilde(3)     = particles(tindx(ti))%results%J(3)

          J(1)          = trslt(ti)%J(1)
          J(2)          = trslt(ti)%J(2)
          J(3)          = trslt(ti)%J(3)

          B(1)          = trslt(ti)%B(1)
          B(2)          = trslt(ti)%B(2)
          B(3)          = trslt(ti)%B(3)

          BTilde(1)     = particles(tindx(ti))%results%B(1)
          BTilde(2)     = particles(tindx(ti))%results%B(2)
          BTilde(3)     = particles(tindx(ti))%results%B(3)


          E_local     = E_local + dot_product( E - Etilde , E - Etilde )
          E_norm_loc  = E_norm_loc + dot_product( E , E )

!          Et_local     = Et_local + dot_product( Et - EtTilde , Et - EtTilde )
!          Et_norm_loc  = Et_norm_loc + dot_product( Et , Et )

          phi_local        = phi_local +  ( phi - phiTilde )**2
          phi_norm_loc     = phi_norm_loc +  phi**2


          A_local     = A_local + dot_product( A - Atilde , A - Atilde )
          A_norm_loc  = A_norm_loc + dot_product( A , A )
          
          dxA_local   = dxA_local    + dot_product( dxA - dxAtilde , dxA - dxAtilde )
          dxA_norm_loc= dxA_norm_loc + dot_product( dxA , dxA )
          
          dyA_local   = dyA_local    + dot_product( dyA - dyAtilde , dyA - dyAtilde )
          dyA_norm_loc= dyA_norm_loc + dot_product( dyA , dyA )

          B_local     = B_local + dot_product( B - BTilde , B - BTilde )
          B_norm_loc  = B_norm_loc + dot_product( B , B )


          J_local     = J_local + dot_product( J - JTilde , J - JTilde )
          J_norm_loc  = J_norm_loc + dot_product( J , J )

          Jirr_local     = Jirr_local + dot_product( Jirr - JirrTilde , Jirr - JirrTilde )
          Jirr_norm_loc  = Jirr_norm_loc + dot_product( Jirr , Jirr )


          rot          = cross_product(v,B)
          rot_tilde    = cross_product(v,Btilde)

          F_dar_loc       = F_dar_loc + (q)**2*( dot_product( E  -  Etilde , E - Etilde )  &
                        +   dot_product( rot - rot_tilde, rot - rot_tilde ) + two*dot_product( E - ETilde , rot - rot_tilde ) )

          F_dar_den_loc  = F_dar_den_loc + (q)**2*( dot_product( E , E)  + two*dot_product( E, rot ) + dot_product( rot, rot ) )

          F_el_loc       = F_el_loc + (q)**2*dot_product( E - Etilde , E - Etilde )

          F_el_den_loc  = F_el_den_loc + (q)**2*dot_product( E , E )


        end do

        call MPI_ALLREDUCE(tn           , tn_global       , 1, MPI_KIND_PARTICLE, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(phi_local    , phi_global      , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(phi_norm_loc , phi_norm_global , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(E_local      , E_global        , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(E_norm_loc   , E_norm_global   , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
!        call MPI_ALLREDUCE(Et_local     , Et_global       , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
!        call MPI_ALLREDUCE(Et_norm_loc  , Et_norm_global  , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(B_norm_loc   , B_norm_global   , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(B_local      , B_global        , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(A_local      , A_global        , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(A_norm_loc   , A_norm_global   , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(dxA_local    , dxA_global      , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(dxA_norm_loc , dxA_norm_global , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(dyA_local    , dyA_global      , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(dyA_norm_loc , dyA_norm_global , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(J_local      , J_global        , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(J_norm_loc   , J_norm_global   , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(Jirr_local   , Jirr_global     , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(Jirr_norm_loc, Jirr_norm_global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(F_dar_loc    , F_dar_glo       , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(F_dar_den_loc, F_dar_den_glo   , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(F_el_loc     , F_el_glo        , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(F_el_den_loc , F_el_den_glo    , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)


        devE                 = sqrt(E_global)/(tn_global-one)
        devPhi               = sqrt(phi_global)/(tn_global-one)

        phi_global           = sqrt(phi_global / phi_norm_global)
        E_global             = sqrt(E_global   / E_norm_global)
!        Et_global            = sqrt(Et_global  / Et_norm_global)
        A_global             = sqrt(A_global   / A_norm_global)
        dxA_global           = sqrt(dxA_global / dxA_norm_global)
        dyA_global           = sqrt(dyA_global / dyA_norm_global)
        B_global             = sqrt(B_global   / B_norm_global)
        J_global             = sqrt(J_global   / J_norm_global)
        Jirr_global          = sqrt(Jirr_global/ Jirr_norm_global)
        F_dar_glo            = sqrt(F_dar_glo  / F_dar_den_glo )
        F_el_glo             = sqrt(F_el_glo   / F_el_den_glo )

        tb = get_time()
        if(root) then
          write(*,'(a,i12)')    " == [direct test] number tested particles         : ", tn
    !      write(*,'(a,es12.4)') " == [direct test] l2 A                            : ", A_norm_global
    !      write(*,'(a,es12.4)') " == [direct test] l2 A error                      : ", A_global
    !      write(*,'(a,es12.4)') " == [direct test] l2 B                            : ", B_norm_global
    !      write(*,'(a,es12.4)') " == [direct test] l2 B error                      : ", B_global
!          write(*,'(a,es12.4)') " == [direct test] MAC                             : ", mac
          write(*,'(a,es12.4)') " == [direct test] Relative error in Pot           : ", phi_global
          write(*,'(a,es12.4)') " == [direct test] Relative error in El            : ", E_global
!          write(*,'(a,es12.4)') " == [direct test] Relative error in Et            : ", Et_global
          write(*,'(a,es12.4)') " == [direct test] Relative error in A             : ", A_global
          write(*,'(a,es12.4)') " == [direct test] Relative error in dxA           : ", dxA_global
          write(*,'(a,es12.4)') " == [direct test] Relative error in dyA           : ", dyA_global
          write(*,'(a,es12.4)') " == [direct test] Relative error in B             : ", B_global
          write(*,'(a,es12.4)') " == [direct test] Relative error in J             : ", J_global
          write(*,'(a,es12.4)') " == [direct test] Relative error in Jirr          : ", Jirr_global
          write(*,'(a,es12.4)') " == [direct test] Relative error in F dar         : ", F_dar_glo
          write(*,'(a,es12.4)') " == [direct test] Relative error in F el          : ", F_el_glo
          write(*,'(a,es12.4)') " == [direct test] L2 error in E                   : ", devE
          write(*,'(a,es12.4)') " == [direct test] L2 error in El Pot              : ", devPhi
          write(*,'(a,es12.4)') " == [direct test] time in test [s]                : ", tb - ta

!          open(unit=rc,file=trim(folder)//trim("monopole2d.dat"),form='formatted',status='unknown',access='append')
!          open(unit=rc,file=trim(folder)//trim("dipole2d.dat"),form='formatted',status='unknown',position='append')
!          open(unit=rc,file=trim(folder)//trim("quadrupole2d.dat"),form='formatted',status='unknown',position='append')
          write(rc,*) sqrt(theta2), phi_global, E_global, A_global, dxA_global, dyA_global,B_global, J_global, Jirr_global, F_dar_glo, F_el_glo
          close(rc)

        end if

        deallocate(tindx)
        deallocate(trnd)
        deallocate(trslt)

      end subroutine test_particles


!      subroutine write_particle_ascii(p,itime)
!        use module_globals, only: woutput
!        implicit none
!        include 'mpif.h'
!
!        type(t_particle), allocatable, intent(in) :: p(:)
!        integer(kind_particle)       , intent(in) :: itime
!        integer(kind_particle)                    :: ip
!        integer(kind_particle),      parameter    :: part   = 100
!        integer(kind_particle),      parameter    :: partp1 = 101
!        integer(kind_particle),      parameter    :: sizeout = 22
!        real(kind_particle)                       :: tmp(sizeout),EE,EE_glo
!
!
!        character(255) str,file0,file1
!
!        write( str, '(i10)' ) itime
!        file0               = "particle/particle_ascii_"
!        file1               = "particle/electric_energy_density"
!!        file0               = file0//str
!!        write(*,*) file0
!!
!
!
!        if (root) then
!            open (unit=part,file=trim(folder)//trim(file0)//trim(adjustl(str))//".dat",action="write",status="replace")
!!            open(unit=partp1,file=trim(file1)//".dat",form='formatted',status='unknown',access='append')
!            write(*,*) "==== Writing particle_ascii - tnp:         ",tnp
!        endif
!
!
!        EE_glo   = 0.0_8
!        EE       = 0.0_8
!
!
!        select case (woutput)
!              case (1)  ! 1D
!
!                do ip = 1,np
!                    tmp(1) = p(ip)%x(1)
!                    tmp(2) = p(ip)%results%E(1)
!                    tmp(3) = p(ip)%results%pot
!                    tmp(4) = p(ip)%results%B(3)
!                    tmp(5) = p(ip)%data%v(1)
!        !            EE     = EE + .5_8*p(ip)%results%E(1)*p(ip)%results%E(1)
!                    write (part,*) tmp(1:5)
!
!                enddo
!
!
!            case (2)  ! 2D
!
!                do ip = 1,np
!                    tmp(1) = p(ip)%x(1)
!                    tmp(2) = p(ip)%x(2)
!                    tmp(3) = p(ip)%results%E(1)
!                    tmp(4) = p(ip)%results%E(2)
!                    tmp(5) = p(ip)%results%pot
!                    tmp(6) = p(ip)%results%A(1)
!                    tmp(7) = p(ip)%results%A(2)
!                    tmp(8)= p(ip)%results%B(3)
!                    tmp(9)= p(ip)%results%J(1)
!                    tmp(10)= p(ip)%results%J(2)
!                    tmp(11)= p(ip)%results%Jirr(1)
!                    tmp(12)= p(ip)%results%Jirr(2)
!                    tmp(13) = p(ip)%data%v(1)
!                    tmp(14) = p(ip)%data%v(2)
!        !            EE     = EE + .5_8*p(ip)%results%E(1)*p(ip)%results%E(1)
!                    write (part,*) tmp(1:sizeout)
!
!                enddo
!
!            case (3)  ! 2D-3V
!
!                do ip = 1,np
!                    tmp(1) = p(ip)%x(1)
!                    tmp(2) = p(ip)%x(2)
!                    tmp(3) = p(ip)%results%E(1)
!                    tmp(4) = p(ip)%results%E(2)
!                    tmp(5) = p(ip)%results%pot
!                    tmp(6) = p(ip)%results%A(1)
!                    tmp(7)= p(ip)%results%A(2)
!                    tmp(8)= p(ip)%results%A(3)
!                    tmp(9)= p(ip)%results%B(1)
!                    tmp(10)= p(ip)%results%B(2)
!                    tmp(11)= p(ip)%results%B(3)
!                    tmp(12) = p(ip)%results%J(1)
!                    tmp(13) = p(ip)%results%J(2)
!                    tmp(14) = p(ip)%results%J(3)
!                    tmp(15) = p(ip)%results%Jirr(1)
!                    tmp(16) = p(ip)%results%Jirr(2)
!                    tmp(17) = p(ip)%data%v(1)
!                    tmp(18) = p(ip)%data%v(2)
!                    tmp(19) = p(ip)%data%v(3)
!        !            EE     = EE + .5_8*p(ip)%results%E(1)*p(ip)%results%E(1)
!                    write (part,*) tmp(1:19)
!!                    write(*,*) ip,p(ip)%results%E(1),p(ip)%results%B(2)
!
!                enddo
!
!            case (4)  ! 3D
!
!                do ip = 1,np
!                    tmp(1) = p(ip)%x(1)
!                    tmp(2) = p(ip)%x(2)
!                    tmp(3) = p(ip)%results%E(1)
!                    tmp(4) = p(ip)%results%E(2)
!                    tmp(5) = p(ip)%results%pot
!                    tmp(6) = p(ip)%results%A(1)
!                    tmp(7) = p(ip)%results%A(2)
!                    tmp(8) = p(ip)%results%A(3)
!                    tmp(9) = p(ip)%results%B(1)
!                    tmp(10)= p(ip)%results%B(2)
!                    tmp(11)= p(ip)%results%B(3)
!                    tmp(12) = p(ip)%results%J(1)
!                    tmp(13) = p(ip)%results%J(2)
!                    tmp(14) = p(ip)%results%J(3)
!                    tmp(15) = p(ip)%results%Jirr(1)
!                    tmp(16) = p(ip)%results%Jirr(2)
!                    tmp(17) = p(ip)%data%v(1)
!                    tmp(18) = p(ip)%data%v(2)
!                    tmp(19) = p(ip)%data%v(3)
!                    write (part,*) tmp(1:sizeout)
!
!                enddo
!
!              case default
!
!                do ip = 1,np
!                    tmp(1) = p(ip)%x(1)
!                    tmp(2) = p(ip)%data%v(1)
!                    tmp(3) = p(ip)%results%E(1)
!                    tmp(4) = p(ip)%results%pot
!                    tmp(5) = p(ip)%results%B(3)
!        !            EE     = EE + .5_8*p(ip)%results%E(1)*p(ip)%results%E(1)
!                    write (part,*) tmp(1:5)
!
!                enddo
!
!        end select
!
!
!!        call MPI_ALLREDUCE(EE, EE_glo, 1, MPI_KIND_PARTICLE, MPI_SUM, MPI_COMM_WORLD, rc)
!
!
!
!        if (root) then
!!            write (partp1,*) dt*itime,  EE_glo
!            close ( unit=part )
!!            close ( unit=partp1 )
!            write(*,*) "==== Writing particle_ascii finished- tnp: ",tnp
!
!        endif
!
!      end subroutine write_particle_ascii


      subroutine write_field_on_grid_ascii(field_grid,itime)    
            use module_pepc_types
            use module_utils
            use encap
            use module_globals, only: folder
            implicit none
!            integer(kind_pe), intent(in) :: my_rank
            integer(kind_default), intent(in)               :: itime
            type(field_grid_t), intent(in)                  :: field_grid
            logical                                         :: firstcall  = .true.
            character(50)                                   :: dir
            character(100)                                  :: filename_i,filename_e,filename_b
            integer(kind_particle)                          :: ip,nl
            
            character(12), parameter                        :: part_dir = 'fields/'
            integer, parameter :: filehandle_i = 40
            integer, parameter :: filehandle_e = 41
            integer, parameter :: filehandle_b = 42

            
            if (firstcall) then
              call create_directory(trim(folder))
              call create_directory(trim(folder)//trim(part_dir))
              firstcall = .false.
            endif
            
!            write(filename_i,'(a,"field_ions_",i6.6,"_",i6.6,".dat")') trim(folder)//trim(part_dir), itime, my_rank
!            write(filename_e,'(a,"field_elec_",i6.6,"_",i6.6,".dat")') trim(folder)//trim(part_dir), itime, my_rank
!            write(filename_b,'(a,"field_beam_",i6.6,"_",i6.6,".dat")') trim(folder)//trim(part_dir), itime, my_rank
            write(filename_i,'(a,"field_",i6.6,"_",i6.6,".dat")') trim(folder)//trim(part_dir), itime, my_rank

            

            open(filehandle_i, file=trim(filename_i), STATUS='REPLACE')
!            open(filehandle_e, file=trim(filename_e), STATUS='REPLACE')
!            open(filehandle_b, file=trim(filename_b), STATUS='REPLACE')
            
            nl = field_grid%nl
            
            do ip=1, nl
!              if (field_grid%p(ip)%label .eq. 1)  then
                write(filehandle_i,'(23(f8.3,x) )') field_grid%p(ip)%x(1:2),field_grid%p(ip)%data%v(1:3), field_grid%p(ip)%results%E(1:2),    &
                field_grid%p(ip)%results%A(1:3), field_grid%p(ip)%results%B(1:3), field_grid%p(ip)%results%J(1:3),  &
                field_grid%p(ip)%results%pot,field_grid%p(ip)%results%dxA(1:3),field_grid%p(ip)%results%dyA(1:3)
!              else if (field_grid%p(ip)%label .eq. -1) then
!                write(filehandle_e,'(14(f8.3,x),i12)') field_grid%p(ip)%x(1:2), field_grid%p(ip)%results%E(1:2),    &
!                field_grid%p(ip)%results%A(1:3), field_grid%p(ip)%results%B(1:3), field_grid%p(ip)%results%J(1:3),  &
!                field_grid%p(ip)%results%pot,field_grid%p(ip)%label
!              else if (field_grid%p(ip)%label .eq. 0) then
!                write(filehandle_b,'(14(f8.3,x),i12)') field_grid%p(ip)%x(1:2), field_grid%p(ip)%results%E(1:2),    &
!                field_grid%p(ip)%results%A(1:3), field_grid%p(ip)%results%B(1:3), field_grid%p(ip)%results%J(1:3),  &
!                field_grid%p(ip)%results%pot,field_grid%p(ip)%label
!              endif
            end do
            close(filehandle_i)
!            close(filehandle_e)
!            close(filehandle_b)

        end subroutine

!      subroutine write_field_on_grid_ascii(field_grid,itime)
!        use module_globals, only: woutput,folder
!        use module_vtk
!        use encap
!        implicit none
!
!        type(field_grid_t), intent(in)         :: field_grid
!        integer           , intent(in)         :: itime
!        integer(kind_particle)                 :: ip,nl
!        integer(kind_particle),      parameter :: ifield = 10
!        integer(kind_particle),      parameter :: sizeout = 19
!        real(kind_particle)                    :: tmp(sizeout)
!        character(255)                         :: str,file0
!
!        write( str, '(i10)' ) itime
!
!        file0 = "fields/field_on_grid_ascii_"
!        if (root ) open (unit=ifield,file=trim(folder)//trim(file0)//trim(adjustl(str))//".dat",action="write",status="replace")
!
!        nl = field_grid%nl
!
!        select case (woutput)
!              case (1)  ! 1D
!
!                do ip = 1,nl
!                    tmp(1) = field_grid%p(ip)%x(1)
!                    tmp(2) = field_grid%p(ip)%results%E(1)
!                    tmp(3) = field_grid%p(ip)%results%pot
!                    tmp(4) = field_grid%p(ip)%results%A(1)
!                    tmp(5) = field_grid%p(ip)%results%B(3)
!
!                    write (ifield,*) tmp(1:4)
!
!                enddo
!              case (2)  !  2D
!
!                do ip = 1,nl
!                    tmp(1) = field_grid%p(ip)%x(1)
!                    tmp(2) = field_grid%p(ip)%x(2)
!                    tmp(3) = field_grid%p(ip)%results%E(1)
!                    tmp(4) = field_grid%p(ip)%results%E(2)
!                    tmp(5) = field_grid%p(ip)%results%pot
!                    tmp(6) = field_grid%p(ip)%results%A(1)
!                    tmp(7) = field_grid%p(ip)%results%A(2)
!                    tmp(8) = field_grid%p(ip)%results%B(3)
!                    tmp(9) = field_grid%p(ip)%results%J(1)
!                    tmp(10)= field_grid%p(ip)%results%J(2)
!                    tmp(11)= field_grid%p(ip)%results%Jirr(1)
!                    tmp(12)= field_grid%p(ip)%results%Jirr(2)
!
!                    write (ifield,*) tmp(1:12)
!
!                enddo
!
!              case (3)  ! 2D-3V
!
!                do ip = 1,nl
!                    tmp(1) = field_grid%p(ip)%x(1)
!                    tmp(2) = field_grid%p(ip)%x(2)
!                    tmp(3) = field_grid%p(ip)%results%E(1)
!                    tmp(4) = field_grid%p(ip)%results%E(2)
!                    tmp(5) = field_grid%p(ip)%results%pot
!                    tmp(6) = field_grid%p(ip)%results%A(1)
!                    tmp(7) = field_grid%p(ip)%results%A(2)
!                    tmp(8) = field_grid%p(ip)%results%A(3)
!                    tmp(9) = field_grid%p(ip)%results%B(1)
!                    tmp(10)= field_grid%p(ip)%results%B(2)
!                    tmp(11)= field_grid%p(ip)%results%B(3)
!                    tmp(12)= field_grid%p(ip)%results%J(1)
!                    tmp(13)= field_grid%p(ip)%results%J(2)
!                    tmp(14)= field_grid%p(ip)%results%J(3)
!                    tmp(15)= field_grid%p(ip)%results%Jirr(1)
!                    tmp(16)= field_grid%p(ip)%results%Jirr(2)
!                    tmp(17)= field_grid%p(ip)%results%Jirr(2)
!
!                    write (ifield,*) tmp(1:17)
!
!                enddo
!
!              case (4)  ! 3D
!
!                do ip = 1,nl
!                    tmp(1) = field_grid%p(ip)%x(1)
!                    tmp(2) = field_grid%p(ip)%x(2)
!                    tmp(3) = field_grid%p(ip)%x(3)
!                    tmp(4) = field_grid%p(ip)%results%E(1)
!                    tmp(5) = field_grid%p(ip)%results%E(2)
!                    tmp(6) = field_grid%p(ip)%results%E(3)
!                    tmp(7) = field_grid%p(ip)%results%pot
!                    tmp(8) = field_grid%p(ip)%results%A(1)
!                    tmp(9) = field_grid%p(ip)%results%A(2)
!                    tmp(10)= field_grid%p(ip)%results%A(3)
!                    tmp(11)= field_grid%p(ip)%results%B(1)
!                    tmp(12)= field_grid%p(ip)%results%B(2)
!                    tmp(13)= field_grid%p(ip)%results%B(3)
!                    tmp(14)= field_grid%p(ip)%results%J(1)
!                    tmp(15)= field_grid%p(ip)%results%J(2)
!                    tmp(16)= field_grid%p(ip)%results%J(3)
!                    tmp(17)= field_grid%p(ip)%results%Jirr(1)
!                    tmp(18)= field_grid%p(ip)%results%Jirr(2)
!                    tmp(19)= field_grid%p(ip)%results%Jirr(2)
!
!                    write (ifield,*) tmp(1:19)
!
!                enddo
!
!              case default
!
!
!                do ip = 1,nl
!                    tmp(1) = field_grid%p(ip)%x(1)
!                    tmp(2) = field_grid%p(ip)%results%E(1)
!                    tmp(3) = field_grid%p(ip)%results%A(1)
!                    tmp(4) = field_grid%p(ip)%results%B(3)
!
!                    write (ifield,*) tmp(1:4)
!
!                enddo
!
!        end select
!
!        if (root )   close ( unit=ifield )
!
!      end subroutine
!
      subroutine write_particles(p)
        use module_vtk
        implicit none

        type(t_particle), allocatable, intent(in) :: p(:)

        integer(kind_particle) :: i
        type(vtkfile_unstructured_grid) :: vtk
        integer :: vtk_step
        real*8 :: time,vect(np)
        real*8 :: ta, tb

        ta = get_time()
        time = dt * step

        if (step .eq. 0) then
          vtk_step = VTK_STEP_FIRST
        else if (step .eq. nt) then
          vtk_step = VTK_STEP_LAST
        else
          vtk_step = VTK_STEP_NORMAL
        endif

        vect = zero

        call vtk%create_parallel("particles", step, my_rank, n_ranks, time, vtk_step)
        call vtk%write_headers(np, 0_kind_particle)
        call vtk%startpoints()
        call vtk%write_data_array("x", p(1:np)%x(1), p(1:np)%x(2), p(1:np)%x(3) )
        call vtk%finishpoints()
        call vtk%startpointdata()
        call vtk%write_data_array("v", p(1:np)%data%v(1), &
                                       p(1:np)%data%v(2), &
                                       p(1:np)%data%v(3) )

        call vtk%write_data_array("E", p(1:np)%results%e(1), &
                                       p(1:np)%results%e(2), &
                                       p(1:np)%results%e(3) )

        call vtk%write_data_array("A", p(1:np)%results%A(1), &
                                       p(1:np)%results%A(2), &
                                       p(1:np)%results%A(3) )
                                       
        call vtk%write_data_array("dxA", p(1:np)%results%dxA(1), &
                                       p(1:np)%results%dxA(2)  , &
                                       p(1:np)%results%dxA(3) )
                                       
        call vtk%write_data_array("dyA", p(1:np)%results%dyA(1), &
                                       p(1:np)%results%dyA(2)  , &
                                       p(1:np)%results%dyA(3) )


        call vtk%write_data_array("B",   p(1:np)%results%B(1),&
                                         p(1:np)%results%B(2),&
                                         p(1:np)%results%B(3) )

        call vtk%write_data_array("J",   p(1:np)%results%J(1),&
                                         p(1:np)%results%J(2),&
                                         p(1:np)%results%J(3) )

        call vtk%write_data_array("Jirr",p(1:np)%results%Jirr(1),&
                                         p(1:np)%results%Jirr(2),&
                                         p(1:np)%results%Jirr(3) )

        call vtk%write_data_array("phi", p(1:np)%results%pot)
        call vtk%write_data_array("q", p(1:np)%data%q)
        call vtk%write_data_array("gamma", p(1:np)%data%g)
        call vtk%write_data_array("m", p(1:np)%data%m)
        call vtk%write_data_array("work", p(1:np)%work)
        call vtk%write_data_array("plabel", p(1:np)%label)
        call vtk%write_data_array("local index", [(i,i=1,np)])
        call vtk%write_data_array("processor", int(np, kind = 4), my_rank)
!        if(particle_test) call vtk%write_data_array("L2 error", direct_L2(1:np))
        call vtk%finishpointdata()
        call vtk%dont_write_cells()
        call vtk%write_final()
        call vtk%close()

        tb = get_time()

        if(root) write(*,'(a,es12.4)') " == [write particles] time in vtk output [s]      : ", tb - ta

      end subroutine write_particles



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

!        subroutine vtk_write_particle_data_results(d, r, vtkf)
!        use module_vtk
!        implicit none
!
!        type(t_particle_data)          , intent(in)    :: d(:)
!        type(t_particle_results)       , intent(in)    :: r(:)
!        type(vtkfile_unstructured_grid), intent(inout) :: vtkf
      
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
        
        
        subroutine write_particles_ascii_nospecies(itime, p)    
            use module_pepc_types
            use module_utils
            use module_globals, only: folder
            implicit none
!            integer(kind_pe), intent(in) :: my_rank
            integer(kind_default), intent(in)               :: itime
            type(t_particle)     , intent(in), dimension(:) :: p
            logical                                         :: firstcall  = .true.
            character(50)                                   :: dir
            character(100)                                  :: filename_i
            integer(kind_particle)                          :: ip
            
            character(12), parameter                        :: part_dir = 'particles/'
            integer, parameter :: filehandle_i = 40
              call vtk_write_particle_data_results(d, r, vtkf)
            end subroutine
            
        end subroutine write_particles_vtk
        
        
        subroutine write_particles_ascii(itime, p)    
            use module_pepc_types
            use module_utils
            use module_globals, only: folder
            implicit none
!            integer(kind_pe), intent(in) :: my_rank
            integer(kind_default), intent(in)               :: itime
            type(t_particle)     , intent(in), dimension(:) :: p
            logical                                         :: firstcall  = .true.
            character(50)                                   :: dir
            character(100)                                  :: filename_i,filename_e,filename_b
            integer(kind_particle)                          :: ip
            
            character(12), parameter                        :: part_dir = 'particles/'
            integer, parameter :: filehandle_i = 40
            integer, parameter :: filehandle_e = 41
            integer, parameter :: filehandle_b = 42

            
            if (firstcall) then
              call create_directory(trim(folder))
              call create_directory(trim(folder)//trim(part_dir))
              firstcall = .false.
            endif
            
            write(filename_i,'(a,"particle_",i6.6,"_",i6.6,".dat")') trim(folder)//trim(part_dir), itime, my_rank
            
            open(filehandle_i, file=trim(filename_i), STATUS='REPLACE')
            
            do ip=1, size(p,kind=kind(ip))
                write(filehandle_i,'(26(f8.3,x))') p(ip)%x(1:3), p(ip)%data%v(1:3), p(ip)%results%E(1:3),&
                p(ip)%results%A(1:3), p(ip)%results%B(1:3), p(ip)%results%J(1:3),p(ip)%results%pot, p(ip)%data%g,&
                p(ip)%results%dxA(1:3),p(ip)%results%dyA(1:3)
            end do
            close(filehandle_i)
            
            if (firstcall) then
              call create_directory(trim(folder))
              call create_directory(trim(folder)//trim(part_dir))
              firstcall = .false.
            endif
            
            write(filename_i,'(a,"particle_ions_",i6.6,"_",i6.6,".dat")') trim(folder)//trim(part_dir), itime, my_rank
            write(filename_e,'(a,"particle_elec_",i6.6,"_",i6.6,".dat")') trim(folder)//trim(part_dir), itime, my_rank
            write(filename_b,'(a,"particle_beam_",i6.6,"_",i6.6,".dat")') trim(folder)//trim(part_dir), itime, my_rank

        end subroutine
        
        
        subroutine write_particles_ascii(itime, p)    
            use module_pepc_types
            use module_utils
            use module_globals, only: folder
            implicit none
!            integer(kind_pe), intent(in) :: my_rank
            integer(kind_default), intent(in)               :: itime
            type(t_particle)     , intent(in), dimension(:) :: p
            logical                                         :: firstcall  = .true.
            character(50)                                   :: dir
            character(100)                                  :: filename_i,filename_e,filename_b
            integer(kind_particle)                          :: ip
            
            character(12), parameter                        :: part_dir = 'particles/'
            integer, parameter :: filehandle_i = 40
            integer, parameter :: filehandle_e = 41
            integer, parameter :: filehandle_b = 42
            

            
            if (firstcall) then
              call create_directory(trim(folder))
              call create_directory(trim(folder)//trim(part_dir))
              firstcall = .false.
            endif
            
            write(filename_i,'(a,"particle_ions_",i6.6,"_",i6.6,".dat")') trim(folder)//trim(part_dir), itime, my_rank
            write(filename_e,'(a,"particle_elec_",i6.6,"_",i6.6,".dat")') trim(folder)//trim(part_dir), itime, my_rank
            write(filename_b,'(a,"particle_beam_",i6.6,"_",i6.6,".dat")') trim(folder)//trim(part_dir), itime, my_rank
            open(filehandle_i, file=trim(filename_i), STATUS='REPLACE')
            open(filehandle_e, file=trim(filename_e), STATUS='REPLACE')
            open(filehandle_b, file=trim(filename_b), STATUS='REPLACE')
            
            do ip=1, size(p,kind=kind(ip))
              if (p(ip)%label .eq. 1)  then
                write(filehandle_i,'(26(f8.3,x))') p(ip)%x(1:3), p(ip)%data%v(1:3), p(ip)%results%E(1:3),&
                p(ip)%results%A(1:3), p(ip)%results%B(1:3), p(ip)%results%J(1:3),p(ip)%results%pot, p(ip)%data%g,&
                p(ip)%results%dxA(1:3),p(ip)%results%dyA(1:3)
              else if (p(ip)%label .eq. -1) then
                write(filehandle_e,'(26(f8.3,x),i12)') p(ip)%x(1:3), p(ip)%data%v(1:3), p(ip)%results%E(1:3),&
                p(ip)%results%A(1:3), p(ip)%results%B(1:3), p(ip)%results%J(1:3),p(ip)%results%pot,  p(ip)%data%g,&
                p(ip)%results%dxA(1:3),p(ip)%results%dyA(1:3)
              else if (p(ip)%label .eq. 0) then
                write(filehandle_b,'(26(f8.3,x),i12)') p(ip)%x(1:3), p(ip)%data%v(1:3), p(ip)%results%E(1:3),&
                p(ip)%results%A(1:3), p(ip)%results%B(1:3), p(ip)%results%J(1:3),p(ip)%results%pot,  p(ip)%data%g,&
                p(ip)%results%dxA(1:3),p(ip)%results%dyA(1:3)
              endif
            end do
            close(filehandle_i)
            close(filehandle_e)
            close(filehandle_b)

            

        end subroutine
  
            open(filehandle_i, file=trim(filename_i), STATUS='REPLACE')
            open(filehandle_e, file=trim(filename_e), STATUS='REPLACE')
            open(filehandle_b, file=trim(filename_b), STATUS='REPLACE')
            
            do ip=1, size(p,kind=kind(ip))
              if (p(ip)%label .eq. 1)  then
                write(filehandle_i,'(26(f8.3,x))') p(ip)%x(1:3), p(ip)%data%v(1:3), p(ip)%results%E(1:3),&
                p(ip)%results%A(1:3), p(ip)%results%B(1:3), p(ip)%results%J(1:3),p(ip)%results%pot, p(ip)%data%g,&
                p(ip)%results%dxA(1:3),p(ip)%results%dyA(1:3)
              else if (p(ip)%label .eq. -1) then
                write(filehandle_e,'(26(f8.3,x),i12)') p(ip)%x(1:3), p(ip)%data%v(1:3), p(ip)%results%E(1:3),&
                p(ip)%results%A(1:3), p(ip)%results%B(1:3), p(ip)%results%J(1:3),p(ip)%results%pot,  p(ip)%data%g,&
                p(ip)%results%dxA(1:3),p(ip)%results%dyA(1:3)
              else if (p(ip)%label .eq. 0) then
                write(filehandle_b,'(26(f8.3,x),i12)') p(ip)%x(1:3), p(ip)%data%v(1:3), p(ip)%results%E(1:3),&
                p(ip)%results%A(1:3), p(ip)%results%B(1:3), p(ip)%results%J(1:3),p(ip)%results%pot,  p(ip)%data%g,&
                p(ip)%results%dxA(1:3),p(ip)%results%dyA(1:3)
              endif
            end do
            close(filehandle_i)
            close(filehandle_e)
            close(filehandle_b)

        end subroutine
  
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
