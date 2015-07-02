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
!  real(kind_particle) :: Lx                       ! box size - x direction
!  real(kind_particle) :: Ly                       ! box size - y direction
!  real(kind_particle) :: Lz                       ! box size - z direction

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
    !use field_helper, only:
    implicit none

    integer, parameter :: fid = 12
    character(255)     :: para_file
    logical            :: read_para_file

    namelist /pepc2dd/ tnp, nt, dt, particle_output, domain_output,  particle_test, &
       diag_interval, periodicity,  spatial_interaction_cutoff,Lx,Ly,Lz,initial_setup,ischeme,restart_file,restart_step!,   &

    ! set default parameter values
    initial_setup   = 2
    tnp             = 1441
    nt              = 20
    dt              = 1e-3
    particle_test   = .true.
    particle_output = .true.
    ischeme         = "leapfrog"
    restart_file    = ""
    restart_step    = 1

    diag_interval   = 2


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

    tnp = tnp

    if(root) then
      write(*,'(a,i12)')    " == total number of particles : ", tnp
      write(*,'(a,i12)')    " == number of time steps      : ", nt
      write(*,'(a,es12.4)') " == time step                 : ", dt
      write(*,*)            "== time integration scheme   : ",   ischeme
      write(*,'(a,es12.4)') " == box size  Lx              : ", Lx
      write(*,'(a,es12.4)') " == box size  Ly              : ", Ly
      write(*,'(a,es12.4)') " == box size  Lz              : ", Lz
      write(*,'(a,es12.4)') " == theta2                    : ", theta2
      write(*,'(a,es12.4)') " == epsilon2                  : ", eps2
      write(*,'(a,l12)')    " == periodicity               : ", periodicity
      write(*,'(a,i12)')    " == diag & IO interval        : ", diag_interval
      write(*,'(a,l12)')    " == particle test             : ", particle_test
      write(*,'(a,l12)')    " == particle output           : ", particle_output
      write(*,'(a,l12)')    " == domain output             : ", domain_output
      write(*,*) ""

    end if

    call pepc_prepare(1_kind_dim)

    ! reste current measurement
    current_q = 0.0_8
    current_I = 0.0_8

  end subroutine

  subroutine load_status_2d(p)
    implicit none
    include 'mpif.h'
    type(t_particle), allocatable, intent(inout) :: p(:)
    integer(kind_particle)                       :: ip,ii,jj,Ncol,ierr,rc,start,stop
    real(kind_particle),allocatable              :: tmp(:,:)

    if (allocated(tmp)) deallocate(tmp)

    allocate(tmp(np,Ncol), stat=rc)

    if(my_rank.eq.0) write(*,'(a)') " == [init] init particles "

    np = tnp / n_ranks
    if(my_rank.eq.(n_ranks-1)) np = np + MOD(tnp, int(n_ranks, kind=kind_particle))

    allocate(particles(np), stat=rc)
    if(rc.ne.0) write(*,*) " === particle allocation error!"

    allocate(direct_L2(np), stat=rc)
    if(rc.ne.0) write(*,*) " === direct_L2 allocation error!"
    direct_L2 = -1.0_8


    if (root)  open(unit=rc, file='test.dat', form='formatted', status='old', action='read')
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    start = 1 + my_rank*np
    stop  = 1 + (my_rank + 1)*np

    do ip = start, stop

      read (rc, *) (tmp (ip, jj), jj=1,Ncol)

      p(ip)%label       = tmp(ip,1)
      p(ip)%data%q      = tmp(ip,2)
      p(ip)%data%m      = 1.0_8

      if(p(ip)%data%q .gt. 0.0) p(ip)%data%m = p(ip)%data%m * 100.0_8

      p(ip)%x(1)        = 0.0_8
      p(ip)%x(2)        = 0.0_8
      p(ip)%x(3)        = 0.0_8

      p(ip)%data%v(1)   = 0.0_8
      p(ip)%data%v(2)   = 0.0_8
      p(ip)%data%v(3)   = 0.0_8

      p(ip)%results%e   = 0.0_8
      p(ip)%results%pot = 0.0_8
      p(ip)%results%B   = 0.0_8
      p(ip)%results%A   = 0.0_8
      p(ip)%results%Jirr= 0.0_8
      p(ip)%results%J   = 0.0_8
      p(ip)%work        = 1.0_8


    end do

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    if (root) close(unit = rc)

    deallocate(tmp)

  end subroutine load_status_2d

  subroutine read_restart_2d(p,filename)
    implicit none
    character(*)                                 :: filename
    type(t_particle), allocatable, intent(out)   :: p(:)
    integer(kind_particle)                       :: ip,rc,jp,size_tmp,io_stat,open_status
    real(kind_particle), allocatable             :: tmp(:,:)
    character(255)                               :: str_proc

    if(my_rank.eq.0) write(*,'(a)') " == [init] init particles "

    np = tnp / n_ranks
    if(my_rank.eq.(n_ranks-1)) np = np + MOD(tnp, int(n_ranks, kind=kind_particle))

    allocate(particles(np), stat=rc)
    if(rc.ne.0) write(*,*) " === particle allocation error!"

    allocate(direct_L2(np), stat=rc)
    if(rc.ne.0) write(*,*) " === direct_L2 allocation error!"
    direct_L2 = -1.0_8

!    if (root) then

    size_tmp = 6
    if (allocated(tmp)) deallocate(tmp)
    allocate(tmp(np,size_tmp), stat = rc)
    write( str_proc , '(i10)' ) my_rank
    open (unit=my_rank,file=trim(filename)//trim(adjustl(str_proc))//".dat",action="read", position='rewind')

    do ip = 1,np
        read(my_rank,*) tmp(ip,1:size_tmp)

        p(ip)%label       = tmp(ip,1)
        p(ip)%data%q      = tmp(ip,2)
        p(ip)%data%m      = 1.0_8

        if(p(ip)%data%q .gt. 0.0) p(ip)%data%m = p(ip)%data%m * 100.0_8

        p(ip)%x(1)        = tmp(ip,3)
        p(ip)%x(2)        = tmp(ip,4)
        p(ip)%x(3)        = 0.0_8

        p(ip)%data%v(1)   = tmp(ip,5)
        p(ip)%data%v(2)   = tmp(ip,6)
        p(ip)%data%v(3)   = 0.0_8

        p(ip)%results%e   = zero
        p(ip)%results%pot = zero
        p(ip)%results%B   = zero
        p(ip)%results%A   = zero
        p(ip)%results%Jirr= zero
        p(ip)%results%J   = zero
        p(ip)%work        = one
        write(*,*) tmp(ip,1:size_tmp)

    end do


    close(unit=my_rank)
    deallocate(tmp)

  end subroutine read_restart_2d


  subroutine write_restart_2d(p,istep)
    implicit none
    type(t_particle), allocatable, intent(in)    :: p(:)
    integer(kind_particle)       , intent(in)    :: istep
    integer(kind_particle)                       :: ip,rc,size_tmp
    real(kind_particle), allocatable             :: tmp(:)
    character(255)                               :: str_istep,str_proc



    if (root)  then

        size_tmp = 6
        if (allocated(tmp)) deallocate(tmp)
        allocate(tmp(size_tmp), stat = rc)
!        open(unit=rc,file="restart.dat",form='formatted',status='unknown',access='append')
        write( str_istep, '(i10)' ) istep
        write( str_proc , '(i10)' ) my_rank
        open (unit=my_rank,file=trim("restart/restart_")//trim("_")//trim(adjustl(str_proc))//".dat",action="write",status="replace")
!        open (unit=my_rank,file=trim("restart/restart_")//trim(adjustl(str_istep))//trim("_")//trim(adjustl(str_proc))//".dat",action="write",status="replace")
!        write (rc,*) " "

        do ip = 1, tnp


              tmp(1)            =   p(ip)%label
              tmp(2)            =   p(ip)%data%q

              tmp(3)            =   p(ip)%x(1)
              tmp(4)            =   p(ip)%x(2)

              tmp(5)            =   p(ip)%data%v(1)
              tmp(6)            =   p(ip)%data%v(2)

              write (rc,*) tmp


        end do

        close(unit = my_rank)
        deallocate(tmp)
    endif



  end subroutine write_restart_2d


  subroutine init_langmuir_waves(p)
    implicit none
    type(t_particle), allocatable, intent(inout) :: p(:)
    integer(kind_particle) :: ip, jp, nppd
    integer(kind_particle) :: rc, iter
    real(kind_particle),parameter ::   prec = 1e-5

    if(my_rank.eq.0) write(*,'(a)') " == [init] init particles "

    ! set initially number of local particles
    np = tnp / n_ranks
    if(my_rank.eq.(n_ranks-1)) np = np + MOD(tnp, int(n_ranks, kind=kind_particle))

    allocate(particles(np), stat=rc)
    if(rc.ne.0) write(*,*) " === particle allocation error!"

    allocate(direct_L2(np), stat=rc)
    if(rc.ne.0) write(*,*) " === direct_L2 allocation error!"
    direct_L2 = -1.0_8

    nppd = sqrt( real(np,8) )

    iter = 1

    do ip=1, nppd

        do jp =1, nppd
          p(iter)%label       = my_rank * (tnp / n_ranks) + iter
          p(iter)%data%q      = (-1.0_8 + 2.0_8*MOD(p(ip)%label,2_kind_particle))/real(tnp, kind=kind_particle)
          p(iter)%data%m      = 1.0_8/real( tnp, kind=kind_particle )

!          p(iter)%x(1)         = (ip - MOD(p(ip)%label,2_kind_particle) )/( nppd - 1.0_8)
!          p(iter)%x(2)         = (jp - MOD(p(ip)%label,2_kind_particle) )/( nppd - 1.0_8)
          p(iter)%x(1)         = (ip - 1.0_8 )/( nppd - 1.0_8)
          p(iter)%x(2)         = (jp - 1.0_8 )/( nppd - 1.0_8)
          p(iter)%x(3)         = 0.0_8

          if(p(iter)%data%q .gt. 0.0)      p(iter)%data%m = p(iter)%data%m * 1846.0_8


          call random_gauss(p(iter)%data%v)

          p(iter)%data%v(3)   = 0.0_8

          p(iter)%results%e   = 0.0_8
          p(iter)%results%pot = 0.0_8
          p(iter)%results%B   = 0.0_8
          p(iter)%results%A   = 0.0_8
          p(iter)%results%Jirr= 0.0_8
          p(iter)%results%J   = 0.0_8
          p(iter)%work        = 1.0_8

          iter       = iter + 1
          next_label = tnp+1

        end do
    end do

  end subroutine init_langmuir_waves


  subroutine init_particles_square(p)
    implicit none
    type(t_particle), allocatable, intent(inout) :: p(:)
    integer(kind_particle) :: ip, jp, nppd
    integer(kind_particle) :: rc, iter
    real(kind_particle),parameter ::   prec = 1e-5

    if(my_rank.eq.0) write(*,'(a)') " == [init] init particles "

    ! set initially number of local particles
    np = tnp / n_ranks
    if(my_rank.eq.(n_ranks-1)) np = np + MOD(tnp, int(n_ranks, kind=kind_particle))

    allocate(particles(np), stat=rc)
    if(rc.ne.0) write(*,*) " === particle allocation error!"

    allocate(direct_L2(np), stat=rc)
    if(rc.ne.0) write(*,*) " === direct_L2 allocation error!"
    direct_L2 = -1.0_8

    nppd = sqrt( real(np,8) )

    iter = 1

    do ip=1, nppd

        do jp =1, nppd
          p(iter)%label       = my_rank * (tnp / n_ranks) + iter
          p(iter)%data%q      = (-1.0_8 + 2.0_8*MOD(p(ip)%label,2_kind_particle))
          p(iter)%data%m      = 1.0_8

!          p(iter)%x(1)         = (ip - MOD(p(ip)%label,2_kind_particle) )/( nppd - 1.0_8)
!          p(iter)%x(2)         = (jp - MOD(p(ip)%label,2_kind_particle) )/( nppd - 1.0_8)
          p(iter)%x(1)         = (ip - 1.0_8 )/( nppd - 1.0_8)
          p(iter)%x(2)         = (jp - 1.0_8 )/( nppd - 1.0_8)
          p(iter)%x(3)         = 0.0_8

          if(p(iter)%data%q .gt. 0.0) then
                p(iter)%data%m = p(iter)%data%m * 100.0_8
!                p(iter)%x(1)         = p(iter)%x(1) + prec
          endif



          call random_gauss(p(iter)%data%v)

          !p(iter)%data%v(1)   = 1.0_8
          !p(iter)%data%v(2)   = 1.0_8
          p(iter)%data%v(3)   = 0.0_8

          p(iter)%results%e   = 0.0_8
          p(iter)%results%pot = 0.0_8
          p(iter)%results%B   = 0.0_8
          p(iter)%results%A   = 0.0_8
          p(iter)%results%Jirr= 0.0_8
          p(iter)%results%J   = 0.0_8
          p(iter)%work        = 1.0_8

          iter       = iter + 1
          next_label = tnp+1

        end do
    end do

  end subroutine init_particles_square


  subroutine init_weibell_instalility(p)
    implicit none
    type(t_particle), allocatable, intent(inout) :: p(:)
    integer(kind_particle) :: ip, jp, nppd
    integer(kind_particle) :: rc, iter
    real(kind_particle),parameter ::   prec = 1e-5

    if(my_rank.eq.0) write(*,'(a)') " == [init] init particles "

    ! set initially number of local particles
    np = tnp / n_ranks
    if(my_rank.eq.(n_ranks-1)) np = np + MOD(tnp, int(n_ranks, kind=kind_particle))

    allocate(particles(np), stat=rc)
    if(rc.ne.0) write(*,*) " === particle allocation error!"

    allocate(direct_L2(np), stat=rc)
    if(rc.ne.0) write(*,*) " === direct_L2 allocation error!"
    direct_L2 = -1.0_8

    nppd = sqrt( real(np,8) )

    iter = 1

    do ip=1, nppd

        do jp =1, nppd
          p(iter)%label       = my_rank * (tnp / n_ranks) + iter
          p(iter)%data%q      = (-1.0_8 + 2.0_8*MOD(p(ip)%label,2_kind_particle))
          p(iter)%data%m      = 1.0_8

          p(iter)%x(1)         = (ip - 1.0_8 )/( nppd - 1.0_8)
          p(iter)%x(2)         = (jp - 1.0_8 )/( nppd - 1.0_8)
          p(iter)%x(3)         = 0.0_8

          if(p(iter)%data%q .gt. 0.0) then
                p(iter)%data%m = p(iter)%data%m * 100.0_8
!                p(iter)%x(1)         = p(iter)%x(1) + prec
          endif



          call random_gauss(p(iter)%data%v)

          !p(iter)%data%v(1)   = 1.0_8
          !p(iter)%data%v(2)   = 1.0_8
          p(iter)%data%v(3)   = 0.0_8

          p(iter)%results%e   = 0.0_8
          p(iter)%results%pot = 0.0_8
          p(iter)%results%B   = 0.0_8
          p(iter)%results%A   = 0.0_8
          p(iter)%results%Jirr= 0.0_8
          p(iter)%results%J   = 0.0_8
          p(iter)%work        = 1.0_8

          iter       = iter + 1
          next_label = tnp+1

        end do
    end do

  end subroutine init_weibell_instalility


subroutine init_particles(p)
    implicit none

    type(t_particle), allocatable, intent(inout) :: p(:)
    integer(kind_particle) :: ip
    integer :: rc
    real(kind_particle)    :: dummy

    real,parameter :: Tkb=1e-1

    if(my_rank.eq.0) write(*,'(a)') " == [init] init particles "

    ! set initially number of local particles
    np = tnp / n_ranks
    if(my_rank.eq.(n_ranks-1)) np = np + MOD(tnp, int(n_ranks, kind=kind_particle))

    allocate(particles(np), stat=rc)
    if(rc.ne.0) write(*,*) " === particle allocation error!"

    allocate(direct_L2(np), stat=rc)
    if(rc.ne.0) write(*,*) " === direct_L2 allocation error!"
    direct_L2 = -1.0_8


    ! set random seed
!    write(*,*) "my_rank/n_ranks", my_rank,n_ranks
    dummy = par_rand(my_rank)

    ! setup random qubic particle cloud
    do ip=1, np

      p(ip)%label       = my_rank * (tnp / n_ranks) + ip
      p(ip)%data%q      = (-1.0_8 + 2.0_8*MOD(p(ip)%label,2_kind_particle))
      p(ip)%data%m      = 1.0_8
!      if(p(ip)%data%q .gt. 0.0) p(ip)%data%m = p(ip)%data%m * 100.0_8

      call random(p(ip)%x)

      call random_gauss(p(ip)%data%v)
      p(ip)%data%v      = p(ip)%data%v * sqrt(Tkb * p(ip)%data%m) / p(ip)%data%m

      p(ip)%x(3)        = 0.0_8
      p(ip)%data%v(3)   = 0.0_8
      !write(*,*) p(ip)%data%v

      p(ip)%results%e   = 0.0_8
      p(ip)%results%pot = 0.0_8
      p(ip)%results%B   = 0.0_8
      p(ip)%results%A   = 0.0_8
      p(ip)%results%Jirr= 0.0_8
      p(ip)%results%J   = 0.0_8
      p(ip)%work        = 1.0_8

    end do

    next_label = tnp+1

  end subroutine init_particles

  subroutine init_landau_damping2d(p)
    use zufall, only: random_gaussian_flux
    implicit none

    type(t_particle), allocatable, intent(inout) :: p(:)
    integer(kind_particle) :: ip,jp,iter
    integer :: rc
    real(kind_particle) :: dummy,dx,L,nppd

    real(kind_particle),parameter :: vtherm = 1e-2
    real(kind_particle),parameter :: v0     = .5_8*1e-2
    real(kind_particle),parameter :: x0     = .5_8*1e-2

    if(my_rank.eq.0) write(*,'(a)') " == [init] init particles "

    ! set initially number of local particles
    np = tnp / n_ranks
    if(my_rank.eq.(n_ranks-1)) np = np + MOD(tnp, int(n_ranks, kind=kind_particle))

    allocate(particles(np), stat=rc)
    if(rc.ne.0) write(*,*) " === particle allocation error!"

    allocate(direct_L2(np), stat=rc)
    if(rc.ne.0) write(*,*) " === direct_L2 allocation error!"
    direct_L2 = -1.0_8



    ! set random seed
    !L   = 2.0_8*pi*sqrt(3.0_8)/real(n_ranks,kind=kind_particle )

    nppd = int( sqrt( real(np,8) ), kind = kind_particle )

    L   = 4.0_8*pi
    dx  = L/real( n_ranks, kind=kind_particle )

    iter = 1

    do ip=1, nppd

        do jp =1, nppd

          p(iter)%label       = my_rank * (tnp / n_ranks) + iter
          p(iter)%data%q      = (-1.0_8 + 2.0_8*MOD(p(ip)%label,2_kind_particle))
          p(iter)%data%m      = 1.0_8

          p(iter)%x(1)         = (ip - 1.0_8 )/( nppd - 1.0_8)
          p(iter)%x(2)         = (jp - 1.0_8 )/( nppd - 1.0_8)
          p(iter)%x(3)         = 0.0_8

          if(p(iter)%data%q .gt. 0.0)  then

                p(iter)%data%m = p(iter)%data%m * 100.0_8
                p(iter)%data%v = 0.0_8

          else

                call random_gaussian_flux(p(ip)%data%v(1),vtherm)
                call random_gaussian_flux(p(ip)%data%v(2),vtherm)
                p(ip)%data%v(3) = 0.0_8

          endif

          p(iter)%results%e   = 0.0_8
          p(iter)%results%pot = 0.0_8
          p(iter)%results%B   = 0.0_8
          p(iter)%results%A   = 0.0_8
          p(iter)%results%Jirr= 0.0_8
          p(iter)%results%J   = 0.0_8
          p(iter)%work        = 1.0_8

          iter       = iter + 1
          next_label = tnp+1

        end do
    end do

    p(:)%x(1)        = p(:)%x(1)*dx + real( my_rank, kind= kind_particle )*dx
    p(:)%x(2)        = p(:)%x(2)*L

    !!! perturbation

    p(:)%x(1)        = p(:)%x(1)      +  x0*cos( 2.0_8*pi/L*p(:)%x(1) )
    p(:)%x(2)        = p(:)%x(2)      +  x0*cos( 2.0_8*pi/L*p(:)%x(2) )

    p(:)%data%v(1)   = p(:)%data%v(1) +  v0*cos( 2.0_8*pi/L*p(:)%x(1) )
    p(:)%data%v(2)   = p(:)%data%v(2) +  v0*cos( 2.0_8*pi/L*p(:)%x(2) )


    next_label = tnp+1

  end subroutine init_landau_damping2d


  subroutine init_landau_damping1d(p)
    use zufall, only: random_gaussian_flux
    implicit none

    type(t_particle), allocatable, intent(inout) :: p(:)
    integer(kind_particle) :: ip,jp
    integer :: rc
    real(kind_particle) :: dummy,dx,L

    real(kind_particle),parameter :: vtherm = 0.20_8
    real(kind_particle),parameter :: v0     = 0.0_8
!    real(kind_particle),parameter :: x0     = 0.0_8
    real(kind_particle),parameter :: x0     = 1e-1

    if(my_rank.eq.0) write(*,'(a)') " == [init] init particles "

    ! set initially number of local particles
    np = tnp / n_ranks
    if(my_rank.eq.(n_ranks-1)) np = np + MOD(tnp, int(n_ranks, kind=kind_particle))

    allocate(particles(np), stat=rc)
    if(rc.ne.0) write(*,*) " === particle allocation error!"

    allocate(direct_L2(np), stat=rc)
    if(rc.ne.0) write(*,*) " === direct_L2 allocation error!"
    direct_L2 = -1.0_8



    ! set random seed
    !L   = 2.0_8*pi*sqrt(3.0_8)/real(n_ranks,kind=kind_particle )

    L   = 4.0_8*pi
    dx  = L/real( n_ranks, kind=kind_particle )

    do ip=1, np

          p(ip)%label       = my_rank * (tnp / n_ranks) + ip
          p(ip)%data%q      = -1.0_8/real( tnp, kind=kind_particle )!(-1.0_8 + 2.0_8*MOD(p(ip)%label,2_kind_particle))
          p(ip)%data%m      = 1.0_8/real( tnp, kind=kind_particle )

          call random(p(ip)%x)

          p(ip)%x(1)         = (real( ip, kind=kind_particle ) - 1.0_8 )/( real( tnp, kind=kind_particle ) - 1.0_8)
          p(ip)%x(2)         = 0.0_8!1e-5*p(ip)%x(2)
          p(ip)%x(3)         = 0.0_8

          if(p(ip)%data%q .gt. 0.0)  p(ip)%data%m = p(ip)%data%m * 100.0_8

           call random_gaussian_flux(p(ip)%data%v(1),vtherm)
!          p(ip)%data%v(1) = vtherm
!          if (mod(ip,2) .eq. 0) p(ip)%data%v(1) = -vtherm


          p(ip)%data%v(2)   = 0.0_8
          p(ip)%data%v(3)   = 0.0_8


          p(ip)%results%e   = 0.0_8
          p(ip)%results%pot = 0.0_8
          p(ip)%results%B   = 0.0_8
          p(ip)%results%A   = 0.0_8
          p(ip)%results%Jirr= 0.0_8
          p(ip)%results%J   = 0.0_8
          p(ip)%work        = 1.0_8

          next_label = tnp+1


    end do

    p(:)%x(1)   = p(:)%x(1)*dx + real( my_rank, kind= kind_particle )*dx
!    p(:)%x(2)   = p(:)%x(2)*L



    !!! perturbation

    p(:)%x(1)       = p(:)%x(1) +  x0*sin( half*p(:)%x(1) )
    p(:)%data%v(1)  = p(:)%data%v(1) +  v0*sin( half*p(:)%x(1) )

    next_label = tnp+1

  end subroutine init_landau_damping1d


  subroutine init_two_stream2d(p)
    use zufall, only: random_gauss_list
    implicit none

    type(t_particle), allocatable, intent(inout) :: p(:)
    integer(kind_particle) :: ip,jp,iter
    integer :: rc
    real(kind_particle) :: dummy,dx,L,nppd

    real(kind_particle),parameter :: utherm = 1e-2  !! thermal velocity - x
    real(kind_particle),parameter :: vtherm = 1e-2  !! thermal velocity - y
    real(kind_particle),parameter :: wtherm = 1e-2  !! thermal velocity - y
    real(kind_particle),parameter :: u0     = 1e-1  !! drift velocity - x
    real(kind_particle),parameter :: v0     = 1e-1  !! drift velocity - y
    real(kind_particle),parameter :: w0     = 1e-1  !! drift velocity - y
!    real(kind_particle),parameter :: x0     = 1e-1

    if(my_rank.eq.0) write(*,'(a)') " == [init] init particles "

    ! set initially number of local particles
    np = tnp / n_ranks
    if(my_rank.eq.(n_ranks-1)) np = np + MOD(tnp, int(n_ranks, kind=kind_particle))

    allocate(particles(np), stat=rc)
    if(rc.ne.0) write(*,*) " === particle allocation error!"

    allocate(direct_L2(np), stat=rc)
    if(rc.ne.0) write(*,*) " === direct_L2 allocation error!"
    direct_L2 = -1.0_8



    ! set random seed
    !L   = 2.0_8*pi*sqrt(3.0_8)/real(n_ranks,kind=kind_particle )

    nppd = int( sqrt( real(np,8) ), kind = kind_particle )

    L   = 4.0_8*pi
    dx  = L/real( n_ranks, kind=kind_particle )

    iter = 1

    do ip=1, nppd

        do jp =1, nppd
          p(iter)%label       = my_rank * (tnp / n_ranks) + iter
          p(iter)%data%q      = (-1.0_8 + 2.0_8*MOD(p(ip)%label,2_kind_particle))
          p(iter)%data%m      = 1.0_8

          p(iter)%x(1)         = (ip - 1.0_8 )/( nppd - 1.0_8)
          p(iter)%x(2)         = (jp - 1.0_8 )/( nppd - 1.0_8)
          p(iter)%x(3)         = 0.0_8

          if(p(iter)%data%q .gt. 0.0)  then

                p(iter)%data%m = p(iter)%data%m * 100.0_8
                p(iter)%data%v = 0.0_8

          else

!                call random_gauss_list( p(iter)%data%v(1),u0,utherm )
!                call random_gaussian_flux(p(iter)%data%v(1),utherm)
!                call random_gaussian_flux(p(iter)%data%v(2),vtherm)
          endif


          p(iter)%data%v(3)   = 0.0_8

          p(iter)%results%e   = 0.0_8
          p(iter)%results%pot = 0.0_8
          p(iter)%results%B   = 0.0_8
          p(iter)%results%A   = 0.0_8
          p(iter)%results%Jirr= 0.0_8
          p(iter)%results%J   = 0.0_8
          p(iter)%work        = 1.0_8

          iter       = iter + 1
          next_label = tnp+1

        end do
    end do

    p(:)%x(1)   = p(:)%x(1)*dx + real( my_rank, kind= kind_particle )*dx
    p(:)%x(2)   = p(:)%x(2)*L

    next_label = tnp+1

  end subroutine init_two_stream2d

  subroutine two_bodies_grav(p)
  implicit none

    type(t_particle), allocatable, intent(inout) :: p(:)
    integer(kind_particle) :: ip
    integer :: rc
    real*8 :: dummy

    real,parameter :: Tkb=1e-9

    if(my_rank.eq.0) write(*,'(a)') " == [init] init particles "

    ! set initially number of local particles
    np = tnp / n_ranks
    if(my_rank.eq.(n_ranks-1)) np = np + MOD(tnp, int(n_ranks, kind=kind_particle))

    allocate(particles(np), stat=rc)
    if(rc.ne.0) write(*,*) " === particle allocation error!"

    allocate(direct_L2(np), stat=rc)
    if(rc.ne.0) write(*,*) " === direct_L2 allocation error!"
    direct_L2 = -1.0_8


    ! set random seed
    dummy = par_rand(my_rank)

    ! setup random qubic particle cloud
    if (tnp == 2) then

        ip                = 1
        p(ip)%x           = 0.0_8
        p(ip)%data%v      = 0.0_8
        p(ip)%data%m      = 1.0_8

        ip                = 2
        p(ip)%x           = 0.0_8
        p(ip)%x(1)        = -1.0_8
        p(ip)%data%v      = 0.0_8
        p(ip)%data%v(2)   = 0.01_8
        p(ip)%data%m      = 1.0_8






        do ip = 1,np
            p(ip)%label       = my_rank * (tnp / n_ranks) + ip
            p(ip)%data%q      = 1.0_8


            p(ip)%results%e   = 0.0_8
            p(ip)%results%pot = 0.0_8
            p(ip)%results%A   = 0.0_8
            p(ip)%results%B   = 0.0_8
            p(ip)%results%J   = 0.0_8
            p(ip)%results%Jirr= 0.0_8
            p(ip)%work        = 1.0_8
        enddo

    endif

    next_label = tnp+1

  end subroutine

  subroutine init_particles3D(p)
    implicit none

    type(t_particle), allocatable, intent(inout) :: p(:)
    integer(kind_particle) :: ip
    integer :: rc
    real*8 :: dummy

    real,parameter :: Tkb=1e-9

    if(my_rank.eq.0) write(*,'(a)') " == [init] init particles "

    ! set initially number of local particles
    np = tnp / n_ranks
    if(my_rank.eq.(n_ranks-1)) np = np + MOD(tnp, int(n_ranks, kind=kind_particle))

    allocate(particles(np), stat=rc)
    if(rc.ne.0) write(*,*) " === particle allocation error!"

    allocate(direct_L2(np), stat=rc)
    if(rc.ne.0) write(*,*) " === direct_L2 allocation error!"
    direct_L2 = -1.0_8


    ! set random seed
    dummy = par_rand(my_rank)

    ! setup random qubic particle cloud
    do ip=1, np

      p(ip)%label       = my_rank * (tnp / n_ranks) + ip
      p(ip)%data%q      = 1.0_8!(-1.0_8 + 2.0_8*MOD(p(ip)%label,2_kind_particle))
      p(ip)%data%m      = 1.0_8
      if(p(ip)%data%q .gt. 0.0) p(ip)%data%m = p(ip)%data%m * 100.0_8

      call random(p(ip)%x)

      call random_gauss(p(ip)%data%v)
!      p(ip)%data%v      = 1.0_8 !p(ip)%data%v * sqrt(Tkb * p(ip)%data%m) / p(ip)%data%m
!      p(ip)%x(1)           = 0.0_8
!      p(ip)%x(3)           = 0.0_8
!      p(ip)%data%v(2)      = 0.0_8
!      p(ip)%data%v(3)      = 0.0_8

      p(ip)%results%e   = 0.0_8
      p(ip)%results%pot = 0.0_8
      p(ip)%results%A   = 0.0_8
      p(ip)%results%B   = 0.0_8
      p(ip)%results%J   = 0.0_8
      p(ip)%results%Jirr= 0.0_8
      p(ip)%work        = 1.0_8

    end do

    next_label = tnp+1

  end subroutine init_particles3D


  subroutine init_particles_ring_2D(p)
    !use pepca_units, only: pi
    implicit none


    type(t_particle), allocatable, intent(inout) :: p(:)
    integer(kind_particle) :: ip
    integer :: rc
    real(kind_particle)                 :: dummy, deltaTheta,theta,eps2,A1,A2,x2,xy,tmp1,tmp2,tmp3
    !real(kind_particle), allocatable    :: theta(:)

    real,parameter :: Tkb=1e-9

    if(my_rank.eq.0) write(*,'(a)') " == [init] init particles "

    ! set initially number of local particles
    np = tnp / n_ranks
    if(my_rank.eq.(n_ranks-1)) np = np + MOD(tnp, int(n_ranks, kind=kind_particle))

    allocate(particles(np), stat=rc)
    if(rc.ne.0) write(*,*) " === particle allocation error!"

    allocate(direct_L2(np), stat=rc)
    if(rc.ne.0) write(*,*) " === direct_L2 allocation error!"
    direct_L2 = -1.0_8

    ! set random seed
    dummy = par_rand(my_rank)

    !call random_seed()
    !call random_number(theta)
    deltaTheta  = 2*pi/real(n_ranks,8)
    eps2 = 1e-4

    A1  = eps2*log( 1.0 + 1.0/eps2  ) -1.0_8
    A2  = eps2*log( 1.0 + 1.0/eps2  ) - log(eps2)

    if ( tnp .eq. 3 ) then


        do ip = 1,2
          p(ip)%label     = ip
          p(ip)%data%q    = 1
          p(ip)%data%m    = 100
          p(ip)%x(3)      = 0.0_8
          p(ip)%data%v(3) = 0.0_8
        enddo

        p(3)%label        = 3
        p(3)%data%q       = -2
        p(3)%data%m       = 1.0_8

        p(1)%x(1)         = -1.0_8
        p(1)%x(2)         = 0.0_8

        p(2)%x(1)         = 1.0_8
        p(2)%x(2)         = 0.0_8

        p(3)%x(1)         = 0.0_8
        p(3)%x(2)         = 1.0_8
        p(3)%x(3)         = 0.0_8

        p(1)%data%v(1)    = ( A2 - 4*A1 )/( .5*A2 - A1 )
        p(2)%data%v(1)    = p(1)%data%v(1)
        p(3)%data%v(1)    = 1.0_8

        p(1)%data%v(2)    = 0.0_8
        p(2)%data%v(2)    = 0.0_8
        p(3)%data%v(2)    = 0.0_8

        p(3)%data%v(3)    = 0.0_8

        do ip = 1,3
            tmp1 = .5*A2*p(2)%data%v(ip)*p(2)%data%q + .5*p(3)%data%q*p(3)%data%v(ip)*A2
            tmp1 = tmp1 - p(2)%data%q*A1*p(2)%data%v(ip)*( p(1)%x(1) - p(2)%x(1) )**2 - p(3)%data%q*A1*p(3)%data%v(ip)*( p(1)%x(1) - p(3)%x(1) )**2
!            tmp1 = tmp1 - A1*p(2)%data%v(ip) + 2*A1*p(3)%data%v(ip)

            tmp2 = .5*A2*p(1)%data%v(ip)*p(1)%data%q + .5*p(3)%data%q*p(3)%data%v(ip)*A2
            tmp2 = tmp2 - p(1)%data%q*A1*p(1)%data%v(ip)*( p(1)%x(1) - p(2)%x(1) )**2 - p(3)%data%q*A1*p(3)%data%v(ip)*( p(2)%x(1) - p(3)%x(1) )**2

            tmp3 = .5*A2*p(1)%data%v(ip)*p(1)%data%q + .5*p(2)%data%q*p(2)%data%v(ip)*A2
            tmp3 = tmp3 - p(1)%data%q*A1*p(1)%data%v(ip)*( p(1)%x(1) - p(3)%x(1) )**2 - p(2)%data%q*A1*p(2)%data%v(ip)*( p(2)%x(1) - p(3)%x(1) )**2


            write(*,*) tmp1, tmp2, tmp3
        enddo

    else
        do ip=1, np

          p(ip)%label       = my_rank * (tnp / n_ranks) + ip
          p(ip)%data%q      = (-1.0_8 + 2.0_8*MOD(p(ip)%label,2_kind_particle))
          p(ip)%data%m      = 1.0_8
          if(p(ip)%data%q .gt. 0.0) p(ip)%data%m = p(ip)%data%m * 100.0_8

          theta              = deltaTheta*( my_rank + ( real(ip,8) -1.0 )/real(np,8) )
    !      theta = 0._8
    !      call init_random_seed()
    !      call random_number(theta)
          p(ip)%x(1)         = cos(theta)
          p(ip)%x(2)         = sin(theta)
          p(ip)%x(3)         = 0.0_8


          p(ip)%data%v(1)      = -sin(theta)!0.05!
          p(ip)%data%v(2)      = cos(theta)
          p(ip)%data%v(3)      = 0.0_8

          p(ip)%results%e   = 0.0_8
          p(ip)%results%pot = 0.0_8
          p(ip)%results%B   = 0.0_8
          p(ip)%results%A   = 0.0_8
          p(ip)%results%Jirr= 0.0_8
          p(ip)%results%J   = 0.0_8
          p(ip)%work        = 1.0_8

        end do

    endif

    next_label = tnp+1

  end subroutine init_particles_ring_2D


subroutine init_particles_wire(p)
    !use pepca_units, only: pi
    implicit none


    type(t_particle), allocatable, intent(inout) :: p(:)
    integer(kind_particle) :: ip,jp
    integer :: rc
    real(kind_particle)                 :: dummy,L,deltaL

    real,parameter :: Tkb=1e-9

    if(my_rank.eq.0) write(*,'(a)') " == [init] init particles "

    ! set initially number of local particles
    np = tnp / n_ranks
    if(my_rank.eq.(n_ranks-1)) np = np + MOD(tnp, int(n_ranks, kind=kind_particle))

    allocate(particles(np), stat=rc)
    if(rc.ne.0) write(*,*) " === particle allocation error!"

    allocate(direct_L2(np), stat=rc)
    if(rc.ne.0) write(*,*) " === direct_L2 allocation error!"
    direct_L2 = -1.0_8


    ! set random seed
    dummy   = par_rand(my_rank)
    L       = 1.0/real(n_ranks,8)
    deltaL  = L/real( np, 8 )


    ! setup random qubic particle cloud
    do ip=1, np

      p(ip)%label       = my_rank * (tnp / n_ranks) + ip
      p(ip)%data%q      = (-1.0_8 + 2.0_8*MOD(p(ip)%label,2_kind_particle))
      p(ip)%data%m      = 1.0_8
      if(p(ip)%data%q .gt. 0.0) p(ip)%data%m = p(ip)%data%m * 100.0_8


      p(ip)%x(1)         =  ( ip - 1.0_8 )*deltaL!(-1.0_8 + 2.0_8*MOD(p(ip)%label,2_kind_particle))*( ip - MOD(p(ip)%label,2_kind_particle) )*deltaL
      p(ip)%x(2)         =  0.0_8
      p(ip)%x(3)         =  0.0_8

!      call random_number(p(ip)%data%v)
      p(ip)%data%v(1)      = 0.05! 2*pi*cos(2*pi*theta(ip))
      p(ip)%data%v(2)      = 0.0_8!-2*pi*sin(2*pi*theta(ip))
      p(ip)%data%v(3)      = 0.0_8

      p(ip)%results%e   = 0.0_8
      p(ip)%results%pot = 0.0_8
      p(ip)%results%B   = 0.0_8
      p(ip)%results%A   = 0.0_8
      p(ip)%results%Jirr= 0.0_8
      p(ip)%results%J   = 0.0_8
      p(ip)%work        = 1.0_8

!      write(*,*) ip, p(ip)%x(1),p(ip)%x(2)

    end do

    next_label = tnp+1

  end subroutine init_particles_wire


  subroutine init_particles_two_wires(p)
    !use pepca_units, only: pi
    implicit none


    type(t_particle), allocatable, intent(inout) :: p(:)
    integer(kind_particle) :: ip,jp
    integer :: rc
    real(kind_particle)                 :: dummy,L,deltaL,solenoidal_radius

    real,parameter :: Tkb=1e-9

    if(my_rank.eq.0) write(*,'(a)') " == [init] init particles "

    ! set initially number of local particles
    np = tnp / n_ranks
    if(my_rank.eq.(n_ranks-1)) np = np + MOD(tnp, int(n_ranks, kind=kind_particle))

    allocate(particles(np), stat=rc)
    if(rc.ne.0) write(*,*) " === particle allocation error!"

    allocate(direct_L2(np), stat=rc)
    if(rc.ne.0) write(*,*) " === direct_L2 allocation error!"
    direct_L2 = -1.0_8


    ! set random seed
    dummy             = par_rand(my_rank)
    L                 = 100.0/real(n_ranks,8)
    deltaL            = 2*L/real( np, 8 )
    solenoidal_radius = 2*1e-18


    ! setup random qubic particle cloud
    do ip=1, np

      p(ip)%label       = my_rank * (tnp / n_ranks) + ip
      p(ip)%data%q      = (-1.0_8 + 2.0_8*MOD(p(ip)%label,2_kind_particle))
      p(ip)%data%m      = 1.0_8
      if(p(ip)%data%q .gt. 0.0) p(ip)%data%m = p(ip)%data%m * 100.0_8


      p(ip)%x(1)         =  -L*( 1 + my_rank  ) + deltaL*real(ip - MOD(p(ip)%label,2_kind_particle) ,8)
      p(ip)%x(2)         =  solenoidal_radius*( MOD(p(ip)%label,2_kind_particle) )
      p(ip)%x(3)         =  0.0_8

!      call random_number(p(ip)%data%v)
      p(ip)%data%v(1)      = 0.05
      p(ip)%data%v(2)      = 0.0_8
      p(ip)%data%v(3)      = 0.0_8

      p(ip)%results%e   = 0.0_8
      p(ip)%results%pot = 0.0_8
      p(ip)%results%B   = 0.0_8
      p(ip)%results%A   = 0.0_8
      p(ip)%results%Jirr= 0.0_8
      p(ip)%results%J   = 0.0_8
      p(ip)%work        = 1.0_8

    end do


    next_label = tnp+1

  end subroutine init_particles_two_wires

  subroutine init_particles_stripe(p)
    !use pepca_units, only: pi
    implicit none


    type(t_particle), allocatable, intent(inout) :: p(:)
    integer(kind_particle) :: ip,jp
    integer :: rc
    real(kind_particle)                 :: dummy,L,wire_section

    real,parameter :: Tkb=1e-9

    if(my_rank.eq.0) write(*,'(a)') " == [init] init particles "

    ! set initially number of local particles
    np = tnp / n_ranks
    if(my_rank.eq.(n_ranks-1)) np = np + MOD(tnp, int(n_ranks, kind=kind_particle))

    allocate(particles(np), stat=rc)
    if(rc.ne.0) write(*,*) " === particle allocation error!"

    allocate(direct_L2(np), stat=rc)
    if(rc.ne.0) write(*,*) " === direct_L2 allocation error!"
    direct_L2 = -1.0_8


    ! set random seed
    dummy            = par_rand(my_rank)
    L                = 0.8/real(n_ranks,8)
    wire_section     = 0.0008

    ! setup random qubic particle cloud
    do ip=1, np

      p(ip)%label       = my_rank * (tnp / n_ranks) + ip
      p(ip)%data%q      = 1.0_8!(-1.0_8 + 2.0_8*MOD(p(ip)%label,2_kind_particle))
      p(ip)%data%m      = 1.0_8
      if(p(ip)%data%q .gt. 0.0) p(ip)%data%m = p(ip)%data%m * 100.0_8

      call random_number(p(ip)%x)

      p(ip)%x(1)         =  2*L*p(ip)%x(1) + L*( real(my_rank,8) - real(n_ranks,8) )
      p(ip)%x(2)         =  2*wire_section*p(ip)%x(2) - wire_section
      p(ip)%x(3)         =  0.0_8

!      call random_number(p(ip)%data%v)
      p(ip)%data%v(1)      = 0.05! 2*pi*cos(2*pi*theta(ip))
      p(ip)%data%v(2)      = 0.0_8!-2*pi*sin(2*pi*theta(ip))
      p(ip)%data%v(3)      = 0.0_8

      p(ip)%results%e   = 0.0_8
      p(ip)%results%pot = 0.0_8
      p(ip)%results%B   = 0.0_8
      p(ip)%results%A   = 0.0_8
      p(ip)%results%Jirr= 0.0_8
      p(ip)%results%J   = 0.0_8
      p(ip)%work        = 1.0_8

!      write(*,*) ip, p(ip)%x(1),p(ip)%x(2)

    end do

    next_label = tnp+1

  end subroutine init_particles_stripe


  subroutine init_particles_disc(p)
!    use module_globals, only: pi
    implicit none


    type(t_particle), allocatable, intent(inout) :: p(:)
    integer(kind_particle) :: ip,jp
    integer :: rc
    real(kind_particle)                 :: dummy,r,theta,rpp

    real,parameter :: Tkb=1e-9

    if(my_rank.eq.0) write(*,'(a)') " == [init] init particles "

    ! set initially number of local particles
    np = tnp / n_ranks
    if(my_rank.eq.(n_ranks-1)) np = np + MOD(tnp, int(n_ranks, kind=kind_particle))

    allocate(particles(np), stat=rc)
    if(rc.ne.0) write(*,*) " === particle allocation error!"

    allocate(direct_L2(np), stat=rc)
    if(rc.ne.0) write(*,*) " === direct_L2 allocation error!"
    direct_L2 = -1.0_8

!
!    !dr               = 2.0_8/tnp !real(n_ranks,8)
!    dr               = 1.0_8/real(n_ranks,8)
!    dtheta           = 0.5_8/(tnp*pi) !real(n_ranks*pi,8)
!
!    ! setup random qubic particle cloud
!    do ip=1, np
!
!      p(ip)%label       = my_rank * (tnp / n_ranks) + ip
!      p(ip)%data%q      = 1.0_8!(-1.0_8 + 2.0_8*MOD(p(ip)%label,2_kind_particle))
!      p(ip)%data%m      = 1.0_8
!      if(p(ip)%data%q .gt. 0.0) p(ip)%data%m = p(ip)%data%m * 100.0_8
!
!      call random_number(p(ip)%x)
!
!      !r                  = 0.5_8*( 1.0_8 + cos(pi*ip*dr) )
!      r                  = real(my_rank + 1,8)*dr*p(ip)%x(1) + real(my_rank,8)*dr
!      theta              = 2*pi*p(ip)%x(2)

    ! set random seed
    dummy            = par_rand(my_rank)
!    L                = 5.0/real(n_ranks,8)
    rpp              = 1/real(n_ranks,8)

    ! setup random qubic particle cloud
    do ip=1, np

      p(ip)%label       = my_rank * (tnp / n_ranks) + ip
      p(ip)%data%q      = 1.0_8!(-1.0_8 + 2.0_8*MOD(p(ip)%label,2_kind_particle))
      p(ip)%data%m      = 1.0_8
      if(p(ip)%data%q .gt. 0.0) p(ip)%data%m = p(ip)%data%m * 100.0_8

      call random_number(p(ip)%x)

      r                  = real(my_rank + 1,8)*rpp*p(ip)%x(1) + real(my_rank,8)*rpp
      r                  = sqrt(r)
      theta              = 2*pi*p(ip)%x(2)

      p(ip)%x(1)         =  r*cos(theta)
      p(ip)%x(2)         =  r*sin(theta)

      p(ip)%x(3)         =  0.0_8

!      call random_number(p(ip)%data%v)
      p(ip)%data%v(1)      =  -r*sin(theta)
      p(ip)%data%v(2)      =  r*cos(theta)
      p(ip)%data%v(3)      =  0.0_8

      p(ip)%results%e   = 0.0_8
      p(ip)%results%pot = 0.0_8
      p(ip)%results%B   = 0.0_8
      p(ip)%results%A   = 0.0_8
      p(ip)%results%Jirr= 0.0_8
      p(ip)%results%J   = 0.0_8
      p(ip)%work        = 1.0_8

!      write(*,*) ip, p(ip)%x(1),p(ip)%x(2)

    end do

    next_label = tnp+1

  end subroutine init_particles_disc


  subroutine init_particles_solenoid(p)
!    use module_globals, only: pi
    implicit none


    type(t_particle), allocatable, intent(inout) :: p(:)
    integer(kind_particle) :: ip,jp
    integer :: rc
    real(kind_particle)                 :: dummy,r,theta,deltaTheta

    real,parameter :: Tkb=1e-9

    if(my_rank.eq.0) write(*,'(a)') " == [init] init particles "

    ! set initially number of local particles
    np = tnp / n_ranks
    if(my_rank.eq.(n_ranks-1)) np = np + MOD(tnp, int(n_ranks, kind=kind_particle))

    allocate(particles(np), stat=rc)
    if(rc.ne.0) write(*,*) " === particle allocation error!"

    allocate(direct_L2(np), stat=rc)
    if(rc.ne.0) write(*,*) " === direct_L2 allocation error!"
    direct_L2 = -1.0_8

    deltaTheta  = 2*pi/real(n_ranks,8)
    r           = 1.0_8

    do ip = 1,np

        p(ip)%label       = my_rank * (tnp / n_ranks) + ip
        p(ip)%data%q      = 1.0_8!(-1.0_8 + 2.0_8*MOD(p(ip)%label,2_kind_particle))
        p(ip)%data%m      = 1.0_8
        if(p(ip)%data%q .gt. 0.0) p(ip)%data%m = p(ip)%data%m * 100.0_8

        theta = deltaTheta*( my_rank + ( real(ip,8) -1.0 )/real(np,8) )
        p(ip)%x(1)  = real( p(ip)%label - 1, kind_particle )/real( tnp , kind_particle )
        p(ip)%x(2)  = r*cos(theta)
        p(ip)%x(3)  = r*sin(theta)

        call random_number(p(ip)%data%v)


        p(ip)%results%e   = 0.0_8
        p(ip)%results%pot = 0.0_8
        p(ip)%results%B   = 0.0_8
        p(ip)%results%A   = 0.0_8
        p(ip)%results%Jirr= 0.0_8
        p(ip)%results%J   = 0.0_8
        p(ip)%work        = 1


    enddo

    next_label = tnp+1

  end subroutine init_particles_solenoid




  subroutine estimation_eps2(p)
  type(t_particle), allocatable, intent(in) :: p(:)
  integer(kind_particle) :: ip,jp
  real(kind_particle)    :: rv

  rv = 0.0_8
!  pi =  2.0_8*acos(0.0_8)

  do ip = 1, np-1
    do jp = ip +1, np
      rv = rv + 1.0_8/sqrt( ( p(ip)%x(1) - p(jp)%x(1) )**2 + ( p(ip)%x(2) - p(jp)%x(2) )**2 + ( p(ip)%x(3) - p(jp)%x(3) )**2  )
    enddo
  enddo

  rv  = 2*rv/real( np, 8 )
  write(*,*) " ==== Estimated smooth paramater eps2 ", ( .0625*3*pi/rv )**2


  end subroutine


  subroutine reallocate_particles(list, oldn, newn)
    implicit none

    type(t_particle), allocatable, intent(inout) :: list(:)
    integer(kind_particle), intent(in) :: oldn, newn

    type(t_particle) :: tmp_p(oldn)

    tmp_p(1:oldn) = list(1:oldn)
    deallocate(list)
    allocate(list(newn))
    list(1:oldn) = tmp_p

  end subroutine


  subroutine random(array)
    implicit none
    real*8 :: array(:)
    integer :: i

    do i = 1,size(array)
       array(i) = par_rand()
    end do

  end subroutine random


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>
  !> portable random number generator, see numerical recipes
  !> check for the random numbers:
  !> the first numbers should be 0.2853809, 0.2533582 and 0.0934685
  !> the parameter iseed is optional
  !>
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function par_rand(iseed)
    implicit none
    real :: par_rand
    integer, intent(in), optional :: iseed

    integer, parameter :: IM1  = 2147483563
    integer, parameter :: IM2  = 2147483399
    real,    parameter :: AM   = 1.0/IM1
    integer, parameter :: IMM1 = IM1-1
    integer, parameter :: IA1  = 40014
    integer, parameter :: IA2  = 40692
    integer, parameter :: IQ1  = 53668
    integer, parameter :: IQ2  = 52774
    integer, parameter :: IR1  = 12211
    integer, parameter :: IR2  = 3791
    integer, parameter :: NTAB = 32
    integer, parameter :: NDIV = 1+IMM1/NTAB
    real,    parameter :: eps_ = 1.2e-7 ! epsilon(eps_)
    real,    parameter :: RNMX = 1.0 - eps_

    integer :: j, k
    integer, volatile, save :: idum  = -1
    integer, volatile, save :: idum2 =  123456789
    integer, volatile, save :: iy    =  0
    integer, volatile, save :: iv(NTAB)


    if (idum <=0 .or. present(iseed)) then
       if (present(iseed)) then
          idum = iseed
       else
          if (-idum < 1) then
             idum = 1
          else
             idum = -idum
          endif
       endif

       idum2 = idum

       do j = NTAB+7,0,-1
          k = idum/IQ1
          idum = IA1 * (idum-k*IQ1) - k*IR1
          if (idum < 0 ) idum = idum + IM1

          if (j<NTAB) iv(j+1) = idum

       end do
       iy = iv(1)
    end if

    k = idum/IQ1
    idum = IA1 * (idum-k*IQ1) - k*IR1
    if (idum < 0) idum = idum + IM1

    k = idum2/IQ2
    idum2 = IA2 * (idum2-k*IQ2) - k*IR2
    if (idum2 < 0) idum2 = idum2 + IM2

    j = iy/NDIV + 1
    iy = iv(j)-idum2
    iv(j) = idum

    if (iy < 1) iy = iy + IMM1
    par_rand = AM*iy
    if (par_rand > RNMX) par_rand = RNMX

  end function par_rand

  real*8 function get_time()
      implicit none
      include 'mpif.h'

      get_time = MPI_WTIME()

  end function get_time

  subroutine random_gauss(list)
    implicit none
    real*8, intent(inout) :: list(:)
    real*8  :: v(2), r, p
    integer :: n, i

!    pi = 2.0_8*acos(0.0_8)
    n  = size(list)

    do i=1, n, 2

       call random(v)

       r = sqrt(-2.0_8 * log(v(1)))
       p = 2.0_8*pi*v(2)

       list(i)                = r * sin(p)
       if((i+1)<=n) list(i+1) = r * cos(p)

    end do

  end subroutine


  subroutine test_particles()

        use module_pepc_types
        use module_directsum
        implicit none
        include 'mpif.h'

        integer(kind_particle), allocatable   :: tindx(:)
        real(kind_particle), allocatable      :: trnd(:)
        type(t_particle_results), allocatable :: trslt(:)
        integer(kind_particle)                :: tn, tn_global, ti
        integer                               :: rc
        real(kind_particle)                   :: Ex,ExTilde,Ey,EyTilde,E_norm_loc,E_norm_global,E_global,E_local
        real(kind_particle)                   :: Bz,BzTilde,B_norm,B_norm_global,B_local,B_global,B_norm_loc
        real(kind_particle)                   :: Ax,AxTilde,Ay,AyTilde,A_norm,A_norm_global,A_local,A_global,A_norm_loc
        real(kind_particle)                   :: phi,phiTilde,phi_global,phi_norm_loc,phi_norm_global,ta,tb,phi_local
        real(kind_particle)                   :: Jirr_norm_global,Jirr_global,Jirr_local,Jirr_norm_loc, Jx,JxTilde,Jy,JyTilde
        real(kind_particle)                   :: J_norm_global,J_global,J_local,J_norm_loc,devE,devPhi, Jirrx,JirrxTilde,Jirry,JirryTilde,Jirrz,JirrzTilde
        real(kind_particle)                   :: F_dar_loc,F_dar_glo,F_dar_den_loc,F_dar_den_glo,F_el_loc,F_el_glo,F_el_den_loc,F_el_den_glo,vx,vy,vz,x,y,z,m,q
        real(kind_particle)                   :: vcrossB(3),vcrossB_tilde(3),Az,AzTilde,Bx,BxTilde,By,ByTilde,Ez,EzTilde,Jz,JzTilde


        ta = get_time()


        E_local       = 0.0_8
        E_norm_loc    = 0.0_8
        E_global      = 0.0_8
        E_norm_global = 0.0_8

        phi_local     = 0.0_8
        phi_norm_loc  = 0.0_8
        phi_global    = 0.0_8

        B_global      = 0.0_8
        B_norm_global = 0.0_8
        B_local       = 0.0_8
        B_norm_loc    = 0.0_8

        A_global      = 0.0_8
        A_norm_global = 0.0_8
        A_local       = 0.0_8
        A_norm_loc    = 0.0_8

        J_global      = 0.0_8
        J_norm_global = 0.0_8
        J_local       = 0.0_8
        J_norm_loc    = 0.0_8

        Jirr_global      = 0.0_8
        Jirr_norm_global = 0.0_8
        Jirr_local       = 0.0_8
        Jirr_norm_loc    = 0.0_8

        F_dar_loc           = 0.0_8
        F_dar_glo           = 0.0_8
        F_dar_den_loc       = 0.0_8
        F_dar_den_glo       = 0.0_8

        F_el_loc           = 0.0_8
        F_el_glo           = 0.0_8
        F_el_den_loc       = 0.0_8
        F_el_den_glo       = 0.0_8

        tn = np!particle_direct / n_ranks
        if(my_rank.eq.(n_ranks-1)) tn = tn + MOD(particle_direct, n_ranks)

        allocate(tindx(tn), trnd(tn), trslt(tn))

        call random(trnd)

        tindx(1:tn) = int(trnd(1:tn) * (np-1)) + 1

        call directforce(particles, tindx, tn, trslt, MPI_COMM_WORLD)

        do ti = 1, tn

          vx          = particles(tindx(ti))%data%v(1)
          vy          = particles(tindx(ti))%data%v(2)
          vz          = particles(tindx(ti))%data%v(3)

          x           = particles(tindx(ti))%x(1)
          y           = particles(tindx(ti))%x(2)
          z           = particles(tindx(ti))%x(3)

          m           = particles(tindx(ti))%data%m
          q           = particles(tindx(ti))%data%q

          Ex          = trslt(ti)%e(1)
          Ey          = trslt(ti)%e(2)
          Ez          = trslt(ti)%e(3)

          ExTilde     = particles(tindx(ti))%results%e(1)
          EyTilde     = particles(tindx(ti))%results%e(2)
          EzTilde     = particles(tindx(ti))%results%e(3)

          Ax          = trslt(ti)%A(1)
          Ay          = trslt(ti)%A(2)
          Az          = trslt(ti)%A(3)

          AxTilde     = particles(tindx(ti))%results%A(1)
          AyTilde     = particles(tindx(ti))%results%A(2)
          AzTilde     = particles(tindx(ti))%results%A(3)


          JirrxTilde  = particles(tindx(ti))%results%Jirr(1)
          JirryTilde  = particles(tindx(ti))%results%Jirr(2)
          JirrzTilde  = particles(tindx(ti))%results%Jirr(3)

          Jirrx       = trslt(ti)%Jirr(1)
          Jirry       = trslt(ti)%Jirr(2)
          Jirrz       = trslt(ti)%Jirr(3)

          JxTilde     = particles(tindx(ti))%results%J(1)
          JyTilde     = particles(tindx(ti))%results%J(2)
          JzTilde     = particles(tindx(ti))%results%J(3)

          Jx          = trslt(ti)%J(1)
          Jy          = trslt(ti)%J(2)
          Jz          = trslt(ti)%J(3)

          Bx          = trslt(ti)%B(1)
          By          = trslt(ti)%B(2)
          Bz          = trslt(ti)%B(3)

          BxTilde     = particles(tindx(ti))%results%B(1)
          ByTilde     = particles(tindx(ti))%results%B(2)
          BzTilde     = particles(tindx(ti))%results%B(3)


          E_local     = E_local + ( ExTilde -  Ex )**2 + ( EyTilde - Ey )**2 + ( EzTilde - Ez )**2
          E_norm_loc  = E_norm_loc + Ex**2 + Ey**2 + EZ**2

          phi         = trslt(ti)%pot
          phiTilde    = particles(tindx(ti))%results%pot

          phi_local        = phi_local +  ( phi - phiTilde )**2
          phi_norm_loc     = phi_norm_loc +  phi**2


          A_local     = A_local + ( AxTilde -  Ax )**2 + ( AyTilde - Ay )**2 + ( AzTilde - Az )**2
          A_norm_loc  = A_norm_loc + Ax**2 + Ay**2  + Az**2

          B_local     = B_local + ( BxTilde -  Bx )**2 + ( ByTilde -  By )**2 + ( BzTilde -  Bz )**2
          B_norm_loc  = B_norm_loc + Bx**2  + By**2 + Bz**2


          J_local     = J_local + ( JxTilde -  Jx )**2 + ( JyTilde - Jy )**2 + ( JzTilde - Jz )**2
          J_norm_loc  = J_norm_loc + Jx**2 + Jy**2 + Jz**2

          Jirr_local     = Jirr_local + ( JirrxTilde -  Jirrx )**2 + ( JirryTilde - Jirry )**2 + ( JirrzTilde - Jirrz )**2
          Jirr_norm_loc  = Jirr_norm_loc + Jirrx**2 + Jirry**2 + Jirrz**2


          vcrossB(1:3)    = (/ vy*( Bz) - vx*(By) ,  vz*(Bx) - vx*(Bz) ,  vx*(By) - vy*(Bx)   /)
          vcrossB_tilde(1:3)    = (/ vy*( BzTilde) - vx*(ByTilde) ,  vz*(BxTilde) - vx*(BzTilde) ,  vx*(ByTilde) - vy*(BxTilde)   /)
          F_dar_loc       = F_dar_loc + (q)**2*( ( ExTilde -  Ex )**2 + ( EyTilde - Ey )**2 + ( EzTilde - Ez )**2 &
                        +   dot_product( vcrossB - vcrossB_tilde, vcrossB - vcrossB_tilde )                       &
                        + 2*dot_product( (/ ExTilde -  Ex, EyTilde -  Ey , EzTilde -  Ez/), vcrossB - vcrossB_tilde ) )

          F_dar_den_loc  = F_dar_den_loc + (q)**2*( Ex**2 +  Ey**2 + Ez**2 + dot_product( vcrossB, vcrossB )        &
                        + 2*dot_product( (/ Ex, Ey , Ez/), vcrossB ) )

          F_el_loc       = F_el_loc + (q)**2*( ( ExTilde -  Ex )**2 + ( EyTilde - Ey )**2  + ( EzTilde - Ez )**2 )

          F_el_den_loc  = F_el_den_loc + (q)**2*( Ex**2 +  Ey**2 +  Ez**2 )


        end do

        call MPI_ALLREDUCE(tn           , tn_global       , 1, MPI_KIND_PARTICLE, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(phi_local    , phi_global      , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(phi_norm_loc , phi_norm_global , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(E_local      , E_global        , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(E_norm_loc   , E_norm_global   , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(B_norm_loc   , B_norm_global   , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(B_local      , B_global        , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(A_local      , A_global        , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(A_norm_loc   , A_norm_global   , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(J_local      , J_global        , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(J_norm_loc   , J_norm_global   , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(Jirr_local   , Jirr_global     , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(Jirr_norm_loc, Jirr_norm_global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(F_dar_loc    , F_dar_glo       , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(F_dar_den_loc, F_dar_den_glo   , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(F_el_loc     , F_el_glo        , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(F_el_den_loc , F_el_den_glo    , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)


        devE                 = sqrt(E_global)/(tn_global-1.0_8)
        devPhi               = sqrt(phi_global)/(tn_global-1.0_8)

        phi_global           = sqrt(phi_global / phi_norm_global)
        E_global             = sqrt(E_global   / E_norm_global)
        A_global             = sqrt(A_global   / A_norm_global)
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
          write(*,'(a,es12.4)') " == [direct test] Relative error in El Pot        : ", phi_global
          write(*,'(a,es12.4)') " == [direct test] Relative error in E             : ", E_global
          write(*,'(a,es12.4)') " == [direct test] Relative error in A             : ", A_global
          write(*,'(a,es12.4)') " == [direct test] Relative error in B             : ", B_global
          write(*,'(a,es12.4)') " == [direct test] Relative error in J             : ", J_global
          write(*,'(a,es12.4)') " == [direct test] Relative error in Jirr          : ", Jirr_global
          write(*,'(a,es12.4)') " == [direct test] Relative error in F dar         : ", F_dar_glo
          write(*,'(a,es12.4)') " == [direct test] Relative error in F el          : ", F_el_glo
          write(*,'(a,es12.4)') " == [direct test] L2 error in E                   : ", devE
          write(*,'(a,es12.4)') " == [direct test] L2 error in El Pot              : ", devPhi
          write(*,'(a,es12.4)') " == [direct test] time in test [s]                : ", tb - ta

!          open(unit=rc,file="monopole2d.dat",form='formatted',status='unknown',access='append')
          open(unit=rc,file="monopole2d.dat",form='formatted',status='unknown',position='append')
!          open (unit=rc,file="quadrupole.dat",action="write",status="replace")
         ! write(rc,*), phi_global, E_global, A_global,B_global, J_global, Jirr_global, F_dar_glo, F_el_glo
          close(rc)

        end if

        deallocate(tindx)
        deallocate(trnd)
        deallocate(trslt)

      end subroutine test_particles

      subroutine energy(p,t)
      implicit none
      include 'mpif.h'
      type(t_particle), allocatable, intent(in) :: p(:)
      real(kind_particle)          , intent(in) :: t
      integer(kind_particle) :: ip,rc=201
      real(kind_particle)    :: upot,ukin,udar,gam,upot_loc,ukin_loc,udar_loc,v2

      upot     = 0.0_8
      ukin     = 0.0_8
      udar     = 0.0_8
      upot_loc = 0.0_8
      ukin_loc = 0.0_8
      udar_loc = 0.0_8

      do ip = 1,np
        v2       = ( p(ip)%data%v(1)**2 + p(ip)%data%v(2)**2 + p(ip)%data%v(3)**2 )
        gam      = 1.0_8/sqrt( 1.0_8 - v2 )
        upot_loc = upot_loc + .5*p(ip)%data%q*p(ip)%results%pot
        udar_loc = udar_loc + .5*p(ip)%data%q/gam*( p(ip)%results%A(1)*p(ip)%data%v(1) + p(ip)%results%A(2)*p(ip)%data%v(2) + p(ip)%results%A(3)*p(ip)%data%v(3) )
        ukin_loc = ukin_loc + .5*p(ip)%data%m*v2

      enddo

      call MPI_ALLREDUCE(upot_loc, upot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
      call MPI_ALLREDUCE(ukin_loc, ukin, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
      call MPI_ALLREDUCE(udar_loc, udar, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)

      if (root) then

!        open(unit=rc,file="energy_midpoint.dat",form='formatted',status='unknown',access='append')
!        open (unit=rc,file=trim("energy_")//trim(adjustl(ischeme))//".dat",action="write",status="append")
!        open(unit=rc,file=trim("energy_")//trim(adjustl(ischeme))//".dat",form='formatted',status='unknown',access='append')
        open(unit=rc,file=trim("energy_")//trim(adjustl(ischeme))//".dat",form='formatted',status='unknown',position='append')
        !write(*,*) " ==== Energy Conservation  ", upot,ukin,udar,upot+udar+ukin
        write(rc,*) t,ukin,upot,upot+ukin!udar,upot+udar+ukin
        close (rc )

      endif

      end subroutine


      subroutine electric_energy(field_grid,t)
      use encap
      implicit none
      include 'mpif.h'
      type(field_grid_t)           , intent(in) :: field_grid
      real(kind_particle)          , intent(in) :: t
      integer(kind_particle)                    :: ip,nl,rc=201
      real(kind_particle)                       :: EE,EE_loc

      EE_loc  = 0.0_8
      EE      = 0.0_8
      nl      = field_grid%nl

      do ip = 1,nl

        EE_loc  = EE_loc + field_grid%p(ip)%results%E(1)**2*field_grid%dx(1)

      enddo

      call MPI_ALLREDUCE(EE_loc, EE, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)

      if (root) then

!        open(unit=rc,file=trim("electric_energy_")//trim(adjustl(ischeme))//".dat",form='formatted',status='unknown',access='append')
        open(unit=rc,file=trim("electric_energy_")//trim(adjustl(ischeme))//".dat",form='formatted',status='unknown',position='append')
        write(rc,*) t,.50_8*EE
        close (rc )

      endif

      end subroutine electric_energy


      subroutine write_field(p)
        type(t_particle), allocatable, intent(in) :: p(:)
        integer(kind_particle) :: indexvec(7), tmp
        character (len=12), dimension(7) :: stringvec

        integer(kind_particle) :: ifield, ip
        indexvec = (/ 1, 2, 3, 4, 5, 6, 7/)
        stringvec = (/ "data/fields/xyz.dat", "data/fields/phi.dat", "data/fields/EEE.dat","data/fields/BBB.dat", "data/fields/rho.dat","data/fields/JJJ.dat","data/fields/AAA.dat"/)

        do ifield = 1,7

            open (unit=indexvec(ifield),file=stringvec(ifield),action="write",status="replace")

        enddo


        write(*,*) " ==  fields written "
        do tmp = 1, np

                do ip = 1, np

                    if ( p(ip)%label .eq. tmp ) then

                        write (1,*) p(ip)%x(1),p(ip)%x(2),p(ip)%x(3)
                        write (2,*) p(ip)%results%pot
                        write (3,*) p(ip)%results%E(1),p(ip)%results%E(2),p(ip)%results%E(3)
                        write (4,*) p(ip)%results%B(2),p(ip)%results%B(2),p(ip)%results%B(3)
                        write (5,*) p(ip)%results%J(1),p(ip)%results%J(2),p(ip)%results%J(3)
                        write (6,*) p(ip)%results%Jirr(1),p(ip)%results%Jirr(2),p(ip)%results%Jirr(3)
                        write (7,*) p(ip)%results%A(1),p(ip)%results%A(2),p(ip)%results%A(3)

                    endif

                enddo

        enddo


        do ifield = 1,7

            close ( unit=indexvec(ifield) )

        enddo

      end subroutine

      subroutine write_field_rectilinear_grid(p)
        use module_vtk
        implicit none
        type(t_particle), allocatable, intent(in) :: p(:)
        integer(kind_particle)  :: i
        integer, dimension(2,3) ::  nnp
        type(vtkfile_rectilinear_grid) :: vtk
        integer :: vtk_step
        real(kind_particle) :: time
        real(kind_particle) :: ta, tb,Ex(np),Ey(np),Ez(np),x(np),y(np),z(np)


        ta = get_time()
        time = dt * step

        if (step .eq. 0) then
          vtk_step = VTK_STEP_FIRST
        else if (step .eq. nt) then
          vtk_step = VTK_STEP_LAST
        else
          vtk_step = VTK_STEP_NORMAL
        endif

        x  = p(:)%x(1)
        y  = p(:)%x(2)
        z  = p(:)%x(3)

        Ex = p(:)%results%E(1)
        Ey = p(:)%results%E(2)
        Ez = p(:)%results%E(3)
    !    nnp(1) = np
    !    nnp(2) = np
    !    nnp(3) = 1


        !nnp = np

        call vtk%create_parallel("fields_on_grid", step, my_rank, n_ranks, time, vtk_step)
        call vtk%write_headers(nnp, nnp)
        call vtk%startcoordinates()
        call vtk%write_data_array("x_coordinate", x )
        call vtk%write_data_array("y_coordinate", y )
        call vtk%write_data_array("z_coordinate", 0 )

        call vtk%finishcoordinates()
        call vtk%startpointdata()

        call vtk%write_data_array("E", Ex,  Ey,  Ez )
              ! no point data here
        call vtk%finishpointdata()
        call vtk%startcelldata()
              ! no cell data here
        call vtk%finishcelldata()
        call vtk%write_final()
        call vtk%close()



        tb = get_time()

        if(root) write(*,'(a,es12.4)') " == [write particles] time in vtk output [s]      : ", tb - ta



           !call vtk%write_data_array(vectorname, vectorvalues(:,:,:,1), vectorvalues(:,:,:,2), vectorvalues(:,:,:,3))

      end subroutine


  subroutine test_fields(p,field_grid_tree,pepc_pars,flag)

        use encap
        use module_pepc_types
        use module_directsum
        use field_helper, only:compute_field
        implicit none
        include 'mpif.h'

        type(t_particle), allocatable, intent(in) :: p(:)
        type(field_grid_t), intent(in)            :: field_grid_tree
        type(pepc_pars_t), intent(in)             :: pepc_pars
        integer(kind_particle), intent(in)        :: flag
!        type(field_grid_t)                        :: field_grid_direct
        integer(kind_particle), allocatable       :: tindx(:)
        real(kind_particle), allocatable          :: trnd(:)
        type(t_particle_results), allocatable     :: trslt(:)
        integer(kind_particle)                    :: ip,rc,nl
        real(kind_particle)                       :: error_loc,error_glo,B_tree,B_an,sum_loc,sum_glo,sumB,sumE,a,r

        error_loc = 0.0_8
        error_glo = 0.0_8
        sum_loc   = 0.0_8
        sumB      = 0.0_8

        allocate(tindx(np), trnd(np), trslt(np))

        call random(trnd)

        tindx(1:np) = int(trnd(1:np) * (np-1)) + 1
        nl = field_grid_tree%nl

!        call directforce(particles, tindx, np, trslt, MPI_COMM_WORLD)
!        call compute_field(pepc_pars, field_grid_direct, particles)

        select case (flag)
            case (4)  !  Ring Current
                              do ip = 1,nl

                                    B_tree      = field_grid_tree%p(ip)%results%B(3)
                                    sumB        = sumB + B_tree**2
                                    B_an        = 2*field_grid_tree%p(ip)%results%E(2)*0.05 - 2*field_grid_tree%p(ip)%results%E(1)*0.05
                                    sumE        = sumE + B_an**2
                                    sum_loc     = sum_loc + B_an**2
!                                    write(*,*) field_grid_tree%p(ip)%data%v(1),B_tree

                                    error_loc   = error_loc + ( B_tree - B_an )**2

                              enddo
            case (5)  !  Wire of Lenght 2L
                              a                 = 5.0
                              do ip = 1,nl

                                    r           = abs( field_grid_tree%p(ip)%x(2) )

                                    if ( r .ne. 0 ) then
                                        B_tree      = field_grid_tree%p(ip)%results%B(3)
                                        sumB        = sumB + B_tree**2 ! field_grid_tree%p(ip)%results%J(1)*
                                        B_an        = 2*a/( r*sqrt( r**2 + a**2 ) )
                                        sumE        = sumE + B_an**2
                                        sum_loc     = sum_loc + B_an**2
                                        error_loc   = error_loc + ( B_tree - B_an )**2
                                    endif

                              enddo


            case default

                              do ip = 1,nl

                                    B_tree      = field_grid_tree%p(ip)%results%B(3)
                                    sumB        = sumB + B_tree**2
                                    B_an        = 2*field_grid_tree%p(ip)%results%E(2)*p(ip)%data%v(1) - 2*field_grid_tree%p(ip)%results%E(1)*p(ip)%data%v(2)
                                    sumE        = sumE + B_an**2
                                    sum_loc     = sum_loc + B_an**2

                                    error_loc   = error_loc + ( B_tree - B_an )**2

                              enddo

        end select


        call MPI_ALLREDUCE(error_loc, error_glo, 1, MPI_KIND_PARTICLE, MPI_SUM, MPI_COMM_WORLD, rc)
        call MPI_ALLREDUCE(sum_loc, sum_glo, 1, MPI_KIND_PARTICLE, MPI_SUM, MPI_COMM_WORLD, rc)

        if(root) then
            select case (flag)
                case(4) ! ring current

                    write(*,'(a,es12.4)') " == Relative error of B [ vs v x E ]      : ",sqrt(error_glo/sum_glo)
!                    write(*,'(a,es12.4)') " == Numerator                             : ",error_loc
!                    write(*,'(a,es12.4)') " == Denominator                           : ",sum_loc
                case(5) ! wire with length 2L
                    write(*,'(a,es12.4)') " == Relative error of B [ Biot-Savart ]   : ",sqrt(error_glo/sum_glo)
    !                write(*,'(a,es12.4)') " == Numerator                             : ",error_loc
    !                write(*,'(a,es12.4)') " == Denominator                           : ",sum_loc
                case default
                    write(*,'(a,es12.4)') " == Relative error of B [ vs v x E ]      : ",sqrt(error_glo/sum_glo)
    !                write(*,'(a,es12.4)') " == Numerator                             : ",error_loc
    !                write(*,'(a,es12.4)') " == Denominator                           : ",sum_loc
            end select
        endif


  end subroutine




      subroutine write_particle_ascii(p,itime)
        implicit none
        include 'mpif.h'

        type(t_particle), allocatable, intent(in) :: p(:)
        integer(kind_particle)       , intent(in) :: itime
        integer(kind_particle)                    :: ip,jp,rc
        integer(kind_particle),      parameter    :: part   = 100
        integer(kind_particle),      parameter    :: partp1 = 101
        integer(kind_particle),      parameter    :: sizeout = 7
        real(kind_particle)                       :: tmp(sizeout),EE,EE_glo


        character(255) str,file0,file1

        write( str, '(i10)' ) itime
        file0               = "langmuir/particle/particle_ascii_"
        file1               = "langmuir/particle/electric_energy_density"
!        file0               = file0//str
!        write(*,*) file0
!


        if (root) then
            open (unit=part,file=trim(file0)//trim(adjustl(str))//".dat",action="write",status="replace")
!            open(unit=partp1,file=trim(file1)//".dat",form='formatted',status='unknown',access='append')
            write(*,*) "==== Writing particle_ascii - tnp:         ",tnp
        endif


        EE_glo   = 0.0_8
        EE       = 0.0_8
        do ip = 1,np
            tmp(1) = p(ip)%x(1)
            tmp(2) = p(ip)%x(2)
            tmp(3) = p(ip)%data%v(1)
            tmp(4) = p(ip)%data%v(2)
!            tmp(3) = field_grid%p(ip)%x(3)
            tmp(5) = p(ip)%results%E(1)
!            tmp(6) = p(ip)%results%E(2)
            tmp(6) = p(ip)%results%pot
!            tmp(6) = field_grid%p(ip)%results%E(3)
!            tmp(7) = p(ip)%results%A(1)
!            tmp(8) = p(ip)%results%A(2)
!            tmp(8) = field_grid%p(ip)%results%A(3)
!            tmp(9) = field_grid%p(ip)%results%B(1)
!            tmp(10)= field_grid%p(ip)%results%B(2)
!            tmp(9) = p(ip)%results%B(3)
!            tmp(10) = p(ip)%results%J(1)
!            tmp(11) = p(ip)%results%J(2)
!            tmp(12) = p(ip)%x(1)**2 + p(ip)%x(2)**2
!            tmp(13) = p(ip)%data%v(1)**2 + p(ip)%data%v(2)**2
!            tmp(14) = p(ip)%data%q
!            tmp(15) = p(ip)%results%Jirr(1)
!            tmp(16) = p(ip)%results%Jirr(2)
!            EE     = EE + .5_8*p(ip)%results%E(1)*p(ip)%results%E(1)
            write (part,*) tmp(1:sizeout)

        enddo

!        call MPI_ALLREDUCE(EE, EE_glo, 1, MPI_KIND_PARTICLE, MPI_SUM, MPI_COMM_WORLD, rc)



        if (root) then
!            write (partp1,*) dt*itime,  EE_glo
            close ( unit=part )
!            close ( unit=partp1 )
            write(*,*) "==== Writing particle_ascii finished- tnp: ",tnp

        endif

      end subroutine write_particle_ascii




      subroutine write_field_on_grid_ascii(field_grid,itime)
        use module_vtk
        use encap
        implicit none

        type(field_grid_t), intent(in)         :: field_grid
        integer(kind_particle)    , intent(in) :: itime
        type(vtkfile_unstructured_grid)        :: vtk
        integer(kind_particle)                 :: ip,nl,jp
        integer(kind_particle),      parameter :: ifield = 10
        integer(kind_particle),      parameter :: sizeout = 11
        real(kind_particle)                    :: tmp(sizeout),r2
        character(255) str,file0,file1

        write( str, '(i10)' ) itime

        file0 = "langmuir/fields/field_on_grid_ascii_"
        if (root ) open (unit=ifield,file=trim(file0)//trim(adjustl(str))//".dat",action="write",status="replace")

        nl = field_grid%nl

        do ip = 1,nl
            tmp(1) = field_grid%p(ip)%x(1)
            tmp(2) = field_grid%p(ip)%x(2)
!            tmp(3) = field_grid%p(ip)%x(3)
            tmp(3) = field_grid%p(ip)%results%E(1)
            tmp(4) = field_grid%p(ip)%results%E(2)
!            tmp(6) = field_grid%p(ip)%results%E(3)
            tmp(5) = field_grid%p(ip)%results%A(1)
            tmp(6) = field_grid%p(ip)%results%A(2)
!            tmp(8) = field_grid%p(ip)%results%A(3)
!            tmp(9) = field_grid%p(ip)%results%B(1)
!            tmp(10)= field_grid%p(ip)%results%B(2)
            tmp(7)= field_grid%p(ip)%results%B(3)
            tmp(8)= field_grid%p(ip)%results%J(1)
            tmp(9)= field_grid%p(ip)%results%J(2)
            tmp(10)= field_grid%p(ip)%results%Jirr(1)
            tmp(11)= field_grid%p(ip)%results%Jirr(2)

            write (ifield,*) tmp

        enddo

        if (root )   close ( unit=ifield )

      end subroutine
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

        vect = 0.0_8

        write(*,*) "particles"
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
        call vtk%write_data_array("m", p(1:np)%data%m)
        call vtk%write_data_array("work", p(1:np)%work)
        call vtk%write_data_array("pelabel", p(1:np)%label)
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



      subroutine write_field_on_grid(field_grid)
        use module_vtk
        use encap
        implicit none

        !type(t_particle), allocatable, intent(in) :: p(:)
        type(field_grid_t), intent(in)  :: field_grid

        integer(kind_particle) :: i
        type(vtkfile_unstructured_grid) :: vtk
        integer :: vtk_step,nl
        real*8 :: time!,vect(np)
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

        nl = field_grid%nl

!        write(*,*) "particles"
        call vtk%create_parallel("fields_on_grid", step, my_rank, n_ranks, time, vtk_step)
        call vtk%write_headers(np, 0_kind_particle)
        call vtk%startpoints()
        call vtk%write_data_array("x", field_grid%p(1:nl)%x(1), field_grid%p(1:nl)%x(2), field_grid%p(1:nl)%x(3) )
        call vtk%finishpoints()
        call vtk%startpointdata()
        call vtk%write_data_array("v", field_grid%p(1:nl)%data%v(1), &
                                       field_grid%p(1:nl)%data%v(2), &
                                       field_grid%p(1:nl)%data%v(3) )

        call vtk%write_data_array("E", field_grid%p(1:nl)%results%e(1), &
                                       field_grid%p(1:nl)%results%e(2), &
                                       field_grid%p(1:nl)%results%e(3) )

        call vtk%write_data_array("A", field_grid%p(1:nl)%results%A(1), &
                                       field_grid%p(1:nl)%results%A(2), &
                                       field_grid%p(1:nl)%results%A(3) )


        call vtk%write_data_array("B",   field_grid%p(1:nl)%results%B(1),&
                                         field_grid%p(1:nl)%results%B(2),&
                                         field_grid%p(1:nl)%results%B(3) )

        call vtk%write_data_array("J",   field_grid%p(1:nl)%results%J(1),&
                                         field_grid%p(1:nl)%results%J(2),&
                                         field_grid%p(1:nl)%results%J(3) )

        call vtk%write_data_array("Jirr",field_grid%p(1:nl)%results%Jirr(1),&
                                         field_grid%p(1:nl)%results%Jirr(2),&
                                         field_grid%p(1:nl)%results%Jirr(3) )

        call vtk%write_data_array("phi", field_grid%p(1:nl)%results%pot)
        call vtk%write_data_array("q", field_grid%p(1:nl)%data%q)
        call vtk%write_data_array("m", field_grid%p(1:nl)%data%m)
!        call vtk%write_data_array("work", p(1:np)%work)
!        call vtk%write_data_array("pelabel", p(1:np)%label)
!        call vtk%write_data_array("local index", [(i,i=1,np)])
!        call vtk%write_data_array("processor", int(np, kind = 4), my_rank)
!        if(particle_test) call vtk%write_data_array("L2 error", direct_L2(1:np))
        call vtk%finishpointdata()
        call vtk%dont_write_cells()
        call vtk%write_final()
        call vtk%close()

        tb = get_time()

        if(root) write(*,'(a,es12.4)') " == [write particles] time in vtk output [s]      : ", tb - ta

      end subroutine write_field_on_grid


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
