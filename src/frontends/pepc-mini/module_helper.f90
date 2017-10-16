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
!

!!!!!!!!!!!!!!!!!!!!
!! helper module
!!!!!!!!!!!!!!!!!!!!

module helper

  use module_pepc_kinds
  use module_pepc_types
  use module_interaction_Specific_types
  implicit none

  ! MPI variables
  integer :: my_rank, n_ranks
  logical :: root

  ! time variables
  real*8 :: dt
  integer :: step

  ! control variables
  integer :: nt               ! number of timesteps
  integer(kind_particle) :: tnp ! total number of particles
  integer(kind_particle) :: np  ! local number of particles
  logical :: particle_output  ! turn vtk output on/off
  logical :: domain_output    ! turn vtk output on/off
  logical :: particle_filter  ! filter particles leaving simulation domain
  logical :: particle_probe   ! turn probin on/off
  logical :: particle_test    ! turn direct summation on/off
  integer :: diag_interval    ! number of timesteps between all diagnostics and IO

  real*8, parameter  :: plasma_width = 5.0_8      ! width of plasma block
  real*8, parameter  :: vessel_width = 5.0_8      ! distance between walls
  real*8             :: vessel_ez                 ! electrical field in vessel in z-direction

  integer(kind_default), parameter :: nprobes    = 200           ! number of probes
  real*8, parameter  :: min_probes = -vessel_width ! min z-pos of probes
  real*8, parameter  :: max_probes =  vessel_width ! max z-pos of probes

  integer, parameter :: particle_direct = 144 ! number of particle for direct summation

  ! current measurement
  real*8 :: current_q, current_I
  real*8, parameter :: current_emission = 1.0e4_8
  integer, parameter :: current_file_id = 14

  ! particle data (position, velocity, mass, charge)
  type(t_particle), allocatable :: particles(:)
  real*8, allocatable           :: direct_L2(:)
  integer(kind_particle)        :: next_label

  ! probe data
  ! first dimension: probe id
  ! second dimension: z-position, number of particles in bin, mass, charge, avg. pot
  real*8, allocatable           :: probes(:,:)


  contains

  subroutine set_parameter()

    use module_pepc
    use module_interaction_specific, only : theta2, eps2, force_law, include_far_field_if_periodic
    use module_mirror_boxes, only: periodicity, spatial_interaction_cutoff
    implicit none

    integer, parameter :: fid = 12
    character(255)     :: para_file
    logical            :: read_para_file

    namelist /pepcmini/ tnp, nt, dt, particle_output, domain_output, particle_filter, particle_test, &
      particle_probe, diag_interval, vessel_ez, periodicity, spatial_interaction_cutoff

    ! set default parameter values
    tnp             = 1441
    nt              = 20
    dt              = 1e-3
    particle_test   = .true.
    particle_output = .true.
    domain_output   = .true.
    particle_filter = .true.
    particle_probe  = .true.
    diag_interval   = 2

    vessel_ez = 100.0_8

    ! read in namelist file
    call pepc_read_parameters_from_first_argument(read_para_file, para_file)

    if (read_para_file) then
      if(root) write(*,'(a)') " == reading parameter file, section pepc-mini: ", para_file
      open(fid,file=para_file)
      read(fid,NML=pepcmini)
      close(fid)
    else
      if(root) write(*,*) " == no param file, using default parameter "
    end if

    tnp = tnp

    if(root) then
      write(*,'(a,i12)')    " == total number of particles : ", tnp
      write(*,'(a,i12)')    " == number of time steps      : ", nt
      write(*,'(a,es12.4)') " == time step                 : ", dt
      write(*,'(a,i12)')    " == diag & IO interval        : ", diag_interval
      write(*,'(a,l12)')    " == particle test             : ", particle_test
      write(*,'(a,l12)')    " == particle output           : ", particle_output
      write(*,'(a,l12)')    " == domain output             : ", domain_output
      write(*,'(a,l12)')    " == filter particles          : ", particle_filter
      write(*,'(a,l12)')    " == field probes              : ", particle_probe
    end if

    call pepc_prepare(3_kind_dim)

    ! reste current measurement
    current_q = 0.0_8
    current_I = 0.0_8

  end subroutine


  subroutine init_particles(p)
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

    ! put probes on z-axis
    if(particle_probe) then

      allocate(probes(nprobes, 5), stat=rc)
      if(rc.ne.0) write(*,*) " === probes allocation error!"

      do ip=1, nprobes
        probes(ip,1) = min_probes + (ip-1)*((max_probes-min_probes)/(nprobes-1))
      end do

    end if

    ! set random seed
    dummy = par_rand(my_rank)

    ! setup random qubic particle cloud
    do ip=1, np

      p(ip)%label       = my_rank * (tnp / n_ranks) + ip
      p(ip)%data%q      = (-1.0_8 + 2.0_8*MOD(p(ip)%label,2_kind_particle))
      p(ip)%data%m      = 1.0_8
      if(p(ip)%data%q .gt. 0.0) p(ip)%data%m = p(ip)%data%m * 100.0_8

      call random(p(ip)%x)
      p(ip)%x(3)         = p(ip)%x(3)*2*plasma_width - plasma_width

      call random_gauss(p(ip)%data%v)
      p(ip)%data%v      = p(ip)%data%v * sqrt(Tkb * p(ip)%data%m) / p(ip)%data%m

      p(ip)%results%e   = 0.0_8
      p(ip)%results%pot = 0.0_8
      p(ip)%work        = 1.0_8

    end do

    next_label = tnp+1

  end subroutine init_particles

  subroutine push_particles(p)
    use module_mirror_boxes
    implicit none

    type(t_particle), allocatable, intent(inout) :: p(:)
    integer(kind_particle) :: ip
    real*8  :: fact

    if(root) write(*,'(a)') " == [pusher] push particles "

    fact = dt

    do ip=1, np
      p(ip)%data%v = p(ip)%data%v + fact * p(ip)%data%q / p(ip)%data%m * p(ip)%results%e
      p(ip)%x      = p(ip)%x      + dt   * p(ip)%data%v
    end do

  end subroutine push_particles

  subroutine filter_particles(p)
    implicit none
    include 'mpif.h'

    type(t_particle), allocatable, intent(inout) :: p(:)
    integer(kind_particle) :: ip, rp
    integer(kind_default) :: rc
    integer :: new_particles
    real*8 :: local_q, global_q

    real,parameter :: Tkb=1e-5

    local_q = 0.0_8
    rp = 1
    do ip=1, np
      if( abs(p(rp)%x(3)) > vessel_width ) then

        local_q = local_q + p(rp)%x(3)/abs(p(rp)%x(3))*p(rp)%data%q

        if(rp .ne. np) p(rp) = p(np)
        np = np - 1
      else

        if(p(rp)%x(1) < 0.0_8) then
          p(rp)%x(1) = p(rp)%x(1) + 1.0_8
        else if(p(rp)%x(1) > 1.0_8) then
          p(rp)%x(1) = p(rp)%x(1) - 1.0_8
        end if

        if(p(rp)%x(2) < 0.0_8) then
          p(rp)%x(2) = p(rp)%x(2) + 1.0_8
        else if(p(rp)%x(2) > 1.0_8) then
          p(rp)%x(2) = p(rp)%x(2) - 1.0_8
        end if

        rp = rp + 1

      end if
    end do

    ! emmit electrons at the catode
    if(root) then

      new_particles = int(dt * current_emission)
      call reallocate_particles(particles, np, np+new_particles)

      do ip=np+1, np+new_particles

        p(ip)%label       = next_label
        next_label        = next_label + 1

        p(ip)%data%q      = -1.0_8
        p(ip)%data%m      =  1.0_8

        call random_gauss(p(ip)%data%v)
        p(ip)%data%v      = p(ip)%data%v * [1.0e-2_8, 1.0e-2_8, 1.0e5_8]*sqrt(Tkb * p(ip)%data%m) / p(ip)%data%m
        if(p(ip)%data%v(3) .gt. 0) p(ip)%data%v(3) = -p(ip)%data%v(3)

        call random(p(ip)%x(1:2))
        p(ip)%x(3)         = vessel_width + dt*p(ip)%data%v(3)

        p(ip)%results%e   = 0.0_8
        p(ip)%results%pot = 0.0_8
        p(ip)%work        = 1.0_8

      end do

      np = np + new_particles

    end if


    call MPI_ALLREDUCE(np, tnp, 1, MPI_KIND_PARTICLE, MPI_SUM, MPI_COMM_WORLD, rc)
    call MPI_ALLREDUCE(local_q, global_q, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)

    current_q = current_q + global_q
    current_I = global_q / dt

    if(root) then

      if(step .eq. 0) then
        open(unit=current_file_id, file="current.dat", status='unknown', position='rewind')
        write(current_file_id,'(a)') "# step(1), time(2), charge(3), current(4), el-field(5)"
      else
        open(unit=current_file_id, file="current.dat", status='unknown', position='append')
      end if

      write(current_file_id,'(i12,4es12.4)') step, step*dt, current_q, current_I, vessel_ez

      close(current_file_id)
    end if

    if(root) write(*,'(a,i12)')    " == [filter] total number of particles            : ", tnp
    if(root) write(*,'(a,es12.4)') " == [current] absorbed charge                     : ", current_q
    if(root) write(*,'(a,es12.4)') " == [current] measured current                    : ", current_I

  end subroutine filter_particles

  subroutine apply_external_field()

    implicit none

    integer(kind_particle) :: ip

    do ip=1, np

      particles(ip)%results%e(3) = particles(ip)%results%e(3) + vessel_ez
      particles(ip)%results%pot  = particles(ip)%results%pot  - particles(ip)%x(3)*vessel_ez

    end do

  end subroutine

  subroutine test_particles()

    use module_pepc_types
    use module_directsum
    implicit none
    include 'mpif.h'

    integer(kind_particle), allocatable   :: tindx(:)
    real*8, allocatable                   :: trnd(:)
    type(t_particle_results), allocatable :: trslt(:)
    integer(kind_particle)                :: tn, tn_global, ti
    integer                               :: rc
    real*8                                :: L2sum_local, L2sum_global, L2
    real*8                                :: ta, tb

    ta = get_time()

    if(allocated(direct_L2)) then
      deallocate(direct_L2)
    end if
    allocate(direct_L2(np))
    direct_L2 = -1.0_8

    tn = particle_direct / n_ranks
    if(my_rank.eq.(n_ranks-1)) tn = tn + MOD(particle_direct, n_ranks)

    allocate(tindx(tn), trnd(tn), trslt(tn))

    call random(trnd)

    tindx(1:tn) = int(trnd(1:tn) * (np-1)) + 1

    call directforce(particles, tindx, tn, trslt, MPI_COMM_WORLD)

    L2sum_local  = 0.0
    L2sum_global = 0.0
    do ti = 1, tn
      L2          = &
                    (particles(tindx(ti))%results%e(1) - trslt(ti)%e(1))**2+ &
                    (particles(tindx(ti))%results%e(2) - trslt(ti)%e(2))**2+ &
                    (particles(tindx(ti))%results%e(3) - trslt(ti)%e(3) - vessel_ez)**2
      L2sum_local = L2sum_local + L2
      direct_L2(tindx(ti)) = L2
    end do

    call MPI_ALLREDUCE(tn, tn_global, 1, MPI_KIND_PARTICLE, MPI_SUM, MPI_COMM_WORLD, rc)
    call MPI_ALLREDUCE(L2sum_local, L2sum_global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)

    L2sum_global = sqrt(L2sum_global) / tn_global

    tb = get_time()
    if(root) then
      write(*,'(a,i12)')    " == [direct test] number tested particles         : ", tn_global
      write(*,'(a,es12.4)') " == [direct test] L2 error in probed particles    : ", L2sum_global
      write(*,'(a,es12.4)') " == [direct test] time in test [s]                : ", tb - ta
    end if

    deallocate(tindx)
    deallocate(trnd)
    deallocate(trslt)

  end subroutine test_particles

  subroutine compute_field()
    use module_pepc
    implicit none

    include 'mpif.h'

    real*8             :: ta, tb
    integer, parameter :: fid = 12
    integer(kind_particle) :: ip
    integer            :: rc
    character(255)     :: fname
    integer            :: pos
    real*8             :: dz, dmin, dsize

    ta = get_time()

    ! gather mean/min/max potential values along the z-axis
    do ip=1, nprobes
      probes(ip, 2:5) = [1.0e-10_8, 0.0_8, 0.0_8, 0.0_8]
    end do
    dmin  = min_probes
    dsize = max_probes - dmin
    dz    = (dsize) / (1.0_8*nprobes)
    do ip=1, np
      pos = int( (particles(ip)%x(3) -dmin -0.5_8*dz)/(dz) ) + 1
      if(pos.ge.1 .and. pos.le.nprobes) then
        probes(pos, 2) = probes(pos, 2) + 1.0_8
        probes(pos, 3) = probes(pos, 3) + particles(ip)%data%m
        probes(pos, 4) = probes(pos, 4) + particles(ip)%data%q
        probes(pos, 5) = probes(pos, 5) + particles(ip)%results%pot
      end if
    end do

    call MPI_ALLREDUCE(MPI_IN_PLACE, probes(:,2), nprobes, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
    call MPI_ALLREDUCE(MPI_IN_PLACE, probes(:,3), nprobes, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
    call MPI_ALLREDUCE(MPI_IN_PLACE, probes(:,4), nprobes, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)
    call MPI_ALLREDUCE(MPI_IN_PLACE, probes(:,5), nprobes, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, rc)

    probes(:,5) = probes(:,5) / probes(:,2)

    ! write results to file
    write(fname, '("probe.",i6.6,".dat")') step
    if(root) write(*,'(2a)') " == [probes] write to file                        : ", trim(fname)
    open(unit=fid, file=fname, status='replace')
    write(fid, '(a)') "z-pos(1), x-y-binned[n(2), m(3), q(4), avg. pot(5)], ext-pot(6)"
    do ip=1, nprobes
      write(fid, '(6es12.4/)') probes(ip,1), probes(ip,2), probes(ip,3), probes(ip,4), probes(ip,5), probes(ip,1)*(-vessel_ez)
    end do
    close(fid)

    tb = get_time()

    if(root) then
      write(*,'(a,es12.4)') " == [probes] time in probes [s]                   : ", tb - ta
    end if

  end subroutine

  subroutine write_particles(p)
    use module_vtk
    use module_vtk_helpers
    implicit none

    include 'mpif.h'

    type(t_particle), intent(in) :: p(:)

    integer :: vtk_step
    real*8 :: ta, tb

    ta = get_time()

    if (step .eq. 0) then
      vtk_step = VTK_STEP_FIRST
    else if (step .eq. nt) then
      vtk_step = VTK_STEP_LAST
    else
      vtk_step = VTK_STEP_NORMAL
    endif

    call vtk_write_particles("particles", MPI_COMM_WORLD, step, dt * step, vtk_step, p, coulomb_and_l2)

    tb = get_time()

    if(root) write(*,'(a,es12.4)') " == [write particles] time in vtk output [s]      : ", tb - ta

    contains

    subroutine coulomb_and_l2(d, r, vtkf)
      use module_vtk
      use module_interaction_specific_types
      implicit none

      type(t_particle_data), intent(in) :: d(:)
      type(t_particle_results), intent(in) :: r(:)
      type(vtkfile_unstructured_grid), intent(inout) :: vtkf

      call vtk_write_particle_data_results(d, r, vtkf)
      if(particle_test) call vtkf%write_data_array("L2 error", direct_L2(:))
    end subroutine
  end subroutine write_particles

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

  real*8 function get_time()
    implicit none
    include 'mpif.h'

    get_time = MPI_WTIME()

  end function get_time

  subroutine random_gauss(list)
    implicit none

    real*8, intent(inout) :: list(:)

    real*8  :: v(2), pi, r, p
    integer :: n, i

    pi = 2.0_8*acos(0.0_8)
    n  = size(list)

    do i=1, n, 2

      call random(v)

      r = sqrt(-2.0_8 * log(v(1)))
      p = 2.0_8*pi*v(2)

      list(i)                = r * sin(p)
      if((i+1)<=n) list(i+1) = r * cos(p)

    end do

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

end module
