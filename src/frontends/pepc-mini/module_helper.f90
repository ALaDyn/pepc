!!!!!!!!!!!!!!!!!!!!
!! helper module
!!!!!!!!!!!!!!!!!!!!

module helper

  use module_pepc_types
  implicit none

  ! MPI variables
  integer :: my_rank, n_ranks
  logical :: root

  ! time variables
  real*8 :: dt
  integer :: step

  ! control variables
  integer :: nt               ! number of timesteps
  integer :: tnp              ! total number of particles
  integer :: np               ! local number of particles
  logical :: particle_output  ! turn vtk output on/off
  logical :: domain_output    ! turn vtk output on/off
  logical :: particle_filter  ! filter particles leaving simulation domain
  logical :: particle_probe   ! turn probin on/off
  logical :: particle_test    ! turn direct summation on/off
  integer :: diag_interval    ! number of timesteps between all diagnostics and IO
  
  real*8, parameter  :: plasma_width = 5.0_8      ! width of plasma block
  real*8, parameter  :: vessel_width = 5.0_8      ! distance between walls
  real*8             :: vessel_ez                 ! electrical field in vessel in z-direction
  
  integer, parameter :: nprobes    = 200           ! number of probes
  real*8, parameter  :: min_probes = -vessel_width ! min z-pos of probes
  real*8, parameter  :: max_probes =  vessel_width ! max z-pos of probes

  integer, parameter :: particle_direct = 144 ! number of particle for direct summation

  ! current measurement
  real*8 :: current_q, current_I
  real*8, parameter :: current_emission = 1.0e6_8
  integer, parameter :: current_file_id = 14
  
  ! particle data (position, velocity, mass, charge)
  type(t_particle), allocatable :: particles(:)
  real*8, allocatable           :: direct_L2(:)
  integer                       :: next_label

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

    namelist /pepcmini/ tnp, nt, dt, particle_output, domain_output, particle_filter, particle_test, particle_probe, diag_interval, vessel_ez
    
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

    call pepc_prepare(3)

    ! reste current measurement
    current_q = 0.0_8
    current_I = 0.0_8

  end subroutine


  subroutine init_particles(p)
    implicit none
    
    type(t_particle), allocatable, intent(inout) :: p(:)
    integer :: i, ip, rc
    integer :: rsize 
    integer, allocatable :: rseed(:)

    real,parameter :: Tkb=1e-9

    if(my_rank.eq.0) write(*,'(a)') " == [init] init particles "
    
    ! set initially number of local particles
    np = tnp / n_ranks
    if(my_rank.eq.(n_ranks-1)) np = np + MOD(tnp, n_ranks)

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

    
    call random_seed(size = rsize)
    allocate(rseed(rsize))
      rseed = my_rank + [(i*144,i=1,rsize)]
    call random_seed(put = rseed)
    deallocate(rseed)
    
    ! setup random qubic particle cloud
    do ip=1, np

      p(ip)%label       = my_rank * (tnp / n_ranks) + ip
      p(ip)%data%q      = (-1.0_8 + 2.0_8*MOD(p(ip)%label,2))
      p(ip)%data%m      = 1.0_8
      if(p(ip)%data%q .gt. 0.0) p(ip)%data%m = p(ip)%data%m * 100.0_8

      call random_number(p(ip)%x)
      p(ip)%x(3)         = p(ip)%x(3)*2*plasma_width - plasma_width

      call random_gauss(p(ip)%data%v) 
      p(ip)%data%v      = p(ip)%data%v * sqrt(Tkb * p(ip)%data%m) / p(ip)%data%m

      p(ip)%results%e   = 0.0_8
      p(ip)%results%pot = 0.0_8
      p(ip)%work        = 0

    end do  

    next_label = tnp+1

  end subroutine init_particles
    
  subroutine push_particles(p)
    use module_mirror_boxes
    implicit none
    
    type(t_particle), allocatable, intent(inout) :: p(:)
    integer :: ip
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
    integer :: ip, rp, rc
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

        call random_number(p(ip)%x(1:2))
        p(ip)%x(3)         = vessel_width + dt*p(ip)%data%v(3)

        p(ip)%results%e   = 0.0_8
        p(ip)%results%pot = 0.0_8
        p(ip)%work        = 0
      
      end do
    
      np = np + new_particles
    
    end if
    

    call MPI_ALLREDUCE(np, tnp, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, rc)
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
    
    integer :: ip
  
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
  
    integer, allocatable                  :: tindx(:)
    real*8, allocatable                   :: trnd(:)
    type(t_particle_results), allocatable :: trslt(:)
    integer                               :: tn, tn_global, ti, rc
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
  
    call random_number(trnd)
  
    tindx = int(trnd * (np-1)) + 1
  
    call directforce(particles, np, tindx, tn, trslt, my_rank, n_ranks, MPI_COMM_WORLD)
  
    L2sum_local  = 0.0
    L2sum_global = 0.0
    do ti = 1, tn
      L2          = &
                    (particles(tindx(ti))%results%e(1) - trslt(ti)%e(1))**2+ &
                    (particles(tindx(ti))%results%e(2) - trslt(ti)%e(2))**2+ &
                    (particles(tindx(ti))%results%e(3) - trslt(ti)%e(3))**2 
      L2sum_local = L2sum_local + L2
      direct_L2(tindx(ti)) = L2
    end do
        
    call MPI_ALLREDUCE(tn, tn_global, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, rc)
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
    integer            :: ip, rc
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
    implicit none
    
    type(t_particle), allocatable, intent(in) :: p(:)

    integer :: i
    type(vtkfile_unstructured_grid) :: vtk
    integer :: vtk_step
    real*8 :: time
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

    call vtk%create_parallel("particles", step, my_rank, n_ranks, time, vtk_step)
    call vtk%write_headers(np, 0)
    call vtk%startpoints()
    call vtk%write_data_array("xyz", np, p(:)%x(1), p(:)%x(2), p(:)%x(3))
    call vtk%finishpoints()
    call vtk%startpointdata()
    call vtk%write_data_array("velocity", np, p(:)%data%v(1), p(:)%data%v(2), p(:)%data%v(3))
    call vtk%write_data_array("el_field", np, p(:)%results%e(1), \
                              p(:)%results%e(2), p(:)%results%e(3))
    call vtk%write_data_array("el_pot", np, p(:)%results%pot)
    call vtk%write_data_array("charge", np, p(:)%data%q)
    call vtk%write_data_array("mass", np, p(:)%data%m)
    call vtk%write_data_array("work", np, p(:)%work)
    call vtk%write_data_array("pelabel", np, p(:)%label)
    call vtk%write_data_array("local index", np, [(i,i=1,np)])
    call vtk%write_data_array("processor", np, p(:)%pid)
    call vtk%write_data_array("L2 error", np, direct_L2(:))
    call vtk%finishpointdata()
    call vtk%dont_write_cells()
    call vtk%write_final()
    call vtk%close()

    tb = get_time()

    if(root) write(*,'(a,es12.4)') " == [write particles] time in vtk output [s]      : ", tb - ta

  end subroutine write_particles
  
  subroutine write_domain(p)
    
    use module_vtk
    use module_treediags
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
    call write_branches_to_vtk(step,  dt * step, vtk_step)
    call write_spacecurve_to_vtk(step, dt * step, vtk_step, p)
    
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
      
      call random_number(v)
      
      r = sqrt(-2.0_8 * log(v(1)))
      p = 2.0_8*pi*v(2)
      
      list(i)                = r * sin(p)
      if((i+1)<=n) list(i+1) = r * cos(p)
    
    end do
  
  end subroutine

  subroutine reallocate_particles(list, oldn, newn)
    implicit none
  
    type(t_particle), allocatable, intent(inout) :: list(:)
    integer, intent(in) :: oldn, newn
    
    type(t_particle) :: tmp_p(oldn)
    
    tmp_p = list
    deallocate(list)
    allocate(list(newn))
    list(1:oldn) = tmp_p
  
  end subroutine

end module