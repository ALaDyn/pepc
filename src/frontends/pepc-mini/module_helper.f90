!!!!!!!!!!!!!!!!!!!!
!! helper module
!!!!!!!!!!!!!!!!!!!!

module helper

  use module_pepc_types
  implicit none

  ! MPI variables
  integer :: my_rank, n_ranks

  ! time variables
  real*8 :: dt
  integer :: step

  ! control variables
  integer :: nt              ! number of timesteps
  integer :: tnp             ! total number of particles
  integer :: np              ! local number of particles
  logical :: particle_output ! turn vtk output on/off
  logical :: domain_output   ! turn vtk output on/off

  ! particle data (position, velocity, mass, charge)
  type(t_particle), allocatable :: particles(:)


  contains

    subroutine set_parameter()
      
      use module_pepc
      use module_interaction_specific, only : theta2, eps2, force_law
      implicit none
      
      character(255) :: para_file
      logical :: read_para_file
      namelist /pepcmini/ tnp, nt, dt, particle_output, domain_output
      
      ! set default parameter values
      tnp             = 1441
      nt              = 20
      dt              = 1e-3
      particle_output = .true.
      domain_output   = .true.
      
      ! read in namelist file
      call pepc_read_parameters_from_first_argument(read_para_file, para_file)
      !call pepc_get_para_file(read_para_file, para_file, my_rank)

      if (read_para_file) then
        if(my_rank .eq. 0) write(*,*) \
          " == reading parameter file, section pepc-mini: ", para_file
        open(10,file=para_file)
        read(10,NML=pepcmini)
        close(10)
      else
        if(my_rank .eq. 0) write(*,*) " == no param file, using default parameter "
      end if
  
      if(my_rank .eq. 0) then
        write(*,*) " == total number of particles : ", tnp
        write(*,*) " == number of time steps      : ", nt
        write(*,*) " == time step                 : ", dt
        write(*,*) " == particle output           : ", particle_output
        write(*,*) " == domain output             : ", domain_output
      end if

      call pepc_prepare(3)

    end subroutine

	subroutine init_particles(p)
	  implicit none
	  
	  type(t_particle), allocatable, intent(inout) :: p(:)
	  integer :: ip, rc
      integer :: rsize 
      integer, allocatable :: rseed(:)

	  if(my_rank.eq.0) write(*,*) " == init particles "
	  
	  ! set initially number of local particles
      np = tnp / n_ranks
      if(my_rank.eq.0) np = np + MOD(tnp, n_ranks)

      allocate(particles(np), stat=rc)
      if(rc.ne.0) write(*,*) " === particle allocation error!"
	  
	  call random_seed(size = rsize)
	  allocate(rseed(rsize))
      rsize = my_rank + 144
	  call random_seed(put = rseed)
	  deallocate(rseed)
	  
	  ! setup random qubic particle cloud
	  do ip=1, np
	    call random_number(p(ip)%x)
	    p(ip)%x      = p(ip)%x*2 - 1
	    p(ip)%data%v = 0.0
	    p(ip)%data%m = 1.0
	    call random_number(p(ip)%data%q)
	    p(ip)%data%q = (2.0*p(ip)%data%q) - 1.0
	    p(ip)%data%q = p(ip)%data%q / abs(p(ip)%data%q)
	    p(ip)%results%e = 0
	    p(ip)%results%pot = 0
	    p(ip)%work = 0
	  end do
	
	  ! setup two massive and charged particles
	  if(my_rank .eq. 0) then
	    p(1)%x      = [0.0, 0.0, -1.8]
	    p(1)%data%v = [0.0, 0.0, 0.0]
	    p(1)%data%m = 1000.0 * tnp
	    p(1)%data%q = 10.0 * tnp
	    p(2)%x      = p(1)%x      * (-1.0)
	    p(2)%data%v = p(1)%data%v
	    p(2)%data%m = p(1)%data%m
	    p(2)%data%q = p(1)%data%q * (-1.0)
	  end if
	
	end subroutine init_particles
	
	subroutine push_particles(p)
	  implicit none
	  
	  type(t_particle), allocatable, intent(inout) :: p(:)
	  integer :: ip
      real*8  :: fact

      if(my_rank.eq.0) write(*,*) " == push particles "

      fact = dt

      do ip=1, np
        p(ip)%data%v = p(ip)%data%v + fact * p(ip)%results%e 
        p(ip)%x      = p(ip)%x      + \
                       fact * p(ip)%data%v * p(ip)%data%q / p(ip)%data%m
      end do

	end subroutine push_particles
	
    subroutine filter_particles(p)
	  implicit none
	  include 'mpif.h'
	  
	  type(t_particle), allocatable, intent(inout) :: p(:)
	  integer :: ip, rc

      real*8 :: dmin(3), dmax(3)

      dmin = [-2.0, -2.0, -2.0]
      dmax = dmin * (-1.0)

      if(my_rank.eq.0) write(*,*) " == filter particles "

      do ip=1, np
        if(any(p(ip)%x .lt. dmin .or. p(ip)%x .gt. dmax)) then
          p(ip) = p(np)
          np = np - 1
        end if
      end do

      call MPI_ALLREDUCE(np, tnp, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, rc)

      !write(*,*) " == particles on rank", my_rank, ": ", np
      if(my_rank.eq.0) write(*,*) " == total number of particles : ", tnp

	end subroutine filter_particles
	
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

	  if(my_rank.eq.0) write(*,*) " == write particles "

      call vtk%create_parallel("particles", step, my_rank, n_ranks, \
           time, vtk_step)
      call vtk%write_headers(np, 0)
	  call vtk%startpoints()
	  call vtk%write_data_array("xyz", np, p(:)%x(1), \
	       p(:)%x(2), p(:)%x(3))
	  call vtk%finishpoints()
	  call vtk%startpointdata()
	  call vtk%write_data_array("velocity", np, p(:)%data%v(1), \
	       p(:)%data%v(2), p(:)%data%v(3))
	  call vtk%write_data_array("el_field", np, \
	       p(:)%results%e(1), \
	       p(:)%results%e(2), p(:)%results%e(3))
	  call vtk%write_data_array("el_pot", np, p(:)%results%pot)
	  call vtk%write_data_array("charge", np, p(:)%data%q)
	  call vtk%write_data_array("mass", np, p(:)%data%m)
	  call vtk%write_data_array("pelabel", np, p(:)%label)
	  call vtk%write_data_array("local index", np, [(i,i=1,np)])
	  call vtk%write_data_array("processor", np, p(:)%pid)
	  call vtk%finishpointdata()
      call vtk%dont_write_cells()
      call vtk%write_final()
      call vtk%close()

      tb = get_time()

      if(my_rank.eq.0) write(*,*) " == time in vtk output [s]: ", tb - ta

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
	
end module
