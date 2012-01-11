!!!!!!!!!!!!!!!!!!!!
!! helper module
!!!!!!!!!!!!!!!!!!!!

module helper

  use module_pepc_types
  implicit none

  ! MPI variables
  integer :: my_rank, n_ranks

  ! force summation variables
  real*8 :: theta, eps, dt

  ! control variables
  integer :: nt         ! number of timesteps
  integer :: tnp        ! total number of particles
  integer :: np         ! local number of particles
  logical :: vtk_output ! turn vtk output on/off

  contains

    subroutine set_parameter()
      
      use module_pepc
      use module_interaction_specific, only : theta2, eps2, force_law
      implicit none
      
      character(255) :: parameterfile
      logical :: read_param_file
      namelist /pepcmini/ theta, eps, tnp, nt, dt, vtk_output
      
      ! set default parameter values
      theta      = 0.6
      eps        = 1e-6
      tnp        = 1441
      nt         = 20
      dt         = 1e-3
      vtk_output = .true.
      
      ! read in namelist file
      call pepc_get_para_file(read_param_file, parameterfile, my_rank)

      if (read_param_file) then
        if(my_rank .eq. 0) write(*,*) \
          " == reading parameter file, section pepc-mini: ", parameterfile
        open(10,file=parameterfile)
        read(10,NML=pepcmini)
        close(10)
      else
        if(my_rank .eq. 0) write(*,*) " == no param file, using default parameter "
      end if
  
      ! initialize calc force params
      theta2      = theta**2
      eps2        = eps**2
      force_law   = 3

      ! set pepc parameter
      call pepc_prepare(3)

      if(my_rank .eq. 0) then
        write(*,*) " == theta                     : ", theta
        write(*,*) " == eps                       : ", eps
        write(*,*) " == total number of particles : ", tnp
        write(*,*) " == number of time steps      : ", nt
        write(*,*) " == time step                 : ", dt
        write(*,*) " == vtk output                : ", vtk_output
      end if

    end subroutine

	subroutine init_particles(p)
	  implicit none
	  
	  type(t_particle), allocatable, intent(inout) :: p(:)
	  integer :: n, ip

	  n = size(p)		
	  if(my_rank.eq.0) write(*,*) " == init particles "
	  
	  call random_seed()
	  
	  ! setup random qubic particle cloud
	  do ip=1, n
	    call random_number(p(ip)%x)
	    p(ip)%x      = p(ip)%x*2 - 1
	    p(ip)%data%v = 0.0
	    p(ip)%data%m = 1.0
	    call random_number(p(ip)%data%q)
	    p(ip)%data%q = (2.0*p(ip)%data%q) - 1.0
	    p(ip)%data%q = p(ip)%data%q / abs(p(ip)%data%q)
	    p(ip)%results%e = 0
	    p(ip)%results%pot = 0
	  end do
	
	  ! setup two massive and charged particles
	  if(my_rank .eq. 0) then
	    p(1)%x      = [-2.0, 0.0, 0.0]
	    p(1)%data%v = [0.0, 0.0, 0.0]
	    p(1)%data%m =  10.0 * tnp
	    p(1)%data%q = -1.0 * tnp
	    p(2)%x      = [2.0, 0.0, 0.0]
	    p(2)%data%v = [0.0, 0.0, 0.0]
	    p(2)%data%m =  10.0 * tnp
	    p(2)%data%q =  1.0 * tnp
	  end if
	
	end subroutine init_particles
	
	subroutine push_particles(p)
	  implicit none
	  
	  type(t_particle), allocatable, intent(inout) :: p(:)
	  integer :: n, ip
      real*8  :: fact

      n = size(p) - 2
      if(my_rank.eq.0) write(*,*) " == push particles "

      fact = dt

      do ip=1, n
        p(ip)%data%v = p(ip)%data%v + fact * p(ip)%results%e 
        p(ip)%x      = p(ip)%x      + \
                       fact * p(ip)%data%v * p(ip)%data%q / p(ip)%data%m
      end do

	end subroutine push_particles
	
	subroutine write_particles(p, step, final)
	  use module_vtk
	  implicit none
	  
	  type(t_particle), allocatable, intent(in) :: p(:)

      integer, intent(in) :: step
      logical, intent(in) :: final
      integer :: i
      type(vtkfile_unstructured_grid) :: vtk
      integer :: vtk_step
	  integer :: n
	  real*8 :: time
	  real*8 :: ta, tb
	  
	  ta = get_time()
	  
	  if(my_rank.eq.0) write(*,*) " == write particles "

      n = size(p) - 2
      time = 0.1_8 * step

      if (step .eq. 0) then
        vtk_step = VTK_STEP_FIRST
      else if (final) then
        vtk_step = VTK_STEP_LAST
      else
        vtk_step = VTK_STEP_NORMAL
      endif

      call vtk%create_parallel("particles", step, my_rank, n_ranks, \
           time, vtk_step)
      call vtk%write_headers(np, 0)
	  call vtk%startpoints()
	  call vtk%write_data_array("xyz", n, p(:)%x(1), \
	       p(:)%x(2), p(:)%x(3))
	  call vtk%finishpoints()
	  call vtk%startpointdata()
	  call vtk%write_data_array("velocity", n, p(:)%data%v(1), \
	       p(:)%data%v(2), p(:)%data%v(3))
	  call vtk%write_data_array("el_field", n, \
	       p(:)%results%e(1), \
	       p(:)%results%e(2), p(:)%results%e(3))
	  call vtk%write_data_array("el_pot", n, p(:)%results%pot)
	  call vtk%write_data_array("charge", np, p(:)%data%q)
	  call vtk%write_data_array("mass", np, p(:)%data%m)
	  call vtk%write_data_array("pelabel", n, p(:)%label)
	  call vtk%write_data_array("local index", n, [(i,i=1,n)])
	  call vtk%write_data_array("processor", n, [(my_rank,i=1,n)])
	  call vtk%finishpointdata()
      call vtk%dont_write_cells()
      call vtk%write_final()
      call vtk%close()

      tb = get_time()

      if(my_rank.eq.0) write(*,*) " == time in vtk output [s]: ", tb - ta

	end subroutine write_particles
	
	real*8 function get_time()
	  implicit none
	  include 'mpif.h'
	  
	  get_time = MPI_WTIME()
	  
	end function get_time
	
end module
