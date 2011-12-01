module physvars

  use treetypes

  type(t_particle), allocatable :: vortex_particles(:)

  real*8, parameter :: pi=3.141592654
    
  !  physics data
  integer  :: n   ! # vortices
  integer  :: np  ! # vortices per PE
  real     :: force_const    ! force constant depending on unit system
  real     :: eps            ! potential/force law cutoff
  real     :: h, m_h         ! initial particle distance and mesh width for remeshing
  integer  :: rem_freq       ! remeshing frequence
  real*8   :: kernel_c       ! mod. remeshing kernel parameter
  real*8   :: thresh         ! vorticity threshold: particles with lower vorticity mag. will be kicked out (mandatory to avoid zero abs_charge)
  real     :: nu             ! viscosity parameter
  real     :: rmax    ! radius of torus segment
  real     :: rl      ! temp radius of current circle
  real     :: r_torus ! radius of torus
  integer  :: nc      ! # circles per torus segment
  integer  :: nphi    ! # torus segments
  integer  :: ns      ! # particles per torus segment
  real     :: g       ! # vorticity amplifier
  real, dimension(3) :: torus_offset  ! shifts coords of both tori (one with +, one with -)

 ! tree stuff
  real :: theta       ! Clumping parameter
  integer :: mac      ! MAC (default=BH)

 ! Variables needing 'copy' for tree routines
  real    :: np_mult
  integer :: my_rank       ! Rank of current task
  integer :: n_cpu         ! # cpus used by program
  integer :: db_level      ! printed o/p debug level

  ! Control stuff
  integer :: ispecial       ! Switch to select special electron configs 
  integer :: weighted       ! load balancing switch
  real    :: dt, ts, te     ! timestep, start-time, end-time
  integer :: nt             ! # timesteps and current timestep
  integer :: rk_stages      ! # Runge-Kutta stages
  integer :: curve_type     !< type of space-filling curve, 0=z-curve, 1=Hilbert-curve
  integer :: ifile_cpu      ! O/P stream

contains


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !>   Initialize physvars variables, allocate frontend arrays, read-in parameters
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine pepc_setup()

        use tree_walk_pthreads
        implicit none
        include 'mpif.h'

        integer :: ierr

        character(50) :: parameterfile
        integer :: read_param_file

        integer*4 IARGC



        namelist /pepcdata/ np_mult, n, num_walk_threads, max_particles_per_thread, &
        mac, theta, eps, ispecial, weighted, curve_type, dt, ts, te, db_level, &
        h, m_h, nu, rem_freq, thresh, &
        rmax, r_torus, nc, nphi, g, torus_offset


        !  Default input set
        db_level        =   0
        np_mult         = -45
        ispecial        =   1
        weighted        =   1
        curve_type      =   0


        ! particles
        np = 0
        n  = 1000 ! total # vortex particles

        ! physics stuff
        force_const  = 0.25D00/pi  ! 3D prefactor for u and af
        mac          = 0
        theta        = 0.6
        eps          = 0.01      ! this is my smoothing radius, will adapt it in special_start
        h            = 0.
        m_h          = 0.
        rem_freq     = 0
        thresh       = 0.
        rmax         = 0.
        r_torus      = 0.
        nc           = 0
        nphi         = 0
        g            = 0
        torus_offset = [0., 0., 0.]

        ! control
        ts           = 0.
        te           = 0.
        dt           = 0.01

        ! rank 0 reads in first command line argument
        read_param_file = 0
        if (my_rank .eq. 0) then
            if( IARGC() .ne. 0 ) then
                call GETARG(1, parameterfile)
                read_param_file = 1
            end if
        end if

        call MPI_BCAST( read_param_file, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )

        ! broadcast file name, read actual inputs from namelist file

        if (read_param_file .eq. 1) then

            call MPI_BCAST( parameterfile, 50, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr )

            if(my_rank .eq. 0) write(*,*) "reading parameter file: ", parameterfile
            open(10,file=parameterfile)
            read(10,NML=pepcdata)
            close(10)
        else
            if(my_rank .eq. 0) write(*,*) "##### using default parameter #####"
        end if

        init: select case(ispecial)
            case(1,2)                         ! Vortex rings
                rl = rmax/(2*nc+1)
                ns = 1+4*nc*(nc+1)
                np = int(2.0*Ns*Nphi/n_cpu)
                n  = 2*Ns*Nphi
                kernel_c = sqrt(nu*rem_freq*dt)/m_h
            case default
                write(*,*) 'ERROR: need to specify setup via ispecial in your .h file'
                call MPI_ABORT(MPI_COMM_WORLD,ierr)
                stop
        end select init

        nt = int((te-ts)/dt) ! Number of timesteps
        rk_stages = 2   ! TODO: inflexible RK time integration scheme, hard-wired so far

        allocate ( vortex_particles(np) )

        if (my_rank == 0) then
            write(*,*) "Starting PEPC-MINI with",n_cpu," Processors, simulating",np, &
            " Particles on each Processor in",nt,"timesteps..."
            write(*,*) "Using",num_walk_threads,"worker-threads and 1 communication thread in treewalk on each processor (i.e. per MPI rank)"
            write(*,*) "Maximum number of particles per work_thread = ", max_particles_per_thread
        end if

    end subroutine pepc_setup


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !>   Deallocate physvars arrays, copy back variables
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine cleanup(my_rank_l,n_cpu_l)

        implicit none
        include 'mpif.h'

        integer, intent(in) :: my_rank_l ! MPI cpu rank
        integer, intent(in) :: n_cpu_l  ! MPI # CPUs

        ! copy call parameters to physvars module

        my_rank = my_rank_l
        n_cpu = n_cpu_l

        ! particle array deallocation in physvars

        deallocate ( vortex_particles )

    end subroutine cleanup


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !>  Set time stamp
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine stamp(istream,ibegin)
        implicit none
        integer, intent(in) :: istream, ibegin
        character :: cdate*8, ctime*10, czone*5

        !      call DATE_AND_TIME(cdate,ctime,czone,vals)
        call DATE_AND_TIME(cdate,ctime,czone)

        if (ibegin.eq.1) then

            write(istream,'(//a20,a12/a20,a12/a20,a12//)') 'PEPC run on ' &
            ,cdate(7:8)//'/'//cdate(5:6)//'/'//cdate(1:4) &
            ,'Time: ',ctime(1:2)//':'//ctime(3:4)//':'//ctime(5:6),' GMT+',czone

        else
            write(istream,'(a,a9)') 'Finished run at time: ',ctime(1:2)//':'//ctime(3:4)//':'//ctime(5:6)
        endif
    end subroutine stamp

end module physvars





