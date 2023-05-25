! This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.
! 
! Copyright (C) 2002-2023 Juelich Supercomputing Centre,
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

module physvars

  use module_pepc_types
  use module_pepc_kinds
  implicit none

  type(t_particle), allocatable :: vortex_particles(:)

  real(kind_physics), parameter :: pi=3.141592654d0
  integer,            parameter :: nDeltar=3     ! ratio btw Rd and m_h
    
  !  physics data
  integer(kind_particle) :: n   ! # vortices
  integer(kind_particle) :: np  ! # vortices per PE
  real(kind_physics)     :: force_const    ! force constant depending on unit system
  real                   :: h, m_h         ! initial particle distance and mesh width for remeshing
  real                   :: Delta_r        ! mesh width for remeshing (m_h and Deltar are the same quantity; m_h should be removed)
  real                   :: Rd             ! diffusive radius
  real                   :: Delta_tavv     ! advective time step
  real                   :: Delta_tdiff    ! diffusive time step
  integer                :: rem_freq       ! remeshing frequence
  integer                :: nv_on_Lref     ! vortices number on Lref
  real(kind_physics)     :: kernel_c       ! mod. remeshing kernel parameter
  real(kind_physics)     :: thresh         ! vorticity threshold: particles with lower vorticity mag. will be kicked out (mandatory to avoid zero abs_charge)
  real                   :: nu             ! viscosity parameter
  real                   :: Uref           ! Reference velocity
  real                   :: Lref           ! Reference length
  real                   :: rmax    ! radius of torus segment (ispecial=1,2)
  real                   :: rl      ! temp radius of current circle (ispecial=1,2)
  real                   :: r_torus ! radius of torus (ispecial=1,2)
  integer                :: nc      ! # circles per torus segment (ispecial=1,2)
  integer                :: nphi    ! # torus segments (ispecial=1,2)
  integer                :: ns      ! # particles per torus segment (ispecial=1,2)
  integer                :: n_in    ! # particles on the sphere (ispecial=3)
  real                   :: Co      ! # Courant number
  real                   :: Re      ! # Reynolds number
  real                   :: g       ! # vorticity/smoothing amplifier (ispecial=1,2,3)
  real, dimension(3)     :: torus_offset  ! shifts coords of both tori, one with +, one with - (ispecial=1,2)

 ! Variables needing 'copy' for tree routines
  integer :: my_rank       ! Rank of current task
  integer :: n_cpu         ! # cpus used by program

  ! Control stuff
  integer :: ispecial       ! Switch to select special electron configs 
  real    :: dt, ts, te     ! timestep, start-time, end-time
  integer :: nt             ! # timesteps and current timestep
  integer :: rk_stages      ! # Runge-Kutta stages

  ! I/O stuff
  integer :: ifile_cpu            ! O/P stream
  integer :: dump_time, cp_time   ! When to dump, when to do a checkpoint (read-in)
  integer :: input_itime          ! Which step should we read in?
  character(50) :: mpifile        ! MPI-IO-file

contains


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !>   Initialize physvars variables, allocate frontend arrays, read-in parameters
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine pepc_setup(itime, trun)

        use module_pepc
        use mpi
        implicit none

        integer, intent(out) :: itime
        real, intent(out) :: trun

        integer :: ierr

        character(255) :: parameterfile
        logical :: read_param_file

        namelist /pepcv/ n, ispecial, ts, te, nu, Co, nv_on_Lref,         &
                         h, Delta_r, thresh, Uref, Lref,                  &
                         rmax, r_torus, nc, nphi, g, torus_offset, n_in,  &
                         dump_time, cp_time, input_itime


        !  Default input set
        ispecial        =   1

        ! physics stuff
        force_const  = 0.25D00/pi  ! 3D prefactor for u and af
        h            = 0.
        Delta_r      = 0.
        Co           = 1.
        rem_freq     = 0
        thresh       = 1.E-7
        rmax         = 0.
        r_torus      = 0.
        nc           = 0
        nphi         = 0
        g            = 2
        torus_offset = [0., 0., 0.]
        Uref         = 1.d0
        Lref         = 1.d0
                
        ! control
        ts           = 0.
        te           = 0.
        dt           = 0.01

        dump_time    = 0
        cp_time      = 0
        

        ! read in first command line argument
        call pepc_get_para_file(read_param_file, parameterfile, my_rank)


        if (read_param_file) then

            if(my_rank .eq. 0) write(*,*) "reading parameter file, section pepcv: ", parameterfile
            open(10,file=parameterfile)
            read(10,NML=pepcv)
            close(10)

        else
            if(my_rank .eq. 0) write(*,*) "##### using default parameter #####"
        end if

        ! New parameters DVH: 0.336 is 0.021 (from DVH) * 16 --> Rd = 4 Dr and fr = 0.021 * Re * Rd^2/L^2 * L/(U Dta)
        ! Dtd = fr * Dta
        if(Delta_r.le.0.) Delta_r = Lref/float(nv_on_Lref)

        m_h = Delta_r
       
        thresh = thresh * Uref/Lref * Delta_r**3

        Re = Uref * Lref/nu
        
        Rd = nDeltar * m_h
 
        Delta_tdiff = 2.1d-2 * Rd*Rd / nu        ! Cfr. formula: 3.1 CiCP Colagrossi 2015
        Delta_tavv  = Co * Delta_r/Uref          ! Cfr. formula: 16 CMAME Rossi 2022
        rem_freq = int(Delta_tdiff/Delta_tavv)   ! Cfr. formula: 17 CMAME Rossi 2022
        
        if(rem_freq.eq.0) then
          if(my_rank .eq. 0) write(*,*) 'Diffusive Dt is dominant over advective'
          if(my_rank .eq. 0) write(*,*) 'Dt_avv/Dt_diff, Dt_diff',Delta_tavv/Delta_tdiff, Delta_tdiff
          Delta_tavv = Delta_tdiff
          rem_freq = 1
        end if
                            
        Delta_tdiff = rem_freq * Delta_tavv    ! Recalculation for convenience: not needed without multi-resolution
        
        dt = Delta_tavv
                        
!       kernel_c = dsqrt(nu * Delta_tdiff) / m_h
        kernel_c = 4.d0 * nu * Delta_tdiff

        if(my_rank.eq.0) then
          write(*,*) 'Particle volume',m_h**3
          write(*,*) 'Remeshing every ',rem_freq
          write(*,*) 'Diffusive radius ',Rd
          write(*,*) 'Delta t ',dt
          write(*,*) 'Circulation threshold ',thresh
        end if
        ! Setup itime to standard value, may change for ispecial=99 (if not: fine)
        itime = 0

        init: select case(ispecial)

            case(1,2,4)                         ! Vortex rings/wakes
 
                rmax = Lref/2.d0
                nc = nv_on_Lref/2
                rl = rmax/(2*nc+1)
                ns = 1+4*nc*(nc+1)
                nphi = nint(pi * r_torus/rl)
                n = 2*ns*nphi
                np = ceiling(1.0*n/n_cpu)+1

            case(3)                           ! Sphere

                n = n_in
                np = ceiling(1.0*n/n_cpu)
                h = sqrt(4.0*pi/n)
                !force_const = force_const*h**3
                !eps = g*h

            case(5)
                                    ! Vortex wake
                h = 2*pi/nc
                n = 2*nc*ceiling(1.0/h)*ceiling(2*pi/h)
                np = ceiling(1.0*n/n_cpu)+1

            case(6)                         ! Single Vortex ring
 
                rmax = Lref/2.d0
                nc = nv_on_Lref/2
                rl = rmax/(2*nc+1)
                ns = 1+4*nc*(nc+1)
                nphi = nint(pi * r_torus/rl)
                n = ns*nphi
                np = ceiling(1.0*n/n_cpu)+1

            case(98)

                n=n_in
                np = ceiling(1.0*n/n_cpu)
                h = 1./n

            case(99)                          ! Setup MPI checkpoint readin

                call setup_MPI_IO_readin(itime)

            case default

                write(*,*) 'ERROR: need to specify setup via ispecial in your .h file'
                call MPI_ABORT(MPI_COMM_WORLD,1,ierr)

        end select init

        ! Setup time variables
        trun = ts+itime*dt
        nt = ceiling((te-ts)/dt) ! Number of timesteps
        rk_stages = 2   ! TODO: inflexible RK time integration scheme, hard-wired so far

        allocate ( vortex_particles(np) )

        if (my_rank == 0) then
            write(*,*) "Starting PEPC-V with",n_cpu," Processors, simulating",np, &
            " Particles on each Processor in",nt,"timesteps..."
            !write(*,*) "Using",num_threads,"worker-threads and 1 communication thread in treewalk on each processor (i.e. per MPI rank)"
            !write(*,*) "Maximum number of particles per work_thread = ", max_particles_per_thread
        end if

    end subroutine pepc_setup



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !>
    !>   Deallocate physvars arrays, copy back variables
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine cleanup(my_rank_l,n_cpu_l)
        use mpi
        implicit none

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
    !>   Get ready for MPI checkpoint readin (should be in module file, but then we`d have a circular dependence...)
    !>   Take care that reading data matches writing it in files.f90:dump()
    !>
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine setup_MPI_IO_readin(itime)
        use mpi
        implicit none

        integer, intent(out) :: itime

        integer :: status(MPI_STATUS_SIZE)
        integer :: ierr, err, fh, tmp_rem_freq, tmp_i
        integer(kind_particle) :: remain
        real ::tmp_dt, tmp_te, tmp_nu, tmp_ts, tmp_h, tmp_m_h
        real(kind_physics) :: tmp_thresh

        write(mpifile,'(a,i6.6,a)') "part_data/particle_", input_itime,".mpi"
        call MPI_FILE_OPEN(MPI_COMM_WORLD,mpifile,IOR(MPI_MODE_RDWR,MPI_MODE_CREATE),MPI_INFO_NULL,fh,ierr)
        if (ierr .ne. MPI_SUCCESS) then
            write(*,*) 'something is wrong here: file open failed',my_rank,ierr,mpifile
            call MPI_ABORT(MPI_COMM_WORLD,err,ierr)
        end if

        ! Set file view to BYTE for header and read
        call MPI_FILE_SET_VIEW(fh,0_MPI_OFFSET_KIND, MPI_BYTE, MPI_BYTE, 'native', MPI_INFO_NULL, ierr)
        call MPI_FILE_READ(fh,n,1,MPI_KIND_PARTICLE,status,ierr)             ! # particles
        call MPI_FILE_READ(fh,tmp_dt,1,MPI_REAL,status,ierr)                 ! timestep
        call MPI_FILE_READ(fh,tmp_ts,1,MPI_REAL,status,ierr)                 ! Starting time
        call MPI_FILE_READ(fh,tmp_i,1,MPI_INTEGER,status,ierr)               ! Last successful timestep (number)
        call MPI_FILE_READ(fh,tmp_te,1,MPI_REAL,status,ierr)                 ! Final time
        call MPI_FILE_READ(fh,tmp_nu,1,MPI_REAL,status,ierr)                 ! Viscousity
        call MPI_FILE_READ(fh,tmp_h,1,MPI_REAL,status,ierr)                  ! Original particle distance
        call MPI_FILE_READ(fh,tmp_m_h,1,MPI_REAL,status,ierr)                ! Remeshing distance
        call MPI_FILE_READ(fh,tmp_rem_freq,1,MPI_INTEGER,status,ierr)        ! Remeshing frequence
        call MPI_FILE_READ(fh,tmp_thresh,1,MPI_KIND_PHYSICS,status,ierr)     ! threshold for pop. control
        call MPI_FILE_CLOSE(fh,ierr)

        dt = tmp_dt
        ts = tmp_ts
        nu = tmp_nu
        h = tmp_h
        m_h = tmp_m_h
        itime = tmp_i
        rem_freq = tmp_rem_freq
        thresh = tmp_thresh

        np = n/n_cpu
        remain = n-np*n_cpu
        if ((remain .gt. 0) .and. (my_rank .lt. remain)) np = np+1

    end subroutine setup_MPI_IO_readin


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





