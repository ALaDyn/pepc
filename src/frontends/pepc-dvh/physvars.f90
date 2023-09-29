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

   real(kind_physics), parameter :: pi = 3.141592654d0

   ! Physics data
   integer(kind_particle) :: n              ! # vortices
   integer(kind_particle) :: np             ! # vortices per PE
   integer(kind_particle) :: nDeltar        ! ratio btw Rd and m_h
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
   real                   :: rmax           ! radius of torus segment (ispecial=1,2)
   real                   :: rl             ! temp radius of current circle (ispecial=1,2)
   real                   :: r_torus        ! radius of torus (ispecial=1,2)
   integer                :: nc             ! # circles per torus segment (ispecial=1,2)
   integer                :: nphi           ! # torus segments (ispecial=1,2)
   integer                :: ns             ! # particles per torus segment (ispecial=1,2)
   integer                :: n_in           ! # particles on the sphere (ispecial=3)
   real                   :: Co             ! # Courant number
   real                   :: Re             ! # Reynolds number
   real                   :: g              ! # vorticity/smoothing amplifier (ispecial=1,2,3)
   real, dimension(3)     :: torus_offset   ! shifts coords of both tori, one with +, one with - (ispecial=1,2)

   ! Variables needing 'copy' for tree routines
   integer :: my_rank              ! Rank of current task
   integer :: n_cpu                ! # cpus used by program

   ! Control stuff
   integer :: ispecial             ! Switch to select special electron configs
   real    :: dt, ts, te           ! timestep, start-time, end-time
   integer :: nt                   ! # timesteps and current timestep
   integer :: rk_stages            ! order of the Runge-Kutta time integration scheme
   logical :: vort_check           ! Control of total vorticity during remeshing

   ! I/O stuff
   integer :: ifile_cpu            ! O/P stream
   integer :: dump_time, cp_time   ! How many outputs in the simulations; when to do a checkpoint (read-in)
   real    :: t_out, dt_out        ! output time and output timestep
   integer :: input_itime, n_out   ! Which step should we read in? Number of printed time instants
   character(50) :: mpifile        ! MPI-IO-file

   ! `Clone` the main datatype t_particle
   ! The aim is to save memory during remeshing. This is _not_ nice or clean and does not
   ! really follow PEPC`s define-centrally-once approach.
   ! This also requires to register the data type with MPI, so make sure init routines take
   ! care of that.
   type t_particle_data_short ! 15*8 byte = 120 byte -> 24 byte
      real(kind_physics) :: alpha(3)    ! vorticity or better: alpha = vorticity * volume
      real(kind_physics) :: vol
      !real(kind_physics) :: x_rk(3)     ! position temp array for Runge-Kutta time integration (required since particles get redistributed between substeps)
      !real(kind_physics) :: alpha_rk(3) ! vorticity temp array for Runge-Kutta time integration (required since particles get redistributed between substeps)
      !real(kind_physics) :: u_rk(3)     ! velocity temp array for Runge-Kutta time integration (required since particles get redistributed between substeps)
      !real(kind_physics) :: af_rk(3)    ! vorticity RHS temp array for Runge-Kutta time integration (required since particles get redistributed between substeps)
   end type t_particle_data_short
   integer, private, parameter :: nprops_particle_data_short = 2 !5
   integer :: MPI_TYPE_PARTICLE_DATA_SHORT_sca
   !type, private t_particle_results_short ! 7*8 byte = 56 byte -> 0 byte
   !   real(kind_physics), dimension(3) :: u   ! velocities
   !   real(kind_physics), dimension(3) :: af  ! RHS for vorticity/alpha ODE
   !   real(kind_physics) :: div               ! divergence
   !end type t_particle_results_short
   !integer, private, parameter :: nprops_particle_results_short = 3
   type t_particle_short ! 7*8 + 56 + 120 byte = 232 byte -> 64 byte
      real(kind_physics) :: x(1:3)      !< coordinates
      real*8 :: work        !< work load from force sum, ATTENTION: the sorting library relies on this being a real*8
      integer(kind_key) :: key      !< particle key, i.e. key on highest tree level
      !integer(kind_node) :: node_leaf !< node index of corresponding leaf (tree node)
      !integer(kind_particle) :: label     !< particle label (only for diagnostic purposes, can be used freely by the frontend
      type(t_particle_data_short) :: data       !< real physics (charge, etc.)
      !type(t_particle_results_short) :: results !< results of calc_force_etc and companions
   end type t_particle_short
   integer, private, parameter :: nprops_particle_short = 4 ! 7
   integer :: MPI_TYPE_PARTICLE_SHORT_sca
   interface assignment(=)
      module procedure assign_particles
   end interface
#ifdef USER_REDUCTION
   interface operator(+)
      module procedure omp_reduction_add
   end interface

!$OMP DECLARE REDUCTION (+ : t_particle_short : omp_out = omp_in + omp_out) &
!$OMP INITIALIZER (omp_priv = omp_orig)
#endif

contains

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !>
   !>   Initialize physvars variables, allocate frontend arrays, read-in parameters
   !>
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine pepc_setup(itime, trun)

      use module_pepc
      use module_interaction_specific, only: sig2
      use mpi
      implicit none

      integer, intent(out) :: itime
      real, intent(out) :: trun
      real*8            :: diff_error, eps_diff, alpha

      integer :: ierr, f_unit

      character(255) :: parameterfile
      logical :: read_param_file

      namelist /pepcv/ n, ispecial, ts, te, nu, Co, nv_on_Lref, rk_stages, &    !&
                       h, Delta_r, thresh, Uref, Lref, &                        !&
                       rmax, r_torus, nc, nphi, g, torus_offset, n_in, &        !&
                       dump_time, cp_time, input_itime, nDeltar, vort_check     !&

      !  Default input set
      ispecial     = 1                                                          !&

      ! Physics stuff
      force_const  = 0.25d0 / pi !& 3D prefactor for u and af                   !&
      h            = 0.                                                         !&
      Delta_r      = 0.                                                         !&
      Co           = 1.                                                         !&
      rem_freq     = 0                                                          !&
      nDeltar      = 3                                                          !&
      rk_stages    = 4                                                          !&
      thresh       = 1.d-7                                                      !&
      rmax         = 0.                                                         !&
      r_torus      = 0.                                                         !&
      nc           = 0                                                          !&
      nphi         = 0                                                          !&
      g            = 2                                                          !&
      torus_offset = [0., 0., 0.]                                               !&
      Uref         = 1.d0                                                       !&
      Lref         = 1.d0                                                       !&
      eps_diff     = 1.d-5                                                      !&

      ! Control
      ts           = 0.                                                         !&
      te           = 0.                                                         !&
      dt           = 0.01                                                       !&

      dump_time    = 1                                                          !&
      cp_time      = 0                                                          !&
       t_out       = 0.                                                         !&
      dt_out       = 0.                                                         !&

      ! Check total vorticity
      vort_check   = .true.                                                     !&

      ! read in first command line argument
      call pepc_get_para_file(read_param_file, parameterfile, my_rank)

      if (read_param_file) then
         if (my_rank .eq. 0) write (*, *) "reading parameter file, section pepcv: ", parameterfile
         open (newunit=f_unit, file=parameterfile)
         read (f_unit, NML=pepcv)
         close (f_unit)
      else
         if (my_rank .eq. 0) write (*, *) "##### using default parameter #####"
      end if

      dt_out = (te - ts)/dfloat(max(dump_time,1))

      if (Delta_r .le. 0.) Delta_r = Lref / float(nv_on_Lref)

      m_h = Delta_r

      sig2 = 0.d0 !&  m_h**2 ! size of the smoothing length set equal to the vortex size
      !           !&  one could also set it equal to zero to represent singular vortex distribution

      thresh = thresh * Uref / Lref * Delta_r**3

      Re = Uref * Lref / nu

      Rd = nDeltar * m_h

! Diffusive Delta_t in 3D has nonlinear expression. Trying the 2D evaluation and correct.
      ! 3D correction
      alpha = 1.d0

100   continue
      diff_error = erfc(1.d0/sqrt(alpha)) + 2.d0 * dexp(-1.d0/alpha)/dsqrt(alpha*pi)

      if(diff_error.gt.eps_diff) then
        alpha = 0.95d0 * alpha
        goto 100
      end if
      if(my_rank.eq.0) write(*,*) 'Diffusion error (csi) ', diff_error, alpha, alpha * Rd*Rd/4.d0/nu

      Delta_tdiff = alpha * Rd * Rd / 4.d0 / nu

      Delta_tavv = Co * Delta_r / Uref              ! Cfr. formula: 16 CMAME Rossi 2022
      if(my_rank.eq.0) write(*,*) 'first evaluation of Dt avv', Delta_tavv

      rem_freq = int(Delta_tdiff / Delta_tavv) + 1  ! Cfr. formula: 17 CMAME Rossi 2022

      if (rem_freq .eq. 1) then
         if (my_rank .eq. 0) write (*, *) 'Diffusive Dt is dominant over advective'
         if (my_rank .eq. 0) write (*, *) 'Dt_avv/Dt_diff, Dt_diff', Delta_tavv / Delta_tdiff, Delta_tdiff
      end if

      Delta_tavv = Delta_tdiff / float(rem_freq)   ! Recalculation for convenience: not needed without multi-resolution
      if(my_rank.eq.0) write(*,*) 'Modified evaluation of Dt avv', Delta_tavv

      dt = Delta_tavv
      if (dt >= dt_out) then
         if(my_rank.eq.0) write(*,*) 'Time step is larger than print time!!!'
         if(my_rank.eq.0) write(*,*) 'dt', dt,'>= dt_out',dt_out
         if(my_rank.eq.0) write(*,*) 'Either decrease the Courant number or decrease dump_time!'
         stop
      endif

      kernel_c = 4.d0 * nu * Delta_tdiff

      if (my_rank .eq. 0) then
         write (*, *) 'Particle volume', m_h**3
         write (*, *) 'Remeshing every ', rem_freq
         write (*, *) 'Diffusive radius ', Rd
         write (*, *) 'Delta t ', dt
         write (*, *) 'Circulation threshold ', thresh
      end if
      ! Setup itime to standard value, may change for ispecial=99 (if not: fine)
      itime = 0

      select case (ispecial)

         ! #### Vortex rings/wakes
      case (1, 2, 4)

         rmax = Lref / 2.d0
         nc = nv_on_Lref / 2
         rl = rmax / (2 * nc + 1)
         ns = 1 + 4 * nc * (nc + 1)
         nphi = nint(pi * r_torus / rl)
         n = 2 * ns * nphi
         np = ceiling(1.0 * n / n_cpu) + 1

         ! #### Sphere
      case (3)

         n = n_in
         np = ceiling(1.0 * n / n_cpu)
         h = sqrt(4.0 * pi / n)
         !force_const = force_const*h**3
         !eps = g*h

         ! #### Vortex wake
      case (5)

         h = 2 * pi / nc
         n = 2 * nc * ceiling(1.0 / h) * ceiling(2 * pi / h)
         np = ceiling(1.0 * n / n_cpu) + 1

         ! #### Single Vortex ring
      case (6)

         rmax = Lref / 2.d0
         nc = nv_on_Lref / 2
         rl = rmax / (2 * nc + 1)
         ns = 1 + 4 * nc * (nc + 1)
         nphi = nint(pi * r_torus / rl)
         n = ns * nphi
         np = ceiling(1.0 * n / n_cpu) + 1

      case (7)

         rmax = Lref / 2.d0
         block
            integer(kind_particle) :: ngrid(3)

            ngrid(1)= NINT((r_torus + rmax)/m_h) + 1
            ngrid(2)= NINT((r_torus + rmax)/m_h) + 1
            ngrid(3)= NINT(           rmax /m_h) + 1

            n = (2*ngrid(1) +1) * (2*ngrid(2) +1) * (2*ngrid(3) +1)
         end block

         np = ceiling(1.0 * n / n_cpu) + 1
         if(norm2(torus_offset).gt.0.d0) np = 2 * np

      case (98)

         n = n_in
         np = ceiling(1.0 * n / n_cpu)
         h = 1./n

         ! #### Setup MPI checkpoint readin
      case (99)

         call setup_MPI_IO_readin(itime)

      case default

         write (*, *) 'ERROR: need to specify setup via ispecial in your .h file'
         call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)

      end select

      ! Setup time variables
      trun = ts + itime * dt
      nt = ceiling((te - ts) / dt) ! Number of timesteps

      allocate (vortex_particles(np))

      if (my_rank .eq. 0) then
         write (*, *) "Starting PEPC-DVH with", n_cpu, " Processors, simulating", np, &
            " Particles on each Processor in", nt, "timesteps..."
         !write(*,*) "Using",num_threads,"worker-threads and 1 communication thread in treewalk on each processor (i.e. per MPI rank)"
         !write(*,*) "Maximum number of particles per work_thread = ", max_particles_per_thread
      end if

      ! Now finish by registering smaller datatype with MPI
      call register_short_mpi_types()

   end subroutine pepc_setup

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !>
   !>   Register shorter datatypes with MPI
   !>
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine register_short_mpi_types()
      use mpi
      integer(kind_default) :: ierr
      ! address calculation
      integer, dimension(1:nprops_particle_short) :: blocklengths, types
      integer(KIND=MPI_ADDRESS_KIND), dimension(1:nprops_particle_short) :: displacements
      integer(KIND=MPI_ADDRESS_KIND), dimension(0:nprops_particle_short) :: address
      ! dummies for address calculation
      type(t_particle_data_short) :: dummy_particle_data
      type(t_particle_short)      :: dummy_particle

      ! first register the interaction-specific short MPI type
      blocklengths(1:nprops_particle_data_short) = [3,1]
      types(1:nprops_particle_data_short) = [MPI_REAL8, MPI_REAL8]
      call MPI_GET_ADDRESS(dummy_particle_data,       address(0), ierr )        !&
      call MPI_GET_ADDRESS(dummy_particle_data%alpha, address(1), ierr )        !&
      call MPI_GET_ADDRESS(dummy_particle_data%vol  , address(2), ierr )        !&
      displacements(1:nprops_particle_data_short) = int(address(1:nprops_particle_data_short) - address(0))
      call MPI_TYPE_CREATE_STRUCT(nprops_particle_data_short, blocklengths, displacements, types, MPI_TYPE_PARTICLE_DATA_SHORT_sca, ierr)
      call MPI_TYPE_COMMIT(MPI_TYPE_PARTICLE_DATA_SHORT_sca, ierr)

      ! register short particle type
      blocklengths(1:nprops_particle_short) = [3, 1, 1, 1]
      types(1:nprops_particle_short) = [MPI_KIND_PHYSICS, MPI_REAL8, MPI_KIND_KEY, MPI_TYPE_PARTICLE_DATA_SHORT_sca]
      call MPI_GET_ADDRESS(dummy_particle,      address(0), ierr)               !&
      call MPI_GET_ADDRESS(dummy_particle%x,    address(1), ierr)               !&
      call MPI_GET_ADDRESS(dummy_particle%work, address(2), ierr)               !&
      call MPI_GET_ADDRESS(dummy_particle%key,  address(3), ierr)               !&
      call MPI_GET_ADDRESS(dummy_particle%data, address(4), ierr)               !&
      displacements(1:nprops_particle_short) = address(1:nprops_particle_short) - address(0)
      call MPI_TYPE_CREATE_STRUCT(nprops_particle_short, blocklengths, displacements, types, MPI_TYPE_PARTICLE_SHORT_sca, ierr)
      call MPI_TYPE_COMMIT(MPI_TYPE_PARTICLE_SHORT_sca, ierr)
   end subroutine register_short_mpi_types

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !>
   !> Subroutine to assign short t_particle to the full version
   !>
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elemental subroutine assign_particles(long_particle, short_particle)
      type(t_particle_short), intent(in) :: short_particle
      type(t_particle), intent(out)      :: long_particle

      long_particle = &                                                         !&
         t_particle( &                                                          !&
            short_particle%x, &                                                 !&
            short_particle%work, &                                              !&
            short_particle%key, &                                               !&
            0_kind_node, &                                                      !&
            0_kind_particle, &                                                  !&
            t_particle_data( &                                                  !&
               short_particle%data%alpha, &                                     !&
               [0._kind_physics, 0._kind_physics, 0._kind_physics], &           !&
               [0._kind_physics, 0._kind_physics, 0._kind_physics], &           !&
               [0._kind_physics, 0._kind_physics, 0._kind_physics], &           !&
               [0._kind_physics, 0._kind_physics, 0._kind_physics], &           !&
               short_particle%data%vol    &                                     !&
            ), &                                                                !&
            t_particle_results( &                                               !&
               [0._kind_physics, 0._kind_physics, 0._kind_physics], &           !&
               [0._kind_physics, 0._kind_physics, 0._kind_physics], &           !&
               [0._kind_physics, 0._kind_physics, 0._kind_physics], &           !&
               0._kind_physics &                                                !&
            ) &                                                                 !&
         )                                                                      !&

   end subroutine assign_particles

#ifdef USER_REDUCTION
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !>
   !> Function to add two short t_particles FOR AN OpenMP REDUCTION CLAUSE
   !>  This assumes/requires the position of the particles to be indentical
   !>
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   elemental function omp_reduction_add(omp_in1, omp_in2)
      type(t_particle_short), intent(in) :: omp_in1, omp_in2
      type(t_particle_short)             :: omp_reduction_add

      omp_reduction_add%x          =                      omp_in1%x             !&
      omp_reduction_add%work       = omp_in2%work       + omp_in1%work          !&
      omp_reduction_add%data%alpha = omp_in2%data%alpha + omp_in1%data%alpha    !&

   end function omp_reduction_add
#endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !>
   !>   Deallocate physvars arrays, copy back variables
   !>
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine cleanup(my_rank_l, n_cpu_l)
      use mpi
      implicit none

      integer, intent(in) :: my_rank_l ! MPI cpu rank
      integer, intent(in) :: n_cpu_l   ! MPI # CPUs

      ! copy call parameters to physvars module

      my_rank = my_rank_l
      n_cpu = n_cpu_l

      ! particle array deallocation in physvars

      deallocate (vortex_particles)

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
      integer :: ierr, err, fh, tmp_rem_freq, tmp_i, tmp_n_out
      integer(kind_particle) :: remain
      real ::tmp_dt, tmp_te, tmp_nu, tmp_ts, tmp_h, tmp_m_h
      real(kind_physics) :: tmp_thresh

      write (mpifile, '(a,i6.6,a)') "part_data/particle_", input_itime, ".mpi"
      call MPI_FILE_OPEN(MPI_COMM_WORLD, mpifile, IOR(MPI_MODE_RDWR, MPI_MODE_CREATE), MPI_INFO_NULL, fh, ierr)
      if (ierr .ne. MPI_SUCCESS) then
         write (*, *) 'something is wrong here: file open failed', my_rank, ierr, mpifile
         call MPI_ABORT(MPI_COMM_WORLD, err, ierr)
      end if

      ! Set file view to BYTE for header and read
      call MPI_FILE_SET_VIEW(fh, 0_MPI_OFFSET_KIND, MPI_BYTE, MPI_BYTE, 'native', MPI_INFO_NULL, ierr)
      call MPI_FILE_READ(fh, n, 1, MPI_KIND_PARTICLE, status, ierr)           ! # particles
      call MPI_FILE_READ(fh, tmp_dt, 1, MPI_REAL, status, ierr)               ! timestep
      call MPI_FILE_READ(fh, tmp_ts, 1, MPI_REAL, status, ierr)               ! Starting time
      call MPI_FILE_READ(fh, tmp_i, 1, MPI_INTEGER, status, ierr)             ! Last successful timestep (number)
      call MPI_FILE_READ(fh, tmp_n_out, 1, MPI_INTEGER, status, ierr)         ! Last successful time output (number)
      call MPI_FILE_READ(fh, tmp_te, 1, MPI_REAL, status, ierr)               ! Final time
      call MPI_FILE_READ(fh, tmp_nu, 1, MPI_REAL, status, ierr)               ! Viscousity
      call MPI_FILE_READ(fh, tmp_h, 1, MPI_REAL, status, ierr)                ! Original particle distance
      call MPI_FILE_READ(fh, tmp_m_h, 1, MPI_REAL, status, ierr)              ! Remeshing distance
      call MPI_FILE_READ(fh, tmp_rem_freq, 1, MPI_INTEGER, status, ierr)      ! Remeshing frequence
      call MPI_FILE_READ(fh, tmp_thresh, 1, MPI_KIND_PHYSICS, status, ierr)   ! threshold for pop. control
      call MPI_FILE_CLOSE(fh, ierr)

      dt = tmp_dt
      ts = tmp_ts
      nu = tmp_nu
      h = tmp_h
      m_h = tmp_m_h
      itime = tmp_i
      n_out = tmp_n_out
      rem_freq = tmp_rem_freq
      thresh = tmp_thresh

      np = n / n_cpu
      remain = n - np * n_cpu
      if ((remain .gt. 0) .and. (my_rank .lt. remain)) np = np + 1

   end subroutine setup_MPI_IO_readin

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !>
   !>  Set time stamp
   !>
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine stamp(istream, ibegin)
      implicit none
      integer, intent(in) :: istream, ibegin
      character :: cdate * 8, ctime * 10, czone * 5

      call DATE_AND_TIME(cdate, ctime, czone)

      if (ibegin .eq. 1) then
         write (istream, '(//a20,a12/a20,a12/a20,a12//)') 'PEPC run on ' &
            , cdate(7:8)//'/'//cdate(5:6)//'/'//cdate(1:4) &
            , 'Time: ', ctime(1:2)//':'//ctime(3:4)//':'//ctime(5:6), ' GMT+', czone
      else
         write (istream, '(a,a9)') 'Finished run at time: ', ctime(1:2)//':'//ctime(3:4)//':'//ctime(5:6)
      end if
   end subroutine stamp

end module physvars
