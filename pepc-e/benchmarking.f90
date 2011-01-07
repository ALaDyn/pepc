module benchmarking

  !  ------------------ PAPI ------------------
#ifdef PEPC_USE_PAPI_GLOBAL
  
#include "f90papi.h"
  
#define PAPI_CNT_NUM 1
  
  integer                      :: papi_events(PAPI_CNT_NUM)
  integer*8                    :: papi_values(PAPI_CNT_NUM)
  character*(PAPI_MAX_STR_LEN) :: papi_name
  integer                      :: papi_check, papi_cnt
  real*8                       :: papi_time_start, papi_time_stop
#endif
  !  ------------------ PAPI ------------------

  real*8 :: time_start, time_pre, time_inner, time_post
  
  integer, private, parameter :: NUM_PARTICLES_FRONT  = 5 !< number of particles from beginning of particle list to use in dump routines
  integer, private, parameter :: NUM_PARTICLES_MID    = 5 !< number of particles from center of particle list to use in dump routines
  integer, private, parameter :: NUM_PARTICLES_BACK   = 5 !< number of particles from end of particle list to use in dump routines
  integer, private, parameter :: NUM_DIAG_PARTICLES = NUM_PARTICLES_FRONT + NUM_PARTICLES_MID + NUM_PARTICLES_BACK !< total number of particles to use in dump routines
  integer, private, parameter :: NUM_DIAG_PROPS     = 12  !< number of properties to be collected for diagnostic purposes (in addition to workload)
  real*8, private :: diag_props(NUM_DIAG_PARTICLES,NUM_DIAG_PROPS)
  integer, private :: diag_work(NUM_DIAG_PARTICLES,1)
  
contains


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>
  !> routine for selecting particles that are used for diagnostic output
  !> for given input idx=1:NUM_DIAG_PARTICLES, it returns a particle number from
  !> the beginning, the center, or the end of the local particle list
  !> the number of particles can be adjusted by modifying the constants
  !> NUM_PARTICLES_FRONT, NUM_PARTICLES_MID, and NUM_PARTICLES_BACK
  !>
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function diagnostic_particle(idx)
    use physvars
    implicit none

    integer, intent(in) :: idx !< index of diagnostic particle, possible input range: (1:NUM_DIAG_PARTICLES)
    integer :: diagnostic_particle

    if (npart_total <= NUM_DIAG_PARTICLES) then
      diagnostic_particle = min(idx, npart_total)
    else
      select case (idx)
        case (1:NUM_PARTICLES_FRONT)
                    diagnostic_particle = idx
        case (NUM_PARTICLES_FRONT+1:NUM_PARTICLES_FRONT+NUM_PARTICLES_MID)
                    diagnostic_particle = npart_total/2 - ((NUM_PARTICLES_MID/2+NUM_PARTICLES_FRONT)-idx)
        case (NUM_PARTICLES_FRONT+NUM_PARTICLES_MID+1:NUM_DIAG_PARTICLES)
                    diagnostic_particle = npart_total-(NUM_DIAG_PARTICLES-idx)
      end select
    end if

  end function

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>
  !>
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine gather_particle_diag()
    use physvars
    implicit none
    include 'mpif.h'

    integer :: r, i, ierr, target_rank, target_particle, target_particle_local
    integer :: fances(0:n_cpu-1)
    real*8  :: diag_props_buf(NUM_DIAG_PROPS)
    integer :: diag_work_buf(1)

    logical :: debug=.false., debug_root

    debug_root = (my_rank.eq.0) .and. debug

    if (debug) write(*,*) "start mpi scan on rank ", my_rank

    call MPI_SCAN(np_local, fances(my_rank), 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER(MPI_IN_PLACE, 0, 0, fances, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

    if(my_rank.eq.0 .and. debug) then
       do i=0, n_cpu-1
          write(*,*) "particle rank fance: rank=", i, " partial sum=", fances(i)
       end do
    end if

    do i=1, NUM_DIAG_PARTICLES
       target_particle = diagnostic_particle(i)
       if(debug_root) write(*,*) "gather diag information for particle ", target_particle
       do r=0, n_cpu-1
          if(target_particle.le.fances(r)) then
             target_rank = r
             exit
          end if
       end do
       if(debug_root) write(*,*) "particle is located on rank ", target_rank

       if(my_rank.eq.target_rank) then
          if(debug) write(*,'(a,4i12)') "particle: ", my_rank, target_particle, fances(my_rank), np_local
          target_particle_local = target_particle-fances(my_rank)+np_local
          if(debug) write(*,'(a,i12,a,i12,a,i12)') "from rank", my_rank, " particle is here, with label ", &
               pelabel(target_particle_local), " local target ", target_particle_local

          diag_props_buf( 1: 3) = [  x(target_particle_local),  y(target_particle_local),  z(target_particle_local)]
          diag_props_buf( 4: 6) = [ ux(target_particle_local), uy(target_particle_local), uz(target_particle_local)]
          diag_props_buf( 7: 8) = [  q(target_particle_local),  m(target_particle_local)]
          diag_props_buf(    9) =  pot(target_particle_local)
          diag_props_buf(10:12) = [ ex(target_particle_local), ey(target_particle_local), ez(target_particle_local)]

          diag_work_buf(1) = int(work(target_particle_local))

          if(debug) write(*,*) "from rank", my_rank, " sending pos_vel ", diag_props_buf

          if(my_rank .ne. 0) then
             call MPI_SEND(diag_props_buf,12, MPI_REAL8,   0, 0, MPI_COMM_WORLD, ierr)
             call MPI_SEND(diag_work_buf,  1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, ierr)
          end if
       end if

       if(my_rank.eq.0 .and. target_rank.ne.0) then
          call MPI_RECV(diag_props_buf,12, MPI_REAL8,  target_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
          call MPI_RECV(diag_work_buf,  1, MPI_INTEGER,target_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
          if(debug) write(*,*) "from rank", my_rank, " received pos_vel ", diag_props_buf
       end if
       
       call MPI_BARRIER(MPI_COMM_WORLD, ierr)

       diag_props(i,:) = diag_props_buf(:)
       diag_work(i,1) = diag_work_buf(1)

    end do


  end subroutine gather_particle_diag


  subroutine benchmarking_dump_diagnostics()
    implicit none

    call dump_trajectory()
    call dump_fields()

  end subroutine benchmarking_dump_diagnostics


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>
  !> output of position of selected particles to file trajectory.dat
  !> for diagnostic purposes.
  !> should only be called once per simulation run
  !>
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine dump_trajectory()
    use physvars
    implicit none
    integer :: p, i

    if(my_rank == 0) write(*,*) "benchmarking: dump_trajectory"
    
    open(91, file="trajectory.dat", STATUS='REPLACE')
    write(91,'(a,i12)') "# npart_total ", npart_total
    write(91,'(a,i12)') "# npart_diag ", NUM_DIAG_PARTICLES
    write(91,'(a, i6, a, i6, a)') "# particle positions for geom ", ispecial, " at timestep ", nt, ": p x y z ux uy uz work"
    do i=1,NUM_DIAG_PARTICLES
       p = diagnostic_particle(i)
       write(91,'(i12,6e20.12,i12)') p, diag_props(i,1), diag_props(i,2), diag_props(i,3), &
            diag_props(i,4), diag_props(i,5), diag_props(i,6), diag_work(i,1)
       
    end do
    close(91)
    
  end subroutine dump_trajectory


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>
  !> output of charge, potential and electric field to file fielddump.dat
  !> for diagnostic purposes.
  !> should only be called once per simulation run
  !>
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine dump_fields()
    use physvars
    use module_fmm_framework
    implicit none
    integer :: p, i

    if(my_rank == 0) write(*,*) "benchmarking: dump_fields"

    open(91, file="fielddump.dat", STATUS='REPLACE')
    write(91,'(a,i12)') "# npart_total ", npart_total
    write(91,'(a,i12)') "# npart_diag ", NUM_DIAG_PARTICLES
    write(91,'(a,e20.12)') "# pot_farfield ", potfarfield
    write(91,'(a,e20.12)') "# pot_nearfield ", potnearfield
    write(91,'(a, i6, a, i6, a)') "# particle charge and field for geom ", ispecial, " at timestep ", nt, ": p q pot ex ey ez"
    do i=1,NUM_DIAG_PARTICLES
       p = diagnostic_particle(i)
       write(91,'(i12,5e20.12)') p, diag_props(i,7), diag_props(i,9:12)
    end do
    close(91)

  end subroutine dump_fields


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>
  !> output of number of interactions of selected particles and in total
  !> into file num_interactions.dat for diagnostic purposes.
  !>
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine dump_num_interactions()
    use physvars
    implicit none
    integer :: i
    logical, save :: print_header = .true.

    if(my_rank == 0) write(*,*) "benchmarking: dump_num_interactions"

    open(91, file="num_interactions.dat", STATUS='UNKNOWN', position='APPEND')

    if (print_header) then
      write(91,*) "# itime, totalnumber, p=", ( diagnostic_particle(i), i=1,NUM_DIAG_PARTICLES )
      print_header = .false.
    end if

    write(91,*) itime, int(sum(work)), ( int(work(diagnostic_particle(i))), i=1,NUM_DIAG_PARTICLES )
    close(91)

  end subroutine dump_num_interactions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write_particles(step)
    
    use physvars
    implicit none

    include 'mpif.h'
    
    integer, intent(in) :: step
    integer :: p, ierr
    real*8 :: writesize, t1, t2, t3, t4
    character(50) :: filename

#ifdef WRITE_PARTICLE_SIONLIB
    !! sionlib related variables
    integer   :: sid, fsblksize
    integer*8 :: chuncksize, left, dumpsize, bwrote
#endif

    if(my_rank == 0) write(*,*) "benchmarking: write particles to file"


!!! write particle data as a binary file
#ifdef WRITE_PARTICLE_BINARY

    if(my_rank == 0) write(*,*) "benchmarking: write particles in binary mode"

    write(filename,'(a,i6.6,a,i6.6,a)') "particle_", step, "_", my_rank, ".dat"

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    t1 = MPI_WTIME()
    open(91, file=filename, STATUS='REPLACE', ACCESS="STREAM")

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    t2 = MPI_WTIME()
    do p=1, np_local
       write(91) x(p), y(p), z(p), ux(p), uy(p), uz(p)
    end do

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    t3 = MPI_WTIME()
    close(91)

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    t4 = MPI_WTIME()

    writesize = np_local*6*kind(x(1))/1024.0/1024.0

#else

!!! write particle data using the sion library
#ifdef WRITE_PARTICLE_SIONLIB

    if(my_rank == 0) write(*,*) "benchmarking: write particles with sionlib"

    write(filename,'(a,i6.6,a,i6.6,a)') "particle_", step, ".dat"

    dumpsize = np_local*kind(x(1))

    chuncksize = 6*dumpsize
    fsblksize  = 2*1024*1024

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    t1=MPI_WTIME()
    call fsion_paropen_mpi(trim(filename), "bw", MPI_COMM_WORLD, chuncksize, fsblksize, my_rank, sid)     

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    t2=MPI_WTIME()
    
    call fsion_write(x,  1, dumpsize, sid, bwrote)
    call fsion_write(y,  1, dumpsize, sid, bwrote)
    call fsion_write(z,  1, dumpsize, sid, bwrote)
    call fsion_write(ux, 1, dumpsize, sid, bwrote)
    call fsion_write(uy, 1, dumpsize, sid, bwrote)
    call fsion_write(uz, 1, dumpsize, sid, bwrote)

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    t3=MPI_WTIME()

    call fsion_parclose_mpi(sid, MPI_COMM_WORLD, ierr)

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    t4=MPI_WTIME()

    writesize=6*dumpsize/1024.0/1024.0

#else
!!! write particle date as a text file    
    
    if(my_rank == 0) write(*,*) "benchmarking: write particles in text mode"

    write(filename,'(a,i6.6,a,i6.6,a)') "particle_", step, "_", my_rank, ".dat"

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    t1=MPI_WTIME()
    open(91, file=filename, STATUS='REPLACE')

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    t2=MPI_WTIME()
    write(91,*) "# particle positions for geom ", ispecial, " at timestep ", nt, ": p x y z vx vy vz"
    do p=1, np_local
       write(91,'(I6.6,6E10.3E2)') p, x(p), y(p), z(p), ux(p), uy(p), uz(p)
    end do

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    t3=MPI_WTIME()
    close(91)

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    t4=MPI_WTIME()

    writesize = np_local*(6*10 + 7)/1024.0/1024.0

#endif
#endif 


    if(my_rank == 0) write(*,*) "benchmarking: rank 0 has written ", writesize, "MB to file"
    if(my_rank == 0) write(*,*) "benchmarking: rank 0 time to open file ", t2-t1, "s"
    if(my_rank == 0) write(*,*) "benchmarking: rank 0 time to write file ", t3-t2, "s"
    if(my_rank == 0) write(*,*) "benchmarking: rank 0 time to close file ", t4-t3, "s"    

    if(my_rank == 0) write(*,*) "benchmarking: rank 0 io performance ", writesize/(t4-t1), "MB/s"

end subroutine write_particles


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  subroutine benchmark_pre

    use physvars
    implicit none

    include 'mpif.h'

    if(my_rank == 0) write(*,*) "benchmarking: pre"

    time_start = MPI_WTIME()
    
  end subroutine benchmark_pre
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine benchmark_inner
    
    use physvars
    implicit none
    
    include 'mpif.h'

    if(my_rank == 0) write(*,*) "benchmarking: inner"


    !  ------------------ PAPI ------------------
#ifdef PEPC_USE_PAPI_GLOBAL
  
    papi_events = (/ PAPI_FP_OPS /)
    
    papi_check = PAPI_VER_CURRENT
    call PAPIF_library_init(papi_check)
    if (papi_check .ne. PAPI_VER_CURRENT) then
       stop "PAPI: error library init"
    endif
    write(*,*) "PAPI: library initialized"

    papi_time_start = MPI_WTIME()
    
    call PAPIF_start_counters(papi_events, PAPI_CNT_NUM, papi_check)
    if (papi_check .ne. PAPI_OK) then
       stop "PAPI: error starting hardware counter"
    endif
    
    write(*,*) "PAPI: start hardware counter measurement"
    
#endif
    !  ------------------ PAPI ------------------
    
    time_pre = MPI_WTIME()
    
  end subroutine benchmark_inner
  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  subroutine benchmark_post
    
    use physvars
    implicit none

    include 'mpif.h'

    if(my_rank == 0) write(*,*) "benchmarking: post"

    time_inner = MPI_WTIME()
    
    !  ------------------ PAPI ------------------
#ifdef PEPC_USE_PAPI_GLOBAL
    
    write(*,*) "PAPI: stop hardware counter"
    call PAPIF_stop_counters(papi_values, PAPI_CNT_NUM, papi_check)
    if (papi_check .ne. PAPI_OK) then
       stop "PAPI: error stopping hardware counter"
    endif
    
    papi_time_stop = MPI_WTIME()
    
    do papi_cnt = 1,PAPI_CNT_NUM
       call PAPIF_event_code_to_name(papi_events(papi_cnt), papi_name, papi_check)
       if (papi_check .ne. PAPI_OK) then
          stop "PAPI: error reading hardware counter name"
       endif
       
       write(*,*) "PAPI: counter ", papi_name, "; value ", papi_values(papi_cnt)
    end do
    
    
    write(*,*) "PAPI: performance [MF/s] ", (papi_values(1))/(papi_time_stop-papi_time_start)/1.0e6
    
#endif
    !  ------------------ PAPI ------------------
    
  end subroutine benchmark_post
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine benchmark_end

    use physvars
    implicit none

    include 'mpif.h'

    if(my_rank == 0) write(*,*) "benchmarking: end"

    time_post = MPI_WTIME()
    
    if(my_rank == 0) then
       write(*,*) "pepc timing - pre: ", time_pre - time_start
       write(*,*) "pepc timing - inner: ", time_inner - time_pre
       write(*,*) "pepc timing - post: ", time_post - time_inner
       write(*,*) "pepc timing - total: ", time_post - time_start
    end if

  end subroutine benchmark_end
  
  
end module benchmarking
