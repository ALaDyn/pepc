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
  
  
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine dump_trajectory()
    
    use physvars
    use utils
    implicit none
    
    integer :: p
    character(50) :: filename

    if(my_rank == 0) write(*,*) "benchmarking: dump_trajectory"
    
    filename = "trajectory.dat"
    open(91, file=filename, STATUS='REPLACE')
    write(91,*) "# particle positions for geom ", ispecial, " at timestep ", nt, ": p x y z"
    do p=1, 4
       write(91,*) p, x(p), y(p), z(p)
    end do
    close(91)
    
  end subroutine dump_trajectory

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write_particles(step)
    
    use physvars
    use utils
    implicit none

    include 'mpif.h'
    
    integer, intent(in) :: step
    integer :: p, ierr
    real*8 :: writesize, t1, t2, t3, t4
    character(50) :: filename

    !! sionlib related variables
    integer   :: sid, fsblksize
    integer*8 :: chuncksize, left, dumpsize, bwrote
    

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
    
    write(*,*) "pepc timing - pre: ", time_pre - time_start
    write(*,*) "pepc timing - inner: ", time_inner - time_pre
    write(*,*) "pepc timing - post: ", time_post - time_inner
    write(*,*) "pepc timing - total: ", time_post - time_start

  end subroutine benchmark_end
  
  
end module benchmarking
