#!/bin/bash

module purge

# set up different module environmets, possibly set EVN variables
#---------------
if [ "$1" == "intel_intel" ]; then
   module load Intel/2017.0.098-GCC-5.4.0 IntelMPI/2017.0.098
   test -f makefile.defs_intel && ln -s makefile.defs_intel makefile.defs
   make allclean
   make pepc-benchmark
   export I_MPI_OFA_TRANSLATION_CACHE=0
   export I_MPI_DAPL_TRANSLATION_CACHE=0
   export I_MPI_FABRIC=shm:ofa
fi

#---------------
if [ "$1" == "gcc_mvapich" ]; then 
   module load GCC/5.4.0 MVAPICH2/2.2-GDR
   test -f makefile.defs_gcc && ln -s makefile.defs_gcc makefile.defs
   make allclean
   make pepc-benchmark
fi

#---------------
if [ "$1" == "intel_mvapich" ]; then 
   module load Intel/2016.4.258-GCC-5.4.0 MVAPICH2/2.2-GDR
   test -f makefile.defs_intel && ln -s makefile.defs_intel makefile.defs
   make allclean
   make pepc-benchmark
fi

#---------------
if [ "$1" == "gcc_ps" ]; then 
   module load GCC/5.4.0 ParaStationMPI/5.1.5-1
   module load pscom/.5.1.2-1
   test -f makefile.defs_gcc && ln -s makefile.defs_gcc makefile.defs
   make allclean
   make pepc-benchmark
fi

#---------------
if [ "$1" == "gcc_pssilly" ]; then 
   module load GCC/5.4.0 ParaStationMPI/5.1.5-1
   test -f makefile.defs_gcc && ln -s makefile.defs_gcc makefile.defs
   make allclean
   make pepc-benchmark
   export PSP_RENDEZVOUS_OPENIB=-1
fi

#---------------
if [ "$1" == "gcc_ompi" ]; then
   module use /usr/local/software/$(cat /etc/FZJ/systemname)/OtherStages
   module load Stages/2016b
   module use ~broemmel/$(cat /etc/FZJ/systemname)/OpenMPI/module/
   module load GCC/5.4.0 OpenMPI/2.1.0
   test -f makefile.defs_gcc && ln -s makefile.defs_gcc makefile.defs
   make allclean
   make pepc-benchmark
fi

#---------------
if [ "$1" == "gcc_ompi_nb" ]; then 
   module use /usr/local/software/$(cat /etc/FZJ/systemname)/OtherStages
   module load Stages/2016b
   module use ~broemmel/$(cat /etc/FZJ/systemname)/OpenMPI/module/
   module load GCC/5.4.0 OpenMPI
   test -f makefile.defs_gcc && ln -s makefile.defs_gcc makefile.defs
   make allclean
   make pepc-benchmark
fi

#---------------
if [ "$1" == "gcc_ps_pin" ]; then 
   module load GCC/5.4.0 ParaStationMPI/5.1.5-1
   module load pscom/.5.1.2-1
   test -f makefile.defs_gcc && ln -s makefile.defs_gcc makefile.defs
   make allclean
   make pepc-benchmark
fi

#---------------
if [ "$1" == "gcc_ps_new" ]; then 
   module load GCC/5.4.0 ParaStationMPI/5.1.5-1
   module load pscom/.5.1.2-1
   test -f makefile.defs_gcc && ln -s makefile.defs_gcc makefile.defs
   make allclean
   make pepc-benchmark
fi

#---------------
if [ "$1" == "gcc_ps_new_malloc" ]; then 
   module load GCC/5.4.0 ParaStationMPI/5.1.5-1
   module load pscom/.5.1.2-1
   test -f makefile.defs_gcc && ln -s makefile.defs_gcc makefile.defs
   make allclean
   make pepc-benchmark
   export MALLOC_ARENA_MAX=1
fi

#---------------
if [ "$1" == "gcc_ps_old" ]; then 
   module load GCC/5.4.0 ParaStationMPI/5.1.5-1
   test -f makefile.defs_gcc && ln -s makefile.defs_gcc makefile.defs
   make allclean
   make pepc-benchmark
fi

#---------------
if [ "$1" == "gcc_ps_old_rdzv" ]; then 
   module load GCC/5.4.0 ParaStationMPI/5.1.5-1
   test -f makefile.defs_gcc && ln -s makefile.defs_gcc makefile.defs
   make allclean
   make pepc-benchmark
   export PSP_RENDEZVOUS_OPENIB=-1
fi
