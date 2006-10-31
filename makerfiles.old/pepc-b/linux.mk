#  --------------------------------------------------------
#
#  Makefile Compiler and Library options for Linux machines
#
#  ---------------------------------------------------------
 
#  Compilers including MPICH wrappers
# full MPI path usually something like /usr/local/mpich/bin

CC          = mpicc
CCC         = mpicxx
F77         = mpif77
CLINKER     = mpicc
CCLINKER    = mpicxx
FLINKER     = mpif77
F90         = mpif90
F90LINKER   = mpif90	  
FC = $(F90)
CPP = /usr/bin/cpp

# KOJAK instrumentation
#FC = mpxlf90_r -qdebug=function_trace
#FC = kinst mpxlf90_r
#FC = kinst-pomp -rcfile $(PREFIX)/pepc-b/opari.rc -- mpxlf90_r

#  Library Archiver commands
RANLIB  = ranlib
AR      = ar

#  Put compiler optimisation in here
QTUNE = -O3 
CFLAGS1= -O3 -g -I/usr/local/include 
FFLAGS1 = $(QTUNE)  
#  Debug options
DB = -g -check bounds -traceback


#   	LIBRARIES
#   ---------------------

#  PEPC library location

LIBPEPC = -L../lpepcsrc -llpepc


#  IBM library for timing routines

#IBMLIB = -lxlf90

#  MPI libraries - should be contained in wrappers 

LIBS_MPI = 
LIBSF_MPI = 

MPITRACE =
#MPITRACE = -L/bgl/local/lib -lmpihpm_f -lbgl_perfctr.rts #for MPI+HPM information


# Visit libraries

# Setup flags for C-preprocessor
# Use VISIT routines with XNBODY visualisation
PREPROC = -DVISIT_NBODY
#PREPROC =

VPREFIX=/tmp/xnbody
APISDIR=$(VPREFIX)/apis

VISITLIBS=-L$(VPREFIX)/xnbody_qt4/apis/nbody3 -llvisit_nbody2 -L$(VPREFIX)/lvisit/lib -llvisit -L$(VPREFIX)/visit_dist/lib -lvisit
#NETCDFLIB = -lnetcdf


