#  --------------------------------------------------------
#
#  Makefile Compiler and Library options for JUMP (IBM-p690)
#
#  ---------------------------------------------------------


#  Compilers

FC = mpxlf90_r
CPP = /usr/lib/cpp 
CC = xlc_r

# For KOJAK instrumentation
#FC = mpxlf90_r -qdebug=function_trace
#FC = kinst mpxlf90_r
#FC = kinst-pomp -rcfile $(PREFIX)/pepc-b/opari.rc -- mpxlf90_r

#  Library Archiver commands
RANLIB  = ranlib
AR      = ar -X64
#AR      = ar 

# Setup flags for C-preprocessor
# Use VISIT routines with XNBODY visualisation

#PREPROC = -WF,-DVISIT_NBODY
PREPROC = -DVISIT_NBODY



#  Profiler information
#PG=-pg -g -qfullpath
#  Debug mode
#DB= -g -qfullpath  
DB= -g -qfullpath  
#DB= -qcheck -qstrict -g -qfullpath
#  Hardware performance counter
HPM = -lhpm
#  Compiler listing
#LISTING=-qreport -qlist -qlistopt -qsource 

FFLAGS1 = -q64 -qsuffix=f=f90:cpp=F90 -qnosave 
#FFLAGS1 = -qsuffix=f=f90:cpp=F90 -qnosave 
IPA= -qipa=inline=key2addr_db -qipa=inline=key2addr -qipa=inline=make_hashentry -qipa=inline=key2node -qipa=inline=next_node
#TUNE= -qarch=pwr4 -qtune=pwr4 -O4 ## $(IPA) 
CFLAGS1 = -q64 -O3 -I/usr/local/include

#  Symbol tables
LMAP = -bnoquiet


#  Auto dependency command for f90
F90DEP=./f90depend -u -I/usr/lpp/ppe.poe/include/thread *.f90 *.F90


#   	LIBRARIES
#   ---------------------


#  Where the PEPC library sits

LIBPEPC = -L../lpepcsrc -llpepc


#  MPI libs contained in mpxlf compiler wrappers

LIBS_MPI =
LIBSF_MPI =

#  IBM library for timing routines

IBMLIB = -lxlf90


#  MPI performance library

#MPITRACE= -L/usr/local/beta/lib -lmpitrace
#MPITRACE= -L/usr/local/beta/lib -lmpiprof
#MPITRACE=  -L/usr/local/beta/lib -lmpihpm -lpmapi
#MPITRACE= -lmpitrace
#MPITRACE= -lmpiprof
#MPITRACE= -lsummary -lpmapi
#MPITRACE=  -lmpihpm -lpmapi


