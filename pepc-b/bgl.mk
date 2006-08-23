#  --------------------------------------------------------
#
#  Makefile Compiler and Library options for JUBL (BlueGene/L)
#
#  ---------------------------------------------------------
 
#  Compilers

BGLSYS = /bgl/BlueLight/ppcfloor/bglsys
BGLLOC = /bgl/local

CC = /usr/local/bin/mpcc
CCP = /usr/local/bin/mpCC
FC = /usr/local/bin/mpxlf90

# KOJAK instrumentation
#FC = mpxlf90_r -qdebug=function_trace
#FC = kinst mpxlf90_r
#FC = kinst-pomp -rcfile $(PREFIX)/pepc-b/opari.rc -- mpxlf90_r

#  Library Archiver commands
RANLIB  = ranlib
AR      = ar


# try -qarch=440 first, then use -qarch=440d for 2nd FPU later on
#  (SIMDization requires at least -O3)
# use -qlist -qsource with 440d and look for Parallel ASM instructions.
#
#IPA=-qipa=inline=key2addr -qipa=inline=make_hashentry -qipa=inline=key2node -qipa=inline=next_node
#QTUNE = -O5 -qtune=440 -qarch=440d ## $(IPA)
QTUNE = -O3 -qtune=440 -qarch=440d ## $(IPA)
CFLAGS1= -O3 -g -I/bgl/local/include -qarch=440 -qtune=440
BGLFLAGS= -I$(BGLSYS)/include -I$(BGLLOC)/include -L$(BGLSYS)/lib  $(QTUNE) 
FFLAGS1 = $(BGLFLAGS) -qsuffix=f=f90 -qsuffix=cpp=F  -qnosave 
DB= -g -qfullpath ##-qcheck
#DB= -g
#

#   	LIBRARIES
#   ---------------------

#  Where the PEPC library sits

LIBPEPC = -L../lpepcsrc -llpepc

#  IBM library for timing routines

IBMLIB = -lxlf90

#  BGL MPI libraries 

LIBS_MPI = -lmpich.rts -lmsglayer.rts -lrts.rts -ldevices.rts
LIBSF_MPI = -lmpich.rts -lfmpich.rts -lmsglayer.rts -lrts.rts -ldevices.rts

#MPITRACE = -L/bgl/local/lib -lmpitrace_f  #for MPI information only
MPITRACE = -L/bgl/local/lib -lmpihpm_f -lbgl_perfctr.rts #for MPI+HPM information
#MPITRACE = -qdebug=function_trace -L/bgl/local/lib -lmpihpm_f -lbgl_perfctr.rts #for per function MPI+HPM information


# Visit libraries

# Setup flags for C-preprocessor
# Use VISIT routines with XNBODY visualisation
PREPROC = -WF,-DVISIT_NBODY,-DNETCDFLIB
#PREPROC = -WF,-DVISIT_NBODY
#PREPROC =

NETCDFLIB = -L/bgl/local/lib -lnetcdf
NCOBJS=ncnbody.o

VISITDIR=/bgl/local/visit.rts
VISITLIBS= -L$(VISITDIR)/lvisit/lib -llvisit -L$(VISITDIR)/lib -lvisit

XNBODYLIBS=-L$(VISITDIR)/apis/nbody3 -llvisit_nbody2


