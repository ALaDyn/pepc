#  --------------------------------------------------------
#
#  Makefile Compiler and Library options for JUBL (BlueGene/L)
#
#  ---------------------------------------------------------
 
#  Compilers

BGLSYS = /bgl/BlueLight/ppcfloor/bglsys

CC = /opt/ibmcmp/vac/7.0/bin/blrts_xlc
CPPC = /opt/ibmcmp/vacpp/7.0/bin/blrts_xlc
FC = /opt/ibmcmp/xlf/9.1/bin/blrts_xlf90


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
#QTUNE = -O3  -qtune=440 -qarch=440d
CFLAGS1= -O3 -g -I/opt/ibmcmp/vac/7.0/include -I/usr/include -I$(BGLSYS)/include -L$(BGLSYS)/lib -qarch=440 -qtune=440
BGLFLAGS= -I$(BGLSYS)/include -L$(BGLSYS)/lib $(QTUNE)
FFLAGS1 = $(BGLFLAGS) -qsuffix=f=f90 -qsuffix=cpp=F90  -qnosave
IPA=-qipa=inline=key2addr -qipa=inline=make_hashentry -qipa=inline=key2node -qipa=inline=next_node
DB= -g -qfullpath
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

MPITRACE = -L/bgl/local/lib -lmpitrace_f  #for MPI information only
#MPITRACE = -L/bgl/local/lib -lmpihpm_f -lbgl_perfctr.rts #for MPI+HPM information


# Visit libraries

# Setup flags for C-preprocessor
# Use VISIT routines with XNBODY visualisation
#PREPROC = -WF,-DVISIT_NBODY
PREPROC =

#NETCDFLIB = -lnetcdf
NETCDFLIB=
#VISITLIBS=-L/usr/local/beta/visit-2.0b/lvisit/apis/spk4 -llvisit_spk -L/usr/local/beta/visit-2.0b/lvisit/lib -llvisit -L/usr/local/beta/visit-2.0b/lib -lvisit
#VISITLIBS=-L/usr/local/beta/visit-2.0b/lvisit/apis/spk5  -llvisit_spk -L/usr/local/beta/visit-2.0b/lvisit/lib -llvisit -L/usr/local/beta/visit-2.0b/lib -lvisit -bloadmap:aa

#XNBODYLIBS=-L/usr/local/beta/visit-2.0b/lvisit/apis/nbody3  -llvisit_nbody2


