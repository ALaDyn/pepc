# 64-bit hardware

#
#  Makefile for treemp and test programs
#
PREFIX = ..


#DBFLAGS = -g -Rabp -ev  -q64 -qrealsize=8  -qsuffix=f=f90 -qcheck
#LDFLAGS = -q64 -qrealsize=8  -g 
# For apprentice
#FFLAGS =   -eA -m4 -O2,inline4
#LDFLAGS = -lpmpi -Xm -lapp

## For IBM AIX compiler
#FFLAGS = -g -q64 -qrealsize=8 -qsuffix=f=f90 -qcheck  
## Optimised
#FFLAGS =  -q64 -qrealsize=8 -qsuffix=f=f90  -O4 -Q -qipa -qarch=pwr4 -qstrict -qreport  -qlist -qsource -I/usr/local/include 


IPA=-qipa=inline=key2addr -qipa=inline=make_hashentry -qipa=inline=key2node -qipa=inline=next_node

#PG=-pg -g -qfullpath
#DB= -g -qfullpath -qcheck 
DB= -g -qfullpath 
HPM = -lhpm
LISTING=-qreport -qlist -qlistopt -qsource 
#MPITRACE= -L/usr/local/beta/lib -lmpitrace
#MPITRACE= -L/usr/local/beta/lib -lmpiprof
#MPITRACE=  -L/usr/local/beta/lib -lmpihpm -lpmapi
#MPITRACE= -lmpitrace
#MPITRACE= -lmpiprof
#MPITRACE= -lsummary -lpmapi
#MPITRACE=  -lmpihpm -lpmapi
#Q64=-q64
FFLAGS1 = -qrealsize=8 -qsuffix=f=f90 -qsuffix=cpp=F90 -qnosave 
TUNE= -qarch=pwr4 -qtune=pwr4 -qhot ## -O5 $(IPA) 
# $(PG) $(HPM)

#FFLAGS =  $(FFLAGS1) $(DB) $(Q64) $(PG) -qhot
#FFLAGS =  $(FFLAGS1) $(DB) $(Q64) $(PG) 
FFLAGS =  $(FFLAGS1) $(DB) $(Q64) $(TUNE)
TFLAGS = $(FFLAGS) $(TUNE)

#ANALYS = -lpmpi -Xm -lpat pat.cld

# Visit libraries
# JUMP
#NETCDFLIB = -lnetcdf
#VISITLIBS=-L/usr/local/beta/visit-2.0b/lvisit/apis/spk4 -llvisit_spk -L/usr/local/beta/visit-2.0b/lvisit/lib -llvisit -L/usr/local/beta/visit-2.0b/lib -lvisit
#VISITLIBS=-L/usr/local/beta/visit-2.0b/lvisit/apis/spk5  -llvisit_spk -L/usr/local/beta/visit-2.0b/lvisit/lib -llvisit -L/usr/local/beta/visit-2.0b/lib -lvisit -bloadmap:aa

#VISITLIBS=-L/home4/jzam04/jzam0401/tree/pepc/src  -llvisit_spk -L/usr/local/beta/visit-2.0b/lvisit/lib -llvisit -L/usr/local/beta/visit-2.0b/lib -lvisit -bloadmap:aa
#XNBODYLIBS=-L/usr/local/beta/visit-2.0b/lvisit/apis/nbody2  -llvisit_nbody2
LMAP = -bnoquiet
LDFLAGS =  $(XNBODYLIBS) $(VISITLIBS) $(NETCDFLIB) 
CFLAGS = -I/usr/local/include
RANLIB  = ranlib
libdir  = .
AR      = ar -X64
LIBPEPC = -L../lpepcsrc -llpepc
IBMLIB = -lxlf90
pepclib = liblpepc.a



#FC = mpxlf90_r
CC = xlc

# KOJAK instrumentation
#FC = mpxlf90_r -qdebug=function_trace
#FC = kinst-pomp -rcfile $(PREFIX)/pepc-b/opari.rc -- mpxlf90_r

# -------------------------------------------------------------------------------------------
# BlueGene Light
# -------------------------------------------------------------------------------------------
BGLSYS = /bgl/BlueLight/ppcfloor/bglsys

CC = /opt/ibmcmp/vac/7.0/bin/blrts_xlc
CPPC = /opt/ibmcmp/vacpp/7.0/bin/blrts_xlc
FC = /opt/ibmcmp/xlf/9.1/bin/blrts_xlf90

# try -qarch=440 first, then use -qarch=440d for 2nd FPU later on
#  (SIMDization requires at least -O3)
# use -qlist -qsource with 440d and look for Parallel ASM instructions.
#
QTUNE = -qarch=440d -qtune=440 
CFLAGS= -O3 -g -I$(BGLSYS)/include -L$(BGLSYS)/lib -qarch=440 -qtune=440
FFLAGS= -O -g $(FFLAGS1) -I$(BGLSYS)/include -L$(BGLSYS)/lib ## $(QTUNE)
#DBFLAGS= -qcheck -qfullpath
DBFLAGS= -qfullpath
#
#MPITRACE=  -L/usr/local/beta/lib -lmpihpm -lpmapi
#MPITRACE= -lmpitrace
#MPITRACE= -lmpiprof
#MPITRACE= -lsummary -lpmapi
TRACE = -L/bgl/local/lib -lmpitrace_f  #for MPI information only
#TRACE = -L/bgl/local/lib -lmpihpm_f -lbgl_perfctr.rts #for MPI+HPM information 
LIBS_MPI = $(TRACE) -lmpich.rts -lmsglayer.rts -lrts.rts -ldevices.rts
LIBSF_MPI = -lmpich.rts -lfmpich.rts -lmsglayer.rts -lrts.rts -ldevices.rts

.SUFFIXES: .f90

# -------------------------------------------------------------------------------------------

# Source files

# -------------------------------------------------------------------------------------------
MODULES.f90 =   physvars.f90 

SOURCES.f90 =   pepcb.f90 \
                setup.f90 setup_arrays.f90 stamp.f90 \
                configure.f90 randion.f90 predef_parts.f90 maxwell1.f90 cold_start.f90 \
                special_start.f90 constrain.f90 cutvector.f90 face.f90 \
                scramble_v.f90 reset_ions.f90 add_electrons.f90 double_target.f90 \
                laser.f90 beam.f90 beam_control.f90 beam_dust.f90 rezone.f90\
                mc_config.f90 \
                velocities.f90 \
                energy_cons.f90 potenergy.f90 kinenergy.f90\
                force_laser.f90 \
                fpond.f90 empond.f90 emplane.f90 laser_bullet.f90 emoblique.f90 \
                push_em.f90 push_full3v.f90 \
                push.f90 reinject.f90 earth_plate.f90 add_ramp.f90 \
                vis_parts_nbody.f90 vis_fields.f90 dump.f90 param_dump.f90 diagnostics.f90 \
                open_files.f90 close_files.f90  \
                densities.f90 sum_fields.f90 dump_fields.f90 \
                track_nc.f90 \
                integrator.f90 abort.f90

# -------------------------------------------------------------------------------------------

# Names of application object files derived from sources

appobjs = $(SOURCES.f90:.f90=.o) ## ncnbody.o
appmods =  $(MODULES.f90:.f90=.o)
mods = $(appmods) $(pepcmods)

default: lpepc pepcb

#  compile mode
.f90.o:   
	$(FC) -c $(FFLAGS) $(VISITLIBS) $(DBFLAGS) $<
.c.o:   
	$(CC) -c $(CFLAGS)  $<

# -----------------------------------------------------------------------------------------
#  targets


pepcb:	$(appmods) $(appobjs) lpepc
	@echo "Creating pepcb ..."
	$(FC)  $(FFLAGS) $(LDFLAGS) $(DBFLAGS) -o $@ $(appobjs) $(appmods) $(pepcmods)  $(LIBPEPC)  $(IBMLIB) $(LIBS_MPI) $(LIBSF_MPI)
	@echo "... done"

# Update PEPC library and headers
lpepc:  
	cd ../lpepcsrc; rm *.mod treevars.o utils.o tree_utils.o; make lpepc; cd ../pepc-b
	/bin/ln -sf $(PREFIX)/lpepcsrc/treevars.o $(PREFIX)/pepc-b/treevars.o
	/bin/ln -sf $(PREFIX)/lpepcsrc/treevars.mod $(PREFIX)/pepc-b/treevars.mod
	/bin/ln -sf $(PREFIX)/lpepcsrc/utils.o $(PREFIX)/pepc-b/utils.o
	/bin/ln -sf $(PREFIX)/lpepcsrc/utils.mod $(PREFIX)/pepc-b/utils.mod
	/bin/ln -sf $(PREFIX)/lpepcsrc/liblpepc.a $(PREFIX)/pepc-b/liblpepc.a
	/bin/ln -sf $(PREFIX)/lpepcsrc/*.opari.inc $(PREFIX)/pepc-b/.



debug:	$(OBJECTS)
	@echo "Loading object files using flags $(DBFLAGS)"
	$(FC) $(LDFLAGS)  $(VISITLIBS) $(DBFLAGS)  -o $@ $(OBJECTS)


rand:	 rantest.o
	$(FC) $(FFLAGS) -o rand rantest.o


clean:
	/bin/rm -f bin/*.o  *% *~ bin/key rand *.o *.pif core *.a *.mod opari.rc *.mod.F90 *.opari.inc

cleanlib: 
	cd ../lpepcsrc;	/bin/rm -f *.o *.mod *.mod.F90 *.opari.inc; cd -

cleanmod:
	rm -f physvars.o physvars.mod

cleancore: 	
	/bin/rm $(pepcobjs) 

cleandata: 
	/bin/rm -f fort.* *.dat *.out *.gle 

depend:
	./f90depend -u -I/usr/lpp/ppe.poe/include/thread *.f90

### DO NOT DELETE OR EDIT THIS LINE
# Everything from here on down is generated by f90depend
# so do NOT add any translation rules below here.

#--- Include-Dependencies
