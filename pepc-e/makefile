#
#  Makefile for treemp and test programs
#

PREFIX = ..

#  Get machine-dependent flags

include ../makefile.defs

FFLAGS += -I../lpepcsrc/$(MACH)

LIBPEPC   = -L../lpepcsrc/$(MACH) -llpepc ../sl_pepc/libsl.a
LIBS      = $(LIBPEPC)

.SUFFIXES: .f90


# -------------------------------------------------------------------------------------------

# Source files

# -------------------------------------------------------------------------------------------

MODULES.f90 =   physvars.f90 utils.f90 benchmarking.f90 pusher.f90 files.f90 energies.f90

SOURCES.f90 =	\
		setup.f90 stamp.f90 configure.f90 face.f90 \
		randion.f90 maxwell1.f90 special_start.f90 cold_start.f90 scramble_v.f90 \
		force_direct.f90 error_test.f90 sum_radial.f90 \
		param_dump.f90 cleanup.f90 pepce.f90

# -------------------------------------------------------------------------------------------

# Names of application object files derived from sources
# Prefix added for multi-arch builds

appobjs = $(addprefix $(MACH)/, $(SOURCES.f90:.f90=.o)) \
	  $(addprefix $(MACH)/, $(NCOBJS) )
appmods =  $(addprefix $(MACH)/, $(MODULES.f90:.f90=.o))
mods = $(appmods)

default: dirs pepce
all: dirs pepce
dirs: 	
	if [ ! -d $(MACH) ]; then mkdir $(MACH); fi

#  compile rules
$(MACH)/%.o: %.f90  
	$(CPP) $(CPPFLAGS) $< $(MACH)/$<
	$(FCPRE) $(FC) -c $(FFLAGS) $(DEFINES) $(MACH)/$< -o $@


#  explicit rules for f90 modules to avoid conflict with built-in ones
$(MACH)/physvars.o: physvars.f90
	$(CPP) $(CPPFLAGS) physvars.f90 $(MACH)/physvars.f90
	$(FCPRE) $(FC) -c $(FFLAGS) $(DEFINES) $(MACH)/physvars.f90 -o $@

$(MACH)/%.o: %.c  
	$(CC) -c $(CFLAGS) $< -o $@


# -----------------------------------------------------------------------------------------
#  targets
pepce:	module.deps lpepc $(appmods) $(appobjs)
	@echo "Creating application binary $(EXECNAME) ..."
	$(LDPRE) $(LD) $(LDFLAGS) -o $(EXECNAME) $(appobjs) $(appmods) $(LIBPEPC) $(LDLIBS)
	@echo "... done"

# Update PEPC library and headers
lpepc:  ../lpepcsrc/$(MACH)/liblpepc.a
#	cd ../lpepcsrc; make lpepc; cd ../pepc-e
	/bin/ln -sf $(PREFIX)/lpepcsrc/$(MACH)/treevars.mod $(PREFIX)/pepc-e/treevars.mod
	/bin/ln -sf $(PREFIX)/lpepcsrc/$(MACH)/tree_utils.mod $(PREFIX)/pepc-e/tree_utils.mod
	/bin/ln -sf $(PREFIX)/lpepcsrc/$(MACH)/utils.mod $(PREFIX)/pepc-e/utils.mod
	/bin/ln -sf $(PREFIX)/lpepcsrc/$(MACH)/timings.mod $(PREFIX)/pepc-e/timings.mod

clean:
	$(RM) $(MACH)/*.o  *% *~ bin/key rand *.o *.pif core *.mod $(MACH)/*.mod module.deps $(EXECNAME)
	$(RM) -r $(MACH)

default: pepce

module.deps:
	../tools/build/f90_mod_deps.py -o $@ --o-prefix "\$$(MACH)/" $(MODULES.f90) $(SOURCES.f90) $(VISIT.f90)

-include module.deps
