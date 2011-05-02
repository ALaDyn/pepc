PREFIX = ..

#  Get machine-dependent flags

include ../makefile.defs

LIBRARIES = ../lpepcsrc/$(MACH)/liblpepc.a ../sl_pepc/libsl.a ../pthreads/$(MACH)/libpthreads.a
PEPCGENERIC = ../pepc-generic

.SUFFIXES: .f90

EXECNAME = ../pepce

# -------------------------------------------------------------------------------------------
# Source files
SOURCES.f90 =   physvars.f90 benchmarking.f90 velocity_setup.f90 \
                param_dump.f90 face.f90 files.f90 pusher.f90 randion.f90 \
                setup.f90 configure.f90 special_start.f90 stamp.f90 \
	 	energies.f90 cleanup.f90 pepce.f90 \
	 	module_vtk.f90 module_diagnostics.f90 module_base64.f90 \
	 	module_checkpoint.f90
# -------------------------------------------------------------------------------------------

# Names of application object files derived from sources
# Prefix added for multi-arch builds
OBJS = $(addprefix $(MACH)/, $(SOURCES.f90:.f90=.o)) \

default: all
all: $(MACH) $(EXECNAME)

$(MACH): 	
	mkdir -p $(MACH)

# preprocess files and put into $(MACH) subdir
# then compile form there
# touch module files after compilation so that the are newer than the .o files
# this prevents make from trying to rebuild them again and again
# we keep the preprocessed files for easier debugging / correct line information etc.
$(MACH)/%.o %.mod: %.f90  
	$(CPP) $(CPPFLAGS) $< $(MACH)/$<
	$(FCPRE) $(FC) -c $(FFLAGS) $(DBFLAGS) $(MACH)/$< -o $(MACH)/$(<:.f90=.o)

# avoid infinite recursion of rule below
../lpepcsrc/%.mod: ../lpepcsrc/%.mod
	@echo "========================================== " $@
# special rules for "external" modules
%.mod: ../lpepcsrc/%.mod
	ln -sf ../lpepcsrc/$@ ./$(@F)


$(EXECNAME): $(LIBRARIES) $(OBJS)
	@echo "Creating application binary $(EXECNAME) ..."
	$(LDPRE) $(LD) $(LDFLAGS) -o $(EXECNAME) $(OBJS) $(LIBPEPC) $(LIBRARIES) $(LDLIBS)
	@echo "... done"

clean:
	$(RM) $(MACH)/*.o  *% *~ bin/key rand *.o *.pif core *.mod $(MACH)/*.mod module.deps $(EXECNAME)
	$(RM) -r $(MACH)

module.deps:
	../tools/build/f90_mod_deps.py -o $@ --dep-template "\1.mod" --mod-template "\1.mod" --o-prefix "\$$(MACH)/" $(SOURCES.f90)

-include module.deps
