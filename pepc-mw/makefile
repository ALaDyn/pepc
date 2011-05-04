PREFIX = ..

#  Get machine-dependent flags

include ../makefile.defs

LIBRARIES = ../lpepcsrc/$(MACH)/liblpepc.a ../sl_pepc/libsl.a ../pthreads/$(MACH)/libpthreads.a
PEPCGENERIC = ../pepc-generic

.SUFFIXES: .f90

EXECNAME  = ../pepcmw

# -------------------------------------------------------------------------------------------
# Source files
SOURCES.f90 =   module_units.f90 physvars.f90 module_io.f90 \
                module_laser.f90 module_pusher.f90 \
                module_workflow.f90 module_param_dump.f90 \
                setup.f90 configure.f90 special_start.f90 density_helmholtz.f90 \
	 	energy_cons.f90 potenergy.f90 kinenergy.f90 cleanup.f90 pepc.f90 \
	 	module_acf.f90 sum_radial.f90 \
	 	module_diagnostics.f90 module_fields.f90 module_icosahedron.f90
SOURCES.GENERIC.f90 = module_vtk.f90 module_base64.f90 module_checkpoint.f90 \
                benchmarking.f90 module_velocity_setup.f90
# -------------------------------------------------------------------------------------------

# Names of application object files derived from sources
# Prefix added for multi-arch builds
SOURCES = $(SOURCES.f90) $(SOURCES.GENERIC.f90)
OBJS = $(addprefix $(MACH)/, $(SOURCES:.f90=.o))
SOURCES.GENERICP.f90 = $(addprefix $(PEPCGENERIC)/, $(SOURCES.GENERIC.f90))

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

$(MACH)/%.o %.mod: $(PEPCGENERIC)/%.f90  
	$(CPP) $(CPPFLAGS) $< $(MACH)/$(<F)
	$(FCPRE) $(FC) -c $(FFLAGS) $(DBFLAGS) $(MACH)/$(<F) -o $(MACH)/$(<F:.f90=.o)

# avoid infinite recursion of rule below
../lpepcsrc/%.mod: ../lpepcsrc/%.mod
	@echo "" $@
# special rules for "external" modules
%.mod: ../lpepcsrc/%.mod
	ln -sf ../lpepcsrc/$@ ./$(@F)

# final rule for executable
$(EXECNAME): $(LIBRARIES) $(OBJS)
	@echo "Creating application binary $(EXECNAME) ..."
	$(LDPRE) $(LD) $(LDFLAGS) -o $(EXECNAME) $(OBJS) $(LIBPEPC) $(LIBRARIES) $(LDLIBS)
	@echo "... done"

clean:
	$(RM) $(MACH)/*.o  *% *~ bin/key rand *.o *.pif core *.mod $(MACH)/*.mod module.deps module.generic.deps $(EXECNAME)
	$(RM) -r $(MACH)

module.deps:
	../tools/build/f90_mod_deps.py -o $@ --dep-template "\1.mod" --mod-template "\1.mod" --o-prefix "\$$(MACH)/" $(SOURCES.f90)
module.generic.deps:
	if [ -z "$(SOURCES.GENERIC.f90)" ]; then \
	  touch module.generic.deps; \
	else \
	../tools/build/f90_mod_deps.py -o $@ --dep-template "\1.mod" --mod-template "\1.mod" --o-prefix "\$$(MACH)/" $(SOURCES.GENERICP.f90); \
	fi

include module.deps
include module.generic.deps
