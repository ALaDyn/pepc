PREFIX = ..

#  Get machine-dependent flags

include ../makefile.defs

LIBRARIES = ../lpepcsrc/$(MACH)/liblpepc.a ../sl_pepc/libsl.a ../pthreads/$(MACH)/libpthreads.a
PEPCGENERIC = ../pepc-generic

.SUFFIXES: .f90

# Override default executable name
EXECNAME = ../pepcb

# -------------------------------------------------------------------------------------------
# Source files
SOURCES.f90 =  module_physvars.f90 \
               module_particle_props.f90 \
               module_utilities.f90 \
               module_velocity_setup.f90 \
               module_field_grid.f90 \
               module_laser.f90 \
               module_io.f90 \
               module_geometry.f90 \
               module_integration_scheme.f90 \
               module_diagnostics.f90 \
               module_particle_beam.f90 \
               module_particle_boundaries.f90 \
               setup_arrays.f90 setup.f90 configure.f90 integrator.f90 diagnostics.f90 pepcb.f90 
#               module_xnbody.f90
SOURCES.GENERIC.f90 = module_gle.f90 

# -------------------------------------------------------------------------------------------


# Names of application object files derived from sources
# Prefix added for multi-arch builds
SOURCES = $(SOURCES.f90) $(SOURCES.GENERIC.f90)
OBJS = $(addprefix $(MACH)/, $(SOURCES:.f90=.o))
SOURCES.GENERICP.f90 = $(addprefix $(PEPCGENERIC)/, $(SOURCES.GENERIC.f90))

default: all
all: $(MACH)
	$(MAKE) $(MFLAGS) $(EXECNAME)

$(MACH): 	
	mkdir -p $(MACH)

# preprocess files and put into $(MACH) subdir
# then compile form there
# touch module files after compilation so that the are newer than the .o files
# this prevents make from trying to rebuild them again and again
# we keep the preprocessed files for easier debugging / correct line information etc.
$(MACH)/%.o: %.f90
	$(CPP) $(CPPFLAGS) $< $(MACH)/$<
	$(FCPRE) $(FC) -c $(FFLAGS) $(DBFLAGS) $(MACH)/$< -o $(MACH)/$(<:.f90=.o)

$(MACH)/%.o: $(PEPCGENERIC)/%.f90
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
	$(LDPRE) $(LD) $(LDFLAGS) -o $(EXECNAME) $(OBJS) $(LIBPEPC) $(LDLIBS) $(LIBRARIES)
	@echo "... done"

clean:
	$(RM) $(MACH)/*.o  *% *~ bin/key rand *.o *.pif core *.mod $(MACH)/*.mod module.deps $(EXECNAME)
	$(RM) -r $(MACH)

module.deps:
	../tools/build/f90_mod_deps.py -o $@ --dep-template "\1.mod" --mod-template "\1.mod" --o-prefix "\$$(MACH)/" $(SOURCES.f90) $(SOURCES.GENERICP.f90)

include module.deps
