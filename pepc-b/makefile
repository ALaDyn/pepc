#
#  Makefile for treemp and test programs
#

PREFIX = ..

#  Get machine-dependent flags

include ../makefile.defs


LIBRARIES   = ../lpepcsrc/$(MACH)/liblpepc.a ../sl_pepc/libsl.a ../pthreads/$(MACH)/libpthreads.a
# ../memwatch/$(MACH)/libmemwatch.a 
# ../pthreads/$(MACH)/libpthreads.a

# Override default executable name
EXECNAME = ../pepcb

.SUFFIXES: .f90

# -------------------------------------------------------------------------------------------
# Source files
# -------------------------------------------------------------------------------------------

SOURCES.f90 = 	module_physvars.f90 \
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
#		module_xnbody.f90 \	



# Names of application object files derived from sources
# Prefix added for multi-arch builds
OBJS = $(addprefix $(MACH)/, $(SOURCES.f90:.f90=.o)) 

default: all
all: $(MACH) $(EXECNAME)

$(MACH): 	
	mkdir -p $(MACH)

# preprocess files and put into $(MACH) subdir
# then compile from there
# touch module files after compilation so that they are newer than the .o files
# this prevents make from trying to rebuild them again and again
# we keep the preprocessed files for easier debugging / correct line information etc.
$(MACH)/%.o: %.f90  
	$(CPP) $(CPPFLAGS) $< $(MACH)/$<
	$(FCPRE) $(FC) -c $(FFLAGS) $(DBFLAGS) $(MACH)/$< -o $@
	touch -c $(@F:.o=.mod)


#  Modules

# override stupid internal rule for .o files to depend on .mod files and call m2c
%.o: %.mod

# for creating the .mod files, we only have to compile the .o-files from the respective .f90 files
# gfortran stores the mod-files in the current directory - just where we need them
# unfortunately, this rule only works if it actually does something, so we echo some useless info
%.mod: $(MACH)/%.o
	@echo 'Making module file' $@

# special rules for "external" modules
treevars.mod:
	ln -sf ../lpepcsrc/treevars.mod ./treevars.mod
tree_utils.mod:
	ln -sf ../lpepcsrc/tree_utils.mod ./tree_utils.mod
module_fmm_framework.mod:
	ln -sf ../lpepcsrc/module_fmm_framework.mod ./module_fmm_framework.mod
module_pepcfields.mod:
	ln -sf ../lpepcsrc/module_pepcfields.mod ./module_pepcfields.mod
#tree_walk_utils.mod:
#	ln -sf ../lpepcsrc/tree_walk_utils.mod ./tree_walk_utils.mod
timings.mod:
	ln -sf ../lpepcsrc/timings.mod ./timings.mod
#	/bin/ln -sf $(PREFIX)/lpepcsrc/$(MACH)/timings.mod $(PREFIX)/pepc-mw/timings.mod

# -------------------------------------------------------------------------------------------


# check for sion usage, add dump_sion or dump_ascii to src-list
ifeq ($(SION),1)
#SOURCES.f90 += dump_sion.f90 predef_parts_sion.f90
LDLIBS += $(LIB_SION)
else
#SOURCES.f90 += dump_ascii.f90 predef_parts_ascii.f90
endif

# Names of application object files derived from sources
# Prefix added for multi-arch builds



$(EXECNAME): $(LIBRARIES) $(OBJS)
	@echo "Creating application binary $(EXECNAME) ..."
	$(LDPRE) $(LD) $(LDFLAGS) -o $(EXECNAME) $(OBJS) $(LIBPEPC) $(LDLIBS) $(LIBRARIES)
	@echo "... done"

clean:
	$(RM) $(MACH)/*.o  *% *~ bin/key rand *.o *.pif core *.mod $(MACH)/*.mod module.deps $(EXECNAME)
	$(RM) -r $(MACH)

module.deps:
	../tools/build/f90_mod_deps.py -o $@ --dep-template "\1.mod" --mod-template "\1.mod" --o-prefix "\$$(MACH)/" $(SOURCES.f90)

-include module.deps





