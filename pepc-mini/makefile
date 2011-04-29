PREFIX = ..

#  Get machine-dependent flags

include ../makefile.defs

LIBRARIES = ../lpepcsrc/$(MACH)/liblpepc.a ../sl_pepc/libsl.a ../pthreads/$(MACH)/libpthreads.a

.SUFFIXES: .f90

EXECNAME = ../pepcmini


# -------------------------------------------------------------------------------------------
# Source files
SOURCES.f90 =   physvars.f90 files.f90 pusher.f90 setup.f90 special_start.f90 stamp.f90 \
	 	cleanup.f90 pepc.f90 
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
$(MACH)/%.o: %.f90  
	$(CPP) $(CPPFLAGS) $< $(MACH)/$<
	$(FCPRE) $(FC) -c $(FFLAGS) $(DBFLAGS) $(MACH)/$< -o $@
	touch -c $(@F:.o=.mod)

# override stupid internal rule for .o files to depend on .mod files and call m2c
%.o: %.mod

# for creating the .mod files, we only have to compile the .o-files from the respective .f90 files
# gfortran stores the mod-files in the current directory - just where we need them
# unfortunately, this rule only works if it actually does somenthing, so we echo some useless info
%.mod: $(MACH)/%.o
	@echo 'Making module file' $@

# special rules for "external" modules
treevars.mod:
	ln -sf ../lpepcsrc/treevars.mod ./treevars.mod
utils.mod:
	ln -sf ../lpepcsrc/utils.mod ./utils.mod
tree_utils.mod:
	ln -sf ../lpepcsrc/tree_utils.mod ./tree_utils.mod
module_fmm_framework.mod:
	ln -sf ../lpepcsrc/module_fmm_framework.mod ./module_fmm_framework.mod
tree_walk_utils.mod:
	ln -sf ../lpepcsrc/tree_walk_utils.mod ./tree_walk_utils.mod
tree_walk_communicator.mod:
	ln -sf ../lpepcsrc/tree_walk_communicator.mod ./tree_walk_communicator.mod
timings.mod:
	ln -sf ../lpepcsrc/timings.mod ./timings.mod
#	/bin/ln -sf $(PREFIX)/lpepcsrc/$(MACH)/timings.mod $(PREFIX)/pepc-mw/timings.mod


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
