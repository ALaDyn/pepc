PREFIX = ..

#  Get machine-dependent flags

include makefile.defs

LIBRARIES = 

.SUFFIXES: .f90


# -------------------------------------------------------------------------------------------
# Source files
SOURCES.f90 =   module_calc_force.f90
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

$(EXECNAME): $(LIBRARIES) $(OBJS)
	@echo "Creating application binary $(EXECNAME) ..."
	$(LDPRE) $(LD) $(LDFLAGS) -o $(EXECNAME) $(OBJS) $(LIBPEPC) $(LIBRARIES) $(LDLIBS)
	@echo "... done"

clean:
	$(RM) $(MACH)/*.o  *% *~ *.o *.pif core *.mod $(MACH)/*.mod module.deps $(EXECNAME)
	$(RM) -r $(MACH)

module.deps:
	../tools/build/f90_mod_deps.py -o $@ --dep-template "\1.mod" --mod-template "\1.mod" --o-prefix "\$$(MACH)/" $(SOURCES.f90)

-include module.deps
