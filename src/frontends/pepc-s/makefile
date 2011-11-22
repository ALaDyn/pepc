MYROOT = ../../../

#  Get machine-dependent flags

include $(MYROOT)makefile.defs

# We set some default backend. TODO: issue a warning message if $BACKEND is not set
ifeq ($(strip $(BACKEND)),)
	BACKEND=coulomb
endif

LIBNAME     = $(MYROOT)$(LIBDIR)lpepcs.a
PEPCGENERIC = $(MYROOT)$(UTILSDIR)
BUILDDIR    = $(MACH)/$(BACKEND)/
EXECNAME    = $(MYROOT)$(BINDIR)pepcs
LIBRARIES   = $(LIBNAME) $(MYROOT)$(LIBDIR)libpepc.$(BACKEND).a $(MYROOT)$(LIBDIR)libpthreads.a $(MYROOT)$(LIBDIR)libsl.a

.SUFFIXES: .f90


# -------------------------------------------------------------------------------------------
# Source files
SOURCES.f90 =   physvars.f90 setup.f90 cleanup.f90 pepc-s.f90 main.f90
SOURCES.GENERIC.f90 = 
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
$(MACH)/%.o: %.f90
	$(CPP) $(CPPFLAGS) $< $(MACH)/$<
	$(FCPRE) $(FC) -c $(FFLAGS) $(DBFLAGS) $(MACH)/$< -o $(MACH)/$(<:.f90=.o)

$(MACH)/%.o: $(PEPCGENERIC)/%.f90
	$(CPP) $(CPPFLAGS) $< $(MACH)/$(<F)
	$(FCPRE) $(FC) -c $(FFLAGS) $(DBFLAGS) $(MACH)/$(<F) -o $(MACH)/$(<F:.f90=.o)

# special rules for "external" modules
%.mod: ../../treecode/lpepcsrc/$(BUILDDIR)%.mod
	ln -sf ../../treecode/lpepcsrc/$(BUILDDIR)$@ ./$(@F)

# final rule for executable
$(LIBNAME): $(OBJS)
	$(RM) $(LIBNAME)
	$(AR) $(ARFLAGS) $(LIBNAME) $^
	$(RANLIB) $(LIBNAME)

$(EXECNAME): $(LIBRARIES)
	$(LDPRE) $(LD) $(LDFLAGS) -o $(EXECNAME) $+ $(LDLIBS)


clean:
	$(RM) $(MACH)/*.o  *% *~ bin/key rand *.o *.pif core *.mod $(MACH)/*.mod module.deps $(LIBNAME) $(EXECNAME)
	$(RM) -rf $(MACH)

module.deps:
	$(MYROOT)$(F90MODDEPS) -o $@ --dep-template "\1.mod" --mod-template "\1.mod" --o-prefix "\$$(MACH)/" $(SOURCES.f90) $(SOURCES.GENERICP.f90)

include module.deps
