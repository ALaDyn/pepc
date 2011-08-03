PREFIX = ..

#  Get machine-dependent flags

include ../makefile.defs

LIBPEPC     = $(PREFIX)/lpepcsrc/$(MACH)/liblpepc.a 
LIBSL       = $(PREFIX)/sl_pepc/libsl.a 
LIBPTHREADS = $(PREFIX)/pthreads/$(MACH)/libpthreads.a
LIBPEPCS    = $(PREFIX)/pepc-s/$(MACH)/libpepcs.a

LIBRARIES   = $(LIBPEPCS) $(LIBPEPC) $(LIBSL) $(LIBPTHREADS)

PEPCGENERIC = ../pepc-generic

.SUFFIXES: .f90

EXECNAME = ../pepcs

# -------------------------------------------------------------------------------------------
# Source files
SOURCES.f90 =   physvars.f90 setup.f90 cleanup.f90 pepc-s.f90 
SOURCES.GENERIC.f90 = 
# -------------------------------------------------------------------------------------------

# Names of application object files derived from sources
# Prefix added for multi-arch builds
SOURCES = $(SOURCES.f90) $(SOURCES.GENERIC.f90)
OBJS = $(addprefix $(MACH)/, $(SOURCES:.f90=.o))
SOURCES.GENERICP.f90 = $(addprefix $(PEPCGENERIC)/, $(SOURCES.GENERIC.f90))
MAINOBJ = $(MACH)/main.o

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
	$(FCPRE) $(FC) -c $(FFLAGS) $(DEFINES) $(MACH)/$< -o $(MACH)/$(<:.f90=.o)

$(MACH)/%.o: $(PEPCGENERIC)/%.f90  
	$(CPP) $(CPPFLAGS) $< $(MACH)/$(<F)
	$(FCPRE) $(FC) -c $(FFLAGS) $(DEFINES) $(MACH)/$(<F) -o $(MACH)/$(<F:.f90=.o)

# avoid infinite recursion of rule below
../lpepcsrc/%.mod: ../lpepcsrc/%.mod
	@echo "" $@
# special rules for "external" modules
%.mod: ../lpepcsrc/%.mod
	ln -sf ../lpepcsrc/$@ ./$(@F)

# final rule for executable
$(EXECNAME): $(LIBRARIES) $(MAINOBJ)
	@echo "Creating application binary $(EXECNAME) ..."
	$(LDPRE) $(LD) $(LDFLAGS) -o $(EXECNAME) $(MAINOBJ) $(LIBRARIES) $(LDLIBS)
	@echo "... done"
	
$(LIBPEPCS): $(OBJS)
	$(AR) $(ARFLAGS) $@ $^ 

clean:
	$(RM) $(MACH)/*.o  *% *~ bin/key rand *.o *.pif core *.mod $(MACH)/*.mod module.deps $(EXECNAME)
	$(RM) -r $(MACH)

module.deps:
	../tools/build/f90_mod_deps.py -o $@ --dep-template "\1.mod" --mod-template "\1.mod" --o-prefix "\$$(MACH)/" $(SOURCES.f90) $(SOURCES.GENERICP.f90)

include module.deps