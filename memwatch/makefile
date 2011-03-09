
include ../makefile.defs

LIBNAME = libmemwatch.a

# -------------------------------------------------------------------------------------------
# Source files
SOURCES.f90 =   memwatchf.f90
SOURCES.c   =   memwatchc.c signalhandler.c
#-------------------------------------------------------------------------------------------

# Names of application object files derived from sources
# Prefix added for multi-arch builds
OBJS = $(addprefix $(MACH)/, $(SOURCES.f90:.f90=_f.o) $(SOURCES.c:.c=_c.o))

# override stupid internal rule for .o files to depend on .mod files and call m2c
%.o: %.mod

# for creating the .mod files, we only have to compile the .o-files from the respective .f90 files
# gfortran stores the mod-files in the current directory - just where we need them
# unfortunately, this rule only works if it actually does somenthing, so we echo some useless info
%.mod: $(MACH)/%.o
	@echo 'Making module file' $@


default: all

clean:
	$(RM) $(MACH)/*.o *.mod $(MACH)/$(LIBNAME) $(MACH)
	
all:    $(MACH) $(MACH)/$(LIBNAME)	

$(MACH):  
	mkdir -p  $(MACH)

$(MACH)/%_f.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

$(MACH)/%_c.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@
	
$(MACH)/$(LIBNAME): $(OBJS)
	$(AR) $(ARFLAGS) $(MACH)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(MACH)/$(LIBNAME)
	
test:   $(MACH)/memwatchtest

$(MACH)/memwatchtest: all
	$(FC) $(FFLAGS) memwatchtest.f90 $(MACH)/$(LIBNAME) -o $(MACH)/memwatchtest
	
	
