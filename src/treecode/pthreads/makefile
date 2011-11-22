MYROOT = ../../../

include $(MYROOT)/makefile.defs

LIBNAME = $(MYROOT)/$(LIBDIR)/libpthreads.a
MODFILE = pthreads_stuff.mod

# -------------------------------------------------------------------------------------------
# Source files
SOURCES.f90 =   pthreads.f90
SOURCES.c   =   pthreads.c bgp_core.c
#-------------------------------------------------------------------------------------------

# Names of application object files derived from sources
# Prefix added for multi-arch builds
OBJS = $(addprefix $(MACH)/, $(SOURCES.f90:.f90=_f.o) $(SOURCES.c:.c=_c.o))

default: all

clean:
	$(RM) $(MACH)/*.o *.mod $(LIBNAME) $(MODNAME) $(MACH)
	
all: $(LIBNAME)

$(MACH):
	mkdir -p $(MACH)

$(MACH)/%_f.o: %.f90
	$(FCPRE) $(FC) $(FFLAGS) -c $< -o $@

$(MACH)/%_c.o: %.c
	$(CCPRE) $(CC) $(CFLAGS) -c $< -o $@
	
$(LIBNAME): $(MACH) $(OBJS)
	$(AR) $(ARFLAGS) $(LIBNAME) $(OBJS)
	$(RANLIB) $(LIBNAME)
