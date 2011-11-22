MYROOT = ../../../

include $(MYROOT)/makefile.defs

# We set some default backend. TODO: issue a warning message if $BACKEND is not set
ifeq ($(strip $(BACKEND)),)
	BACKEND=coulomb
endif

LIBNAME = $(MYROOT)$(LIBDIR)libpepc.$(BACKEND).a
BUILDDIR = $(MACH)/$(BACKEND)/
FULLBACKENDDIR = $(MYROOT)$(INTSPECDIR)$(BACKEND)/

.SUFFIXES: .f90
# -------------------------------------------------------------------------------------------
# Source files
PEPCSRC.f90 =   treevars.f90 tree_utils.f90 timings.f90 treetypes.f90 \
                module_debug.f90 module_math_tools.f90 module_fmm_framework.f90 \
                tree_aswalk_pthreads.f90 libpepc_setup.f90 module_tree_domains.f90 \
                module_allocation.f90 libpepc_finalize.f90 tree_stats.f90 \
                module_tree.f90 module_spacefilling.f90 module_htable.f90 \
                module_pepcfields.f90 module_branching.f90 tree_walk_communicator.f90
PEPCBACKEND.f90 = module_calc_force.f90 module_interaction_specific.f90 module_pepc_wrappers.f90
#-------------------------------------------------------------------------------------------

# Names of application object files derived from sources
# Prefix added for multi-arch builds
OBJS        = $(addprefix $(BUILDDIR), $(PEPCSRC.f90:.f90=.o))
OBJSBACKEND = $(addprefix $(BUILDDIR), $(PEPCBACKEND.f90:.f90=.o))

PEPCBACKSRC = $(addprefix $(FULLBACKENDDIR), $(PEPCBACKEND.f90))

default: $(BUILDDIR) libpepc

libpepc: $(OBJSBACKEND) $(OBJS)
	$(RM) $(LIBNAME)
	$(AR) $(ARFLAGS) $(LIBNAME) $^
	$(RANLIB) $(LIBNAME)

$(BUILDDIR):
	mkdir -p  $(BUILDDIR)

# preprocess files and put into $(BUILDDIR) subdir
# then compile form there
$(BUILDDIR)%.o: %.f90
	$(CPP) $(CPPFLAGS) $< $(BUILDDIR)/$(<F)
	cd $(BUILDDIR) ; $(FCPRE) $(FC) -c $(FFLAGS) $(DBFLAGS) $(<F) -o $(<:.f90=.o)

$(BUILDDIR)%.o: $(FULLBACKENDDIR)%.f90
	$(CPP) $(CPPFLAGS) $< $(BUILDDIR)/$(<F)
	cd $(BUILDDIR) ; $(FCPRE) $(FC) -c $(FFLAGS) $(DBFLAGS) $(<F) -o $(<F:.f90=.o)


# special rules for "external" modules
pthreads_stuff.mod:
	ln -sf ../../../pthreads/pthreads_stuff.mod $(BUILDDIR)pthreads_stuff.mod


clean:
	$(RM) $(BUILDDIR)liblpepc.a $(OBJS) $(BUILDDIR)/*.mod *.mod *% *~ bin/key rand *.o *.pif core module.deps
	$(RM) -rf $(BUILDDIR) $(MACH)

module.deps:
	$(MYROOT)$(F90MODDEPS) -o $@ --dep-template "\1.mod" --mod-template "\1.mod" --o-prefix "\$$(BUILDDIR)" $(PEPCSRC.f90) $(PEPCBACKSRC)

include module.deps
