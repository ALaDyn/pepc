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
                module_debug.f90 module_math_tools.f90 module_mirror_boxes.f90\
                module_setup.f90 module_tree_domains.f90 \
                module_allocation.f90 tree_stats.f90 \
                module_tree.f90 module_spacefilling.f90 module_htable.f90 \
                module_pepcfields.f90 module_branching.f90 module_tree_walk_communicator.f90
ifeq ($(SMPSS), 1)
PEPCSMPSS.f90   = module_tree_walk_smpss_utils.f90 module_tree_walk_smpss.f90
else
PEPCSRC.f90    += module_tree_walk_pthreads.f90
endif
PEPCBACKEND.f90 = module_calc_force.f90 module_interaction_specific.f90
#-------------------------------------------------------------------------------------------

-include $(FULLBACKENDDIR)makefile.backend

# Names of application object files derived from sources
# Prefix added for multi-arch builds
OBJS        = $(addprefix $(BUILDDIR), $(PEPCSRC.f90:.f90=.o))
OBJSBACKEND = $(addprefix $(BUILDDIR), $(PEPCBACKEND.f90:.f90=.o))

PEPCBACKSRC = $(addprefix $(FULLBACKENDDIR), $(PEPCBACKEND.f90))

default: $(BUILDDIR) $(LIBNAME)

$(LIBNAME): $(OBJSBACKEND) $(OBJS)
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
$(BUILDDIR)pthreads_stuff.mod:
	ln -sf ../../../pthreads/pthreads_stuff.mod $(BUILDDIR)pthreads_stuff.mod


clean:
	$(RM) $(MYROOT)$(LIBDIR)libpepc.*.a $(OBJS) $(BUILDDIR)*.mod $(BUILDDIR)*.o $(BUILDDIR)module.deps
	$(RM) -rf $(BUILDDIR) $(MACH)

$(BUILDDIR)module.deps: $(BUILDDIR)
	$(MYROOT)$(F90MODDEPS) -o $@ --dep-template "\1.mod" --mod-template "\1.mod" --o-prefix "\$$(BUILDDIR)" $(PEPCSRC.f90) $(PEPCBACKSRC)

include $(BUILDDIR)module.deps
