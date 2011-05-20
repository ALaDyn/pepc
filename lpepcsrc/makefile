include ../makefile.defs

pepclib = $(MACH)/liblpepc.a

.SUFFIXES: .f90
# -------------------------------------------------------------------------------------------
# Source files
PEPCSRC.f90 =   treevars.f90 tree_utils.f90 timings.f90 treetypes.f90 \
                module_debug.f90 module_math_tools.f90 module_fmm_framework.f90 \
                module_calc_force.f90 tree_aswalk_pthreads.f90 libpepc_setup.f90 \
                tree_domains.f90 restore.f90 module_allocation.f90 libpepc_finalize.f90 \
		tree_stats.f90 tree_local.f90 tree_exchange.f90 tree_global.f90 \
		make_hashentry.f90 get_address.f90 get_next_node.f90 \
	        diagnose_tree.f90 check_table.f90 module_pepcfields.f90 mac.f90
#-------------------------------------------------------------------------------------------

# Names of application object files derived from sources
# Prefix added for multi-arch builds
OBJS = $(addprefix $(MACH)/, $(PEPCSRC.f90:.f90=.o))

default: all
all: $(MACH) module.deps $(pepclib)

$(MACH):  
	mkdir -p  $(MACH)

# preprocess files and put into $(MACH) subdir
# then compile form there
# touch module files after compilation so that the are newer than the .o files
# this prevents make from trying to rebuild them again and again
# we keep the preprocessed files for easier debugging / correct line information etc.
$(MACH)/%.o %.mod: %.f90  
	$(CPP) $(CPPFLAGS) $< $(MACH)/$<
	$(FCPRE) $(FC) -c $(FFLAGS) $(DBFLAGS) $(MACH)/$< -o $(MACH)/$(<:.f90=.o)

# special rules for "external" modules
pthreads_stuff.mod:
	ln -sf ../pthreads/pthreads_stuff.mod ./pthreads_stuff.mod

$(pepclib):   $(OBJS)
	@echo "Creating PEPC library ..."
	$(RM) $(MACH)/liblpepc.a
	$(AR) $(ARFLAGS) $(MACH)/liblpepc.a $^
	$(RANLIB) $(MACH)/liblpepc.a
	@echo "... done"


clean:
	$(RM) $(MACH)/liblpepc.a $(OBJS) $(MACH)/*.mod *.mod *% *~ bin/key rand *.o *.pif core module.deps
	$(RM) -r $(MACH)

module.deps:
	../tools/build/f90_mod_deps.py -o $@ --dep-template "\1.mod" --mod-template "\1.mod" --o-prefix "\$$(MACH)/" $(PEPCSRC.f90)

-include module.deps
