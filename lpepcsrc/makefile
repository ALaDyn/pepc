#
#  Makefile for treemp and test programs
#

include ../makefile.defs

pepclib = liblpepc.a

.SUFFIXES: .f90
# -------------------------------------------------------------------------------------------

# Source files

#-------------------------------------------------------------------------------------------
# Pepc library

PEPCMOD.f90 =   treevars.f90 utils.f90 tree_utils.f90 timings.f90 treetypes.f90 draw_gle.f90
PEPCSRC.f90 =   tree_domains.f90 setup_treearrays.f90 restore.f90 tree_allocate.f90 tree_deallocate.f90 \
		tree_build.f90  tree_branches.f90 tree_fill.f90 tree_properties.f90 \
		tree_aswalk.f90 tree_prefetch.f90 sum_force.f90 sum_bfield.f90 \
	        tree_update.f90 tree_stats.f90 tree_local.f90 tree_exchange.f90 tree_global.f90 \
		make_hashentry.f90 get_address.f90 get_next_node.f90 \
	        diagnose_tree.f90 check_table.f90 sum_lennardjones.f90 sum_bond.f90 \
		fields.f90 fields_p.f90 cleanup_treearrays.f90 mac.f90

LIBOBJS = $(addprefix $(MACH)/, $(PEPCSRC.f90:.f90=.o))
LIBMODS = $(addprefix $(MACH)/, $(PEPCMOD.f90:.f90=.o))
POBJS = $(PEPCSRC.f90:.f90=.o)
PMODS =  $(PEPCMOD.f90:.f90=.o)


all: dirs mods lpepc

dirs:  
	if [ ! -d $(MACH) ]; then mkdir $(MACH); fi

#  compile mode for PEPC library routines
$(MACH)/%.o: %.f90  
	$(CPP) $(CPPFLAGS) $< $(MACH)/$<
	$(FCPRE) $(FC) -c $(FFLAGS) $(DBFLAGS) $(MACH)/$< -o $@

$(pepclib): mods lpepc
lpepc:   $(LIBOBJS) $(LIBMODS) 
	@echo "Creating PEPC library ..."
	$(RM) $(MACH)/liblpepc.a
	$(AR) $(ARFLAGS) $(MACH)/liblpepc.a $^
	$(RANLIB) $(MACH)/liblpepc.a
	if [ -e treevars.mod ]; then cp *.mod $(MACH)/.; fi
	@echo "... done"


mods:   $(LIBMODS) 

debug:	$(OBJECTS)
	@echo "Loading object files using flags $(DBFLAGS)"
	$(FC) $(LDFLAGS)  $(VISITLIBS) $(DBFLAGS)  -o $@ $(OBJECTS)

clean:
	$(RM) $(MACH)/liblpepc.a $(LIBOBJS) $(LIBMODS) $(MACH)/*.mod *.mod *% *~ bin/key rand *.o *.pif core *.ps module.deps
	$(RM) -r $(MACH)

cleanps: 
	$(RM) *.ps 

%.mod: $(MACH)/%.o
	@true

module.deps:
	../tools/build/f90_mod_deps.py -o $@ --o-prefix "\$$(MACH)/" $(PEPCMOD.f90) $(PEPCSRC.f90)

-include module.deps
