#
#  Makefile for treemp and test programs
#

PREFIX = ..

#  Get machine-dependent flags

include ../makefile.defs

FFLAGS += -I../lpepcsrc/$(MACH)

LIBPEPC   = -L../lpepcsrc/$(MACH) -llpepc ../sl_pepc/libsl.a
LIBS      = $(LIBPEPC)

.SUFFIXES: .f90


# -------------------------------------------------------------------------------------------

# Source files

# -------------------------------------------------------------------------------------------

MODULES.f90 =   physvars.f90 

SOURCES.f90 =	\
                setup.f90 setup_arrays.f90 stamp.f90 \
                plasma_start.f90 predef_parts.f90 maxwell1.f90 cold_start.f90 \
                special_start.f90 constrain.f90 cutvector.f90 face.f90 \
                scramble_v.f90 reset_ions.f90 add_electrons.f90 stretch_sphere.f90 \
                cluster_sa.f90 laser.f90 beam.f90 beam_dust.f90 rezone.f90\
                mc_config.f90 \
                velocities.f90 perturb_temp.f90 \
                energy_cons.f90 potenergy.f90 kinenergy.f90\
                force_laser.f90 emplane_lin.f90 \
                fpond.f90 fpond_lin.f90 empond.f90 emplane.f90 laser_bullet.f90 emoblique.f90\
		fpond_helm.f90 fpond_helmc.f90\
		density_helmholtz.f90 em_helmholtz.f90 trisolve.f90 \
                push_em.f90 push_full3v.f90 \
                push.f90 reinject.f90 earth_plate.f90 add_ramp.f90 \
                dump.f90 param_dump.f90 \
                open_files.f90 close_files.f90  \
                densities.f90 sum_fields.f90 dump_fields.f90 sum_radial.f90 \
                track_nc.f90 sum_fieldave.f90 \
                integrator.f90 abort.f90 push_nonrel.f90 \
	        error_test.f90 force_direct.f90 laser_hist.f90 \
		pepcb.f90 configure.f90 diagnostics.f90 	

VISIT.f90   =   vis_parts_nbody.f90 vis_fields_nbody.f90 vis_control.f90 \
		vis_vecfields_nbody.f90 vis_config.f90


# -------------------------------------------------------------------------------------------

# Names of application object files derived from sources
# Prefix added for multi-arch builds

appobjs = $(addprefix $(MACH)/, $(SOURCES.f90:.f90=.o)) \
	  $(addprefix $(MACH)/, $(NCOBJS) )
appmods =  $(addprefix $(MACH)/, $(MODULES.f90:.f90=.o))
visitobjs = $(addprefix $(MACH)/, $(VISIT.f90:.f90=.o) )
mods = $(appmods)

default: dirs module.deps pepcb
all: dirs module.deps pepcb
dirs: 	
	if [ ! -d $(MACH) ]; then mkdir $(MACH); fi

#  compile rules
$(MACH)/%.o: %.f90  
	$(CPP) $(CPPFLAGS) $< $(MACH)/$<
	$(FCPRE) $(FC) -c $(FFLAGS) $(DEFINES) $(MACH)/$< -o $@


#  explicit rules for f90 modules to avoid conflict with built-in ones
$(MACH)/physvars.o: physvars.f90
	$(CPP) $(CPPFLAGS) physvars.f90 $(MACH)/physvars.f90
	$(FCPRE) $(FC) -c $(FFLAGS) $(DEFINES) $(MACH)/physvars.f90 -o $@

$(MACH)/%.o: %.c  
	$(CC) -c $(CFLAGS) $< -o $@


# -----------------------------------------------------------------------------------------
#  targets
pepcb:	lpepc $(appmods) $(appobjs)
	@echo "Creating application binary $(EXECNAME) ..."
	$(LDPRE) $(LD) $(LDFLAGS) -o $(EXECNAME) $(appobjs) $(appmods) $(LIBPEPC) $(LDLIBS)
	@echo "... done"

# Update PEPC library and headers
lpepc:  ../lpepcsrc/$(MACH)/liblpepc.a
#	cd ../lpepcsrc; make lpepc; cd ../pepc-b
	ln -sf $(PREFIX)/lpepcsrc/$(MACH)/treevars.mod $(PREFIX)/pepc-b/treevars.mod
	ln -sf $(PREFIX)/lpepcsrc/$(MACH)/tree_utils.mod $(PREFIX)/pepc-b/tree_utils.mod
	ln -sf $(PREFIX)/lpepcsrc/$(MACH)/utils.mod $(PREFIX)/pepc-b/utils.mod
	ln -sf $(PREFIX)/lpepcsrc/$(MACH)/timings.mod $(PREFIX)/pepc-b/timings.mod

clean:
	$(RM) $(MACH)/*.o  *% *~ bin/key rand *.o *.pif core *.mod $(MACH)/*.mod $(EXECNAME) module.deps
	$(RM) -r $(MACH)

default: pepcb

module.deps:
	./f90_mod_deps.py -o $@ --o-prefix "\$$(MACH)/" $(MODULES.f90) $(SOURCES.f90) $(VISIT.f90)

-include module.deps




