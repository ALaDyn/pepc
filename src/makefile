# 64-bit hardware

#
#  Makefile for treemp and test programs
#


#DBFLAGS = -g -Rabp -ev  -q64 -qrealsize=8  -qsuffix=f=f90 -qcheck
#LDFLAGS = -q64 -qrealsize=8  -g 
# For apprentice
#FFLAGS =   -eA -m4 -O2,inline4
#LDFLAGS = -lpmpi -Xm -lapp

## For IBM AIX compiler
#FFLAGS = -g -q64 -qrealsize=8 -qsuffix=f=f90 -qcheck  
## Optimised
#FFLAGS =  -q64 -qrealsize=8 -qsuffix=f=f90  -O4 -Q -qipa -qarch=pwr4 -qstrict -qreport  -qlist -qsource -I/usr/local/include 
FZJ = -lllapi -lfzj

IPA=-qipa=inline=key2addr -qipa=inline=make_hashentry -qipa=inline=key2node -qipa=inline=next_node

PG=-pg
DB= -g -qfullpath -qcheck
HPM = -lhpm
LISTING=-qreport -qlist -qlistopt -qsource 
#MPITRACE= -L/usr/local/beta/lib -lmpitrace
#MPITRACE= -L/usr/local/beta/lib -lmpiprof
#MPITRACE=  -L/usr/local/beta/lib -lmpihpm -lpmapi
#MPITRACE= -lmpitrace
#MPITRACE= -lmpiprof
#MPITRACE= -lsummary -lpmapi
#MPITRACE=  -lmpihpm -lpmapi
Q64=-q64
FFLAGS1 = -qrealsize=8 -qsuffix=f=f90 -qsuffix=cpp=F90  
TUNE= -qarch=pwr4 -qtune=pwr4 -O5 $(IPA) 
# $(PG) $(HPM)

#FFLAGS =  $(FFLAGS1) $(DB) $(Q64) $(PG)
#FFLAGS =  $(FFLAGS1) $(DB) $(Q64) ##$(PG) 
FFLAGS =   $(FFLAGS1) $(DB) $(Q64) $(TUNE)
#TFLAGS = $(FFLAGS) $(TUNE)

#ANALYS = -lpmpi -Xm -lpat pat.cld

# Visit libraries
# JUMP

VISITLIBS=-L/usr/local/beta/visit-2.0b/lvisit/apis/spk4 -llvisit_spk -L/usr/local/beta/visit-2.0b/lvisit/lib -llvisit -L/usr/local/beta/visit-2.0b/lib -lvisit
#VISITLIBS=-L/usr/local/beta/visit-2.0b/lvisit/apis/spk5  -llvisit_spk -L/usr/local/beta/visit-2.0b/lvisit/lib -llvisit -L/usr/local/beta/visit-2.0b/lib -lvisit -bloadmap:aa

#VISITLIBS=-L/home4/jzam04/jzam0401/tree/pepc/src  -llvisit_spk -L/usr/local/beta/visit-2.0b/lvisit/lib -llvisit -L/usr/local/beta/visit-2.0b/lib -lvisit -bloadmap:aa

LDFLAGS = $(VISITLIBS) $(MPITRACE)
RANLIB  = ranlib
libdir  = .
AR      = ar -X64
LIBPEPC = -L$(libdir) -llpepc
pepclib = liblpepc.a

.SUFFIXES: .f90


FC = mpxlf90_r

# KOJAK instrumentation
#FC = mpxlf90_r -qdebug=function_trace
#FC = kinst mpxlf90_r

# -------------------------------------------------------------------------------------------

# Source files

#-------------------------------------------------------------------------------------------
# Pepc library

PEPCMOD.f90 =  treevars.f90 utils.f90 
PEPCSRC.f90 =  make_domains.f90 setup_arrays.f90 \
		tree_build.f90  make_branches.f90 tree_fill.f90 tree_properties.f90 \
		tree_walk.f90	tree_prefetch.f90 sum_force.f90 sum_bfield.f90 \
		make_hashentry.f90 get_address.f90  get_node.f90 key_test.f90 get_next_node.f90 \
	        diagnose_tree.f90 check_table.f90 \
		draw_domains.f90 draw_tree2d.f90 draw_lists.f90

pepcobjs = $(PEPCSRC.f90:.f90=.o)
pepcmods =  $(PEPCMOD.f90:.f90=.o)

#  compile mode for PEPC routines 
#
$(pepcmods):   
	$(FC) -c $(FFLAGS)  $<
$(pepcobjs):   
	$(FC) -c $(FFLAGS)  $<

# -------------------------------------------------------------------------------------------
MODULES.f90 =   physvars.f90 

SOURCES.f90 =	treemp.f90 \
		setup.f90 stamp.f90 \
		configure.f90 randion.f90 predef_parts.f90 maxwell1.f90 cold_start.f90 \
		special_start.f90 constrain.f90 cutvector.f90 face.f90 \
		scramble_v.f90 reset_ions.f90 add_electrons.f90 double_target.f90 \
		beam.f90 beam_control.f90 beam_dust.f90 \
		mc_config.f90 \
		velocities.f90 \
	 	energy_cons.f90 potenergy.f90 kinenergy.f90\
		forces.f90 sum_lennardjones.f90\
		fpond.f90 empond.f90 emplane.f90 laser_bullet.f90 laser.f90 \
		push_em.f90 push_full3v.f90 \
		push.f90 constrain.f90 reinject.f90 earth_plate.f90 add_ramp.f90 \
		vis_parts.f90 vis_fields.f90 dump.f90 param_dump.f90 diagnostics.f90 \
		open_files.f90 close_files.f90	\
		densities.f90 sum_fields.f90 dump_fields.f90 \
		track_nc.f90 rezone.f90		


# -------------------------------------------------------------------------------------------

# Names of application object files derived from sources

appobjs = $(SOURCES.f90:.f90=.o)
appmods =  $(MODULES.f90:.f90=.o)
mods = $(appmods) $(pepcmods)
all: mods pepc

#  compile mode
.f90.o:   
	$(FC) -c $(FFLAGS) $(VISITLIBS) $(DBFLAGS) $<

# -----------------------------------------------------------------------------------------
#  targets

pepc:	$(appobjs) $(appmods) $(pepcmods) $(pepclib)
	@echo "\n Creating pepc application binary ...\n"
	$(FC)  $(FFLAGS) $(LDFLAGS) $(LIBPEPC) -o $@ $(appobjs) $(appmods) $(pepcmods)
	@echo "... done\n"

pepc2:	$(appobjs) $(appmods) $(pepcmods) $(pepclib)
	@echo "\n Creating pepc application binary ...\n"
	$(FC)  $(FFLAGS) $(LDFLAGS) $(LIBPEPC) -o $@ $(appobjs) $(appmods) $(pepcmods)
	@echo "... done\n"

$(pepclib): lpepc
lpepc:  $(pepcobjs) $(pepcmods)
	@echo "\n Creating PEPC library ...\n"
	rm -f liblpepc.a
	$(AR) q liblpepc.a $(pepcobjs)
	$(RANLIB) liblpepc.a
	@echo "... done\n"

mods:   $(MODS) 


pepck:	$(OBJECTS) 
	$(FC)  $(FFLAGS) $(LDFLAGS) -o $@ $(OBJECTS) $(MODS) 

debug:	$(OBJECTS)
	@echo "Loading object files using flags $(DBFLAGS)"
	$(FC) $(LDFLAGS)  $(VISITLIBS) $(DBFLAGS)  -o $@ $(OBJECTS)


rand:	 rantest.o
	$(FC) $(FFLAGS) -o rand rantest.o


clean:
	/bin/rm -f bin/*.o  *% *~ bin/key rand *.o *.pif core *.mod 

cleanbin: /bin/rm -r pepc debug core

cleancore: 	
	/bin/rm $(pepcobjs) 

cleandata: 
	/bin/rm -f fort.* *.dat *.out *.gle 

depend:
	./f90depend -u -I/usr/lpp/ppe.poe/include/thread *.f90
default: treemp

### DO NOT DELETE OR EDIT THIS LINE
# Everything from here on down is generated by f90depend
# so do NOT add any translation rules below here.

#--- Include-Dependencies

treevars.o: /usr/lpp/ppe.poe/include/thread/mpif.h


#--- Module-Dependencies

add_electrons.o: physvars.o treevars.o utils.o
add_ramp.o: physvars.o treevars.o utils.o
beam.o: physvars.o treevars.o utils.o
beam_control.o: physvars.o treevars.o utils.o
beam_dust.o: physvars.o treevars.o utils.o
boundaries.o: physvars.o treevars.o
check_table.o: treevars.o utils.o
close_files.o: treevars.o
cold_start.o: treevars.o utils.o
configure.o: physvars.o treevars.o utils.o
constrain.o: physvars.o treevars.o utils.o
cutvector.o: physvars.o treevars.o
densities.o: physvars.o treevars.o
diagnose_tree.o: treevars.o
diagnostics.o: physvars.o treevars.o utils.o
double_target.o: physvars.o treevars.o utils.o
draw2d.o: treevars.o
draw2d_hash.o: treevars.o
draw_domains.o: treevars.o
draw_lists.o: treevars.o
draw_map.o: treevars.o
draw_tree2d.o: treevars.o
driver.o: utils.o my_mpidefs.o
dump.o: physvars.o treevars.o
dump_fields.o: physvars.o treevars.o
earth_plate.o: physvars.o treevars.o utils.o
energy_cons.o: physvars.o treevars.o utils.o
face.o: physvars.o treevars.o utils.o
forces.o: physvars.o treevars.o utils.o
get_address.o: treevars.o utils.o
get_next_node.o: treevars.o utils.o
get_node.o: treevars.o utils.o
key_test.o: treevars.o utils.o
kinenergy.o: physvars.o treevars.o
make_branches.o: treevars.o utils.o
make_domains.o: treevars.o utils.o
make_hashentry.o: treevars.o utils.o
mc_config.o: physvars.o treevars.o utils.o
open_files.o: treevars.o
param_dump.o: physvars.o treevars.o utils.o
potenergy.o: physvars.o treevars.o utils.o
predef_parts.o: physvars.o treevars.o utils.o
push.o: physvars.o treevars.o
push_em.o: physvars.o treevars.o
push_full3v.o: physvars.o treevars.o
randion.o: physvars.o treevars.o utils.o
reinject.o: utils.o
reset_ions.o: physvars.o treevars.o
scramble_v.o: physvars.o treevars.o utils.o
setup.o: physvars.o treevars.o utils.o
setup_arrays.o: utils.o treevars.o physvars.o
special_start.o: physvars.o treevars.o utils.o
sum_bond.o: physvars.o treevars.o
sum_fields.o: physvars.o treevars.o
sum_force.o: treevars.o
sum_lennardjones.o: physvars.o treevars.o
track_nc.o: physvars.o treevars.o
tree_build.o: treevars.o utils.o
tree_fill.o: treevars.o utils.o
tree_prefetch.o: treevars.o utils.o
tree_properties.o: treevars.o utils.o
tree_walk.o: treevars.o utils.o
tree_walk_async.o: treevars.o utils.o
treemp.o: physvars.o treevars.o utils.o
treevars.o: my_mpidefs.o
utils.o: my_mpidefs.o
velocities.o: physvars.o treevars.o utils.o
vis_fields.o: physvars.o treevars.o
vis_parts.o: physvars.o treevars.o


