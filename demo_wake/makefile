# 64-bit hardware

#
#  Makefile for treemp and test programs
#
# Generated automatically from Makefile.in by configure.
ALL: default
##### User configurable options #####
PAL_ROOT = /usr/local/vt


#ARCHFLAGS = -xtarget=ultra2i -xcache=16/32/1:2048/64/1 -xarch=v9a -xtypemap=real:64,double:128,integer:64  
#DBFLAGS = -g -Rabp
LDFLAGS = -lpmpi -Xm
# For apprentice
#FFLAGS =   -eA -m4 -O3,inline4
#LDFLAGS = -lpmpi -Xm -lapp

# For PAT analysis
FFLAGS = -m3 -O3 -Rapb
#ANALYS = -lpmpi -Xm -lpat pat.cld

# Visit libraries
VISITLIBS=-L /usr_FZJ/jzam04/jzam0401/lvisit_cray -llvisit_spk -L/usr/local/visit/lvisit/lib -llvisit -L/usr/local/visit/visit20/lib -lvisit 


FC = f90

# Source files

MODULES.f90 =  treevars.f90 utils.f90


SOURCES.f90 =	treemp.f90 \
		setup.f90 \
		configure.f90 randion.f90 predef_parts.f90 maxwell1.f90 scramble_v.f90 \
		make_keys.f90 make_domains.f90 \
		treebuild.f90  make_branches.f90 tree_fill.f90 tree_properties.f90 \
		tree_walk.f90 \
		potenergy.f90 sumpot.f90 \
		forces.f90 sum_force.f90 \
		make_hashentry.f90 get_address.f90  get_node.f90 \
	        diagnose_tree.f90 \
		draw_domains.f90 draw_tree2d.f90 draw_lists.f90\
		visualize.f90 \
		open_files.f90 close_files.f90		

# -----------------------
# Source combinations

SOURCES = $(MODULES.f90) \
	  $(SOURCES.f90) 

SRCKEY = $(MODULES.f90) \
	 $(SRCKRY.f90) 

# Dependencies

$(SOURCES): $(MODULES)

# Names of object files derived from sources

OBJECTS = $(SOURCES:.f90=.o)


#  compile mode
.f90.o:	treevars.f90
	$(FC) -c $(FFLAGS) $(DBFLAGS) $<

#  targets

rantest: 	randtest.o
	$(CF90)  $(FFLAGS) -o $@  randtest.o

treemp:	$(OBJECTS)
	@echo "Loading object files using flags $(FFLAGS)"
	$(CF90)  $(ANALYS) $(VISITLIBS) -o $@ $(OBJECTS) 

debug:	$(OBJECTS)
	@echo "Loading object files using flags $(DBFLAGS)"
	f90 $(LDFLAGS) -o $@ $(OBJECTS)


rand:	 rantest.o
	f90 $(FFLAGS) -o rand rantest.o


clean:
	/bin/rm -f bin/*.o  *% *~ bin/key rand *.o *.pif core debug treemp

cleandata: 
	/bin/rm -f fort.* *.dat *.out *.gle *.gz wf_*

default: treemp
