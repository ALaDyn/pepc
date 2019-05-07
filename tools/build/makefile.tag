#
#  makefile for pepc to tag source files according to chosen frontend
#
SUBFILES=makefile.paths makefile.defs makefile.frontend

include $(SUBFILES)
include makefile.envs

CPPFLAGS += $(CPPFLAGS_BACKEND)
CPPFLAGS += -DFFLAGS="\"$(FFLAGS)\""     \
            -DCFLAGS="\"$(CFLAGS)\""     \
            -DLDFLAGS="\"$(LDFLAGS)\""   \
            -DCOMPILER="\"$(COMPILER)\"" \
            -DSVNREVISION="\"$(SVNREVISION)\"" \
            -DMACH="\"$(MACH)\"" \
            -DWALKALGORITHM="\"$(WALK)\""

FSRC = $(filter %.f90, $(SRC))
FPRC = $(FSRC:%.f90=pp_%.f90)
FOBJ = $(FSRC:%.f90=%.o)

CSRC = $(filter %.c, $(SRC))
COBJ = $(CSRC:%.c=%.o)

RED=\e[0;31m
GREEN=\e[0;32m
BC=\e[1m# bold
UL=\e[4m# underline
NC=\e[0m# No Color, default font

both_tags: $(ROOTDIR)/src/frontends/$(FRONTEND)/tags tags

$(ROOTDIR)/src/frontends/$(FRONTEND)/tags: $(FOBJ) $(COBJ)
	@printf "==== tagging source files for frontend : $(BC)$(FRONTEND)$(NC)\n"
	@cd $(ROOTDIR)/src/frontends/$(FRONTEND); ctags -R ../../frontends/$(FRONTEND)/* ../../interaction_specific/$(BACKEND)/* ../../treecode/* ../../utils/*

tags: $(FOBJ) $(COBJ)
	@printf "==== tagging source files in build directory\n"
	@ctags -R $(FPRC) $(CSRC)

