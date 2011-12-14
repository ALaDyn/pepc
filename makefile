#
#  top-level makefile for pepc
#

include makefiles/makefile.paths
include makefile.defs

ROOTDIR     = $(shell "pwd")
SVNREVISION = $(shell "svnversion")

help: info
	@echo -e $(HELP)

info:
	@echo -e "======== make info"
	@echo -e "==== target architecture: $(MACH)"
	@echo -e "==== code version: $(SVNREVISION)"
	@echo -e "==== pepc directory: $(ROOTDIR)"
	@echo -e ""

buildenv:
	@echo -e "======== build environment"
	@echo -e "== CPP      : $(CPP)"
	@echo -e "== CPPFLAGS : $(CPPFLAGS)"
	@echo -e "== FCPRE    : $(FCPRE)"
	@echo -e "== FC       : $(FC)"
	@echo -e "== FCFLAGS  : $(FCFLAGS)"
	@echo -e "== CCPRE    : $(CCPRE)"
	@echo -e "== CC       : $(CC)"
	@echo -e "== CCFLAGS  : $(CCFLAGS)"
	@echo -e "== LDPRE    : $(LDPRE)"
	@echo -e "== LD       : $(LD)"
	@echo -e "== LDFLAGS  : $(LDFLAGS)"
	@echo -e ""

readme:
	cat README | less

libsl: 
$(LIBDIR)/libsl.a: $(LIBDIR)
	@echo -e "==== building libsl"
	@ln -sf $(ROOTDIR)/makefile.defs $(SLPEPCDIR)/makefile.defs
	@$(MAKE) -C $(SLPEPCDIR) $(MFLAGS)
	@cp -p $(SLPEPCDIR)/libsl.a $(LIBDIR)/libsl.a

clean:
	@echo "==== cleaning build and bin"
	@$(RM) makefile.envs
	@$(RM) $(BUILDDIR) $(BINDIR)
	@echo -e ""

cleanlib:
	@echo "==== cleaning libraries"
	@$(RM) $(LIBDIR)
	@ln -sf $(ROOTDIR)/makefile.defs $(SLPEPCDIR)/makefile.defs
	@cd src/treecode/sl_pepc && $(MAKE) $(MAKEFLAGS) clean 
	@echo -e ""

cleanall: cleanlib clean
	@echo "==== all cleaned"

allclean: cleanall

pepc-%: pepclogo info buildenv $(LIBDIR)/libsl.a
	@echo "======== start building frontend ** $@ **"
	@echo "==== date: " $(shell "date")
	@echo "==== make target: " $@
	@$(RM) makefile.envs
	@echo "FRONTEND=$@" >> makefile.envs
	@echo "ROOTDIR=$(ROOTDIR)" >> makefile.envs
	@echo "SVNREVISION=$(SVNREVISION)" >> makefile.envs
	@echo "WORKDIR=$(BUILDDIR)/$(MACH)/$@" >> makefile.envs
	@mkdir -p $(BUILDDIR)/$(MACH)/$@
	@$(MAKE) $(MFLAGS) -f $(MAKEDIR)/makefile.prepare
	@mkdir -p $(BINDIR)
	@cp -p $(BUILDDIR)/$(MACH)/$@/$@ $(BINDIR)
	@echo -e ""
	@echo "======== successfully build frontend ** $@ **"
	@echo -e ""

$(LIBDIR):
	@mkdir $(LIBDIR)


MAKEFILEDEFSINFO = "\n\n !!! To create be able to build pepc, you first have to create a file called makefile.defs\n\
inside the pepc root directory. The makefiles directory contains a number of \n\
samples, which you usually can use via\n\n\
 >  ln -sf makefiles/makefile.defs.extension ./makefile.defs\n\n\
After creating this link to a certain file, a call of\n\n\
 >  make help\n\n\
might give you further information. Additionally, consider calling\n\n\
 >  make readme\n\n\
for a detailed description of what has to be done for getting pepc running.\n\
"

makefile.defs:
	@echo -e $(MAKEFILEDEFSINFO)
	@exit 1

pepclogo:
	@echo -e "    ____    ____    ____    ____                                       "
	@echo -e "   /\  _\`\ /\  _\`\ /\  _\`\ /\  _\`\                       "
	@echo -e "   \ \ \L\ \ \ \L\_\ \ \L\ \ \ \/\_\      The Pretty Efficient      "
	@echo -e "    \ \ ,__/\ \  _\L\ \ ,__/\ \ \/_/_           Parallel Coulomb Solver "
	@echo -e "     \ \ \/  \ \ \L\ \ \ \/  \ \ \L\ \                                 "
	@echo -e "      \ \_\   \ \____/\ \_\   \ \____/           p.gibbon@fz-juelich.de"
	@echo -e "       \/_/    \/___/  \/_/    \/___/                                  "
	@echo -e ""
