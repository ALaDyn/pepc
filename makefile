#
#  top-level makefile for pepc
#

ROOTDIR      = $(shell pwd)
SVNREVISION  = $(shell svnversion)

include tools/build/makefile.paths
include makefile.defs

ALLFRONTENDS = $(shell ls $(FRONTENDDIR))

help: info
	@printf $(HELP)
	@echo ""

info:
	@echo "======== make info"
	@echo "==== target architecture : $(MACH)"
	@echo "==== code version        : $(SVNREVISION)"
	@echo "==== pepc directory      : $(ROOTDIR)"
	@echo "==== available frontends : $(ALLFRONTENDS)"
	@echo ""

buildenv:
	@echo "======== build environment"
	@echo "== CPP      : $(CPP)"
	@echo "== CPPFLAGS : $(CPPFLAGS)"
	@echo "== FCPRE    : $(FCPRE)"
	@echo "== FC       : $(FC)"
	@echo "== FCFLAGS  : $(FCFLAGS)"
	@echo "== CCPRE    : $(CCPRE)"
	@echo "== CC       : $(CC)"
	@echo "== CCFLAGS  : $(CCFLAGS)"
	@echo "== LDPRE    : $(LDPRE)"
	@echo "== LD       : $(LD)"
	@echo "== LDFLAGS  : $(LDFLAGS)"
	@echo ""

readme:
	cat README | less

all:
	-$(MAKE) $(MFLAGS) $(ALLFRONTENDS)
	@echo ""
	@echo "======== build all results:"
	@echo "== available: " $(ALLFRONTENDS)
	@echo "== build    : " $(shell cd $(BINDIR) ; ls $(ALLFRONTENDS))
	@echo ""

libsl: $(LIBDIR)/libsl.a

$(LIBDIR)/libsl.a: $(LIBDIR)
	@echo "==== building libsl"
	@ln -sf $(ROOTDIR)/makefile.defs $(SLPEPCDIR)/makefile.defs
	@$(MAKE) -C $(SLPEPCDIR) $(MFLAGS)
	@cp -p $(SLPEPCDIR)/libsl.a $(LIBDIR)/libsl.a

clean:
	@echo "==== cleaning build and bin"
	@$(RM) makefile.envs
	@$(RM) $(BUILDDIR) $(BINDIR)
	@echo ""

cleanlib:
	@echo "==== cleaning libraries"
	@$(RM) $(LIBDIR)
	@ln -sf $(ROOTDIR)/makefile.defs $(SLPEPCDIR)/makefile.defs
	@cd src/treecode/sl_pepc && $(MAKE) $(MFLAGS) clean 
	@echo ""

cleanall: cleanlib clean
	@echo "==== all cleaned"

allclean: cleanall

pepc-%: pepclogo info buildenv $(LIBDIR)/libsl.a
	@echo "======== start building frontend { $@ }"
	@echo "==== date: " $(shell "date")
	@echo "==== make target: " $@
	@mkdir -p $(BUILDDIR)/$(MACH)/$@
	@mkdir -p $(BINDIR)
	@-$(RM) $(BINDIR)/$@
	@FRONTEND=$@ ROOTDIR=$(ROOTDIR) SVNREVISION=$(SVNREVISION) WORKDIR=$(BUILDDIR)/$(MACH)/$@ $(MAKE) $(MFLAGS) -f $(MAKEDIR)/makefile.prepare
	@cp -p $(BUILDDIR)/$(MACH)/$@/$@ $(BINDIR)
	@echo ""
	@echo "======== successfully built frontend { $@ } :-)"
	@echo ""

$(LIBDIR):
	@mkdir $(LIBDIR)


MAKEFILEDEFSINFO = "\n\n\
!!! To be able to build pepc, you first have to create a \n\
file called 'makefile.defs' inside the pepc root directory. \n\
A number of example files is contained in the\n\
'makefiles' directory. They usually can be used use by creating a symbolic link\n\n\
 >  ln -sf makefiles/makefile.defs.extension ./makefile.defs\n\n\
After creating this link to a certain file, a call of\n\n\
 >  make help\n\n\
might give you further information. Additionally, consider calling\n\n\
 >  make readme\n\n\
for a detailed description of what has to be done for getting pepc running.\n\
For some special configuration you might have to adapt\n\
the makefile.defs file to refelct your particular setup\n\
(copiler names, pathes, etc.)\n\
"

makefile.defs:
	@printf $(MAKEFILEDEFSINFO)
	@echo ""
	@exit 1

pepclogo:
	@echo "    ____    ____    ____    ____                                       "
	@echo "   /\  _\`\ /\  _\`\ /\  _\`\ /\  _\`\                       "
	@echo "   \ \ \L\ \ \ \L\_\ \ \L\ \ \ \/\_\      The Pretty Efficient      "
	@echo "    \ \ ,__/\ \  _\L\ \ ,__/\ \ \/_/_           Parallel Coulomb Solver "
	@echo "     \ \ \/  \ \ \L\ \ \ \/  \ \ \L\ \                                 "
	@echo "      \ \_\   \ \____/\ \_\   \ \____/           p.gibbon@fz-juelich.de"
	@echo "       \/_/    \/___/  \/_/    \/___/                                  "
	@echo ""
