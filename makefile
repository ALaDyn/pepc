#
#  top-level makefile for pepc
#

ROOTDIR      = $(shell pwd)
SVNREVISION  = $(shell svnversion)

include tools/build/makefile.paths
include makefile.defs

ALLFRONTENDS = $(shell ls $(FRONTENDDIR))

RED=\e[0;31m
GREEN=\e[0;32m
BC=\e[1m# bold
UL=\e[4m# underline
NC=\e[0m# No Color, default font

help: info
	@printf "$(HELP)\n"
	@echo ""

info:
	@printf "======== $(UL)make info$(NC)\n"
	@printf "==== target architecture : $(BC)$(MACH)$(NC)\n"
	@printf "==== code version        : $(BC)$(SVNREVISION)$(NC)\n"
	@printf "==== pepc directory      : $(BC)$(ROOTDIR)$(NC)\n"
	@printf "==== available frontends : $(BC)$(ALLFRONTENDS)$(NC)\n"
	@echo ""

buildenv:
	@printf "======== $(UL)build environment$(NC)\n"
	@printf "== CPP      : $(BC)$(CPP)$(NC)\n"
	@printf "== CPPFLAGS : $(BC)$(CPPFLAGS)$(NC)\n"
	@printf "== FCPRE    : $(BC)$(FCPRE)$(NC)\n"
	@printf "== FC       : $(BC)$(FC)$(NC)\n"
	@printf "== FFLAGS   : $(BC)$(FFLAGS)$(NC)\n"
	@printf "== CCPRE    : $(BC)$(CCPRE)$(NC)\n"
	@printf "== CC       : $(BC)$(CC)$(NC)\n"
	@printf "== CFLAGS   : $(BC)$(CFLAGS)$(NC)\n"
	@printf "== LDPRE    : $(BC)$(LDPRE)$(NC)\n"
	@printf "== LD       : $(BC)$(LD)$(NC)\n"
	@printf "== LDFLAGS  : $(BC)$(LDFLAGS)$(NC)\n"
	@echo ""

readme:
	cat README | less

all:
	-$(MAKE) $(MFLAGS) -k $(ALLFRONTENDS)
	@echo ""
	-$(MAKE) $(MFLAGS) allresult

allresult:
	@printf "======== $(UL)build all results$(NC)\n"
	@for f in $(ALLFRONTENDS); do if [ -e ${BINDIR}/$$f ]; then printf "== %-20s $(GREEN)OK$(NC)\n" $$f ; else printf "== %-20s $(RED)FAILED$(NC)\n" $$f; fi; done
	@echo ""

libsl: $(LIBDIR)/libsl.a

$(LIBDIR)/libsl.a: $(LIBDIR)
	@printf "==== $(UL)building libsl$(NC)\n"
	@ln -sf "$(ROOTDIR)/makefile.defs" "$(SLPEPCDIR)/makefile.defs"
	@$(MAKE) -C "$(SLPEPCDIR)" $(MFLAGS)
	@cp -p "$(SLPEPCDIR)/libsl.a" "$(LIBDIR)/libsl.a"

libopa: $(LIBDIR)/libopa.a

$(LIBDIR)/libopa.a: $(LIBDIR)
	@printf "==== $(UL)building openpa$(NC)\n"
	@ROOTDIR="$(ROOTDIR)" $(MAKE) -C "$(OPADIR)" $(MFLAGS)

clean:
	@printf "==== $(UL)cleaning build directory and binaries in bin directory$(NC)\n"
	@$(RM) makefile.envs
	@$(RM) "$(BUILDDIR)" $(addprefix $(BINDIR)/, $(ALLFRONTENDS))
	@echo ""

cleanlib:
	@printf "==== $(UL)cleaning libraries$(NC)\n"
	@$(RM) "$(LIBDIR)"
	@ln -sf "$(ROOTDIR)/makefile.defs" "$(SLPEPCDIR)/makefile.defs"
	@cd src/treecode/sl_pepc && $(MAKE) $(MFLAGS) clean 
	@ROOTDIR="$(ROOTDIR)" $(MAKE) -C "$(OPADIR)" $(MFLAGS) clean
	@echo ""

cleandoc:
	@printf "==== $(UL)cleaning Doxygen documentation$(NC)\n"
	@$(RM) "$(DOCDIR)"

cleanall: cleanlib cleandoc clean
	@-$(RM) "$(BINDIR)"
	@printf "==== $(UL)all cleaned$(NC)\n"

allclean: cleanall

pepc-%: pepclogo info buildenv $(LIBDIR)/libsl.a $(LIBDIR)/libopa.a
	@if [ ! -d "$(FRONTENDDIR)/$@" ]; then printf "======== $(RED)Frontend $@ does not exist. Aborting.$(NC)\n"; false; fi
	@printf "======== start building frontend $(BC){ $@ }$(NC)\n"
	@printf "==== date: $(BC)$(shell "date")$(NC)\n"
	@printf "==== make target: $(BC)$@$(NC)\n"
	@mkdir -p "$(BUILDDIR)/$(MACH)/$@"
	@mkdir -p "$(BINDIR)"
	@-$(RM) "$(BINDIR)/$@"
	@FRONTEND="$@" ROOTDIR="$(ROOTDIR)" SVNREVISION="$(SVNREVISION)" WORKDIR="$(BUILDDIR)/$(MACH)/$@" $(MAKE) $(MFLAGS) -f "$(MAKEDIR)/makefile.prepare"
	@cp -p "$(BUILDDIR)/$(MACH)/$@/$@" "$(BINDIR)"
	@echo ""
	@printf "======== $(GREEN)successfully built frontend { $@ } :-)$(NC)\n"
	@echo ""

$(LIBDIR):
	@mkdir "$(LIBDIR)"

$(DOCDIR):
	@mkdir -p "$(DOCDIR)"

doc: $(DOCDIR) $(TOOLSDIR)/Doxyfile
	@printf "======== start building Doxygen documentation\n"
	@doxygen "$(TOOLSDIR)/Doxyfile"
	@printf "=== you can view the source code documentation by opening $(BOLD)$(DOCDIR)/index.html$(NC) with your favourite web browser\n"

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
	@echo "    _____   ____   _____   _____                              "
	@echo "   /\  _ \`\/\  __\/\  _ \`\/\  __\`\       The                  "
	@echo "   \ \ \L\ \ \ \_L\ \ \L\ \ \ \/\_\        Pretty Efficient   "
	@echo "    \ \ ,__/\ \  _\\ \ ,__/\ \ \/_/_        Parallel Coulomb   "
	@echo "     \ \ \/  \ \ \_L\ \ \/  \ \ \_\ \    Solver               "
	@echo "      \ \_\   \ \____\ \_\   \ \____/                         "
	@echo "       \/_/    \/____/\/_/    \/___/     pepc@fz-juelich.de   "
	@echo ""

