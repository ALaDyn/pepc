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
	@printf $(HELP)
	@echo -e ""

info:
	@echo -e "======== $(UL)make info$(NC)"
	@echo -e "==== target architecture : $(BC)$(MACH)$(NC)"
	@echo -e "==== code version        : $(BC)$(SVNREVISION)$(NC)"
	@echo -e "==== pepc directory      : $(BC)$(ROOTDIR)$(NC)"
	@echo -e "==== available frontends : $(BC)$(ALLFRONTENDS)$(NC)"
	@echo -e ""

buildenv:
	@echo -e "======== $(UL)build environment$(NC)"
	@echo -e "== CPP      : $(BC)$(CPP)$(NC)"
	@echo -e "== CPPFLAGS : $(BC)$(CPPFLAGS)$(NC)"
	@echo -e "== FCPRE    : $(BC)$(FCPRE)$(NC)"
	@echo -e "== FC       : $(BC)$(FC)$(NC)"
	@echo -e "== FFLAGS   : $(BC)$(FFLAGS)$(NC)"
	@echo -e "== CCPRE    : $(BC)$(CCPRE)$(NC)"
	@echo -e "== CC       : $(BC)$(CC)$(NC)"
	@echo -e "== CFLAGS   : $(BC)$(CFLAGS)$(NC)"
	@echo -e "== LDPRE    : $(BC)$(LDPRE)$(NC)"
	@echo -e "== LD       : $(BC)$(LD)$(NC)"
	@echo -e "== LDFLAGS  : $(BC)$(LDFLAGS)$(NC)"
	@echo -e ""

readme:
	cat README | less

all:
	-$(MAKE) $(MFLAGS) -k $(ALLFRONTENDS)
	@echo -e ""
	-$(MAKE) $(MFLAGS) allresult

allresult:
	@echo -e "======== $(UL)build all results$(NC)"
	@for f in $(ALLFRONTENDS); do if [ -e ${BINDIR}/$$f ]; then printf "== %-20s $(GREEN)OK$(NC)\n" $$f ; else printf "== %-20s $(RED)FAILED$(NC)\n" $$f; fi; done
	@echo -e ""

libsl: $(LIBDIR)/libsl.a

$(LIBDIR)/libsl.a: $(LIBDIR)
	@echo -e "==== $(UL)building libsl$(NC)"
	@ln -sf "$(ROOTDIR)/makefile.defs" "$(SLPEPCDIR)/makefile.defs"
	@$(MAKE) -C "$(SLPEPCDIR)" $(MFLAGS)
	@cp -p "$(SLPEPCDIR)/libsl.a" "$(LIBDIR)/libsl.a"

libopa: $(LIBDIR)
	@echo -e "==== $(UL)building openpa$(NC)"
	@ROOTDIR="$(ROOTDIR)" $(MAKE) -C "$(OPADIR)" $(MFLAGS)

clean:
	@echo -e "==== $(UL)cleaning build directory and binaries in bin directory$(NC)"
	@$(RM) makefile.envs
	@$(RM) "$(BUILDDIR)" $(addprefix $(BINDIR)/, $(ALLFRONTENDS))
	@echo -e ""

cleanlib:
	@echo -e "==== $(UL)cleaning libraries$(NC)"
	@$(RM) "$(LIBDIR)"
	@ln -sf "$(ROOTDIR)/makefile.defs" "$(SLPEPCDIR)/makefile.defs"
	@cd src/treecode/sl_pepc && $(MAKE) $(MFLAGS) clean 
	@ROOTDIR="$(ROOTDIR)" $(MAKE) -C "$(OPADIR)" $(MFLAGS) clean
	@echo -e ""

cleandoc:
	@echo -e "==== $(UL)cleaning Doxygen documentation$(NC)"
	@$(RM) "$(DOCDIR)"

cleanall: cleanlib cleandoc clean
	@-$(RM) "$(BINDIR)"
	@echo -e "==== $(UL)all cleaned$(NC)"

allclean: cleanall

pepc-%: pepclogo info buildenv $(LIBDIR)/libsl.a libopa
	@echo -e "======== start building frontend $(BC){ $@ }$(NC)"
	@echo -e "==== date: $(BC)$(shell "date")$(NC)"
	@echo -e "==== make target: $(BC)$@$(NC)"
	@mkdir -p "$(BUILDDIR)/$(MACH)/$@"
	@mkdir -p "$(BINDIR)"
	@-$(RM) "$(BINDIR)/$@"
	@FRONTEND="$@" ROOTDIR="$(ROOTDIR)" SVNREVISION="$(SVNREVISION)" WORKDIR="$(BUILDDIR)/$(MACH)/$@" $(MAKE) $(MFLAGS) -f "$(MAKEDIR)/makefile.prepare"
	@cp -p "$(BUILDDIR)/$(MACH)/$@/$@" "$(BINDIR)"
	@echo -e ""
	@echo -e "======== $(GREEN)successfully built frontend { $@ } :-)$(NC)"
	@echo -e ""

$(LIBDIR):
	@mkdir "$(LIBDIR)"

$(DOCDIR):
	@mkdir -p "$(DOCDIR)"

doc: $(DOCDIR) $(TOOLSDIR)/Doxyfile
	@echo -e "======== start building Doxygen documentation"
	@doxygen "$(TOOLSDIR)/Doxyfile"
	@echo -e "=== you can view the source code documentation by opening $(BOLD)$(DOCDIR)/index.html$(NC) with your favourite web browser"

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
	@echo "      \ \_\   \ \____/\ \_\   \ \____/           pepc@fz-juelich.de"
	@echo "       \/_/    \/___/  \/_/    \/___/                                  "
	@echo ""
