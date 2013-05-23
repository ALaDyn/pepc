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
NC=\e[0m # No Color

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
	@echo "== FFLAGS   : $(FFLAGS)"
	@echo "== CCPRE    : $(CCPRE)"
	@echo "== CC       : $(CC)"
	@echo "== CFLAGS   : $(CFLAGS)"
	@echo "== LDPRE    : $(LDPRE)"
	@echo "== LD       : $(LD)"
	@echo "== LDFLAGS  : $(LDFLAGS)"
	@echo ""

readme:
	cat README | less

all:
	-$(MAKE) $(MFLAGS) $(ALLFRONTENDS)
	@echo ""
	-$(MAKE) $(MFLAGS) allresult

allresult:
	@echo "======== build all results:"
	@for f in $(ALLFRONTENDS); do if [ -e ${BINDIR}/$$f ]; then printf "== %-20s $(GREEN)OK$(NC)\n" $$f ; else printf "== %-20s $(RED)FAILED$(NC)\n" $$f; fi; done
	@echo ""

libsl: $(LIBDIR)/libsl.a

$(LIBDIR)/libsl.a: $(LIBDIR)
	@echo "==== building libsl"
	@ln -sf $(ROOTDIR)/makefile.defs $(SLPEPCDIR)/makefile.defs
	@$(MAKE) -C $(SLPEPCDIR) $(MFLAGS)
	@cp -p $(SLPEPCDIR)/libsl.a $(LIBDIR)/libsl.a

libopa: $(LIBDIR)
	@echo "==== building openpa"
	@ROOTDIR=$(ROOTDIR) $(MAKE) -C $(OPADIR) $(MFLAGS)

clean:
	@echo "==== cleaning build directory and binaries in bin directory"
	@$(RM) makefile.envs
	@$(RM) $(BUILDDIR) $(addprefix $(BINDIR)/, $(ALLFRONTENDS))
	@echo ""

cleanlib:
	@echo "==== cleaning libraries"
	@$(RM) $(LIBDIR)
	@ln -sf $(ROOTDIR)/makefile.defs $(SLPEPCDIR)/makefile.defs
	@cd src/treecode/sl_pepc && $(MAKE) $(MFLAGS) clean 
	@ROOTDIR=$(ROOTDIR) $(MAKE) -C $(OPADIR) $(MFLAGS) clean
	@echo ""

cleandoc:
	@echo "==== cleaning Doxygen documentation"
	@$(RM) $(DOCDIR)

cleanall: cleanlib cleandoc clean
	@-$(RM) $(BINDIR)
	@echo "==== all cleaned"

allclean: cleanall

pepc-%: pepclogo info buildenv $(LIBDIR)/libsl.a libopa
	@echo "======== start building frontend { $@ }"
	@echo "==== date: " $(shell "date")
	@echo "==== make target: " $@
	@mkdir -p $(BUILDDIR)/$(MACH)/$@
	@mkdir -p $(BINDIR)
	@-$(RM) $(BINDIR)/$@
	@FRONTEND=$@ ROOTDIR=$(ROOTDIR) SVNREVISION=$(SVNREVISION) WORKDIR=$(BUILDDIR)/$(MACH)/$@ $(MAKE) $(MFLAGS) -f $(MAKEDIR)/makefile.prepare
	@ln -s $(BUILDDIR)/$(MACH)/$@/$@ $(BINDIR)
	@echo ""
	@echo "======== successfully built frontend { $@ } :-)"
	@echo ""

$(LIBDIR):
	@mkdir $(LIBDIR)

$(DOCDIR):
	@mkdir -p $(DOCDIR)

doc: $(DOCDIR) $(TOOLSDIR)/Doxyfile
	@echo "======== start building Doxygen documentation"
	@doxygen $(TOOLSDIR)/Doxyfile
	@echo "=== you can view the source code documentation by opening $(DOCDIR)/index.html with your favourite web browser"

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
	@echo '    _____   ____   _____   _____                              '
	@echo '   /\  _ `\/\  __\/\  _ `\/\  __`\       The                  '
	@echo '   \ \ \L\ \ \ \_L\ \ \L\ \ \ \/\_\        Pretty Efficient   '
	@echo '    \ \ ,__/\ \  _\\ \ ,__/\ \ \/_/_       Parallel Coulomb   '
	@echo '     \ \ \/  \ \ \_L\ \ \/  \ \ \_\ \    Solver               '
	@echo '      \ \_\   \ \____\ \_\   \ \____/                         '
	@echo '       \/_/    \/____/\/_/    \/___/     pepc@fz-juelich.de   '
	@echo ''

