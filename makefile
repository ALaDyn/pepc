#
#  Makefile for pepc
#

include makefile.defs

export BACKEND

default: pepce

benchmark: pepce

help:
	@echo -e "## target architecture: $(MACH)"
	@echo -e "## code version:" `$(SVNVERSION)`
	@echo -e $(HELP)
	
readme:
	cat README | less
	
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

all: pepce pepcmini pepcmw pepcs pepcb pepcnn pepcv pepcsph


$(LIBDIR)libpthreads.a:
	@echo "============  Making PThreads Fortan wrapper library  ============="
	cd $(PTHREADSDIR) && $(MAKE) $(MFLAGS)

libsl:
	@echo "============  Making PEPC Sorting library  ============="
	cd $(SLPEPCDIR) && $(MAKE) $(MFLAGS)

libpepc.%: $(LIBDIR)libpthreads.a libsl
	$(eval BACKEND:=$*)
	@echo "============  Making PEPC Library - $(BACKEND) version ============="
	cd $(LPEPCDIR) && $(MAKE) $(MFLAGS)

pepcmw: libpepc.coulomb
	@echo "============  Making Frontend PEPC-MW (Mathias Winkel version)  ============="
	cd $(FRONTENDDIR)pepc-mw && $(MAKE) $(MFLAGS)

pepcb: libpepc.coulomb
	@echo "============  Making Frontend PEPC-B (Laser/beam-plasma with magnetic fields)  ============="
	cd $(FRONTENDDIR)pepc-b && $(MAKE) $(MFLAGS)

pepcmini: libpepc.coulomb
	@echo "============  Making Frontend PEPC-MINI (minial version)  ============="
	cd $(FRONTENDDIR)pepc-mini && $(MAKE) $(MFLAGS)

pepcnn: libpepc.nearestneighbour
	@echo "============  Making Frontend PEPC-NN (nearest neighbour search)  ============="
	cd $(FRONTENDDIR)pepc-nn && $(MAKE) $(MFLAGS)

pepcsph: libpepc.nearestneighbour
	@echo "============  Making Frontend PEPC-SPH (next neighbour + SPH)  ============="
	cd $(FRONTENDDIR)pepc-sph && $(MAKE) $(MFLAGS)

pepce: libpepc.coulomb
	@echo "============  Making Frontend PEPC-E (Benchmark version)  ============="
	cd $(FRONTENDDIR)pepc-e && $(MAKE) $(MFLAGS)

pepcs: libpepc.coulomb
	@echo "============  Making Frontend PEPC-S (ScaFaCoS-library version + minimal frontend)  ============="
	cd $(FRONTENDDIR)pepc-s && $(MAKE) $(MFLAGS)

pepcv: libpepc.vortex
	@echo "============  Making Frontend PEPC-V (Vortex version)  ============="
	cd $(FRONTENDDIR)pepc-v && $(MAKE) $(MFLAGS) 

clean: clean-doc
	cd $(SLPEPCDIR)   && $(MAKE) $(MFLAGS) clean
	cd $(PTHREADSDIR) && $(MAKE) $(MFLAGS) clean
	cd $(LPEPCDIR)    && $(MAKE) $(MFLAGS) clean
	cd $(FRONTENDDIR)pepc-s     && $(MAKE) $(MFLAGS) clean
	cd $(FRONTENDDIR)pepc-e     && $(MAKE) $(MFLAGS) clean
	cd $(FRONTENDDIR)pepc-b     && $(MAKE) $(MFLAGS) clean
	cd $(FRONTENDDIR)pepc-mw    && $(MAKE) $(MFLAGS) clean
	cd $(FRONTENDDIR)pepc-mini  && $(MAKE) $(MFLAGS) clean
	cd $(FRONTENDDIR)pepc-nn    && $(MAKE) $(MFLAGS) clean
	cd $(FRONTENDDIR)pepc-sph   && $(MAKE) $(MFLAGS) clean
	cd $(FRONTENDDIR)pepc-v     && $(MAKE) $(MFLAGS) clean

clean-doc:
	rm -rf ./doc

doc: clean-doc
	mkdir ./doc
	doxygen ./tools/Doxyfile
	@echo "--- you can view the source code documentation by opening ./doc/index.html with your favourite web browser ---"

clean-dist:
	rm -rf ./benchmark
	rm -f benchmark-VER.tgz

dist: clean-dist	
	@echo "--- exporting svn directory structure ---"
	svn export ./ ./benchmark
	rm -rf ./benchmark/jube
	@echo "--- creating tarball ---"
	tar -czvf ./benchmark-VER.tgz ./benchmark/
	@echo "--- removing temporary files ---"
	rm -rf ./benchmark
	@echo "--- before publishing do not forget to update revision number in filename ---"

.PHONY: readme
