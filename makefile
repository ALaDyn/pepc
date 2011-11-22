#
#  Makefile for pepc
#

include makefile.defs

export BACKEND

default: pepce

benchmark: pepce

info:
	@echo $(LIBDIR)

all: pepce pepcmini pepcmw pepcs pepcb


libpthreads:
	@echo "============  Making PThreads Fortan wrapper library  ============="
	cd $(PTHREADSDIR) && $(MAKE) $(MFLAGS)

libsl:
	@echo "============  Making PEPC Sorting library  ============="
	cd $(SLPEPCDIR) && $(MAKE) $(MFLAGS)

libpepc.%: libpthreads libsl
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

pepce: libpepc.coulomb
	@echo "============  Making Frontend PEPC-E (Benchmark version)  ============="
	cd $(FRONTENDDIR)pepc-e && $(MAKE) $(MFLAGS)

pepcs: libpepc.coulomb
	@echo "============  Making Frontend PEPC-S (ScaFaCoS-library version + minimal frontend)  ============="
	cd $(FRONTENDDIR)pepc-s && $(MAKE) $(MFLAGS)

clean: clean-doc
	cd $(SLPEPCDIR)   && $(MAKE) $(MFLAGS) clean
	cd $(PTHREADSDIR) && $(MAKE) $(MFLAGS) clean
	cd $(LPEPCDIR)    && $(MAKE) $(MFLAGS) clean
	cd $(FRONTENDDIR)pepc-s     && $(MAKE) $(MFLAGS) clean
	cd $(FRONTENDDIR)pepc-e     && $(MAKE) $(MFLAGS) clean
	cd $(FRONTENDDIR)pepc-b     && $(MAKE) $(MFLAGS) clean
	cd $(FRONTENDDIR)pepc-mw    && $(MAKE) $(MFLAGS) clean
	cd $(FRONTENDDIR)pepc-mini  && $(MAKE) $(MFLAGS) clean

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
	
