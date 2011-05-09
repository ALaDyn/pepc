#
#  Makefile for pepc
#

include makefile.defs

default: pepce

benchmark: pepce

all: pepce pepcmini pepcmw pepcs

pepcmw: pepcbasics
	@echo "============  Making Frontend PEPC-MW (Mathias Winkel version)  ============="
	cd pepc-mw && $(MAKE) $(MFLAGS)
	
pepcmini: pepcbasics
	@echo "============  Making Frontend PEPC-MINI (minial version)  ============="
	cd pepc-mini && $(MAKE) $(MFLAGS)
	
pepce:  pepcbasics
	@echo "============  Making Frontend PEPC-E (Benchmark version)  ============="
	cd pepc-e && $(MAKE) $(MFLAGS)

pepcs:  pepcbasics
	@echo "============  Making Frontend PEPC-S (ScaFaCoS-library version + minimal frontend)  ============="
	cd pepc-s && $(MAKE) $(MFLAGS)

pepcbasics:
	@echo "============  Making PEPC Sorting Library  ============="
	cd sl_pepc && $(MAKE) $(MFLAGS)
	@echo "============  Making PEPC Pthreads Interface  ============="
	cd pthreads && $(MAKE) $(MFLAGS)
	@echo "============  Making PEPC Library  ============="
	cd lpepcsrc && $(MAKE) $(MFLAGS) -j1

clean: clean-doc
	cd sl_pepc  && $(MAKE) $(MFLAGS) clean
	cd pthreads && $(MAKE) $(MFLAGS) clean
	cd lpepcsrc && $(MAKE) $(MFLAGS) clean
	cd pepc-e   && $(MAKE) $(MFLAGS) clean
	cd pepc-mw  && $(MAKE) $(MFLAGS) clean
	cd pepc-mini  && $(MAKE) $(MFLAGS) clean

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
	
