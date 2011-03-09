#
#  Makefile for pepc
#

include makefile.defs

default: pepce

benchmark: pepce

all: pepce

pepce:  pepcbasics
	@echo "============  Making Frontend PEPC-E (Benchmark version)  ============="
	cd pepc-e && $(MAKE) $(MFLAGS)

pepcbasics:
	@echo "============  Making PEPC Sorting Library  ============="
	cd sl_pepc && $(MAKE) $(MFLAGS)
	@echo "============  Making PEPC Memory Bookkeeper =============" 
	cd memwatch && $(MAKE) $(MFLAGS) 
	@echo "============  Making PEPC Pthreads Interface  ============="
	cd pthreads && $(MAKE) $(MFLAGS)
	@echo "============  Making PEPC Library  ============="
	cd lpepcsrc && $(MAKE) $(MFLAGS)

clean: clean-doc
	cd sl_pepc  && $(MAKE) $(MFLAGS) clean
	cd memwatch && $(MAKE) $(MFLAGS) clean
	cd pthreads && $(MAKE) $(MFLAGS) clean
	cd lpepcsrc && $(MAKE) $(MFLAGS) clean
	cd pepc-e   && $(MAKE) $(MFLAGS) clean

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
	