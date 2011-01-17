#
#  Makefile for pepc
#

include makefile.defs

pepce:
	cd sl_pepc && $(MAKE)
	cd lpepcsrc && $(MAKE)
	cd pepc-e && $(MAKE)

clean: clean-doc
	cd sl_pepc && $(MAKE) clean
	cd lpepcsrc && $(MAKE) clean && cd ..
	cd pepc-e && $(MAKE) clean && cd ..

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
	