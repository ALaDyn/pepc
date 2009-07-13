#
#  Makefile for pepc
#

include makefile.defs

pepce:
	cd lpepcsrc && $(MAKE)
	cd pepc-e && $(MAKE)

clean:
	cd lpepcsrc && $(MAKE) clean && cd ..
	cd pepc-e && $(MAKE) clean && cd ..
