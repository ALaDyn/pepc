#
#  Makefile for pepc
#

include makefile.defs

pepce:
	cd sl_pepc && $(MAKE)
	cd lpepcsrc && $(MAKE)
	cd pepc-e && $(MAKE)

clean:
#	cd sl_pepc && $(MAKE) clean
	cd lpepcsrc && $(MAKE) clean && cd ..
	cd pepc-e && $(MAKE) clean && cd ..
