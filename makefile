#
#  Makefile for pepc
#

include makefile.defs

info:
	echo "info kram"

pepce:
	cd sl_pepc && $(MAKE)
	cd lpepcsrc && $(MAKE)
	cd pepc-e && $(MAKE)

pepcb:
	cd sl_pepc && $(MAKE)
	cd lpepcsrc && $(MAKE)
	cd pepc-b && $(MAKE)

clean:
	cd sl_pepc && $(MAKE) clean
	cd lpepcsrc && $(MAKE) clean
	cd pepc-e && $(MAKE) clean
	cd pepc-b && $(MAKE) clean
