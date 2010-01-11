#
#  Makefile for pepc
#

include makefile.defs

usage:
	@echo "PEPC makefile usage"
	@echo " "
	@echo "to build pepc-b, use: make pepcb"
	@echo "to build pepc-b with sionlib, use: make pepcb SION=1, or set SION=1 in makefile.defs"
	@echo " "
	@echo "to build pepc-e, use: make pepce"
	@echo "to build pepc-e with sionlib, use: make pepce SION=1, or set SION=1 in makefile.defs"
	@echo " "
	@echo "to clean up, use: make clean"

info:
	@echo ""
	@echo "PEPC makefile info"
	@echo "TARGET=$(TARGET)"
	@echo "SION=$(SION)"
	@echo "MACH=$(MACH)"

pepce:
	cd sl_pepc && $(MAKE)
	cd lpepcsrc && $(MAKE)
	cd pepc-e && $(MAKE)
	make info TARGET=$@

pepcb:
	cd sl_pepc && $(MAKE)
	cd lpepcsrc && $(MAKE)
	cd pepc-b && $(MAKE)
	make info

clean:
	cd sl_pepc && $(MAKE) clean
	cd lpepcsrc && $(MAKE) clean
	cd pepc-e && $(MAKE) clean
	cd pepc-b && $(MAKE) clean
