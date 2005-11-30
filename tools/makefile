#
#  Makefile for PEPC tools
#
ALL: default
##### User configurable options #####

#LDFLAGS = -lpmpi -Xm


.SUFFIXES: .f90


## Intel compiler
# FC = ifort
#FFLAGS = -g -O3

## IBM compiler 
FC = /usr/bin/xlf90
FFLAGS = -O3 -qsuffix=f=f90

# Source files

SOURCES.f90 =	slicer.f90 


# Names of object files derived from sources

OBJECTS = $(SOURCES:.f90=.o)


#  compile mode
.f90.o: 	
	$(FC) -c $(FFLAGS) $(DBFLAGS) $<

#  targets

slicer: slicer.o
	$(FC)  $(FFLAGS) -o $@ slicer.o 
	cp slicer ../bin/.

clean:
	/bin/rm -f *.o  ../bin/slicer *% *~ rand *.o *.pif core 






