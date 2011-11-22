
##############################################################
#          General Makefile for BGP_hpccount
##############################################################

.PHONY: all clean realclean dummy_hpc

IHPCT_BASE ?= .

#############
### Targets #
#############

DEST		= ${PWD}

LIB_DIR		= lib

LIBRARY		= $(LIB_DIR)/libhpccnt.a

LIB_SO		= $(LIB_DIR)/libhpccnt.so

DUMMY_HPC	= dummy_hpc

TARGETS		= $(DUMMY_HPC) $(LIBRARY) $(LIB_SO)

# PRESET VARIABLES
#################################
### Compilers, Programs and Flags
#################################

### C Compiler
CC		= mpicc
CFLAGS		= -g -Iinclude
CSOFLAGS	= -fPIC -I$(IHPCT_BASE)/include $(CFLAGS)

### Linker/Loader
LD		= mpicc
LDFLAGS		= 
LDSOFLAGS	= -shared -fPIC -L$(IHPCT_BASE)/lib -lhpc


#################################
### Sources and Objects         #
#################################

LIB_SRC		= src/libhpccnt_c.c src/libhpccnt_f.c

LIB_OBJ		= $(LIB_SRC:%.c=%.o)

LIB_SO_OBJ	= $(LIB_SRC:%.c=%.os)

MAKEFILE	= Makefile

.SUFFIXES: .os
.PHONY: clean realclean clobber

########################
### Main rules         #
########################

all:		$(TARGETS)

$(DUMMY_HPC):;	(cd dummy_hpc; $(MAKE))

$(LIBRARY):	$(LIB_OBJ) $(MAKEFILE)
		@echo "Making $(LIBRARY) ..."
		@mkdir -p $(LIB_DIR)
		@$(AR) $(ARFLAGS) $(LIBRARY) $(LIB_OBJ)
		@echo "done"

$(LIB_SO):	$(LIB_SO_OBJ) $(MAKEFILE)
		@echo "Making $(LIB_SO) ..."
		@mkdir -p $(LIB_DIR)
		@$(LD) $(LDSOFLAGS) -o $(LIB_SO) $(LIB_SO_OBJ)
		@echo "done"

clean:;		(cd dummy_hpc; $(MAKE) clean)
		rm -f $(LIB_OBJ) $(LIB_SO_OBJ) core

realclean:      clobber

clobber:	clean
		(cd dummy_hpc; $(MAKE) realclean)
		rm -f $(LIBRARY) $(LIB_SO) core
		rm -rf $(LIB_DIR)

$(LIB_OBJ):	Makefile
$(LIB_SO_OBJ):	Makefile

# These are the general compilation rules (for unusual suffixes)

.c.o:
	$(CC) $(CFLAGS) -c -o $@ $<

.c.os:
	$(CC) $(CSOFLAGS) -c -o $@ $<
