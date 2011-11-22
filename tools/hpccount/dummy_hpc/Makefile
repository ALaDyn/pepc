#############
### Targets #
#############

DEST		= .

LIB_DIR		= ../lib

LIBRARY		= $(LIB_DIR)/libhpc.a

LIB_SO		= $(LIB_DIR)/libhpc.so

TARGET		= $(LIBRARY) $(LIB_SO)

# PRESET VARIABLES
#################################
### Compilers, Programs and Flags
#################################

### C Compiler
CC		= mpicc
CFLAGS		= -g -fPIC -I../include

### Linker/Loader
LD		= mpicc
LDFLAGS		= 
LDSOFLAGS	= -shared -fPIC $(LDFLAGS)


#################################
### Sources Headers and Objects #
#################################

HDRS		= ../include/libhpc.h

EXTHDRS		=

SRCS		= libhpc.c

SYSHDRS		=

OBJS		= libhpc.o

LDL		=

LIBEXP_LIST	=

MAKEFILE	= ./Makefile

#######################
### Program Libraries #
#######################

LIBS		=

########################
### Additional Targets #
########################

all:		$(TARGET)

$(LIBRARY):	$(OBJS) $(MAKEFILE)
		@echo "Making $(LIBRARY) ..."
		mkdir -p $(LIB_DIR)
		@$(AR) $(ARFLAGS) $(LIBRARY) $(OBJS)
		@echo "done"

$(LIB_SO):	$(OBJS) $(MAKEFILE)
		@echo "Making $(LIB_SO) ..."
		@mkdir -p $(LIB_DIR)
		@$(LD) $(LDSOFLAGS) $(LIBEXP_LIST) -o $(LIB_SO) \
			$(OBJS) $(LDL)
		@echo "done"

clean:;		rm -f $(OBJS) core

realclean:      clobber

clobber:;	rm -rf $(OBJS) $(TARGET) core tags

$(OBJS): Makefile
