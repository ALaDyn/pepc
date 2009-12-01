# This is a Makefile!

#  
#  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
#  
#  file: Makefile
#  timestamp: 2009-11-13 12:31:36 +0100
#  


# only on top level
ifeq ($(mod),)

# initial config
MAKEFILE_IN:=

SL_USE_MPI:=
ALL_AS_MPI:=
DI_VERSIONS:=
FAR:=perl
TUNEABLE:=#TUNEABLE#
TARGET:=

WRAPPER_PREFIX:=
WRAPPER_MPI_PREFIX:=mpi_

INCDIR:=./
LIBDIR:=

ALL:=libsl


or2 = $(if $(1),$(1),$(2))
or3 = $(if $(1),$(1),$(if $(2),$(2),$(3)))
or3 = $(if $(1),$(1),$(if $(2),$(2),$(if $(3),$(3),$(4))))

and2 = $(if $(1),$(2),)
and3 = $(if $(1),$(if $(2),$(3),),)
and4 = $(if $(1),$(if $(2),$(if $(3),$(4),),),)


mod_list:=

MODULES:=$(call or2,$(MODULES),$(dir $(wildcard sl_*/.)))
#$(info modules: $(MODULES))

MAKEFILE_IN:=$(call or2,$(MAKEFILE_IN),Makefile.in)

-include $(MAKEFILE_IN)

RM:=$(call or2,$(RM),rm -f)
ECHO:=$(call or2,$(ECHO),echo)


ifneq ($(DEBUG),)
 CPPFLAGS+=-DDEBUG=$(DEBUG)
endif


LIBSL_SEQ:=$(call or2,$(LIBSL_SEQ),$(LIBDIR)libsl_seq.a)
LIBSL_MPI:=$(call or2,$(LIBSL_MPI),$(LIBDIR)libsl_mpi.a)
LIBSL:=$(call or2,$(LIBSL),$(LIBDIR)libsl.a)

TUNE:=sltune
DEVL:=sldevel
DEMO:=sldemo

PREREQ:=Makefile $(MAKEFILE_IN)

M_SL_USE_MPI:=

M_INFO:=

M_CORE_SEQ:=
M_CORE_MPI:=
M_CORE:=

M_LIBSL_SEQ:=
M_LIBSL_MPI:=
M_LIBSL:=

M_TUNE_SEQ:=
M_TUNE_MPI:=
M_TUNE:=

M_DEVL_SEQ:=
M_DEVL_MPI:=
M_DEVL:=

M_DEMO_SEQ:=
M_DEMO_MPI:=
M_DEMO:=

M_CLEAN:=

M_INCDIR:=
M_CODE_INCLUDE:=

M_OBJS_SL_CORE_SEQ:=
M_OBJS_SL_CORE_MPI:=
M_OBJS_SL_CORE:=

M_OBJS_SL_WRAP_SEQ:=
M_OBJS_SL_WRAP_MPI:=
M_OBJS_SL_WRAP:=

all: $(call or2,$(ALL),demo)

else

MODULES:=#MODULES#

endif


# push mod on the list
mod_list:=$(mod) $(mod_list)

# include modules
$(foreach mod,$(MODULES),$(eval include Makefile))

# pop mod from the list
mod:=$(firstword $(mod_list))
mod_list:=$(wordlist 2,$(words $(mod_list)),$(mod_list))


# current module
cmod:=$(patsubst %/,%,$(mod))
#cmod_:=$(subst /,_,$(subst .,_,$(cmod)))
cmod_:=$(cmod)
cmod_dir:=$(call and2,$(cmod),$(cmod)/)

#$(warning mod: $(mod))
#$(warning cmod: $(cmod))
#$(warning cmod_: $(cmod_))
#$(warning cmod_dir: $(cmod_dir))


# base directories
$(cmod_)_SRCDIR:=$(cmod_dir)src/
$(cmod_)_TMPDIR:=$(cmod_dir)tmp/
#$(cmod_)_LIBDIR:=$(cmod_dir)
$(cmod_)_INCDIR:=$(cmod_dir)include/

M_INCDIR+=$($(cmod_)_INCDIR)

# include files
$(cmod_)_SRC_INCLUDE:=$($(cmod_)_SRCDIR)include/

$(cmod_)_CODE_SL_INCLUDE:=$(wildcard $($(cmod_)_SRC_INCLUDE)*.h)

$(cmod_)_CODE_INCLUDE:=$(wildcard $($(cmod_)_INCDIR)*.h $(call and2,$(INCDIR),$(INCDIR)*.h))

M_CODE_INCLUDE+=$($(cmod_)_CODE_INCLUDE)


# CC, CPP, LD and flags
$(cmod_)_CC:=$(call or2,$($(cmod_)_CC),$(CC))
$(cmod_)_CFLAGS:=$(CFLAGS) $($(cmod_)_CFLAGS)
$(cmod_)_CFLAGS_SL:=$(CFLAGS_SL) $($(cmod_)_CFLAGS_SL)
$(cmod_)_CPPFLAGS:=$(CPPFLAGS) $($(cmod_)_CPPFLAGS) $(addprefix -I,$(INCDIR))
$(cmod_)_CPPFLAGS_SL:=$(CPPFLAGS_SL) $($(cmod_)_CPPFLAGS_SL) -I$($(cmod_)_SRC_INCLUDE)

$(cmod_)_MPICC:=$(call or2,$($(cmod_)_MPICC),$(MPICC))
$(cmod_)_MPICFLAGS:=$(MPICFLAGS) $($(cmod_)_MPICFLAGS)
$(cmod_)_MPICFLAGS_SL:=$(MPICFLAGS_SL) $($(cmod_)_MPICFLAGS_SL)
$(cmod_)_MPICPPFLAGS:=$(MPICPPFLAGS) $($(cmod_)_MPICPPFLAGS)
$(cmod_)_MPICPPFLAGS_SL:=$(MPICPPFLAGS_SL) $($(cmod_)_MPICPPFLAGS_SL)

$(cmod_)_LD:=$(call or2,$($(cmod_)_LD),$(LD))
$(cmod_)_LDFLAGS:=$(LDFLAGS) $($(cmod_)_LDFLAGS)

$(cmod_)_MPILD:=$(call or2,$($(cmod_)_MPILD),$(MPILD))
$(cmod_)_MPILDFLAGS:=$(MPILDFLAGS) $($(cmod_)_MPILDFLAGS)


# options
$(cmod_)_SL_USE_MPI:=$(call or2,$($(cmod_)_SL_USE_MPI),$(SL_USE_MPI))

# SL_USE_MPI = 0 | 1 | detect (default)
ifneq ($($(cmod_)_SL_USE_MPI),0) # not 0
 ifneq ($($(cmod_)_SL_USE_MPI),1) # not 1
  # else: use configuration file
  ifneq ($(wildcard $($(cmod_)_SRC_INCLUDE)sl_config.h),)
   ifneq ($(shell $(PRINTF) "\#include \"sl_config.h\"\nSL_USE_MPI\n" | $(CPP) $($(cmod_)_CPPFLAGS_SL) 2>/dev/null | $(GREP) SL_USE_MPI),SL_USE_MPI)
    $(cmod_)_SL_USE_MPI:=1
   else
    $(cmod_)_SL_USE_MPI:=0
   endif
  endif
 endif
endif

ifeq ($($(cmod_)_SL_USE_MPI),1)
  M_SL_USE_MPI:=1
endif

$(cmod_)_ALL_AS_MPI:=$(call or2,$($(cmod_)_ALL_AS_MPI),$(ALL_AS_MPI))

# ALL_AS_MPI = 0 | 1 | SL_USE_MPI (default)
ifneq ($($(cmod_)_ALL_AS_MPI),0) # not 0
 ifneq ($($(cmod_)_ALL_AS_MPI),1) # not 1
  # else: use SL_USE_MPI or not
#  $(cmod_)_ALL_AS_MPI:=$($(cmod_)_SL_USE_MPI)
  $(cmod_)_ALL_AS_MPI:=
 endif
endif

$(cmod_)_DI_VERSION:=$(call or3,$($(cmod_)_DI_VERSION),$(DI_VERSIONS),false)

$(cmod_)_TUNEABLE:=$(call or2,$($(cmod_)_TUNEABLE),$(TUNEABLE))

$(cmod_)_TARGET:=$(call or3,$($(cmod_)_TARGET),$(TARGET),NOT)

$(cmod_)_PREFIX:=$(call and2,$(patsubst sl_%,%,$(cmod)),$(patsubst sl_%,%,$(cmod))_)


# option dependent flags
$(cmod_)_CFLAGS_SL+=$(call and2,$($(cmod_)_TUNEABLE),-DSL_TUNEABLE) -D$($(cmod_)_TARGET) -DPREFIX=$($(cmod_)_PREFIX)


# target files
$(cmod_)_LIBSL_SEQ:=$(call or2,$($(cmod_)_LIBSL_SEQ),$(LIBSL_SEQ))
$(cmod_)_LIBSL_MPI:=$(call or2,$($(cmod_)_LIBSL_MPI),$(LIBSL_MPI))
$(cmod_)_LIBSL:=$(call or2,$($(cmod_)_LIBSL),$(LIBSL))

$(cmod_)_TUNE:=$(call or2,$($(cmod)_TUNE),$(if $(cmod),$(cmod)_$(TUNE),$(TUNE)))

$(cmod_)_DEVL:=$(call or2,$($(cmod)_DEVL),$(if $(cmod),$(cmod)_$(DEVL),$(DEVL)))

$(cmod_)_DEMO:=$(call or2,$($(cmod)_DEMO),$(if $(cmod),$(cmod)_$(DEMO),$(DEMO)))


# core files
$(cmod_)_SRC_CORE_SEQ:=$($(cmod_)_SRCDIR)core/
$(cmod_)_SRC_CORE_MPI:=$($(cmod_)_SRCDIR)core_mpi/
$(cmod_)_TMP_CORE_SEQ:=$($(cmod_)_TMPDIR)core/
$(cmod_)_TMP_CORE_MPI:=$($(cmod_)_TMPDIR)core_mpi/

$(cmod_)_CODE_SL_CORE_SEQ_INCLUDE:=$(wildcard $($(cmod_)_SRC_CORE_SEQ)*.h)
$(cmod_)_CODE_SL_CORE_MPI_INCLUDE:=$(wildcard $($(cmod_)_SRC_CORE_MPI)*.h)

# filter-out existing _di-versions (to get a clean list of core-files)
$(cmod_)_CODE_SL_CORE_SEQ:=$(filter-out $(wildcard $($(cmod_)_SRC_CORE_SEQ)*_di.c),$(wildcard $($(cmod_)_SRC_CORE_SEQ)*.c))
$(cmod_)_CODE_SL_CORE_MPI:=$(filter-out $(wildcard $($(cmod_)_SRC_CORE_MPI)*_di.c),$(wildcard $($(cmod_)_SRC_CORE_MPI)*.c))
# add _di-versions if not disabled
ifneq ($($(cmod_)_DI_VERSION),false)
 $(cmod_)_CODE_SL_CORE_SEQ+=$(addsuffix _di.c,$(basename $($(cmod_)_CODE_SL_CORE_SEQ)))
 $(cmod_)_CODE_SL_CORE_MPI+=$(addsuffix _di.c,$(basename $($(cmod_)_CODE_SL_CORE_MPI)))
endif

$(cmod_)_OBJS_SL_CORE_SEQ:=$(addsuffix .o, $(basename $(notdir $($(cmod_)_CODE_SL_CORE_SEQ))))
$(cmod_)_OBJS_SL_CORE_MPI:=$(addsuffix .o, $(basename $(notdir $($(cmod_)_CODE_SL_CORE_MPI))))

$(cmod_)_TMP_OBJS_SL_CORE_SEQ:=$(addprefix $($(cmod_)_TMP_CORE_SEQ),$($(cmod_)_OBJS_SL_CORE_SEQ))
$(cmod_)_TMP_OBJS_SL_CORE_MPI:=$(addprefix $($(cmod_)_TMP_CORE_MPI),$($(cmod_)_OBJS_SL_CORE_MPI))


# tune files
$(cmod_)_SRC_TUNE_SEQ:=$($(cmod_)_SRCDIR)tune/
$(cmod_)_SRC_TUNE_MPI:=$($(cmod_)_SRCDIR)tune_mpi/
$(cmod_)_TMP_TUNE_SEQ:=$($(cmod_)_TMPDIR)tune/
$(cmod_)_TMP_TUNE_MPI:=$($(cmod_)_TMPDIR)tune_mpi/

$(cmod_)_CODE_SL_TUNE_SEQ_INCLUDE:=$(wildcard $($(cmod_)_SRC_TUNE_SEQ)*.h)
$(cmod_)_CODE_SL_TUNE_MPI_INCLUDE:=$(wildcard $($(cmod_)_SRC_TUNE_MPI)*.h)

$(cmod_)_CODE_SL_TUNE_SEQ:=$(wildcard $($(cmod_)_SRC_TUNE_SEQ)*.c)
$(cmod_)_CODE_SL_TUNE_MPI:=$(wildcard $($(cmod_)_SRC_TUNE_MPI)*.c)

$(cmod_)_OBJS_SL_TUNE_SEQ:=$(addsuffix .o, $(basename $(notdir $($(cmod_)_CODE_SL_TUNE_SEQ))))
$(cmod_)_OBJS_SL_TUNE_MPI:=$(addsuffix .o, $(basename $(notdir $($(cmod_)_CODE_SL_TUNE_MPI))))

$(cmod_)_TMP_OBJS_SL_TUNE_SEQ:=$(addprefix $($(cmod_)_TMP_TUNE_SEQ),$($(cmod_)_OBJS_SL_TUNE_SEQ))
$(cmod_)_TMP_OBJS_SL_TUNE_MPI:=$(addprefix $($(cmod_)_TMP_TUNE_MPI),$($(cmod_)_OBJS_SL_TUNE_MPI))


# devel files
$(cmod_)_SRC_DEVL_SEQ:=$($(cmod_)_SRCDIR)devel/
$(cmod_)_SRC_DEVL_MPI:=$($(cmod_)_SRCDIR)devel_mpi/
$(cmod_)_TMP_DEVL_SEQ:=$($(cmod_)_TMPDIR)devel/
$(cmod_)_TMP_DEVL_MPI:=$($(cmod_)_TMPDIR)devel_mpi/

$(cmod_)_CODE_SL_DEVL_SEQ_INCLUDE:=$(wildcard $($(cmod_)_SRC_DEVL_SEQ)*.h)
$(cmod_)_CODE_SL_DEVL_MPI_INCLUDE:=$(wildcard $($(cmod_)_SRC_DEVL_MPI)*.h)

$(cmod_)_CODE_SL_DEVL_SEQ:=$(wildcard $($(cmod_)_SRC_DEVL_SEQ)*.c)
$(cmod_)_CODE_SL_DEVL_MPI:=$(wildcard $($(cmod_)_SRC_DEVL_MPI)*.c)

$(cmod_)_OBJS_SL_DEVL_SEQ:=$(addsuffix .o, $(basename $(notdir $($(cmod_)_CODE_SL_DEVL_SEQ))))
$(cmod_)_OBJS_SL_DEVL_MPI:=$(addsuffix .o, $(basename $(notdir $($(cmod_)_CODE_SL_DEVL_MPI))))

$(cmod_)_TMP_OBJS_SL_DEVL_SEQ:=$(addprefix $($(cmod_)_TMP_DEVL_SEQ),$($(cmod_)_OBJS_SL_DEVL_SEQ))
$(cmod_)_TMP_OBJS_SL_DEVL_MPI:=$(addprefix $($(cmod_)_TMP_DEVL_MPI),$($(cmod_)_OBJS_SL_DEVL_MPI))


# demo files
$(cmod_)_SRC_DEMO_SEQ:=$($(cmod_)_SRCDIR)demo/
$(cmod_)_SRC_DEMO_MPI:=$($(cmod_)_SRCDIR)demo_mpi/
$(cmod_)_TMP_DEMO_SEQ:=$($(cmod_)_TMPDIR)demo/
$(cmod_)_TMP_DEMO_MPI:=$($(cmod_)_TMPDIR)demo_mpi/

$(cmod_)_CODE_SL_DEMO_SEQ_INCLUDE:=$(wildcard $($(cmod_)_SRC_DEMO_SEQ)*.h)
$(cmod_)_CODE_SL_DEMO_MPI_INCLUDE:=$(wildcard $($(cmod_)_SRC_DEMO_MPI)*.h)

$(cmod_)_CODE_SL_DEMO_SEQ:=$(wildcard $($(cmod_)_SRC_DEMO_SEQ)*.c)
$(cmod_)_CODE_SL_DEMO_MPI:=$(wildcard $($(cmod_)_SRC_DEMO_MPI)*.c)

$(cmod_)_OBJS_SL_DEMO_SEQ:=$(addsuffix .o, $(basename $(notdir $($(cmod_)_CODE_SL_DEMO_SEQ))))
$(cmod_)_OBJS_SL_DEMO_MPI:=$(addsuffix .o, $(basename $(notdir $($(cmod_)_CODE_SL_DEMO_MPI))))

$(cmod_)_TMP_OBJS_SL_DEMO_SEQ:=$(addprefix $($(cmod_)_TMP_DEMO_SEQ),$($(cmod_)_OBJS_SL_DEMO_SEQ))
$(cmod_)_TMP_OBJS_SL_DEMO_MPI:=$(addprefix $($(cmod_)_TMP_DEMO_MPI),$($(cmod_)_OBJS_SL_DEMO_MPI))


# wrapper files
$(cmod_)_SRC_WRAP_SEQ:=$(cmod_dir)
$(cmod_)_SRC_WRAP_MPI:=$(cmod_dir)
$(cmod_)_TMP_WRAP_SEQ:=$(cmod_dir)
$(cmod_)_TMP_WRAP_MPI:=$(cmod_dir)

$(cmod_)_CODE_SL_WRAP_MPI_INCLUDE:=$(wildcard $($(cmod_)_SRC_WRAP_MPI)$(WRAPPER_MPI_PREFIX)*.h)
$(cmod_)_CODE_SL_WRAP_SEQ_INCLUDE:=$(filter-out $($(cmod_)_CODE_SL_WRAP_MPI_INCLUDE),$(wildcard $($(cmod_)_SRC_WRAP_SEQ)*.h))

$(cmod_)_CODE_SL_WRAP_MPI:=$(wildcard $($(cmod_)_SRC_WRAP_MPI)$(WRAPPER_MPI_PREFIX)*.c)
$(cmod_)_CODE_SL_WRAP_SEQ:=$(filter-out $($(cmod_)_CODE_SL_WRAP_MPI),$(wildcard $($(cmod_)_SRC_WRAP_SEQ)*.c))

$(cmod_)_OBJS_SL_WRAP_SEQ:=$(addsuffix .o, $(basename $(notdir $($(cmod_)_CODE_SL_WRAP_SEQ))))
$(cmod_)_OBJS_SL_WRAP_MPI:=$(addsuffix .o, $(basename $(notdir $($(cmod_)_CODE_SL_WRAP_MPI))))

$(cmod_)_TMP_OBJS_SL_WRAP_SEQ:=$(addprefix $($(cmod_)_TMP_WRAP_SEQ),$($(cmod_)_OBJS_SL_WRAP_SEQ))
$(cmod_)_TMP_OBJS_SL_WRAP_MPI:=$(addprefix $($(cmod_)_TMP_WRAP_MPI),$($(cmod_)_OBJS_SL_WRAP_MPI))


# find-and-replace scripts
$(cmod_)_NAMES_OXL_ALL_FUNCTIONS:=$($(cmod_)_SRCDIR)names_oxl_all_functions
$(cmod_)_FAR_SCRIPT_DI:=$($(cmod_)_TMPDIR)far_script_di
$(cmod_)_FIND_AND_REPLACE:=./$(cmod_dir)scripts/find_and_replace.sh "$(FAR)"

# list of required object-files
$(cmod_)_OBJS_SL_CORE:=$($(cmod_)_OBJS_SL_CORE_SEQ)
$(cmod_)_OBJS_SL_TUNE:=$($(cmod_)_OBJS_SL_TUNE_SEQ)
$(cmod_)_OBJS_SL_DEVL:=$($(cmod_)_OBJS_SL_DEVL_SEQ)
$(cmod_)_OBJS_SL_DEMO:=$($(cmod_)_OBJS_SL_DEMP_SEQ)
$(cmod_)_OBJS_SL_WRAP:=$($(cmod_)_OBJS_SL_WRAP_SEQ)

$(cmod_)_TMP_OBJS_SL_CORE:=$($(cmod_)_TMP_OBJS_SL_CORE_SEQ)
$(cmod_)_TMP_OBJS_SL_TUNE:=$($(cmod_)_TMP_OBJS_SL_TUNE_SEQ)
$(cmod_)_TMP_OBJS_SL_DEVL:=$($(cmod_)_TMP_OBJS_SL_DEVL_SEQ)
$(cmod_)_TMP_OBJS_SL_DEMO:=$($(cmod_)_TMP_OBJS_SL_DEMO_SEQ)
$(cmod_)_TMP_OBJS_SL_WRAP:=$($(cmod_)_TMP_OBJS_SL_WRAP_SEQ)

ifeq ($($(cmod_)_SL_USE_MPI),1)
 $(cmod_)_OBJS_SL_CORE+=$($(cmod_)_OBJS_SL_CORE_MPI)
 $(cmod_)_OBJS_SL_TUNE+=$($(cmod_)_OBJS_SL_TUNE_MPI)
 $(cmod_)_OBJS_SL_DEVL:=$($(cmod_)_OBJS_SL_DEVL_MPI)
 $(cmod_)_OBJS_SL_DEMO:=$($(cmod_)_OBJS_SL_DEMO_MPI)
 $(cmod_)_OBJS_SL_WRAP+=$($(cmod_)_OBJS_SL_WRAP_MPI)

 $(cmod_)_TMP_OBJS_SL_CORE+=$($(cmod_)_TMP_OBJS_SL_CORE_MPI)
 $(cmod_)_TMP_OBJS_SL_TUNE+=$($(cmod_)_TMP_OBJS_SL_TUNE_MPI)
 $(cmod_)_TMP_OBJS_SL_DEVL:=$($(cmod_)_TMP_OBJS_SL_DEVL_MPI)
 $(cmod_)_TMP_OBJS_SL_DEMO:=$($(cmod_)_TMP_OBJS_SL_DEMO_MPI)
 $(cmod_)_TMP_OBJS_SL_WRAP+=$($(cmod_)_TMP_OBJS_SL_WRAP_MPI)
else
 ifeq ($(M_SL_USE_MPI),1)
  $(cmod_)_OBJS_SL_WRAP+=$($(cmod_)_OBJS_SL_WRAP_MPI)

  $(cmod_)_TMP_OBJS_SL_WRAP+=$($(cmod_)_TMP_OBJS_SL_WRAP_MPI)
 endif
endif


# final CC, LD and flags
$(cmod_)_CFLAGS_:=$($(cmod_)_CFLAGS) $($(cmod_)_CPPFLAGS)
$(cmod_)_CFLAGS_SL_:=$($(cmod_)_CFLAGS_SL) -DSL_USE_MPI_IGNORE $($(cmod_)_CPPFLAGS_SL)
$(cmod_)_MPICFLAGS_:=$($(cmod_)_CFLAGS) $($(cmod_)_MPICFLAGS) $($(cmod_)_CPPFLAGS) $($(cmod_)_MPICPPFLAGS)
$(cmod_)_MPICFLAGS_SL_:=$($(cmod_)_CFLAGS_SL) $($(cmod_)_MPICFLAGS_SL) -DSL_USE_MPI_FORCE $($(cmod_)_CPPFLAGS_SL) $($(cmod_)_MPICPPFLAGS_SL)

ifeq ($($(cmod_)_ALL_AS_MPI),1)
 $(cmod_)_CC:=$($(cmod_)_MPICC)
 $(cmod_)_CFLAGS_:=$($(cmod_)_MPICFLAGS_)
 $(cmod_)_CFLAGS_SL_:=$($(cmod_)_MPICFLAGS_SL_)
endif

$(cmod_)_LDFLAGS_:=$($(cmod_)_LDFLAGS)
$(cmod_)_MPILDFLAGS_:=$($(cmod_)_LDFLAGS) $($(cmod_)_MPILDFLAGS)

ifeq ($($(cmod_)_SL_USE_MPI),1)
 $(cmod_)_LD:=$($(cmod_)_MPILD)
 $(cmod_)_LDFLAGS_:=$($(cmod_)_MPILDFLAGS_)
endif


# prerequisites
$(cmod_)_PREREQ:=$(PREREQ) $($(cmod_)_PREREQ)


# commands
define $(cmod_)_comp
@echo "compiling $(TYPE) '$<' to '$@'"
-@$(MKDIR_P) $(dir $@)
$(CC) -c $(CFLAGS) -o $@ $<
endef

#@$(ECHO) "adding $(TYPE) '$<' to '$@'"
define $(cmod_)_add
-@$(MKDIR_P) $(dir $@)
$(AR) -r $@ $<
endef

define $(cmod_)_make
@$(ECHO) "making $(TYPE) '$@'"
-@$(MKDIR_P) $(dir $@)
$(LD) $(LDFLAGS) -o $@ $^
endef


# info
$(cmod_)_info: cmod:=$(cmod)
$(cmod_)_info: cmod_:=$(cmod_)
$(cmod_)_info:
ifneq ($(cmod),)
	@$(ECHO) "$(cmod)"
endif
	@$(ECHO) " SL_USE_MPI: $($(cmod_)_SL_USE_MPI)"
	@$(ECHO) " ALL_AS_MPI: $($(cmod_)_ALL_AS_MPI)"
	@$(ECHO) " PREFIX: $($(cmod_)_PREFIX)"

	@$(ECHO) " CC: $($(cmod_)_CC)"
	@$(ECHO) " CFLAGS: $($(cmod_)_CFLAGS)"
	@$(ECHO) " CPPFLAGS: $($(cmod_)_CPPFLAGS)"
	@$(ECHO) " CFLAGS_SL: $($(cmod_)_CFLAGS_SL)"
	@$(ECHO) " CPPFLAGS_SL: $($(cmod_)_CPPFLAGS_SL)"
	@$(ECHO) " CFLAGS_: $($(cmod_)_CFLAGS_)"
	@$(ECHO) " CFLAGS_SL_: $($(cmod_)_CFLAGS_SL_)"
	@$(ECHO) " MPICC: $($(cmod_)_MPICC)"
	@$(ECHO) " MPICFLAGS: $($(cmod_)_MPICFLAGS)"
	@$(ECHO) " MPICPPFLAGS: $($(cmod_)_MPICPPFLAGS)"
	@$(ECHO) " MPICFLAGS_SL: $($(cmod_)_MPICFLAGS_SL)"
	@$(ECHO) " MPICPPFLAGS_SL: $($(cmod_)_MPICPPFLAGS_SL)"
	@$(ECHO) " MPICFLAGS_: $($(cmod_)_MPICFLAGS_)"
	@$(ECHO) " MPICFLAGS_SL_: $($(cmod_)_MPICFLAGS_SL_)"
	@$(ECHO) " LD: $($(cmod_)_LD)"
	@$(ECHO) " LDFLAGS_: $($(cmod_)_LDFLAGS_)"
	@$(ECHO) " MPILD: $($(cmod_)_MPILD)"
	@$(ECHO) " MPILDFLAGS_: $($(cmod_)_MPILDFLAGS_)"
	@$(ECHO) " PREREQ: $($(cmod_)_PREREQ)"

	@$(ECHO) " src-dir: $($(cmod_)_SRCDIR)"
	@$(ECHO) " src-core-dir: $($(cmod_)_SRC_CORE_SEQ)"
	@$(ECHO) " src-core-mpi-dir: $($(cmod_)_SRC_CORE_MPI)"
	@$(ECHO) " src-tune-dir: $($(cmod_)_SRC_TUNE_SEQ)"
	@$(ECHO) " src-tune-mpi-dir: $($(cmod_)_SRC_TUNE_MPI)"
	@$(ECHO) " src-devel-dir: $($(cmod_)_SRC_DEVL_SEQ)"
	@$(ECHO) " src-devel-mpi-dir: $($(cmod_)_SRC_DEVL_MPI)"
	@$(ECHO) " src-demo-dir: $($(cmod_)_SRC_DEMO_SEQ)"
	@$(ECHO) " src-demo-mpi-dir: $($(cmod_)_SRC_DEMO_MPI)"

	@$(ECHO) " tmp-dir: $($(cmod_)_TMPDIR)"
	@$(ECHO) " tmp-core-dir: $($(cmod_)_TMP_CORE_SEQ)"
	@$(ECHO) " tmp-core-mpi-dir: $($(cmod_)_TMP_CORE_MPI)"
	@$(ECHO) " tmp-tune-dir: $($(cmod_)_TMP_TUNE_SEQ)"
	@$(ECHO) " tmp-tune-mpi-dir: $($(cmod_)_TMP_TUNE_MPI)"
	@$(ECHO) " tmp-devel-dir: $($(cmod_)_TMP_DEVL_SEQ)"
	@$(ECHO) " tmp-devel-mpi-dir: $($(cmod_)_TMP_DEVL_MPI)"
	@$(ECHO) " tmp-demo-dir: $($(cmod_)_TMP_DEMO_SEQ)"
	@$(ECHO) " tmp-demo-mpi-dir: $($(cmod_)_TMP_DEMO_MPI)"

	@$(ECHO) " code-include: $($(cmod_)_CODE_INCLUDE)"
	@$(ECHO) " code-sl-include: $($(cmod_)_CODE_SL_INCLUDE)"

	@$(ECHO) " code-core-seq: $($(cmod_)_CODE_SL_CORE_SEQ)"
	@$(ECHO) " code-core-mpi: $($(cmod_)_CODE_SL_CORE_MPI)"
	@$(ECHO) " objs-core-seq: $($(cmod_)_OBJS_SL_CORE_SEQ)"
	@$(ECHO) " objs-core-mpi: $($(cmod_)_OBJS_SL_CORE_MPI)"
	@$(ECHO) " tmp-objs-core-seq: $($(cmod_)_TMP_OBJS_SL_CORE_SEQ)"
	@$(ECHO) " tmp-objs-core-mpi: $($(cmod_)_TMP_OBJS_SL_CORE_MPI)"

	@$(ECHO) " code-tune-seq: $($(cmod_)_CODE_SL_TUNE_SEQ)"
	@$(ECHO) " code-tune-mpi: $($(cmod_)_CODE_SL_TUNE_MPI)"
	@$(ECHO) " objs-tune-seq: $($(cmod_)_OBJS_SL_TUNE_SEQ)"
	@$(ECHO) " objs-tune-mpi: $($(cmod_)_OBJS_SL_TUNE_MPI)"
	@$(ECHO) " tmp-objs-tune-seq: $($(cmod_)_TMP_OBJS_SL_TUNE_SEQ)"
	@$(ECHO) " tmp-objs-tune-mpi: $($(cmod_)_TMP_OBJS_SL_TUNE_MPI)"

	@$(ECHO) " code-devel-seq: $($(cmod_)_CODE_SL_DEVL_SEQ)"
	@$(ECHO) " code-devel-mpi: $($(cmod_)_CODE_SL_DEVL_MPI)"
	@$(ECHO) " objs-devel-seq: $($(cmod_)_OBJS_SL_DEVL_SEQ)"
	@$(ECHO) " objs-devel-mpi: $($(cmod_)_OBJS_SL_DEVL_MPI)"
	@$(ECHO) " tmp-objs-devel-seq: $($(cmod_)_TMP_OBJS_SL_DEVL_SEQ)"
	@$(ECHO) " tmp-objs-devel-mpi: $($(cmod_)_TMP_OBJS_SL_DEVL_MPI)"

	@$(ECHO) " code-demo-seq: $($(cmod_)_CODE_SL_DEMO_SEQ)"
	@$(ECHO) " code-demo-mpi: $($(cmod_)_CODE_SL_DEMO_MPI)"
	@$(ECHO) " objs-demo-seq: $($(cmod_)_OBJS_SL_DEMO_SEQ)"
	@$(ECHO) " objs-demo-mpi: $($(cmod_)_OBJS_SL_DEMO_MPI)"
	@$(ECHO) " tmp-objs-demo-seq: $($(cmod_)_TMP_OBJS_SL_DEMO_SEQ)"
	@$(ECHO) " tmp-objs-demo-mpi: $($(cmod_)_TMP_OBJS_SL_DEMO_MPI)"

	@$(ECHO) " code-wrapper-seq: $($(cmod_)_CODE_SL_WRAP_SEQ)"
	@$(ECHO) " code-wrapper-mpi: $($(cmod_)_CODE_SL_WRAP_MPI)"
	@$(ECHO) " objs-wrapper-seq: $($(cmod_)_OBJS_SL_WRAP_SEQ)"
	@$(ECHO) " objs-wrapper-mpi: $($(cmod_)_OBJS_SL_WRAP_MPI)"
	@$(ECHO) " tmp-objs-wrapper-seq: $($(cmod_)_TMP_OBJS_SL_WRAP_SEQ)"
	@$(ECHO) " tmp-objs-wrapper-mpi: $($(cmod_)_TMP_OBJS_SL_WRAP_MPI)"

	@$(ECHO) " libsl_seq: $($(cmod_)_LIBSL_SEQ)"
	@$(ECHO) " libsl_seq(*): $($(cmod_)_LIBSL_SEQ)($($(cmod_)_OBJS_SL_CORE_SEQ))"
	@$(ECHO) " libsl_mpi: $($(cmod_)_LIBSL_MPI)"
	@$(ECHO) " libsl_mpi(*): $($(cmod_)_LIBSL_MPI)($($(cmod_)_OBJS_SL_CORE_MPI))"
	@$(ECHO) " libsl: $($(cmod_)_LIBSL)"
	@$(ECHO) " libsl(*): $($(cmod_)_LIBSL)($($(cmod_)_OBJS_SL_CORE))"
	@$(ECHO) " tune: $($(cmod_)_TUNE)"
	@$(ECHO) " devel: $($(cmod_)_DEVL)"
	@$(ECHO) " demo: $($(cmod_)_DEMO)"

	@$(ECHO) " NAMES_OXL_ALL_FUNCTIONS: $($(cmod_)_NAMES_OXL_ALL_FUNCTIONS)"
	@$(ECHO) " FAR_SCRIPT_DI: $($(cmod_)_FAR_SCRIPT_DI)"
	@$(ECHO) " FIND_AND_REPLACE: $($(cmod_)_FIND_AND_REPLACE)"

M_INFO+=$(cmod_)_info


# core
$(cmod_)_core_seq: $($(cmod_)_TMP_OBJS_SL_CORE_SEQ)
$(cmod_)_core_mpi: $($(cmod_)_TMP_OBJS_SL_CORE_MPI)
$(cmod_)_core: $($(cmod_)_TMP_OBJS_SL_CORE)

M_CORE_SEQ+=$(cmod_)_core_seq
M_CORE_MPI+=$(cmod_)_core_mpi
M_CORE+=$(cmod_)_core

# core files
$($(cmod_)_TMP_OBJS_SL_CORE_SEQ): CC:=$($(cmod_)_CC)
$($(cmod_)_TMP_OBJS_SL_CORE_SEQ): CFLAGS:=$(strip $($(cmod_)_CFLAGS_) $($(cmod_)_CFLAGS_SL_))
$($(cmod_)_TMP_OBJS_SL_CORE_SEQ): TYPE:=core

$($(cmod_)_TMP_OBJS_SL_CORE_MPI): CC:=$($(cmod_)_MPICC)
$($(cmod_)_TMP_OBJS_SL_CORE_MPI): CFLAGS:=$(strip $($(cmod_)_MPICFLAGS_) $($(cmod_)_MPICFLAGS_SL_))
$($(cmod_)_TMP_OBJS_SL_CORE_MPI): TYPE:=core-mpi

$($(cmod_)_TMP_CORE_SEQ)%.o: $($(cmod_)_SRC_CORE_SEQ)%.c $($(cmod_)_CODE_SL_INCLUDE) $($(cmod_)_PREREQ)
	$($(cmod_)_comp)

$($(cmod_)_TMP_CORE_MPI)%.o: $($(cmod_)_SRC_CORE_MPI)%.c $($(cmod_)_CODE_SL_INCLUDE) $($(cmod_)_PREREQ)
	$($(cmod_)_comp)

ifneq ($($(cmod_)_DI_VERSIONS),config)
 ifneq ($($(cmod_)_DI_VERSIONS),false)
$($(cmod_)_FAR_SCRIPT_DI): cmod:=$(cmod)
$($(cmod_)_FAR_SCRIPT_DI): cmod_:=$(cmod_)
$($(cmod_)_FAR_SCRIPT_DI): $($(cmod_)_NAMES_OXL_ALL_FUNCTIONS) $($(cmod_)_PREREQ)
	-@$(MKDIR_P) $(dir $@)
	@$($(cmod_)_FIND_AND_REPLACE) begin > $@
	@$($(cmod_)_FIND_AND_REPLACE) mod_list 1 $< $($(cmod)_PREFIX): $($(cmod)_PREFIX):_di >> $@
	@$($(cmod_)_FIND_AND_REPLACE) end >> $@

$($(cmod_)_SRC_CORE_SEQ)%_di.c: cmod:=$(cmod)
$($(cmod_)_SRC_CORE_SEQ)%_di.c: cmod_:=$(cmod_)
$($(cmod_)_SRC_CORE_SEQ)%_di.c: $($(cmod_)_SRC_CORE_SEQ)%.c $($(cmod_)_FAR_SCRIPT_DI)
	@$(ECHO) "making di-version '$@' from '$<'"
	@$(ECHO) "#define SL_DATA_IGNORE" > $@
	@$($(cmod_)_FIND_AND_REPLACE) exec $(word 2,$^) $< >> $@

$($(cmod_)_SRC_CORE_MPI)%_di.c: cmod:=$(cmod)
$($(cmod_)_SRC_CORE_MPI)%_di.c: cmod_:=$(cmod_)
$($(cmod_)_SRC_CORE_MPI)%_di.c: $($(cmod_)_SRC_CORE_MPI)%.c $($(cmod_)_FAR_SCRIPT_DI)
	@$(ECHO) "making mpi-di-version '$@' from '$<'"
	@$(ECHO) "#define SL_DATA_IGNORE" > $@
	@$($(cmod_)_FIND_AND_REPLACE) exec $(word 2,$^) $< >> $@

.PRECIOUS: $($(cmod_)_SRC_CORE_SEQ)%_di.c $($(cmod_)_SRC_CORE_MPI)%_di.c
 endif
endif


# tune
#$(cmod_)_tune_seq: $(cmod_)_libsl_seq $($(cmod_)_TUNE_SEQ)
#$(cmod_)_tune_mpi: $(cmod_)_libsl_mpi $($(cmod_)_TUNE_MPI)
$(cmod_)_tune: $(cmod_)_libsl $($(cmod_)_TUNE)

#M_TUNE_SEQ+=$(cmod_)_tune_seq
#M_TUNE_MPI+=$(cmod_)_tune_mpi
M_TUNE+=$(cmod_)_tune

# tune files
$($(cmod_)_TMP_OBJS_SL_TUNE_SEQ): CC:=$($(cmod_)_CC)
$($(cmod_)_TMP_OBJS_SL_TUNE_SEQ): CFLAGS:=$(strip $($(cmod_)_CFLAGS_) $($(cmod_)_CFLAGS_SL_) $(call or2,$(findstring -DSL_TUNEABLE,$($(cmod_)_CFLAGS_SL_)),-DSL_TUNEABLE))
$($(cmod_)_TMP_OBJS_SL_TUNE_SEQ): TYPE=tune

$($(cmod_)_TMP_OBJS_SL_TUNE_MPI): CC:=$($(cmod_)_MPICC)
$($(cmod_)_TMP_OBJS_SL_TUNE_MPI): CFLAGS:=$(strip $($(cmod_)_MPICFLAGS_) $($(cmod_)_MPICFLAGS_SL_) $(call or2,$(findstring -DSL_TUNEABLE,$($(cmod_)_MPICFLAGS_SL_)),-DSL_TUNEABLE))
$($(cmod_)_TMP_OBJS_SL_TUNE_MPI): TYPE=tune-mpi

$($(cmod_)_TMP_TUNE_SEQ)%.o: $($(cmod_)_SRC_TUNE_SEQ)%.c $($(cmod_)_SRC_TUNE_SEQ_INCLUDE) $($(cmod_)_CODE_SL_INCLUDE) $($(cmod_)_PREREQ)
	$($(cmod_)_comp)

$($(cmod_)_TMP_TUNE_MPI)%.o: $($(cmod_)_SRC_TUNE_MPI)%.c $($(cmod_)_SRC_TUNE_MPI_INCLUDE) $($(cmod_)_CODE_SL_INCLUDE) $($(cmod_)_PREREQ)
	$($(cmod_)_comp)

$($(cmod_)_TUNE): LD:=$($(cmod_)_LD)
$($(cmod_)_TUNE): LDFLAGS:=$(strip $($(cmod_)_LDFLAGS_))
$($(cmod_)_TUNE): TYPE:=tune

$($(cmod_)_TUNE): $($(cmod_)_TMP_OBJS_SL_TUNE) $($(cmod_)_LIBSL)
	$($(cmod_)_make)


# devel
#$(cmod_)_devel_seq: $(cmod_)_libsl_seq $($(cmod_)_DEVL_SEQ)
#$(cmod_)_devel_mpi: $(cmod_)_libsl_mpi $($(cmod_)_DEVL_MPI)
$(cmod_)_devel: $(cmod_)_libsl $($(cmod_)_DEVL)

#M_DEVL_SEQ+=$(cmod_)_devel_seq
#M_DEVL_MPI+=$(cmod_)_devel_mpi
M_DEVL+=$(cmod_)_devel

# devel files
$($(cmod_)_TMP_OBJS_SL_DEVL_SEQ): CC:=$($(cmod_)_CC)
$($(cmod_)_TMP_OBJS_SL_DEVL_SEQ): CFLAGS:=$(strip $($(cmod_)_CFLAGS_) $($(cmod_)_CFLAGS_SL_))
$($(cmod_)_TMP_OBJS_SL_DEVL_SEQ): TYPE=devel

$($(cmod_)_TMP_OBJS_SL_DEVL_MPI): CC:=$($(cmod_)_MPICC)
$($(cmod_)_TMP_OBJS_SL_DEVL_MPI): CFLAGS:=$(strip $($(cmod_)_MPICFLAGS_) $($(cmod_)_MPICFLAGS_SL_))
$($(cmod_)_TMP_OBJS_SL_DEVL_MPI): TYPE=devel-mpi

$($(cmod_)_TMP_DEVL_SEQ)%.o: $($(cmod_)_SRC_DEVL_SEQ)%.c $($(cmod_)_SRC_DEVL_SEQ_INCLUDE) $($(cmod_)_CODE_SL_INCLUDE) $($(cmod_)_PREREQ)
	$($(cmod_)_comp)

$($(cmod_)_TMP_DEVL_MPI)%.o: $($(cmod_)_SRC_DEVL_MPI)%.c $($(cmod_)_SRC_DEVL_MPI_INCLUDE) $($(cmod_)_CODE_SL_INCLUDE) $($(cmod_)_PREREQ)
	$($(cmod_)_comp)

$($(cmod_)_DEVL): LD:=$($(cmod_)_LD)
$($(cmod_)_DEVL): LDFLAGS:=$(strip $($(cmod_)_LDFLAGS_))
$($(cmod_)_DEVL): TYPE:=devel

$($(cmod_)_DEVL): $($(cmod_)_TMP_OBJS_SL_DEVL) $($(cmod_)_LIBSL)
	$($(cmod_)_make)


# demo
#$(cmod_)_demo_seq: $(cmod_)_libsl_seq $($(cmod_)_DEMO_SEQ)
#$(cmod_)_demo_mpi: $(cmod_)_libsl_mpi $($(cmod_)_DEMO_MPI)
$(cmod_)_demo: $(cmod_)_libsl $($(cmod_)_DEMO)

#M_DEMO_SEQ+=$(cmod_)_demo_seq
#M_DEMO_MPI+=$(cmod_)_demo_mpi
M_DEMO+=$(cmod_)_demo

# demo files
$($(cmod_)_TMP_OBJS_SL_DEMO_SEQ): CC:=$($(cmod_)_CC)
$($(cmod_)_TMP_OBJS_SL_DEMO_SEQ): CFLAGS:=$(strip $($(cmod_)_CFLAGS_) $($(cmod_)_CFLAGS_SL_))
$($(cmod_)_TMP_OBJS_SL_DEMO_SEQ): TYPE=demo

$($(cmod_)_TMP_OBJS_SL_DEMO_MPI): CC:=$($(cmod_)_MPICC)
$($(cmod_)_TMP_OBJS_SL_DEMO_MPI): CFLAGS:=$(strip $($(cmod_)_MPICFLAGS_) $($(cmod_)_MPICFLAGS_SL_))
$($(cmod_)_TMP_OBJS_SL_DEMO_MPI): TYPE=demo-mpi

$($(cmod_)_TMP_DEMO_SEQ)%.o: $($(cmod_)_SRC_DEMO_SEQ)%.c $($(cmod_)_SRC_DEMO_SEQ_INCLUDE) $($(cmod_)_CODE_SL_INCLUDE) $($(cmod_)_PREREQ)
	$($(cmod_)_comp)

$($(cmod_)_TMP_DEMO_MPI)%.o: $($(cmod_)_SRC_DEMO_MPI)%.c $($(cmod_)_SRC_DEMO_MPI_INCLUDE) $($(cmod_)_CODE_SL_INCLUDE) $($(cmod_)_PREREQ)
	$($(cmod_)_comp)

$($(cmod_)_DEMO): LD:=$($(cmod_)_LD)
$($(cmod_)_DEMO): LDFLAGS:=$(strip $($(cmod_)_LDFLAGS_))
$($(cmod_)_DEMO): TYPE:=demo

$($(cmod_)_DEMO): $($(cmod_)_TMP_OBJS_SL_DEMO) $($(cmod_)_LIBSL)
	$($(cmod_)_make)


# wrapper
$(cmod_)_wrapper_seq: $($(cmod_)_TMP_OBJS_SL_WRAP_SEQ)
$(cmod_)_wrapper_mpi: $($(cmod_)_TMP_OBJS_SL_WRAP_MPI)
$(cmod_)_wrapper: $($(cmod_)_TMP_OBJS_SL_WRAP)

M_WRAP_SEQ+=$(cmod_)_wrapper_seq
M_WRAP_MPI+=$(cmod_)_wrapper_mpi
M_WRAP+=$(cmod_)_wrapper

# wrapper files
$($(cmod_)_TMP_OBJS_SL_WRAP_SEQ): CC:=$($(cmod_)_CC)
$($(cmod_)_TMP_OBJS_SL_WRAP_SEQ): CFLAGS:=$(strip $($(cmod_)_CFLAGS_) $(addprefix -I,$(M_INCDIR) $(INCDIR)))
$($(cmod_)_TMP_OBJS_SL_WRAP_SEQ): TYPE:=wrapper

$($(cmod_)_TMP_OBJS_SL_WRAP_MPI): CC:=$($(cmod_)_MPICC)
$($(cmod_)_TMP_OBJS_SL_WRAP_MPI): CFLAGS:=$(strip $($(cmod_)_MPICFLAGS_) $(addprefix -I,$(M_INCDIR) $(INCDIR)))
$($(cmod_)_TMP_OBJS_SL_WRAP_MPI): TYPE:=wrapper-mpi

$($(cmod_)_TMP_WRAP_SEQ)%.o: $($(cmod_)_SRC_WRAP_SEQ)%.c $($(cmod_)_CODE_SL_WRAP_SEQ_INCLUDE) $(M_CODE_INCLUDE) $($(cmod_)_PREREQ)
	$($(cmod_)_comp)

$($(cmod_)_TMP_WRAP_MPI)%.o: $($(cmod_)_SRC_WRAP_MPI)%.c $($(cmod_)_CODE_SL_WRAP_MPI_INCLUDE) $(M_CODE_INCLUDE) $($(cmod_)_PREREQ)
	$($(cmod_)_comp)


# libsl
$(cmod_)_libsl_seq: $(cmod_)_core_seq $(cmod_)_wrapper_seq $($(cmod_)_LIBSL_SEQ)($($(cmod_)_OBJS_SL_CORE_SEQ) $($(cmod_)_OBJS_SL_WRAP_SEQ))
$(cmod_)_libsl_mpi: $(cmod_)_core_mpi $(cmod_)_wrapper_mpi $($(cmod_)_LIBSL_MPI)($($(cmod_)_OBJS_SL_CORE_MPI) $($(cmod_)_OBJS_SL_WRAP_MPI))
$(cmod_)_libsl: $(cmod_)_core $(cmod_)_wrapper $($(cmod_)_LIBSL)($($(cmod_)_OBJS_SL_CORE) $($(cmod_)_OBJS_SL_WRAP))

M_LIBSL_SEQ+=$(cmod_)_libsl_seq
M_LIBSL_MPI+=$(cmod_)_libsl_mpi
M_LIBSL+=$(cmod_)_libsl

ifneq ($($(cmod_)_LIBSL_SEQ),$(LIBSL_SEQ))
$($(cmod_)_LIBSL_SEQ): $($(cmod_)_LIBSL_SEQ)($($(cmod_)_OBJS_SL_CORE_SEQ) $($(cmod_)_OBJS_SL_WRAP_SEQ))
else
M_OBJS_SL_CORE_SEQ+=$($(cmod_)_OBJS_SL_CORE_SEQ)
M_OBJS_SL_WRAP_SEQ+=$($(cmod_)_OBJS_SL_WRAP_SEQ)
endif

ifneq ($($(cmod_)_LIBSL_MPI),$(LIBSL_MPI))
$($(cmod_)_LIBSL_MPI): $($(cmod_)_LIBSL_MPI)($($(cmod_)_OBJS_SL_CORE_MPI) $($(cmod_)_OBJS_SL_WRAP_MPI))
else
M_OBJS_SL_CORE_MPI+=$($(cmod_)_OBJS_SL_CORE_MPI)
M_OBJS_SL_WRAP_MPI+=$($(cmod_)_OBJS_SL_WRAP_MPI)
endif

ifneq ($($(cmod_)_LIBSL),$(LIBSL))
$($(cmod_)_LIBSL): $($(cmod_)_LIBSL)($($(cmod_)_OBJS_SL_CORE) $($(cmod_)_OBJS_SL_WRAP)
else
M_OBJS_SL_CORE+=$($(cmod_)_OBJS_SL_CORE)
M_OBJS_SL_WRAP+=$($(cmod_)_OBJS_SL_WRAP)
endif

%($($(cmod_)_OBJS_SL_CORE_SEQ)): TYPE:=core
%($($(cmod_)_OBJS_SL_CORE_MPI)): TYPE:=core-mpi
%($($(cmod_)_OBJS_SL_WRAP_SEQ)): TYPE:=wrapper
%($($(cmod_)_OBJS_SL_WRAP_MPI)): TYPE:=wrapper-mpi

$($(cmod_)_LIBSL_SEQ)(%): $($(cmod_)_TMP_CORE_SEQ)%
	$($(cmod_)_add)

$($(cmod_)_LIBSL_SEQ)(%): $($(cmod_)_TMP_WRAP_SEQ)%
	$($(cmod_)_add)

$($(cmod_)_LIBSL_MPI)(%): $($(cmod_)_TMP_CORE_MPI)%
	$($(cmod_)_add)

$($(cmod_)_LIBSL_MPI)(%): $($(cmod_)_TMP_WRAP_MPI)%
	$($(cmod_)_add)

$($(cmod_)_LIBSL)(%): $($(cmod_)_TMP_CORE_SEQ)%
	$($(cmod_)_add)

$($(cmod_)_LIBSL)(%): $($(cmod_)_TMP_CORE_MPI)%
	$($(cmod_)_add)

$($(cmod_)_LIBSL)(%): $($(cmod_)_TMP_WRAP_SEQ)%
	$($(cmod_)_add)

$($(cmod_)_LIBSL)(%): $($(cmod_)_TMP_WRAP_MPI)%
	$($(cmod_)_add)


# clean
$(cmod_)_TMPDIRS:=
ifneq ($($(cmod_)_TMP_CORE_SEQ),$($(cmod_)_SRC_CORE_SEQ))
$(cmod_)_TMPDIRS+=$($(cmod_)_TMP_CORE_SEQ)
endif
ifneq ($($(cmod_)_TMP_CORE_MPI),$($(cmod_)_SRC_CORE_MPI))
$(cmod_)_TMPDIRS+=$($(cmod_)_TMP_CORE_MPI)
endif
ifneq ($($(cmod_)_TMP_TUNE_SEQ),$($(cmod_)_SRC_TUNE_SEQ))
$(cmod_)_TMPDIRS+=$($(cmod_)_TMP_TUNE_SEQ)
endif
ifneq ($($(cmod_)_TMP_TUNE_MPI),$($(cmod_)_SRC_TUNE_MPI))
$(cmod_)_TMPDIRS+=$($(cmod_)_TMP_TUNE_MPI)
endif
ifneq ($($(cmod_)_TMP_DEVL_SEQ),$($(cmod_)_SRC_DEVL_SEQ))
$(cmod_)_TMPDIRS+=$($(cmod_)_TMP_DEVL_SEQ)
endif
ifneq ($($(cmod_)_TMP_DEVL_MPI),$($(cmod_)_SRC_DEVL_MPI))
$(cmod_)_TMPDIRS+=$($(cmod_)_TMP_DEVL_MPI)
endif
ifneq ($($(cmod_)_TMP_DEMO_SEQ),$($(cmod_)_SRC_DEMO_SEQ))
$(cmod_)_TMPDIRS+=$($(cmod_)_TMP_DEMO_SEQ)
endif
ifneq ($($(cmod_)_TMP_DEMO_MPI),$($(cmod_)_SRC_DEMO_MPI))
$(cmod_)_TMPDIRS+=$($(cmod_)_TMP_DEMO_MPI)
endif
ifneq ($($(cmod_)_TMPDIR),$($(cmod_)_SRCDIR))
$(cmod_)_TMPDIRS+=$($(cmod_)_TMPDIR)
endif


$(cmod_)_clean: cmod:=$(cmod)
$(cmod_)_clean: cmod_:=$(cmod_)
$(cmod_)_clean:
ifneq ($(cmod),)
	@$(ECHO) "$(cmod):"
endif
	@$(ECHO) " cleaning temporary, object & archive files"
	$(RM) $($(cmod_)_FAR_SCRIPT_DI)
ifneq ($($(cmod_)_DI_VERSIONS),config)
	$(RM) $($(cmod_)_SRC_CORE_SEQ)*_di.c
	$(RM) $($(cmod_)_SRC_CORE_MPI)*_di.c
endif
	$(RM) $(strip $($(cmod_)_TMP_OBJS_SL_CORE_SEQ) $($(cmod_)_TMP_OBJS_SL_CORE_MPI))
	$(RM) $(strip $($(cmod_)_TMP_OBJS_SL_TUNE_SEQ) $($(cmod_)_TMP_OBJS_SL_TUNE_MPI) $($(cmod_)_TUNE))
	$(RM) $(strip $($(cmod_)_TMP_OBJS_SL_DEVL_SEQ) $($(cmod_)_TMP_OBJS_SL_DEVL_MPI)	$($(cmod_)_DEVL))
	$(RM) $(strip $($(cmod_)_TMP_OBJS_SL_DEMO_SEQ) $($(cmod_)_TMP_OBJS_SL_DEMO_MPI)	$($(cmod_)_DEMO))
	$(RM) $(strip $($(cmod_)_TMP_OBJS_SL_WRAP_SEQ) $($(cmod_)_TMP_OBJS_SL_WRAP_MPI))
ifneq ($($(cmod_)_LIBSL_SEQ),$(LIBSL_SEQ))
	$(RM) $($(cmod_)_LIBSL_SEQ)
endif
ifneq ($($(cmod_)_LIBSL_MPI),$(LIBSL_MPI))
	$(RM) $($(cmod_)_LIBSL_MPI)
endif
ifneq ($($(cmod_)_LIBSL),$(LIBSL))
	$(RM) $($(cmod_)_LIBSL)
endif
	$(RM) -r $($(cmod_)_TMPDIRS)

M_CLEAN+=$(cmod_)_clean


# top level (no module)
ifeq ($(mod),)


info: $(M_INFO)
	@$(ECHO) "objs-core-seq: $(M_OBJS_SL_CORE_SEQ)"
	@$(ECHO) "objs-core-mpi: $(M_OBJS_SL_CORE_MPI)"
	@$(ECHO) "objs-core: $(M_OBJS_SL_CORE)"
#	@$(ECHO) "tmp-objs-core: $(TMP_OBJS_SL_CORE_SEQ)"
#	@$(ECHO) "tmp-objs-core-mpi: $(TMP_OBJS_SL_CORE_MPI)"

#	@$(ECHO) "objs-tune: $(M_OBJS_SL_TUNE_SEQ)"
#	@$(ECHO) "objs-tune-mpi: $(M_OBJS_SL_TUNE_MPI)"
#	@$(ECHO) "tmp-objs-tune: $(TMP_OBJS_SL_TUNE_SEQ)"
#	@$(ECHO) "tmp-objs-tune-mpi: $(TMP_OBJS_SL_TUNE_MPI)"

#	@$(ECHO) "objs-devel: $(M_OBJS_SL_DEVL_SEQ)"
#	@$(ECHO) "objs-devel-mpi: $(M_OBJS_SL_DEVL_MPI)"
#	@$(ECHO) "tmp-objs-devel: $(TMP_OBJS_SL_DEVL_SEQ)"
#	@$(ECHO) "tmp-objs-devel-mpi: $(TMP_OBJS_SL_DEVL_MPI)"

#	@$(ECHO) "objs-demo: $(M_OBJS_SL_DEMO_SEQ)"
#	@$(ECHO) "objs-demo-mpi: $(M_OBJS_SL_DEMO_MPI)"
#	@$(ECHO) "tmp-objs-demo: $(TMP_OBJS_SL_DEMO_SEQ)"
#	@$(ECHO) "tmp-objs-demo-mpi: $(TMP_OBJS_SL_DEMO_MPI)"

	@$(ECHO) "objs-wrapper-seq: $(M_OBJS_SL_WRAP_SEQ)"
	@$(ECHO) "objs-wrapper-mpi: $(M_OBJS_SL_WRAP_MPI)"
	@$(ECHO) "objs-wrapper: $(M_OBJS_SL_WRAP)"
#	@$(ECHO) "tmp-objs-wrapper: $(TMP_OBJS_SL_WRAP_SEQ)"
#	@$(ECHO) "tmp-objs-wrapper-mpi: $(TMP_OBJS_SL_WRAP_MPI)"

core_seq: $(M_CORE_SEQ)
core_mpi: $(M_CORE_MPI)
core: $(M_CORE)

wrapper_seq: $(M_WRAP_SEQ)
wrapper_mpi: $(M_WRAP_MPI)
wrapper: $(M_WRAP)

libsl_seq: $(M_LIBSL_SEQ)
libsl_mpi: $(M_LIBSL_MPI)
libsl: $(M_LIBSL)

$(LIBSL_SEQ): $(LIBSL_SEQ)($(M_OBJS_SL_CORE_SEQ) $(M_OBJS_SL_WRAP_SEQ))
$(LIBSL_MPI): $(LIBSL_MPI)($(M_OBJS_SL_CORE_MPI) $(M_OBJS_SL_WRAP_MPI))
$(LIBSL): $(LIBSL)($(M_OBJS_SL_CORE) $(M_OBJS_SL_WRAP))

tune: $(M_TUNE)

devel: $(M_DEVL)

demo: $(M_DEMO)

clean: $(M_CLEAN)
	$(RM) $(LIBSL_SEQ) $(LIBSL_MPI) $(LIBSL) $(TUNE) $(DEVL) $(DEMO)

endif
