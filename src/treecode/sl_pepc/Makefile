# This is a Makefile!

#  
#  SL - Sorting Library, v0.1, (michael.hofmann@informatik.tu-chemnitz.de)
#  
#  file: Makefile
#  


# only on top level
ifeq ($(mod),)

all: libsl.a
	@echo "==== building the sorting library" 

# initial config
MAKEFILE_IN:=Makefile.pepc

SL_USE_MPI:=
ALL_AS_MPI:=
DI_VERSIONS:=
FAR:=perl
TUNEABLE:=#TUNEABLE#
TARGET:=
QUIET:=

#SLDEBUG:=5
#SLDEBUG_OUTPUT:=5
#SLDEBUG_OUTPUT_NOT:=true
#SLDEBUG_ALLOC:=1
#SLDEBUG_ALLOC_NOT:=true

WRAPPER_SRC:=
WRAPPER_PREFIX:=
WRAPPER_MPI_SRC:=
WRAPPER_MPI_PREFIX:=mpi_

EXEC:=a.out

EXEC_SRC:=DUMMY
EXEC_PREFIX:=
EXEC_MPI_SRC:=DUMMY
EXEC_MPI_PREFIX:=mpi_

SRCBASE:=

INCDIR:=.
LIBDIR:=
INTDIR:=#INTDIR#
TMPDIR:=#TMPDIR#

ALL:=libsl

DEPRECATED:=#DEPRECATED#
DEVEL:=#DEVEL#

BULK_AR:=true

# defaults
DEFAULT_SL_USE_MPI:=true
DEFAULT_ALL_AS_MPI:=

SL_CONFIG_INCLUDE=sl_config.h sl_tune.h sl_environment.h


or2 = $(if $(1),$(1),$(2))
or3 = $(if $(1),$(1),$(if $(2),$(2),$(3)))
or3 = $(if $(1),$(1),$(if $(2),$(2),$(if $(3),$(3),$(4))))

and2 = $(if $(1),$(2),)
and3 = $(if $(1),$(if $(2),$(3),),)
and4 = $(if $(1),$(if $(2),$(if $(3),$(4),),),)


mod_list:=

ifeq ($(MODULES),)
 MODULES:=$(call or2,$(MODULES),$(dir $(wildcard sl_*/.)))
endif
#$(info modules: $(MODULES))

MAKEFILE_IN:=$(call or2,$(MAKEFILE_IN),Makefile.in)

# workaround problems with 'include' and non-existing files
MAKEFILE_IN:=$(wildcard $(MAKEFILE_IN))

-include $(MAKEFILE_IN)

RM:=$(call or2,$(RM),rm -f)
ECHO:=$(call or2,$(ECHO),echo)
MKDIR_P:=$(call or2,$(MKDIR_P),mkdir -p)
RANLIB:=$(call or2,$(RANLIB),ranlib)

SL_CFLAGS:=
SL_CPPFLAGS:=
SL_MPICFLAGS:=
SL_MPICPPFLAGS:=
SL_LDFLAGS:=
SL_MPILDFLAGS:=
SL_LDLIBS:=
SL_MPILDLIBS:=

ifneq ($(DEBUG),)
 SLDEBUG:=$(DEBUG)
endif

ifneq ($(SLDEBUG),)
 SL_CPPFLAGS+=-DSLDEBUG=$(SLDEBUG)
endif

ifneq ($(SLDEBUG_OUTPUT),)
 SL_CPPFLAGS+=-DSLDEBUG_OUTPUT=$(SLDEBUG_OUTPUT)
endif

ifneq ($(SLDEBUG_ALLOC),)
 SL_CPPFLAGS+=-DSLDEBUG_ALLOC=$(SLDEBUG_ALLOC)
endif

ifneq ($(SLDEBUG_ALLOC_NOT),)
 SL_CPPFLAGS+=-DSLDEBUG_ALLOC_NOT=$(SLDEBUG_ALLOC_NOT)
endif


LIBSL_A_SEQ:=$(call or2,$(LIBSL_A_SEQ),$(LIBDIR)libsl_seq.a)
LIBSL_A_MPI:=$(call or2,$(LIBSL_A_MPI),$(LIBDIR)libsl_mpi.a)
LIBSL_A:=$(call or2,$(LIBSL_A),$(LIBDIR)libsl.a)

EXEC:=$(call or2,$(EXEC),$(LIBDIR)a.out)

TUNE:=sltune
DEVL:=sldevel
DEMO:=sldemo

PREREQ+=Makefile $(MAKEFILE_IN)

M_SL_USE_MPI:=

M_INFO:=

M_CORE_SEQ:=
M_CORE_MPI:=
M_CORE:=

M_LIBSL_SEQ:=
M_LIBSL_MPI:=
M_LIBSL:=

M_EXEC:=

M_TUNE_SEQ:=
M_TUNE_MPI:=
M_TUNE:=

M_DEVL_SEQ:=
M_DEVL_MPI:=
M_DEVL:=

M_DEMO_SEQ:=
M_DEMO_MPI:=
M_DEMO:=

M_WRAP_SEQ:=
M_WRAP_MPI:=
M_WRAP:=

M_CLEAN:=

M_INCDIR:=
M_CODE_INCLUDE:=

#M_OBJS_SL_CORE_SEQ:=
#M_OBJS_SL_CORE_MPI:=
#M_OBJS_SL_CORE:=
#
#M_TMP_OBJS_SL_CORE_SEQ:=
#M_TMP_OBJS_SL_CORE_MPI:=
#M_TMP_OBJS_SL_CORE:=
#
#M_OBJS_SL_WRAP_SEQ:=
#M_OBJS_SL_WRAP_MPI:=
#M_OBJS_SL_WRAP:=
#
#M_TMP_OBJS_SL_WRAP_SEQ:=
#M_TMP_OBJS_SL_WRAP_MPI:=
#M_TMP_OBJS_SL_WRAP:=

M_TMP_OBJS_SL_EXEC_SEQ:=
M_TMP_OBJS_SL_EXEC_MPI:=
M_TMP_OBJS_SL_EXEC:=


sl_all: $(call or2,$(ALL),demo)

else

MODULES:=#MODULES#
MODULES_NOT:=#MODULES#

mod:=$(patsubst %/,%,$(mod))

endif


# push mod on the list
mod_list:=$(mod) $(mod_list)

$(mod)_MODULES:=$(patsubst %/,%,$(MODULES))
$(mod)_MODULES_NOT:=$(patsubst %/,%,$(MODULES_NOT))

# include modules
ifneq ($(MODULES),)
 ifeq ($(findstring sl_pepckeys,$($(mod)_MODULES_NOT)),)
  mod:=sl_pepckeys/
  include Makefile
  mod:=$(firstword $(mod_list))
 endif
 ifeq ($(findstring sl_pepcparts,$($(mod)_MODULES_NOT)),)
  mod:=sl_pepcparts/
  include Makefile
  mod:=$(firstword $(mod_list))
 endif
endif


#$(warning $(mod)_MODULES: $($(mod)_MODULES))
#$(warning $(mod)_MODULES_NOT: $($(mod)_MODULES_NOT))

# pop mod from the list
mod:=$(firstword $(mod_list))
mod_list:=$(wordlist 2,$(words $(mod_list)),$(mod_list))


# current module
cmod:=$(patsubst %/,%,$(mod))
#cmod_:=$(subst /,_,$(subst .,_,$(cmod)))
cmod_:=$(cmod)
cmod_dir:=$(call and2,$(cmod),$(cmod)/)
cmod_prefix:=$(call and2,$(cmod),$(cmod)_)

#$(warning mod: $(mod))
#$(warning cmod: $(cmod))
#$(warning cmod_: $(cmod_))
#$(warning cmod_dir: $(cmod_dir))
#$(warning cmod_prefix: $(cmod_prefix))

# source directory
ifeq ($($(cmod_)_SRCDIR),)
 ifneq ($(wildcard $(cmod_dir)src/.),)
  $(cmod_)_SRCDIR:=$(cmod_dir)src/
 else
  ifneq ($(wildcard $(SRCBASE)src/.),)
   $(cmod_)_SRCDIR:=$(SRCBASE)src/
  endif
 endif
endif

# config directory
ifeq ($($(cmod_)_CFGDIR),)
 ifneq ($(wildcard $(cmod_dir)src/include/.),)
  $(cmod_)_CFGDIR:=$(cmod_dir)src/include/
 else
  $(cmod_)_CFGDIR:=$(cmod_dir)
 endif
endif

# temp directory
$(cmod_)_TMPDIR:=$(call or3,$($(cmod_)_TMPDIR),$(TMPDIR),$(cmod_dir)tmp/)

# library directory
$(cmod_)_LIBDIR:=$(call or3,$($(cmod_)_LIBDIR),$(LIBDIR))

# interface directory
ifeq ($($(cmod_)_INTDIR),)
 ifneq ($(wildcard $(cmod_dir)include/.),)
  $(cmod_)_INTDIR:=$(cmod_dir)include/
 endif
endif

# include files
$(cmod_)_INCLUDE:=$($(cmod_)_INCDIR) $(INCDIR)
$(cmod_)_SL_INCLUDE:=$(call and2,$($(cmod_)_SRCDIR),$($(cmod_)_SRCDIR)include/)
$(cmod_)_SL_CFG_INCLUDE:=$($(cmod_)_CFGDIR)

$(cmod_)_CODE_INCLUDE:=$(foreach d,$($(cmod_)_INCLUDE),$(wildcard $(d)*.h))
$(cmod_)_CODE_SL_INCLUDE:=$(filter-out $(addprefix %,$(SL_CONFIG_INCLUDE)),$(foreach d,$($(cmod_)_SL_INCLUDE),$(wildcard $(d)*.h)))
$(cmod_)_CODE_SL_CFG_INCLUDE:=$(filter $(addprefix %,$(SL_CONFIG_INCLUDE)),$(foreach d,$($(cmod_)_SL_CFG_INCLUDE),$(wildcard $(d)*.h)))

ifneq ($($(cmod_)_SL_INCLUDE),$($(cmod_)_SL_CFG_INCLUDE))
 DOUBLE_CFG_INCLUDE:=$(filter $(addprefix %,$(SL_CONFIG_INCLUDE)),$(foreach d,$($(cmod_)_SL_INCLUDE),$(wildcard $(d)*.h)))
 ifneq ($(DOUBLE_CFG_INCLUDE),)
  $(error found config files '$(notdir $(DOUBLE_CFG_INCLUDE))' in includes of '$($(cmod_)_SL_INCLUDE)', this is bad!)
 endif
endif

ifneq ($(findstring sl_config.h,$($(cmod_)_CODE_SL_CFG_INCLUDE)),)
 cmod_cfg:=1
else
 cmod_cfg:=
endif

#$(warning cmod_cfg: $(cmod_cfg))

ifeq ($(cmod_cfg),)

 ifeq ($(wildcard sl_*/.),)
   $(error sl_config.h not found in '$($(cmod_)_SL_CFG_INCLUDE)', the sorting library does not seem to be properly configured)
 endif
endif

M_INCDIR+=$($(cmod_)_INCDIR) $($(cmod_)_INTDIR)

M_CODE_INCLUDE+=$($(cmod_)_CODE_INCLUDE)


# CC, CPP, LD and flags
$(cmod_)_CC:=$(call or2,$($(cmod_)_CC),$(CC))
$(cmod_)_CFLAGS:=$(CFLAGS) $(SL_CFLAGS) $($(cmod_)_CFLAGS)
$(cmod_)_CFLAGS_SL:=$(CFLAGS_SL) $($(cmod_)_CFLAGS_SL)
$(cmod_)_CFLAGS_NSL:=$(CFLAGS_NSL) $($(cmod_)_CFLAGS_NSL)
$(cmod_)_CPPFLAGS:=$(CPPFLAGS) $(SL_CPPFLAGS) $($(cmod_)_CPPFLAGS) $(addprefix -I,$($(cmod_)_INCLUDE))
$(cmod_)_CPPFLAGS_SL:=$(CPPFLAGS_SL) $($(cmod_)_CPPFLAGS_SL) $(addprefix -I,$($(cmod_)_SL_CFG_INCLUDE) $(filter-out $($(cmod_)_SL_CFG_INCLUDE),$($(cmod_)_SL_INCLUDE)))
$(cmod_)_CPPFLAGS_NSL:=$(CPPFLAGS_NSL) $($(cmod_)_CPPFLAGS_NSL) $(addprefix -I,$(M_INCDIR)) $(addprefix -DNOT_,$($(mod)_MODULES_NOT))

$(cmod_)_MPICC:=$(call or2,$($(cmod_)_MPICC),$(MPICC))
$(cmod_)_MPICFLAGS:=$(MPICFLAGS) $(SL_MPICFLAGS) $($(cmod_)_MPICFLAGS)
$(cmod_)_MPICFLAGS_SL:=$(MPICFLAGS_SL) $($(cmod_)_MPICFLAGS_SL)
$(cmod_)_MPICFLAGS_NSL:=$(MPICFLAGS_NSL) $($(cmod_)_MPICFLAGS_NSL)
$(cmod_)_MPICPPFLAGS:=$(MPICPPFLAGS) $(SL_MPICPPFLAGS) $($(cmod_)_MPICPPFLAGS)
$(cmod_)_MPICPPFLAGS_SL:=$(MPICPPFLAGS_SL) $($(cmod_)_MPICPPFLAGS_SL)
$(cmod_)_MPICPPFLAGS_NSL:=$(MPICPPFLAGS_NSL) $($(cmod_)_MPICPPFLAGS_NSL) -DSL_USE_MPI

$(cmod_)_LD:=$(call or2,$($(cmod_)_LD),$(LD))
$(cmod_)_LDFLAGS:=$(LDFLAGS) $(SL_LDFLAGS) $($(cmod_)_LDFLAGS)
$(cmod_)_LDLIBS:=$(LDLIBS) $(SL_LDLIBS) $($(cmod_)_LDLIBS)

$(cmod_)_MPILD:=$(call or2,$($(cmod_)_MPILD),$(MPILD))
$(cmod_)_MPILDFLAGS:=$(MPILDFLAGS) $(SL_MPILDFLAGS) $($(cmod_)_MPILDFLAGS)
$(cmod_)_MPILDLIBS:=$(MPILDLIBS) $(SL_MPILDLIBS) $($(cmod_)_MPILDLIBS)


# options
$(cmod_)_SL_USE_MPI:=$(call or3,$($(cmod_)_SL_USE_MPI),$(SL_USE_MPI),$(DEFAULT_SL_USE_MPI))

# SL_USE_MPI = true | false | detect (default)
ifneq ($($(cmod_)_SL_USE_MPI),true) # not true
 ifneq ($($(cmod_)_SL_USE_MPI),false) # not false
  # else: use configuration file
  ifneq ($(cmod_cfg),)
   ifneq ($(shell $(PRINTF) "\#include \"sl_config.h\"\nSL_USE_MPI\n" | $(CPP) $($(cmod_)_CPPFLAGS_SL) 2>/dev/null | $(GREP) SL_USE_MPI),SL_USE_MPI)
    $(cmod_)_SL_USE_MPI:=true
   else
    $(cmod_)_SL_USE_MPI:=false
   endif
  endif
 endif
endif

ifneq ($(cmod_cfg),)
 ifeq ($($(cmod_)_SL_USE_MPI),true)
  M_SL_USE_MPI:=true
 endif
endif

$(cmod_)_ALL_AS_MPI:=$(call or3,$($(cmod_)_ALL_AS_MPI),$(ALL_AS_MPI),$(DEFAULT_ALL_AS_MPI))

# ALL_AS_MPI = true | false | SL_USE_MPI (default)
ifneq ($($(cmod_)_ALL_AS_MPI),true) # not true
 ifneq ($($(cmod_)_ALL_AS_MPI),false) # not false
  # else: use SL_USE_MPI or not
  $(cmod_)_ALL_AS_MPI:=$($(cmod_)_SL_USE_MPI)
#  $(cmod_)_ALL_AS_MPI:=
 endif
endif

$(cmod_)_DI_VERSION:=$(call or3,$($(cmod_)_DI_VERSION),$(DI_VERSIONS),false)

$(cmod_)_TUNEABLE:=$(call or2,$($(cmod_)_TUNEABLE),$(TUNEABLE))

#$(cmod_)_TARGET:=$(call or3,$($(cmod_)_TARGET),$(TARGET),NOT)

$(cmod_)_PREFIX:=$(call and2,$(patsubst sl_%,%,$(cmod)),$(patsubst sl_%,%,$(cmod))_)

$(cmod_)_DEPRECATED:=$(call or2,$($(cmod_)_DEPRECATED),$(DEPRECATED))
$(cmod_)_DEVEL:=$(call or2,$($(cmod_)_DEVEL),$(DEVEL))


# option dependent flags
$(cmod_)_CPPFLAGS_SL+=$(call and2,$($(cmod_)_TUNEABLE),-DSL_TUNEABLE) -DSL_PREFIX=$($(cmod_)_PREFIX)
$(cmod_)_CPPFLAGS_SL+=$(call and2,$($(cmod_)_DEPRECATED),-DSL_DEPRECATED)
$(cmod_)_CPPFLAGS_SL+=$(call and2,$($(cmod_)_DEVEL),-DSL_DEVEL)


# target files
$(cmod_)_LIBSL_A_SEQ:=$(call or2,$($(cmod_)_LIBSL_A_SEQ),$(LIBSL_A_SEQ))
$(cmod_)_LIBSL_A_MPI:=$(call or2,$($(cmod_)_LIBSL_A_MPI),$(LIBSL_A_MPI))
$(cmod_)_LIBSL_A:=$(call or2,$($(cmod_)_LIBSL_A),$(LIBSL_A))

$(cmod_)_EXEC:=$(call or2,$($(cmod_)_EXEC),$(EXEC))

# FIXME: why cmod and not cmod_?
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
ifneq ($($(cmod_)_DEPRECATED),)
 $(cmod_)_CODE_SL_CORE_SEQ_INCLUDE+=$(wildcard $($(cmod_)_SRC_CORE_SEQ)deprecated/*.h)
 $(cmod_)_CODE_SL_CORE_MPI_INCLUDE+=$(wildcard $($(cmod_)_SRC_CORE_MPI)deprecated/*.h)
endif
ifneq ($($(cmod_)_DEVEL),)
 $(cmod_)_CODE_SL_CORE_SEQ_INCLUDE+=$(wildcard $($(cmod_)_SRC_CORE_SEQ)devel/*.h)
 $(cmod_)_CODE_SL_CORE_MPI_INCLUDE+=$(wildcard $($(cmod_)_SRC_CORE_MPI)devel/*.h)
endif

# filter-out existing _di-versions (to get a clean list of core-files)
$(cmod_)_CODE_SL_CORE_SEQ:=$(filter-out $(wildcard $($(cmod_)_SRC_CORE_SEQ)*_di.c),$(wildcard $($(cmod_)_SRC_CORE_SEQ)*.c))
$(cmod_)_CODE_SL_CORE_MPI:=$(filter-out $(wildcard $($(cmod_)_SRC_CORE_MPI)*_di.c),$(wildcard $($(cmod_)_SRC_CORE_MPI)*.c))
ifneq ($($(cmod_)_DEPRECATED),)
 $(cmod_)_CODE_SL_CORE_SEQ+=$(filter-out $(wildcard $($(cmod_)_SRC_CORE_SEQ)deprecated/*_di.c),$(wildcard $($(cmod_)_SRC_CORE_SEQ)deprecated/*.c))
 $(cmod_)_CODE_SL_CORE_MPI+=$(filter-out $(wildcard $($(cmod_)_SRC_CORE_MPI)deprecated/*_di.c),$(wildcard $($(cmod_)_SRC_CORE_MPI)deprecated/*.c))
endif
ifneq ($($(cmod_)_DEVEL),)
 $(cmod_)_CODE_SL_CORE_SEQ+=$(filter-out $(wildcard $($(cmod_)_SRC_CORE_SEQ)devel/*_di.c),$(wildcard $($(cmod_)_SRC_CORE_SEQ)devel/*.c))
 $(cmod_)_CODE_SL_CORE_MPI+=$(filter-out $(wildcard $($(cmod_)_SRC_CORE_MPI)devel/*_di.c),$(wildcard $($(cmod_)_SRC_CORE_MPI)devel/*.c))
endif
# add _di-versions if not disabled
ifneq ($($(cmod_)_DI_VERSION),false)
 $(cmod_)_CODE_SL_CORE_SEQ+=$(addsuffix _di.c,$(basename $($(cmod_)_CODE_SL_CORE_SEQ)))
 $(cmod_)_CODE_SL_CORE_MPI+=$(addsuffix _di.c,$(basename $($(cmod_)_CODE_SL_CORE_MPI)))
endif

ifneq ($(cmod_cfg),)
 $(cmod_)_OBJS_SL_CORE_SEQ:=$(addprefix $(cmod_prefix), $(addsuffix .o, $(basename $(notdir $($(cmod_)_CODE_SL_CORE_SEQ)))))
 $(cmod_)_OBJS_SL_CORE_MPI:=$(addprefix $(cmod_prefix), $(addsuffix .o, $(basename $(notdir $($(cmod_)_CODE_SL_CORE_MPI)))))
endif

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

ifneq ($(cmod_cfg),)
 $(cmod_)_OBJS_SL_TUNE_SEQ:=$(addsuffix .o, $(basename $(notdir $($(cmod_)_CODE_SL_TUNE_SEQ))))
 $(cmod_)_OBJS_SL_TUNE_MPI:=$(addsuffix .o, $(basename $(notdir $($(cmod_)_CODE_SL_TUNE_MPI))))
endif

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


# exec files
$(cmod_)_SRC_EXEC_SEQ:=$(cmod_dir)$(EXEC_SRC)
$(cmod_)_SRC_EXEC_MPI:=$(cmod_dir)$(EXEC_MPI_SRC)
$(cmod_)_TMP_EXEC_SEQ:=$(cmod_dir)$(EXEC_SRC)
$(cmod_)_TMP_EXEC_MPI:=$(cmod_dir)$(EXEC_MPI_SRC)

$(cmod_)_CODE_SL_EXEC_MPI_INCLUDE:=$(wildcard $($(cmod_)_SRC_EXEC_MPI)$(EXEC_MPI_PREFIX)*.h)
$(cmod_)_CODE_SL_EXEC_SEQ_INCLUDE:=$(filter-out $($(cmod_)_CODE_SL_EXEC_MPI_INCLUDE),$(wildcard $($(cmod_)_SRC_EXEC_SEQ)$(EXEC_PREFIX)*.h))
$(cmod_)_CODE_SL_EXEC_INCLUDE:=$($(cmod_)_CODE_SL_EXEC_SEQ_INCLUDE) $($(cmod_)_CODE_SL_EXEC_MPI_INCLUDE)

$(cmod_)_CODE_SL_EXEC_MPI:=$(wildcard $($(cmod_)_SRC_EXEC_MPI)$(EXEC_MPI_PREFIX)*.c)
$(cmod_)_CODE_SL_EXEC_SEQ:=$(filter-out $($(cmod_)_CODE_SL_EXEC_MPI),$(wildcard $($(cmod_)_SRC_EXEC_SEQ)$(EXEC_PREFIX)*.c))
$(cmod_)_CODE_SL_EXEC:=$($(cmod_)_CODE_SL_EXEC_SEQ) $($(cmod_)_CODE_SL_EXEC_MPI)

$(cmod_)_OBJS_SL_EXEC_SEQ:=$(addsuffix .o, $(basename $(notdir $($(cmod_)_CODE_SL_EXEC_SEQ))))
$(cmod_)_OBJS_SL_EXEC_MPI:=$(addsuffix .o, $(basename $(notdir $($(cmod_)_CODE_SL_EXEC_MPI))))

$(cmod_)_TMP_OBJS_SL_EXEC_SEQ:=$(addprefix $($(cmod_)_TMP_EXEC_SEQ),$($(cmod_)_OBJS_SL_EXEC_SEQ))
$(cmod_)_TMP_OBJS_SL_EXEC_MPI:=$(addprefix $($(cmod_)_TMP_EXEC_MPI),$($(cmod_)_OBJS_SL_EXEC_MPI))


# wrapper files
$(cmod_)_SRC_WRAP_SEQ:=$(cmod_dir)
$(cmod_)_SRC_WRAP_MPI:=$(cmod_dir)
$(cmod_)_TMP_WRAP_SEQ:=$(cmod_dir)
$(cmod_)_TMP_WRAP_MPI:=$(cmod_dir)

$(cmod_)_CODE_SL_WRAP_MPI_INCLUDE:=$(filter-out $($(cmod_)_CODE_SL_EXEC_INCLUDE),$(wildcard $($(cmod_)_SRC_WRAP_MPI)$(WRAPPER_MPI_PREFIX)*.h))
$(cmod_)_CODE_SL_WRAP_SEQ_INCLUDE:=$(filter-out $($(cmod_)_CODE_SL_EXEC_INCLUDE) $($(cmod_)_CODE_SL_WRAP_MPI_INCLUDE),$(wildcard $($(cmod_)_SRC_WRAP_SEQ)$(WRAPPER_PREFIX)*.h))
$(cmod_)_CODE_SL_WRAP_INCLUDE:=$($(cmod_)_CODE_SL_WRAP_SEQ_INCLUDE) $($(cmod_)_CODE_SL_WRAP_MPI_INCLUDE)

$(cmod_)_CODE_SL_WRAP_MPI:=$(filter-out $($(cmod_)_CODE_SL_EXEC),$(wildcard $($(cmod_)_SRC_WRAP_MPI)$(WRAPPER_MPI_PREFIX)*.c))
$(cmod_)_CODE_SL_WRAP_SEQ:=$(filter-out $($(cmod_)_CODE_SL_EXEC) $($(cmod_)_CODE_SL_WRAP_MPI),$(wildcard $($(cmod_)_SRC_WRAP_SEQ)$(WRAPPER_PREFIX)*.c))

$(cmod_)_OBJS_SL_WRAP_SEQ:=$(addsuffix .o, $(basename $(notdir $($(cmod_)_CODE_SL_WRAP_SEQ))))
$(cmod_)_OBJS_SL_WRAP_MPI:=$(addsuffix .o, $(basename $(notdir $($(cmod_)_CODE_SL_WRAP_MPI))))

$(cmod_)_TMP_OBJS_SL_WRAP_SEQ:=$(addprefix $($(cmod_)_TMP_WRAP_SEQ),$($(cmod_)_OBJS_SL_WRAP_SEQ))
$(cmod_)_TMP_OBJS_SL_WRAP_MPI:=$(addprefix $($(cmod_)_TMP_WRAP_MPI),$($(cmod_)_OBJS_SL_WRAP_MPI))


# find-and-replace scripts
$(cmod_)_NAMES_OXL_ALL_FUNCTIONS:=$($(cmod_)_SRCDIR)names_oxl_all_functions
$(cmod_)_FAR_SCRIPT_DI:=$($(cmod_)_TMPDIR)$(cmod_prefix)far_script_di
$(cmod_)_FIND_AND_REPLACE:=./$(cmod_dir)scripts/find_and_replace.sh "$(FAR)"

# list of required object-files
$(cmod_)_OBJS_SL_CORE:=$($(cmod_)_OBJS_SL_CORE_SEQ)
$(cmod_)_OBJS_SL_TUNE:=$($(cmod_)_OBJS_SL_TUNE_SEQ)
$(cmod_)_OBJS_SL_DEVL:=$($(cmod_)_OBJS_SL_DEVL_SEQ)
$(cmod_)_OBJS_SL_DEMO:=$($(cmod_)_OBJS_SL_DEMP_SEQ)
$(cmod_)_OBJS_SL_WRAP:=$($(cmod_)_OBJS_SL_WRAP_SEQ)
$(cmod_)_OBJS_SL_EXEC:=$($(cmod_)_OBJS_SL_EXEC_SEQ)

$(cmod_)_TMP_OBJS_SL_CORE:=$($(cmod_)_TMP_OBJS_SL_CORE_SEQ)
$(cmod_)_TMP_OBJS_SL_TUNE:=$($(cmod_)_TMP_OBJS_SL_TUNE_SEQ)
$(cmod_)_TMP_OBJS_SL_DEVL:=$($(cmod_)_TMP_OBJS_SL_DEVL_SEQ)
$(cmod_)_TMP_OBJS_SL_DEMO:=$($(cmod_)_TMP_OBJS_SL_DEMO_SEQ)
$(cmod_)_TMP_OBJS_SL_WRAP:=$($(cmod_)_TMP_OBJS_SL_WRAP_SEQ)
$(cmod_)_TMP_OBJS_SL_EXEC:=$($(cmod_)_TMP_OBJS_SL_EXEC_SEQ)

ifeq ($($(cmod_)_SL_USE_MPI),true)
 $(cmod_)_OBJS_SL_CORE+=$($(cmod_)_OBJS_SL_CORE_MPI)
 $(cmod_)_OBJS_SL_TUNE+=$($(cmod_)_OBJS_SL_TUNE_MPI)
 $(cmod_)_OBJS_SL_DEVL:=$($(cmod_)_OBJS_SL_DEVL_MPI)
 $(cmod_)_OBJS_SL_DEMO:=$($(cmod_)_OBJS_SL_DEMO_MPI)
 $(cmod_)_OBJS_SL_WRAP+=$($(cmod_)_OBJS_SL_WRAP_MPI)
 $(cmod_)_OBJS_SL_EXEC+=$($(cmod_)_OBJS_SL_EXEC_MPI)

 $(cmod_)_TMP_OBJS_SL_CORE+=$($(cmod_)_TMP_OBJS_SL_CORE_MPI)
 $(cmod_)_TMP_OBJS_SL_TUNE+=$($(cmod_)_TMP_OBJS_SL_TUNE_MPI)
 $(cmod_)_TMP_OBJS_SL_DEVL:=$($(cmod_)_TMP_OBJS_SL_DEVL_MPI)
 $(cmod_)_TMP_OBJS_SL_DEMO:=$($(cmod_)_TMP_OBJS_SL_DEMO_MPI)
 $(cmod_)_TMP_OBJS_SL_WRAP+=$($(cmod_)_TMP_OBJS_SL_WRAP_MPI)
 $(cmod_)_TMP_OBJS_SL_EXEC+=$($(cmod_)_TMP_OBJS_SL_EXEC_MPI)
else
 ifeq ($(M_SL_USE_MPI),true)
  $(cmod_)_OBJS_SL_WRAP+=$($(cmod_)_OBJS_SL_WRAP_MPI)
  $(cmod_)_OBJS_SL_EXEC+=$($(cmod_)_OBJS_SL_EXEC_MPI)

  $(cmod_)_TMP_OBJS_SL_WRAP+=$($(cmod_)_TMP_OBJS_SL_WRAP_MPI)
  $(cmod_)_TMP_OBJS_SL_EXEC+=$($(cmod_)_TMP_OBJS_SL_EXEC_MPI)
 endif
endif


# final CC, LD and flags
$(cmod_)_CFLAGS_:=$($(cmod_)_CFLAGS) $($(cmod_)_CPPFLAGS)
$(cmod_)_CFLAGS_SL_:=$($(cmod_)_CFLAGS_SL) -DSL_USE_MPI_IGNORE $($(cmod_)_CPPFLAGS_SL)
$(cmod_)_CFLAGS_NSL_:=$($(cmod_)_CFLAGS_NSL) $($(cmod_)_CPPFLAGS_NSL)
$(cmod_)_MPICFLAGS_:=$($(cmod_)_CFLAGS) $($(cmod_)_MPICFLAGS) $($(cmod_)_CPPFLAGS) $($(cmod_)_MPICPPFLAGS)
$(cmod_)_MPICFLAGS_SL_:=$($(cmod_)_CFLAGS_SL) $($(cmod_)_MPICFLAGS_SL) -DSL_USE_MPI_FORCE $($(cmod_)_CPPFLAGS_SL) $($(cmod_)_MPICPPFLAGS_SL)
$(cmod_)_MPICFLAGS_NSL_:=$($(cmod_)_CFLAGS_NSL) $($(cmod_)_MPICFLAGS_NSL) $($(cmod_)_CPPFLAGS_NSL) $($(cmod_)_MPICPPFLAGS_NSL)

ifeq ($($(cmod_)_ALL_AS_MPI),true)
 $(cmod_)_CC:=$($(cmod_)_MPICC)
 $(cmod_)_CFLAGS_:=$($(cmod_)_MPICFLAGS_)
 $(cmod_)_CFLAGS_SL_:=$($(cmod_)_MPICFLAGS_SL_)
 $(cmod_)_CFLAGS_NSL_:=$($(cmod_)_MPICFLAGS_NSL_)
endif

$(cmod_)_LDFLAGS_:=$($(cmod_)_LDFLAGS)
$(cmod_)_MPILDFLAGS_:=$($(cmod_)_LDFLAGS) $($(cmod_)_MPILDFLAGS)

$(cmod_)_LDLIBS_:=$($(cmod_)_LDLIBS)
$(cmod_)_MPILDLIBS_:=$($(cmod_)_LDLIBS) $($(cmod_)_MPILDLIBS)

ifeq ($($(cmod_)_SL_USE_MPI),true)
 $(cmod_)_LD:=$($(cmod_)_MPILD)
 $(cmod_)_LDFLAGS_:=$($(cmod_)_MPILDFLAGS_)
 $(cmod_)_LDLIBS_:=$($(cmod_)_MPILDLIBS_)
endif


# prerequisites
$(cmod_)_PREREQ:=$(PREREQ) $($(cmod_)_PREREQ)


# commands
ifeq ($(QUIET),)
define $(cmod_)_comp
@echo "compiling $(TYPE) '$<' to '$@'"
-@$(MKDIR_P) $(dir $@)
$(CC) -c $(CFLAGS) -o $@ $<
endef
else
define $(cmod_)_comp
-@$(MKDIR_P) $(dir $@)
$(CC) -c $(CFLAGS) -o $@ $<
endef
endif

#@$(ECHO) "adding $(TYPE) '$^' to '$@'"
define $(cmod_)_add
-@$(MKDIR_P) $(dir $@)
$(AR) $(ARFLAGS) $@ $^
endef

define $(cmod_)_make
@$(ECHO) "making $(TYPE) '$@'"
-@$(MKDIR_P) $(dir $@)
$(LD) $(LDFLAGS) -o $@ $^ $(LDLIBS)
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
	@$(ECHO) " CFLAGS_NSL: $($(cmod_)_CFLAGS_NSL)"
	@$(ECHO) " CPPFLAGS_NSL: $($(cmod_)_CPPFLAGS_NSL)"
	@$(ECHO) " CFLAGS_: $($(cmod_)_CFLAGS_)"
	@$(ECHO) " CFLAGS_SL_: $($(cmod_)_CFLAGS_SL_)"
	@$(ECHO) " CFLAGS_NSL_: $($(cmod_)_CFLAGS_NSL_)"
	@$(ECHO) " MPICC: $($(cmod_)_MPICC)"
	@$(ECHO) " MPICFLAGS: $($(cmod_)_MPICFLAGS)"
	@$(ECHO) " MPICPPFLAGS: $($(cmod_)_MPICPPFLAGS)"
	@$(ECHO) " MPICFLAGS_SL: $($(cmod_)_MPICFLAGS_SL)"
	@$(ECHO) " MPICPPFLAGS_SL: $($(cmod_)_MPICPPFLAGS_SL)"
	@$(ECHO) " MPICFLAGS_NSL: $($(cmod_)_MPICFLAGS_NSL)"
	@$(ECHO) " MPICPPFLAGS_NSL: $($(cmod_)_MPICPPFLAGS_NSL)"
	@$(ECHO) " MPICFLAGS_: $($(cmod_)_MPICFLAGS_)"
	@$(ECHO) " MPICFLAGS_SL_: $($(cmod_)_MPICFLAGS_SL_)"
	@$(ECHO) " MPICFLAGS_NSL_: $($(cmod_)_MPICFLAGS_NSL_)"
	@$(ECHO) " LD: $($(cmod_)_LD)"
	@$(ECHO) " LDFLAGS_: $($(cmod_)_LDFLAGS_)"
	@$(ECHO) " LDLIBS_: $($(cmod_)_LDLIBS_)"
	@$(ECHO) " MPILD: $($(cmod_)_MPILD)"
	@$(ECHO) " MPILDFLAGS_: $($(cmod_)_MPILDFLAGS_)"
	@$(ECHO) " MPILDLIBS_: $($(cmod_)_MPILDLIBS_)"
	@$(ECHO) " PREREQ: $($(cmod_)_PREREQ)"

	@$(ECHO) " SRCDIR: $($(cmod_)_SRCDIR)"
	@$(ECHO) " CFGDIR: $($(cmod_)_CFGDIR)"
	@$(ECHO) " TMPDIR: $($(cmod_)_TMPDIR)"
	@$(ECHO) " LIBDIR: $($(cmod_)_LIBDIR)"
	@$(ECHO) " INTDIR: $($(cmod_)_INTDIR)"

	@$(ECHO) " include-dir: $($(cmod_)_INCLUDE)"
	@$(ECHO) " sl-include-dir: $($(cmod_)_SL_INCLUDE)"
	@$(ECHO) " sl-cfg-include-dir: $($(cmod_)_SL_CFG_INCLUDE)"

	@$(ECHO) " src-dir: $($(cmod_)_SRCDIR)"
	@$(ECHO) " src-core-dir: $($(cmod_)_SRC_CORE_SEQ)"
	@$(ECHO) " src-core-mpi-dir: $($(cmod_)_SRC_CORE_MPI)"
	@$(ECHO) " src-tune-dir: $($(cmod_)_SRC_TUNE_SEQ)"
	@$(ECHO) " src-tune-mpi-dir: $($(cmod_)_SRC_TUNE_MPI)"
	@$(ECHO) " src-devel-dir: $($(cmod_)_SRC_DEVL_SEQ)"
	@$(ECHO) " src-devel-mpi-dir: $($(cmod_)_SRC_DEVL_MPI)"
	@$(ECHO) " src-demo-dir: $($(cmod_)_SRC_DEMO_SEQ)"
	@$(ECHO) " src-demo-mpi-dir: $($(cmod_)_SRC_DEMO_MPI)"
	@$(ECHO) " src-wrapper-dir: $($(cmod_)_SRC_WRAP_SEQ)"
	@$(ECHO) " src-wrapper-mpi-dir: $($(cmod_)_SRC_WRAP_MPI)"
	@$(ECHO) " src-exec-dir: $($(cmod_)_SRC_EXEC_SEQ)"
	@$(ECHO) " src-exec-mpi-dir: $($(cmod_)_SRC_EXEC_MPI)"

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
	@$(ECHO) " code-sl-cfg-include: $($(cmod_)_CODE_SL_CFG_INCLUDE)"

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

	@$(ECHO) " code-wrapper-include-seq: $($(cmod_)_CODE_SL_WRAP_SEQ_INCLUDE)"
	@$(ECHO) " code-wrapper-include-mpi: $($(cmod_)_CODE_SL_WRAP_MPI_INCLUDE)"
	@$(ECHO) " code-wrapper-include: $($(cmod_)_CODE_SL_WRAP_INCLUDE)"
	@$(ECHO) " code-wrapper-seq: $($(cmod_)_CODE_SL_WRAP_SEQ)"
	@$(ECHO) " code-wrapper-mpi: $($(cmod_)_CODE_SL_WRAP_MPI)"
	@$(ECHO) " objs-wrapper-seq: $($(cmod_)_OBJS_SL_WRAP_SEQ)"
	@$(ECHO) " objs-wrapper-mpi: $($(cmod_)_OBJS_SL_WRAP_MPI)"
	@$(ECHO) " tmp-objs-wrapper-seq: $($(cmod_)_TMP_OBJS_SL_WRAP_SEQ)"
	@$(ECHO) " tmp-objs-wrapper-mpi: $($(cmod_)_TMP_OBJS_SL_WRAP_MPI)"

	@$(ECHO) " code-exec-include-seq: $($(cmod_)_CODE_SL_EXEC_SEQ_INCLUDE)"
	@$(ECHO) " code-exec-include-mpi: $($(cmod_)_CODE_SL_EXEC_MPI_INCLUDE)"
	@$(ECHO) " code-exec-include: $($(cmod_)_CODE_SL_EXEC_INCLUDE)"
	@$(ECHO) " code-exec-seq: $($(cmod_)_CODE_SL_EXEC_SEQ)"
	@$(ECHO) " code-exec-mpi: $($(cmod_)_CODE_SL_EXEC_MPI)"
	@$(ECHO) " objs-exec-seq: $($(cmod_)_OBJS_SL_EXEC_SEQ)"
	@$(ECHO) " objs-exec-mpi: $($(cmod_)_OBJS_SL_EXEC_MPI)"
	@$(ECHO) " tmp-objs-exec-seq: $($(cmod_)_TMP_OBJS_SL_EXEC_SEQ)"
	@$(ECHO) " tmp-objs-exec-mpi: $($(cmod_)_TMP_OBJS_SL_EXEC_MPI)"

	@$(ECHO) " libsl_seq: $($(cmod_)_LIBSL_A_SEQ)"
	@$(ECHO) " libsl_seq(*): $($(cmod_)_LIBSL_A_SEQ)($($(cmod_)_OBJS_SL_CORE_SEQ))"
	@$(ECHO) " libsl_mpi: $($(cmod_)_LIBSL_A_MPI)"
	@$(ECHO) " libsl_mpi(*): $($(cmod_)_LIBSL_A_MPI)($($(cmod_)_OBJS_SL_CORE_MPI))"
	@$(ECHO) " libsl: $($(cmod_)_LIBSL_A)"
	@$(ECHO) " libsl(*): $($(cmod_)_LIBSL_A)($($(cmod_)_OBJS_SL_CORE))"
	@$(ECHO) " tune: $($(cmod_)_TUNE)"
	@$(ECHO) " devel: $($(cmod_)_DEVL)"
	@$(ECHO) " demo: $($(cmod_)_DEMO)"

	@$(ECHO) " NAMES_OXL_ALL_FUNCTIONS: $($(cmod_)_NAMES_OXL_ALL_FUNCTIONS)"
	@$(ECHO) " FAR_SCRIPT_DI: $($(cmod_)_FAR_SCRIPT_DI)"
	@$(ECHO) " FIND_AND_REPLACE: $($(cmod_)_FIND_AND_REPLACE)"

M_INFO+=$(cmod_)_info


# core
$(cmod_)_core_seq:
$(cmod_)_core_mpi:
$(cmod_)_core:

ifneq ($(cmod_cfg),)
$(cmod_)_core_seq: $($(cmod_)_TMP_OBJS_SL_CORE_SEQ)
$(cmod_)_core_mpi: $($(cmod_)_TMP_OBJS_SL_CORE_MPI)
$(cmod_)_core: $($(cmod_)_TMP_OBJS_SL_CORE)
endif

M_CORE_SEQ+=$(cmod_)_core_seq
M_CORE_MPI+=$(cmod_)_core_mpi
M_CORE+=$(cmod_)_core

# core files
$($(cmod_)_TMP_OBJS_SL_CORE_SEQ): override CC:=$($(cmod_)_CC)
$($(cmod_)_TMP_OBJS_SL_CORE_SEQ): override CFLAGS:=$(strip $($(cmod_)_CFLAGS_) $($(cmod_)_CFLAGS_SL_))
$($(cmod_)_TMP_OBJS_SL_CORE_SEQ): TYPE:=core

$($(cmod_)_TMP_OBJS_SL_CORE_MPI): override CC:=$($(cmod_)_MPICC)
$($(cmod_)_TMP_OBJS_SL_CORE_MPI): override CFLAGS:=$(strip $($(cmod_)_MPICFLAGS_) $($(cmod_)_MPICFLAGS_SL_))
$($(cmod_)_TMP_OBJS_SL_CORE_MPI): TYPE:=core-mpi

$($(cmod_)_TMP_CORE_SEQ)$(cmod_prefix)%.o: $($(cmod_)_SRC_CORE_SEQ)%.c $($(cmod_)_CODE_INCLUDE) $($(cmod_)_CODE_SL_CFG_INCLUDE) $($(cmod_)_CODE_SL_INCLUDE) $($(cmod_)_PREREQ)
	$($(cmod_)_comp)

$($(cmod_)_TMP_CORE_MPI)$(cmod_prefix)%.o: $($(cmod_)_SRC_CORE_MPI)%.c $($(cmod_)_CODE_INCLUDE) $($(cmod_)_CODE_SL_CFG_INCLUDE) $($(cmod_)_CODE_SL_INCLUDE) $($(cmod_)_PREREQ)
	$($(cmod_)_comp)

ifneq ($($(cmod_)_DEPRECATED),)
$($(cmod_)_TMP_CORE_SEQ)$(cmod_prefix)%.o: $($(cmod_)_SRC_CORE_SEQ)deprecated/%.c $($(cmod_)_CODE_INCLUDE) $($(cmod_)_CODE_SL_CFG_INCLUDE) $($(cmod_)_CODE_SL_INCLUDE) $($(cmod_)_PREREQ)
	$($(cmod_)_comp)

$($(cmod_)_TMP_CORE_MPI)$(cmod_prefix)%.o: $($(cmod_)_SRC_CORE_MPI)deprecated/%.c $($(cmod_)_CODE_INCLUDE) $($(cmod_)_CODE_SL_CFG_INCLUDE) $($(cmod_)_CODE_SL_INCLUDE) $($(cmod_)_PREREQ)
	$($(cmod_)_comp)
endif

ifneq ($($(cmod_)_DEVEL),)
$($(cmod_)_TMP_CORE_SEQ)$(cmod_prefix)%.o: $($(cmod_)_SRC_CORE_SEQ)devel/%.c $($(cmod_)_CODE_INCLUDE) $($(cmod_)_CODE_SL_CFG_INCLUDE) $($(cmod_)_CODE_SL_INCLUDE) $($(cmod_)_PREREQ)
	$($(cmod_)_comp)

$($(cmod_)_TMP_CORE_MPI)$(cmod_prefix)%.o: $($(cmod_)_SRC_CORE_MPI)devel/%.c $($(cmod_)_CODE_INCLUDE) $($(cmod_)_CODE_SL_CFG_INCLUDE) $($(cmod_)_CODE_SL_INCLUDE) $($(cmod_)_PREREQ)
	$($(cmod_)_comp)
endif

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
#$(cmod_)_tune_seq:
#$(cmod_)_tune_mpi:
$(cmod_)_tune:

ifneq ($(cmod_cfg),)
#$(cmod_)_tune_seq: $(cmod_)_libsl_seq $($(cmod_)_TUNE_SEQ)
#$(cmod_)_tune_mpi: $(cmod_)_libsl_mpi $($(cmod_)_TUNE_MPI)
$(cmod_)_tune: $(cmod_)_libsl $($(cmod_)_TUNE)
endif

#M_TUNE_SEQ+=$(cmod_)_tune_seq
#M_TUNE_MPI+=$(cmod_)_tune_mpi
M_TUNE+=$(cmod_)_tune

# tune files
$($(cmod_)_TMP_OBJS_SL_TUNE_SEQ): override CC:=$($(cmod_)_CC)
$($(cmod_)_TMP_OBJS_SL_TUNE_SEQ): override CFLAGS:=$(strip $($(cmod_)_CFLAGS_) $($(cmod_)_CFLAGS_SL_) $(call or2,$(findstring -DSL_TUNEABLE,$($(cmod_)_CFLAGS_SL_)),-DSL_TUNEABLE))
$($(cmod_)_TMP_OBJS_SL_TUNE_SEQ): TYPE=tune

$($(cmod_)_TMP_OBJS_SL_TUNE_MPI): CC:=$($(cmod_)_MPICC)
$($(cmod_)_TMP_OBJS_SL_TUNE_MPI): CFLAGS:=$(strip $($(cmod_)_MPICFLAGS_) $($(cmod_)_MPICFLAGS_SL_) $(call or2,$(findstring -DSL_TUNEABLE,$($(cmod_)_MPICFLAGS_SL_)),-DSL_TUNEABLE))
$($(cmod_)_TMP_OBJS_SL_TUNE_MPI): TYPE=tune-mpi

$($(cmod_)_TMP_TUNE_SEQ)%.o: $($(cmod_)_SRC_TUNE_SEQ)%.c $($(cmod_)_SRC_TUNE_SEQ_INCLUDE) $($(cmod_)_CODE_SL_INCLUDE) $($(cmod_)_PREREQ)
	$($(cmod_)_comp)

$($(cmod_)_TMP_TUNE_MPI)%.o: $($(cmod_)_SRC_TUNE_MPI)%.c $($(cmod_)_SRC_TUNE_MPI_INCLUDE) $($(cmod_)_CODE_SL_INCLUDE) $($(cmod_)_PREREQ)
	$($(cmod_)_comp)

$($(cmod_)_TUNE): override LD:=$($(cmod_)_LD)
$($(cmod_)_TUNE): override LDFLAGS:=$(strip $($(cmod_)_LDFLAGS_))
$($(cmod_)_TUNE): TYPE:=tune

$($(cmod_)_TUNE): $($(cmod_)_TMP_OBJS_SL_TUNE) $($(cmod_)_LIBSL_A)
	$($(cmod_)_make)


# devel
#$(cmod_)_devel_seq:
#$(cmod_)_devel_mpi:
$(cmod_)_devel:

ifneq ($(cmod_cfg),)
#$(cmod_)_devel_seq: $(cmod_)_libsl_seq $($(cmod_)_DEVL_SEQ)
#$(cmod_)_devel_mpi: $(cmod_)_libsl_mpi $($(cmod_)_DEVL_MPI)
$(cmod_)_devel: $(cmod_)_libsl $($(cmod_)_DEVL)
endif

#M_DEVL_SEQ+=$(cmod_)_devel_seq
#M_DEVL_MPI+=$(cmod_)_devel_mpi
M_DEVL+=$(cmod_)_devel

# devel files
$($(cmod_)_TMP_OBJS_SL_DEVL_SEQ): override CC:=$($(cmod_)_CC)
$($(cmod_)_TMP_OBJS_SL_DEVL_SEQ): override CFLAGS:=$(strip $($(cmod_)_CFLAGS_) $($(cmod_)_CFLAGS_SL_))
$($(cmod_)_TMP_OBJS_SL_DEVL_SEQ): TYPE=devel

$($(cmod_)_TMP_OBJS_SL_DEVL_MPI): override CC:=$($(cmod_)_MPICC)
$($(cmod_)_TMP_OBJS_SL_DEVL_MPI): override CFLAGS:=$(strip $($(cmod_)_MPICFLAGS_) $($(cmod_)_MPICFLAGS_SL_))
$($(cmod_)_TMP_OBJS_SL_DEVL_MPI): TYPE=devel-mpi

$($(cmod_)_TMP_DEVL_SEQ)%.o: $($(cmod_)_SRC_DEVL_SEQ)%.c $($(cmod_)_SRC_DEVL_SEQ_INCLUDE) $($(cmod_)_CODE_SL_INCLUDE) $($(cmod_)_PREREQ)
	$($(cmod_)_comp)

$($(cmod_)_TMP_DEVL_MPI)%.o: $($(cmod_)_SRC_DEVL_MPI)%.c $($(cmod_)_SRC_DEVL_MPI_INCLUDE) $($(cmod_)_CODE_SL_INCLUDE) $($(cmod_)_PREREQ)
	$($(cmod_)_comp)

$($(cmod_)_DEVL): override LD:=$($(cmod_)_LD)
$($(cmod_)_DEVL): override LDFLAGS:=$(strip $($(cmod_)_LDFLAGS_))
$($(cmod_)_DEVL): TYPE:=devel

$($(cmod_)_DEVL): $($(cmod_)_TMP_OBJS_SL_DEVL) $($(cmod_)_LIBSL_A)
	$($(cmod_)_make)


# demo
#$(cmod_)_demo_seq:
#$(cmod_)_demo_mpi:
$(cmod_)_demo:

ifneq ($(cmod_cfg),)
#$(cmod_)_demo_seq: $(cmod_)_libsl_seq $($(cmod_)_DEMO_SEQ)
#$(cmod_)_demo_mpi: $(cmod_)_libsl_mpi $($(cmod_)_DEMO_MPI)
$(cmod_)_demo: $(cmod_)_libsl $($(cmod_)_DEMO)
endif

#M_DEMO_SEQ+=$(cmod_)_demo_seq
#M_DEMO_MPI+=$(cmod_)_demo_mpi
M_DEMO+=$(cmod_)_demo

# demo files
$($(cmod_)_TMP_OBJS_SL_DEMO_SEQ): override CC:=$($(cmod_)_CC)
$($(cmod_)_TMP_OBJS_SL_DEMO_SEQ): override CFLAGS:=$(strip $($(cmod_)_CFLAGS_) $($(cmod_)_CFLAGS_SL_))
$($(cmod_)_TMP_OBJS_SL_DEMO_SEQ): TYPE=demo

$($(cmod_)_TMP_OBJS_SL_DEMO_MPI): override CC:=$($(cmod_)_MPICC)
$($(cmod_)_TMP_OBJS_SL_DEMO_MPI): override CFLAGS:=$(strip $($(cmod_)_MPICFLAGS_) $($(cmod_)_MPICFLAGS_SL_))
$($(cmod_)_TMP_OBJS_SL_DEMO_MPI): TYPE=demo-mpi

$($(cmod_)_TMP_DEMO_SEQ)%.o: $($(cmod_)_SRC_DEMO_SEQ)%.c $($(cmod_)_SRC_DEMO_SEQ_INCLUDE) $($(cmod_)_CODE_SL_INCLUDE) $($(cmod_)_PREREQ)
	$($(cmod_)_comp)

$($(cmod_)_TMP_DEMO_MPI)%.o: $($(cmod_)_SRC_DEMO_MPI)%.c $($(cmod_)_SRC_DEMO_MPI_INCLUDE) $($(cmod_)_CODE_SL_INCLUDE) $($(cmod_)_PREREQ)
	$($(cmod_)_comp)

$($(cmod_)_DEMO): override LD:=$($(cmod_)_LD)
$($(cmod_)_DEMO): override LDFLAGS:=$(strip $($(cmod_)_LDFLAGS_))
$($(cmod_)_DEMO): TYPE:=demo

$($(cmod_)_DEMO): $($(cmod_)_TMP_OBJS_SL_DEMO) $($(cmod_)_LIBSL_A)
	$($(cmod_)_make)


# wrapper
$(cmod_)_wrapper_seq: $($(cmod_)_TMP_OBJS_SL_WRAP_SEQ)
$(cmod_)_wrapper_mpi: $($(cmod_)_TMP_OBJS_SL_WRAP_MPI)
$(cmod_)_wrapper: $($(cmod_)_TMP_OBJS_SL_WRAP)

M_WRAP_SEQ+=$(cmod_)_wrapper_seq
M_WRAP_MPI+=$(cmod_)_wrapper_mpi
M_WRAP+=$(cmod_)_wrapper

# wrapper files
$($(cmod_)_TMP_OBJS_SL_WRAP_SEQ): override CC:=$($(cmod_)_CC)
$($(cmod_)_TMP_OBJS_SL_WRAP_SEQ): override CFLAGS:=$(strip $($(cmod_)_CFLAGS_) $($(cmod_)_CFLAGS_NSL_))
$($(cmod_)_TMP_OBJS_SL_WRAP_SEQ): TYPE:=wrapper
$($(cmod_)_TMP_OBJS_SL_WRAP_SEQ): $($(cmod_)_CODE_SL_WRAP_SEQ_INCLUDE) $(M_CODE_INCLUDE) $($(cmod_)_PREREQ)

$($(cmod_)_TMP_OBJS_SL_WRAP_MPI): override CC:=$($(cmod_)_MPICC)
$($(cmod_)_TMP_OBJS_SL_WRAP_MPI): override CFLAGS:=$(strip $($(cmod_)_MPICFLAGS_) $($(cmod_)_MPICFLAGS_NSL_))
$($(cmod_)_TMP_OBJS_SL_WRAP_MPI): TYPE:=wrapper-mpi
$($(cmod_)_TMP_OBJS_SL_WRAP_MPI): $($(cmod_)_CODE_SL_WRAP_MPI_INCLUDE) $(M_CODE_INCLUDE) $($(cmod_)_PREREQ)

$($(cmod_)_TMP_WRAP_SEQ)%.o: $($(cmod_)_SRC_WRAP_SEQ)%.c
	$($(cmod_)_comp)

$($(cmod_)_TMP_WRAP_MPI)%.o: $($(cmod_)_SRC_WRAP_MPI)%.c
	$($(cmod_)_comp)


# exec files
M_TMP_OBJS_SL_EXEC_SEQ+=$($(cmod_)_TMP_OBJS_SL_EXEC_SEQ)
M_TMP_OBJS_SL_EXEC_MPI+=$($(cmod_)_TMP_OBJS_SL_EXEC_MPI)
M_TMP_OBJS_SL_EXEC+=$($(cmod_)_TMP_OBJS_SL_EXEC)


# exec files
$($(cmod_)_TMP_OBJS_SL_EXEC_SEQ): override CC:=$($(cmod_)_CC)
$($(cmod_)_TMP_OBJS_SL_EXEC_SEQ): override CFLAGS:=$(strip $($(cmod_)_CFLAGS_) $($(cmod_)_CFLAGS_NSL_))
$($(cmod_)_TMP_OBJS_SL_EXEC_SEQ): TYPE:=exec
$($(cmod_)_TMP_OBJS_SL_EXEC_SEQ): $($(cmod_)_CODE_SL_EXEC_SEQ_INCLUDE) $($(cmod_)_CODE_SL_WRAP_SEQ_INCLUDE) $(M_CODE_INCLUDE) $($(cmod_)_PREREQ)

$($(cmod_)_TMP_OBJS_SL_EXEC_MPI): override CC:=$($(cmod_)_MPICC)
$($(cmod_)_TMP_OBJS_SL_EXEC_MPI): override CFLAGS:=$(strip $($(cmod_)_MPICFLAGS_) $($(cmod_)_MPICFLAGS_NSL_))
$($(cmod_)_TMP_OBJS_SL_EXEC_MPI): TYPE:=exec-mpi
$($(cmod_)_TMP_OBJS_SL_EXEC_MPI): $($(cmod_)_CODE_SL_EXEC_INCLUDE) $($(cmod_)_CODE_SL_WRAP_INCLUDE) $(M_CODE_INCLUDE) $($(cmod_)_PREREQ)

$($(cmod_)_TMP_EXEC_SEQ)%.o: $($(cmod_)_SRC_EXEC_SEQ)%.c
	$($(cmod_)_comp)

$($(cmod_)_TMP_EXEC_MPI)%.o: $($(cmod_)_SRC_EXEC_MPI)%.c
	$($(cmod_)_comp)


# libsl files
ifeq ($(cmod),)
 $(cmod_)_MY_LIBSL_A_SEQ:=1
 $(cmod_)_MY_LIBSL_A_MPI:=1
 $(cmod_)_MY_LIBSL_A:=1
else
 $(cmod_)_MY_LIBSL_A_SEQ:=$(if $(subst $(LIBSL_A_SEQ),,$($(cmod_)_LIBSL_A_SEQ)),1)
 $(cmod_)_MY_LIBSL_A_MPI:=$(if $(subst $(LIBSL_A_MPI),,$($(cmod_)_LIBSL_A_MPI)),1)
 $(cmod_)_MY_LIBSL_A:=$(if $(subst $(LIBSL_A),,$($(cmod_)_LIBSL_A)),1)
endif

ifneq ($(cmod_cfg),)
$($(cmod_)_LIBSL_A_SEQ): $($(cmod_)_TMP_OBJS_SL_CORE_SEQ)
$($(cmod_)_LIBSL_A_MPI): $($(cmod_)_TMP_OBJS_SL_CORE_MPI)
$($(cmod_)_LIBSL_A): $($(cmod_)_TMP_OBJS_SL_CORE)
endif
$($(cmod_)_LIBSL_A_SEQ): $($(cmod_)_TMP_OBJS_SL_WRAP_SEQ)
$($(cmod_)_LIBSL_A_MPI): $($(cmod_)_TMP_OBJS_SL_WRAP_MPI)
$($(cmod_)_LIBSL_A): $($(cmod_)_TMP_OBJS_SL_WRAP)

ifeq ($(BULK_AR),)

 ifneq ($(cmod_cfg),)
$($(cmod_)_LIBSL_A_SEQ): $(if $($(cmod_)_OBJS_SL_CORE_SEQ),$($(cmod_)_LIBSL_A_SEQ)($($(cmod_)_OBJS_SL_CORE_SEQ)))
$($(cmod_)_LIBSL_A_MPI): $(if $($(cmod_)_OBJS_SL_CORE_MPI),$($(cmod_)_LIBSL_A_MPI)($($(cmod_)_OBJS_SL_CORE_MPI)))
$($(cmod_)_LIBSL_A): $(if $($(cmod_)_OBJS_SL_CORE),$($(cmod_)_LIBSL_A)($($(cmod_)_OBJS_SL_CORE)))
 endif

$($(cmod_)_LIBSL_A_SEQ): $(if $($(cmod_)_OBJS_SL_WRAP_SEQ),$($(cmod_)_LIBSL_A_SEQ)($($(cmod_)_OBJS_SL_WRAP_SEQ)))
$($(cmod_)_LIBSL_A_MPI): $(if $($(cmod_)_OBJS_SL_WRAP_MPI),$($(cmod_)_LIBSL_A_MPI)($($(cmod_)_OBJS_SL_WRAP_MPI)))
$($(cmod_)_LIBSL_A): $(if $($(cmod_)_OBJS_SL_WRAP),$($(cmod_)_LIBSL_A)($($(cmod_)_OBJS_SL_WRAP)))

%($($(cmod_)_OBJS_SL_CORE_SEQ)): TYPE:=core
%($($(cmod_)_OBJS_SL_CORE_MPI)): TYPE:=core-mpi
%($($(cmod_)_OBJS_SL_WRAP_SEQ)): TYPE:=wrapper
%($($(cmod_)_OBJS_SL_WRAP_MPI)): TYPE:=wrapper-mpi


$($(cmod_)_LIBSL_A_SEQ)(%): $($(cmod_)_TMP_CORE_SEQ)%
	$($(cmod_)_add)

$($(cmod_)_LIBSL_A_SEQ)(%): $($(cmod_)_TMP_WRAP_SEQ)%
	$($(cmod_)_add)

$($(cmod_)_LIBSL_A_MPI)(%): $($(cmod_)_TMP_CORE_MPI)%
	$($(cmod_)_add)

$($(cmod_)_LIBSL_A_MPI)(%): $($(cmod_)_TMP_WRAP_MPI)%
	$($(cmod_)_add)

$($(cmod_)_LIBSL_A)(%): $($(cmod_)_TMP_CORE_SEQ)%
	$($(cmod_)_add)

$($(cmod_)_LIBSL_A)(%): $($(cmod_)_TMP_CORE_MPI)%
	$($(cmod_)_add)

$($(cmod_)_LIBSL_A)(%): $($(cmod_)_TMP_WRAP_SEQ)%
	$($(cmod_)_add)

$($(cmod_)_LIBSL_A)(%): $($(cmod_)_TMP_WRAP_MPI)%
	$($(cmod_)_add)

else

 ifneq ($($(cmod_)_MY_LIBSL_A_SEQ),)
$($(cmod_)_LIBSL_A_SEQ):
	$($(cmod_)_add)
 endif

 ifneq ($($(cmod_)_MY_LIBSL_A_MPI),)
$($(cmod_)_LIBSL_A_MPI):
	$($(cmod_)_add)
 endif

 ifneq ($($(cmod_)_MY_LIBSL_A),)
$($(cmod_)_LIBSL_A):
	$($(cmod_)_add)
 endif

endif


# libsl
$(cmod_)_libsl_seq:
$(cmod_)_libsl_mpi:
$(cmod_)_libsl:

ifneq ($($(cmod_)_MY_LIBSL_A_SEQ),)
M_LIBSL_SEQ+=$(cmod_)_libsl_seq
$(cmod_)_libsl_seq: $($(cmod_)_LIBSL_A_SEQ)
endif

ifneq ($($(cmod_)_MY_LIBSL_A_MPI),)
M_LIBSL_MPI+=$(cmod_)_libsl_mpi
$(cmod_)_libsl_mpi: $($(cmod_)_LIBSL_A_MPI)
endif

ifneq ($($(cmod_)_MY_LIBSL_A),)
M_LIBSL+=$(cmod_)_libsl
$(cmod_)_libsl: $($(cmod_)_LIBSL_A)
endif


# exec
ifeq ($(cmod_),)
$(cmod_)_exec: libsl $($(cmod_)_EXEC)

M_EXEC+=$(cmod_)_exec

$($(cmod_)_EXEC): override LD:=$($(cmod_)_LD)
$($(cmod_)_EXEC): override LDFLAGS:=$(strip $($(cmod_)_LDFLAGS_))
$($(cmod_)_EXEC): TYPE:=exec
$($(cmod_)_EXEC): $(M_TMP_OBJS_SL_EXEC) $($(cmod_)_LIBSL_A)
	$($(cmod_)_make)
endif


# clean
$(cmod_)_TMPFILES:=
ifneq ($(cmod_cfg),)
 $(cmod_)_TMPFILES+=$($(cmod_)_FAR_SCRIPT_DI)
 ifneq ($($(cmod_)_DI_VERSIONS),config)
  $(cmod_)_TMPFILES+=$($(cmod_)_SRC_CORE_SEQ)*_di.c
  $(cmod_)_TMPFILES+=$($(cmod_)_SRC_CORE_MPI)*_di.c
 endif
 $(cmod_)_TMPFILES+=$($(cmod_)_TMP_OBJS_SL_CORE_SEQ) $($(cmod_)_TMP_OBJS_SL_CORE_MPI)
 $(cmod_)_TMPFILES+=$($(cmod_)_TMP_OBJS_SL_TUNE_SEQ) $($(cmod_)_TMP_OBJS_SL_TUNE_MPI) $($(cmod_)_TUNE)
 $(cmod_)_TMPFILES+=$($(cmod_)_TMP_OBJS_SL_DEVL_SEQ) $($(cmod_)_TMP_OBJS_SL_DEVL_MPI) $($(cmod_)_DEVL)
 $(cmod_)_TMPFILES+=$($(cmod_)_TMP_OBJS_SL_DEMO_SEQ) $($(cmod_)_TMP_OBJS_SL_DEMO_MPI) $($(cmod_)_DEMO)
endif
$(cmod_)_TMPFILES+=$($(cmod_)_TMP_OBJS_SL_WRAP_SEQ) $($(cmod_)_TMP_OBJS_SL_WRAP_MPI)
$(cmod_)_TMPFILES+=$($(cmod_)_TMP_OBJS_SL_EXEC_SEQ) $($(cmod_)_TMP_OBJS_SL_EXEC_MPI)
ifneq ($($(cmod_)_LIBSL_A_SEQ),$(LIBSL_A_SEQ))
 $(cmod_)_TMPFILES+=$($(cmod_)_LIBSL_A_SEQ)
endif
ifneq ($($(cmod_)_LIBSL_A_MPI),$(LIBSL_A_MPI))
 $(cmod_)_TMPFILES+=$($(cmod_)_LIBSL_A_MPI)
endif
ifneq ($($(cmod_)_LIBSL_A),$(LIBSL_A))
 $(cmod_)_TMPFILES+=+=$($(cmod_)_LIBSL_A)
endif
$(cmod_)_TMPFILES:=$(strip $($(cmod_)_TMPFILES))

$(cmod_)_TMPDIRS:=
ifneq ($(cmod_cfg),)
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
endif
$(cmod_)_TMPDIRS:=$(strip $($(cmod_)_TMPDIRS))


$(cmod_)_clean: cmod:=$(cmod)
$(cmod_)_clean: cmod_:=$(cmod_)
$(cmod_)_clean:
ifneq ($(cmod),)
	@$(ECHO) "$(cmod):"
endif
	@$(ECHO) " cleaning temporary, object & archive files"
ifneq ($($(cmod_)_TMPFILES),)
	$(RM) $($(cmod_)_TMPFILES)
endif
ifneq ($($(cmod_)_TMPDIRS),)
	$(RM) -r $($(cmod_)_TMPDIRS)
endif

M_CLEAN+=$(cmod_)_clean


# top level (no module)
ifeq ($(mod),)

%.a:
#	@$(ECHO) "global ranlib $@"
	$(RANLIB) $@

info: $(M_INFO)
	@$(ECHO) "core-seq: $(M_CORE_SEQ)"
	@$(ECHO) "core-mpi: $(M_CORE_MPI)"
	@$(ECHO) "core: $(M_CORE)"

#	@$(ECHO) "objs-core-seq: $(M_OBJS_SL_CORE_SEQ)"
#	@$(ECHO) "objs-core-mpi: $(M_OBJS_SL_CORE_MPI)"
#	@$(ECHO) "objs-core: $(M_OBJS_SL_CORE)"
#	@$(ECHO) "tmp-objs-core-seq: $(M_TMP_OBJS_SL_CORE_SEQ)"
#	@$(ECHO) "tmp-objs-core-mpi: $(M_TMP_OBJS_SL_CORE_MPI)"
#	@$(ECHO) "tmp-objs-core: $(M_TMP_OBJS_SL_CORE)"

	@$(ECHO) "wrap-seq: $(M_WRAP_SEQ)"
	@$(ECHO) "wrap-mpi: $(M_WRAP_MPI)"
	@$(ECHO) "wrap: $(M_WRAP)"

#	@$(ECHO) "objs-wrap-seq: $(M_OBJS_SL_WRAP_SEQ)"
#	@$(ECHO) "objs-wrap-mpi: $(M_OBJS_SL_WRAP_MPI)"
#	@$(ECHO) "objs-wrap: $(M_OBJS_SL_WRAP)"
#	@$(ECHO) "tmp-objs-wrap-seq: $(M_TMP_OBJS_SL_WRAP_SEQ)"
#	@$(ECHO) "tmp-objs-wrap-mpi: $(M_TMP_OBJS_SL_WRAP_MPI)"
#	@$(ECHO) "tmp-objs-wrap: $(M_TMP_OBJS_SL_WRAP)"

	@$(ECHO) "exec-seq: $(M_EXEC_SEQ)"
	@$(ECHO) "exec-mpi: $(M_EXEC_MPI)"
	@$(ECHO) "exec: $(M_EXEC)"

	@$(ECHO) "tmp-objs-exec-seq: $(M_TMP_OBJS_SL_EXEC_SEQ)"
	@$(ECHO) "tmp-objs-exec-mpi: $(M_TMP_OBJS_SL_EXEC_MPI)"
	@$(ECHO) "tmp-objs-exec: $(M_TMP_OBJS_SL_EXEC)"

	@$(ECHO) "libsl-seq: $(M_LIBSL_SEQ)"
	@$(ECHO) "libsl-mpi: $(M_LIBSL_MPI)"
	@$(ECHO) "libsl: $(M_LIBSL)"


core_seq: $(M_CORE_SEQ)
core_mpi: $(M_CORE_MPI)
core: $(M_CORE)

wrapper_seq: $(M_WRAP_SEQ)
wrapper_mpi: $(M_WRAP_MPI)
wrapper: $(M_WRAP)

libsl_seq: $(M_LIBSL_SEQ)
libsl_mpi: $(M_LIBSL_MPI)
libsl: $(M_LIBSL)

exec: $(M_EXEC)

tune: $(M_TUNE)

devel: $(M_DEVL)

demo: $(M_DEMO)

clean: $(M_CLEAN)
	$(RM) $(LIBSL_A_SEQ) $(LIBSL_A_MPI) $(LIBSL_A) $(EXEC) $(TUNE) $(DEVL) $(DEMO) $(CLEAN_FILES)

endif


# adding object files to a single archive can not be done in parallel (see "Dangers When Using Archives" in the GNU Make Manual)
ifeq ($(BULK_AR),)
.NOTPARALLEL:
endif
