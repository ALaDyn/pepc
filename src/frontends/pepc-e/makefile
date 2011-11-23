
MYROOT = ../../../
include $(MYROOT)makefile.defs


FRONTENDNAME = pepce
BACKEND      = coulomb

# -------------------------------------------------------------------------------------------
# Source files
SOURCES.f90 =   physvars.f90 param_dump.f90 face.f90 files.f90 pusher.f90 randion.f90 \
                setup.f90 configure.f90 special_start.f90 \
                energies.f90 cleanup.f90 pepce.f90 module_diagnostics.f90
SOURCES.UTILS.f90 = module_vtk.f90 module_base64.f90 module_checkpoint.f90 \
                benchmarking.f90 module_velocity_setup.f90
# -------------------------------------------------------------------------------------------

include $(MYROOT)$(FRONTENDDIR)makefile.frontend