# This is a pepc-frontend makefile.
# In principle, it only informs the buildsystem about the frontends (and thus
# the final binaries) name, the backend to be included and all source files.
# 
# Everything between the two dashed lines can and should be adapted to your
# needs:
#   FRONTENDNAME = name of your final binary (ideally similar to your frontend name)
#   BACKEND      = directory name of the backend to be used
#   SOURCES.f90  = file names of all source files inside your frontend directory
#   SOURCES.UTILS.f90 = names of all files from src/utils directory that must be included
#
# Please do *not* modify any other files, especially do not touch makefile.frontend unless
# you really know what you are doing and you stay compatible to *all* frontends.
# Also, do not modify anything else outside the dashed lines.
#
# If your frontend needs specific modifications to the things set up in makefile frontend, 
# you can insert them below the 'include'-statement at the end of this file
#
#
MYROOT = ../../../
include $(MYROOT)makefile.defs


# -------------------------------------------------------------------------------------------
FRONTENDNAME = pepcmw
BACKEND      = coulomb
# Source files
SOURCES.f90 =   module_units.f90 physvars.f90 module_io.f90 module_laser.f90 module_pusher.f90 \
                module_workflow.f90 module_param_dump.f90 setup.f90 density_helmholtz.f90 \
                module_energy.f90 cleanup.f90 pepc.f90 \
                module_acf.f90 sum_radial.f90 module_setup.f90 \
                module_diagnostics.f90 module_fields.f90 module_icosahedron.f90
SOURCES.UTILS.f90 = module_vtk.f90 module_base64.f90 module_checkpoint.f90 \
                benchmarking.f90 module_velocity_setup.f90 module_treediags.f90 \
                module_vtk_helpers.f90 module_directsum.f90
# -------------------------------------------------------------------------------------------

include $(MYROOT)$(FRONTENDDIR)makefile.frontend
