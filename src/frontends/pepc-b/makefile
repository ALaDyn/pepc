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
FRONTENDNAME = pepcb
BACKEND      = coulomb
# Source files
SOURCES.f90 =  module_physvars.f90 \
               module_particle_props.f90 \
               module_utilities.f90 \
               module_velocity_setup.f90 \
               module_field_grid.f90 \
               module_laser.f90 \
               module_io.f90 \
               module_geometry.f90 \
               module_integration_scheme.f90 \
               module_diagnostics.f90 \
               module_particle_beam.f90 \
               module_particle_boundaries.f90 \
               setup_arrays.f90 setup.f90 configure.f90 integrator.f90 diagnostics.f90 pepcb.f90 
#               module_xnbody.f90
SOURCES.UTILS.f90 = module_gle.f90 
# -------------------------------------------------------------------------------------------

include $(MYROOT)$(FRONTENDDIR)makefile.frontend
