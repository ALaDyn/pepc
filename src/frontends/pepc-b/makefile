
MYROOT = ../../../
include $(MYROOT)makefile.defs


FRONTENDNAME = pepcb
BACKEND      = coulomb

# -------------------------------------------------------------------------------------------
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
