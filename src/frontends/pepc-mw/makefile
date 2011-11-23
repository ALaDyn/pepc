MYROOT = ../../../
include $(MYROOT)makefile.defs


FRONTENDNAME = pepcmw
BACKEND      = coulomb

# -------------------------------------------------------------------------------------------
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
