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
FRONTENDNAME = pepcs
BACKEND      = coulomb
# Source files
SOURCES.f90 = pepc-s.f90 main.f90
SOURCES.UTILS.f90 = 
# -------------------------------------------------------------------------------------------

include $(MYROOT)$(FRONTENDDIR)makefile.frontend