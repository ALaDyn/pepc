# PEPC -  Pretty Efficient Parallel Coulomb-solver

##### Authors:  
Paul Gibbon, Mathias Winkel, Robert Speck, Lukas Arnold  
Forschungszentrum Juelich GmbH  
Juelich Supercomputing Centre  

#### Webpage:  
[http://www.fz-juelich.de/ias/jsc/pepc](http://www.fz-juelich.de/ias/jsc/pepc)  
#### E-Mail:  
[pepc@fz-juelich.de](mail-to:pepc@fz-juelich.de)  
#### CI:  
[![pipeline status](https://gitlab.jsc.fz-juelich.de/SLPP/pepc/pepc/badges/master/pipeline.svg)](https://gitlab.jsc.fz-juelich.de/SLPP/pepc/pepc/-/commits/master) 


# 0. LICENSE

This file is part of PEPC - The Pretty Efficient Parallel Coulomb Solver.

Copyright (C) 2002-2021  
Juelich Supercomputing Centre,   
Forschungszentrum Juelich GmbH,  
Germany

PEPC is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

PEPC is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with PEPC. If not, see [http://www.gnu.org/licenses/](http://www.gnu.org/licenses/).


# 1. REQUIREMENTS

 - A reasonably modern Fortran compiler with support for
   Fortran 2003 object orientation, e.g.:
    - GCC >= 4.6
    - Intel >= 12.1
    - IBM XL Fortran >= 12
 - A C compiler and a C preprocessor with support for
   variadic macros.
 - An MPI library with support for MPI_THREAD_MULTIPLE.
   PEPC will complain about missing support at run time.
 - Support for POSIX threads (pthreads).


# 2. COMPILATION

For compiling this program, the machine-dependent compiler 
names and flags have to be set in a definition file called
`makefile.defs` in the root directory.

Several exemplary definitions are available in the subdirectory
`./makefiles`. For example, on JUQUEEN you can simply call

```sh
ln -sf ./makefiles/makefile.defs.juqueen_xlf ./makefile.defs
```

from the root directory to create a symbolic link the the 
appropriate definitions. Some makefile.defs contain
comments on correct MPI implementations etc. - do not
forget to take a look there for hints.

After providing the makefile.defs file, you can simply
call 

```
make help
```

to show information/recommendation on machine specifics.

After setting the proper environment (modules, paths) 
call 

```sh
make pepc-mini
```

to build the `pepc-mini` frontend into the `./bin/` directory. 
Parallel make (i.e. `make -j`) should work.

There are several different frontends available at the moment:

* `pepc-b`:  
Laser/beam-plasma with magnetic fields  
* `pepc-e`: 
pure electrostatics,  
simple molecular dynamics,  
some diagnostics for energy, performance, vtk output etc.
* `pepc-mini`:  
pure electrostatics  
simple molecular dynamics  
no diagnostics  
minimum requirements to get pepc running  
* `pepc-mw`:  
frontend for laser-plasma interaction (experimental)  
* `pepc-s`:  
library version for inclusion in ScaFaCoS project  
* `pepc-nn`:  
tree-based nearest neighbour search  
* `pepc-sph`:  
gravitation and smoothed particle hydrodynamics (sph) for astrophysics research (experimental)  
`pepc-v`:  
vortex dynamics using the vortex particle method  

To build an alternative frontend, just call

```sh
make pepc-e
```
or
```sh
make pepc-mw
```

All frontends can be built using 

```sh
make (-j) All
```

At the current stage, there is no real documentation available for the
different frontends. However, you might simply want to take a look at 
the respective sourcecode to find out what they are doing.


# 3. RUNNING THE PROGRAM

Usually, the frontend's source directories contain a sample
input deck for the frontends. For running pepc-essential
it is for example called `params`. It contains user-adjustable
parameters and can be fed to the executable as first command
line parameter:
```sh
cd bin
mpirun -np 32 pepc-essential ../src/frontneds/pepc/essential/params
```


# 4. DOCUMENTATION

Rudimentary doxygen documentation is available by calling

```sh
make doc
```

from the root directory. A users guide is in preparation. 


# 5. DIRECTORY STRUCTURE / ADDING OWN FUNCTIONALITY

Inside the `./src/` directory, you will find four subdirectories:
- treecode: PEPC kernel, everything that is necessary for the pure algorithmic part of the treecode
- interaction_specific: interaction specific backends. The different
  subdirectories herein (currently: coulomb, vortex, and
  nearestneighbour) provide data structures and functions
  for the different applications. See inline documentation
  in the sourcecode (especially inside the coulomb-subdir,
  which should be well documented) to finc out about what the
  functions should do and which of them are necessary.
  The only files that must be provided in this directory are
  (names may not be changed, public functions and datastructures
  inside these files are mandatory):
  * `module_interaction_specific.f90`: data structures and functions for manipulating them
  * `module_calc_force.f90`: functions for actual force-law and multipole acceptance criterion etc.
  * `makefile.backend`: backend specific modifications to treecode makefile, may be empty
- utils: source code of utilities (mainly for treecode diagnostics, vtk-output etc.)
- frontends: different applications that utilize the treecode for their respective very specific purpose
            
In case you want to use pepc for developing an own treecode-based N-body code,
you might start by copying and modifying the pepc-mini frontend, which is a 
very simple coulomb-MD programme. It uses the coulomb backend, that implements
an expansion of the plummer potential 1/sqrt(r^2+eps^2) up to quadrupole order.

Take care that your frontend-directory is called "pepc-something" with no 
further minus sign ("-") to be automatically recognized in the build system.

If you want to provide an own interaction-specific backend (for using other
multipole orders and/or force laws), just copy and modify the coulomb subdirectory
there. The backends do not have to be registered in some makefile, but are selected 
inside the main makefile by creating a dependency on lippepc.backend_directory_name
instead of libpepc.coulomb. In case you only need to modify the interaction specific
types but not the functions that are dealing with them (for example for adding 
velocity, mass, etc.), take a look at the coulomb-backend. There, this is done for the
pepc-mini frontend, that uses other types than for example pepc-b while still
using the same force expression. See also variable BACKENDTYPE in pepc-mini's 
`makefile.include`.