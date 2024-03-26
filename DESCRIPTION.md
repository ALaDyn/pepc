<!---
This file is automatically included in the Helmholtz Research Software Directory: https://helmholtz.software/software/pepc
-->

## Algorithm

The oct-tree method was originally introduced by Josh Barnes & Piet Hut in the mid 1980s to speed up astrophysical N-body simulations with long range interactions, see [Nature 324, 446 (1986)](http://dx.doi.org/10.1038/324446a0). Their idea was to use successively larger multipole-groupings of distant particles to reduce the computational effort in the force calculation from the usual O(N2) operations needed for brute-force summation to a more amenable O(N log N). Though mathematically less elegant than the Fast Multipole Method, the Barnes-Hut algorithm is well suited to dynamic, nonlinear problems, can be quickly adapted to numerous interaction kernels and combined with multiple-timestep integrators.

![image](https://www.fz-juelich.de/en/ias/jsc/about-us/structure/simulation-and-data-labs/sdl-plasma-physics/pepc/pepc_overview.png/@@images/image/preview)

*Simulation of interaction of laser radiation with neutral nanocluster with 3871 ions and electrons. Bottom left: spacefilling Hilbert curve used in the code for domain decomposition; Bottom right: Oct-tree domain decomposition represented by the so-called branch nodes.*

The PEPC project (Pretty Efficient Parallel Coulomb Solver) is a public tree code that has been developed at [Jülich Supercomputing Centre](https://www.fz-juelich.de/en/ias/jsc) since the early 2000s. Our tree code is a non-recursive version of the Barnes-Hut algorithm, using a level-by-level approach to both tree construction and traversals. The parallel version is a hybrid MPI/pthreads implementation of the Warren-Salmon 'Hashed Oct-Tree' scheme, including several refinements of the tree traversal routine - the most challenging component in terms of scalability.

The code is structurally divided into three parts:

  1. **kernel routines** that handle all tree code specific data structures and communication as well as the actual tree traversal.
  2. **interaction-specific modules**, i.e. routines that apply to specific interaction kernels and multipole expansions. Currently, the following interaction kernels are available:
     - Coulomb-interaction/gravitation,
     - algebraic kernels for vortex methods,
     - Darwin for magnetoinductive plasmas (no EM wave propagation),
     - nearest-neighbour interactions for smooth particle hydrodynamics (SPH).
  3. **'front-end'** applications. For example
     - PEPC-mini, a skeleton molecular dynamics program with different diagnostics including VTK output for convenient visualization,
     - PEPC-b, a code for laser- or particle beam-plasma interactions as well as plasma-wall interactions,
     - PEPC-s, a library version for the ScaFaCoS project,
     - PEPC-v, an application for simulating vortex dynamics using the vortex particle method,
     - PEPC-dvh, vortex dynamics using the diffused vortex hydrodynamics method,
     - PEPC-g, gravitational interaction and optional smooth particle hydrodynamics frontend (SPH) for simulating stellar discs consisting of gas and dust, developed together with Max Planck Institute for Radio Astronomy (MPIfR) Bonn,
     - several internal experimental frontends.

Due to this structure, the adaption to further interaction kernels as well as additional applications and experimental variations to the tree code algorithm can be implemented conveniently.

![image](https://www.fz-juelich.de/en/ias/jsc/about-us/structure/simulation-and-data-labs/sdl-plasma-physics/pepc/pepc_structure-1.png/@@images/image/preview)

*Structure of the tree code framework. Currently, PEPC supports three interaction-specific modules and several respective frontends for different applications ranging from plasma physics through vortex-particle methods to smooth particle hydrodynamics. Due to well-defined interfaces, symbolized by the grey blocks, the implementation of further interaction models as well as additional frontends is straightforward.*

## Implementation details and scaling

PEPC itself is written in Fortran 2003 with some C wrappers for POSIX functions. In addition to the main pthreads-parallelised tree traversal, a version based on OpenMP is available and there has also been a version using SMPSs. Besides the hybrid parallelisation, a number of improvements to the original Warren-Sallmon ‘hashed oct-tree’ scheme have been included and allow for an excellent scaling of the code on different architectures with up to 2,048,000,000 particles across up to 294,912 processors (at the then largest machine at JSC: JUGENE).

In addition, it has been adapted to > 32 parallel threads per MPI rank to take full advantage of the capabilities of the BlueGene/Q installation JUQUEEN that was available at JSC. As a result it showed great parallel scalability across the full machine and qualified for the [High-Q Club](https://www.fz-juelich.de/ias/jsc/high-q-club).

![image](https://www.fz-juelich.de/en/ias/jsc/about-us/structure/simulation-and-data-labs/sdl-plasma-physics/pepc/pepc_scaling_juqueen.png/@@images/image/preview)

*Weak scaling and parallel efficiency of the parallel Barnes-Hut tree code PEPC across the full Blue Gene/Q installation JUQUEEN at JSC.*
