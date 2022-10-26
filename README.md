# SedimentModeling 

### About
SedimentModeling is a C++ code developed to model the visco-elastic compaction of sediment observed in nature. This code solves coupled partial differential equations and phase-field order parameter equations using the finite-element method to model sediment compaction. The master branch has a one-dimensional implementation. The other branches consist of 2D extensions of the sediment compaction model. The salient features of the computational implementation are : high-fidelity finite-element implementation with support for direct solvers and fully implicit time-stepping schemes. 

### Installation:

SedimentModeling code builds on top of the deal.II library.

1) Install CMake, PETSc, Trilinos, SLEPc, p4est, and deal.II (version 9.3.0 recommended)<br>

2) Clone the SedimentModel GitHub repository <br>
```
$ git clone git clone https://github.com/cmmg/SedimentModeling.git
$ cd SedimentModeling
$ git checkout master
$ cmake .
$ make -j nprocs
  ```
[here nprocs denotes the number of processors]

### Visualization:

  Output of the primary fields is in the vtk format (parallel:*.pvtu, serial:*.vtu files). These can be visualized with the following open-source applications:
  1. VisIt (https://visit-dav.github.io/visit-website/releases-as-tables/)
  2. Paraview (http://www.paraview.org/download/)


License
-------
GNU Lesser General Public License (LGPL). Please see the file LICENSE for details.
