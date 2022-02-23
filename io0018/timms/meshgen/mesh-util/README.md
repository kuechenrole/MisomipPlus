## `MESH2D: A MATLAB-based mesh generator`

`MESH2D` is a `MATLAB` / `OCTAVE`-based unstructured mesh-generator for two-dimensional polygonal geometries. `MESH2D` provides a range of simple, yet effective two-dimensional meshing algorithms, including a variation on the "classical" Delaunay refinement technique, a new "Frontal"-Delaunay refinement scheme, a non-linear mesh optimisation method, and auxiliary mesh and geometry pre- and post-processing facilities. 

Algorithms implemented in `MESH2D` are "provably-good" - ensuring convergence, geometrical and topological correctness, and providing guarantees on algorithm termination and worst-case element quality bounds. Support for user-defined "mesh-spacing" functions and "multi-part" geometry definitions is provided, allowing `MESH2D` to handle a wide range of complex domain types and user-defined constraints. `MESH2D` typically generates very high-quality output, appropriate for a variety of finite-volume and/or finite-element type applications.

<p align="center">
  <img src = "../master/poly-data/lake-1-small.png"> &nbsp &nbsp &nbsp &nbsp
  <img src = "../master/poly-data/lake-2-small.png">
</p>

The `MESH2D` package has been around now for a while (!), having first been released around 2006 as an extension of the <a href="http://persson.berkeley.edu/distmesh/">`DISTMESH`</a> algorithm of Persson and Strang. 

As of `MESH2D release 3.0.0`, the code has been completely rewritten, and is now essentially a simplified version of my <a href="https://github.com/dengwirda/jigsaw-matlab/">`JIGSAW`</a> mesh-generation algorithm (a `C++` code), providing a more straightforward and accessible implementation of Delaunay-based triangulation and mesh optimisation techniques. 


## `Starting Out`

After downloading and unzipping the current <a href="https://github.com/dengwirda/mesh2d/archive/master.zip">repository</a>, navigate to the installation directory within <a href="http://www.mathworks.com">`MATLAB`</a> / <a href="https://www.gnu.org/software/octave">`OCTAVE`</a> and run the set of examples contained in `tridemo.m`:
```
tridemo( 0); % a very simple example to get everything started.
tridemo( 1); % investigate the impact of the "radius-edge" threshold.
tridemo( 2); % Frontal-Delaunay vs. Delaunay-refinement refinement.
tridemo( 3); % explore impact of user-defined mesh-size constraints.
tridemo( 4); % explore impact of "hill-climbing" mesh optimisations.
tridemo( 5); % assemble triangulations for "multi-part" geometries.
tridemo( 6); % assemble triangulations with "internal" constraints.
tridemo( 7); % investigate the use of "quadtree"-type refinement.
tridemo( 8); % explore use of custom, user-defined mesh-size functions.
tridemo( 9); % larger-scale problem, mesh refinement + optimisation. 
tridemo(10); % medium-scale problem, mesh refinement + optimisation. 
```

For <a href="https://www.gnu.org/software/octave">`OCTAVE`</a> users, performance can be improved by compiling elements of the `MESH2D` library. Running `compile.m` within the `MESH2D` installation directory will complete the build process (note: requires a `-dev` installation of <a href="https://www.gnu.org/software/octave">`OCTAVE`</a>).


## `Attribution!`

If you make use of `MESH2D` please include a reference to the following! `MESH2D` is designed to provide a simple and easy-to-understand implementation of Delaunay-based mesh-generation techniques. For a much more advanced, and fully three-dimensional mesh-generation library, see the <a href="https://github.com/dengwirda/jigsaw-matlab/">`JIGSAW`</a> package. `MESH2D` makes use of the <a href="https://github.com/dengwirda/aabb-tree">`AABBTREE`</a> and <a href="https://github.com/dengwirda/find-tria">`FINDTRIA`</a> packages to compute efficient spatial queries and intersection tests. 

`[1]` - Darren Engwirda, <a href="http://hdl.handle.net/2123/13148">Locally-optimal Delaunay-refinement and optimisation-based mesh generation</a>, Ph.D. Thesis, School of Mathematics and Statistics, The University of Sydney, September 2014.
