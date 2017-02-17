# Geometry Processing – Parameterization

> **To get started:** Fork this repository then issue
> 
>     git clone --recursive http://github.com/[username]/geometry-processing-parameterization.git
>

## Installation, Layout, and Compilation

See
[introduction](http://github.com/alecjacobson/geometry-processing-introduction).

## Execution

Once built, you can execute the assignment from inside the `build/` by running
on a given mesh:

    ./parameterization [path to mesh.obj]

## Background

In this assignment we will explore how to _flatten_ a surface
[embedded](https://en.wikipedia.org/wiki/Embedding) (or even just
[immersed](https://en.wikipedia.org/wiki/Immersion)) in $\R^3$ to the flat
plane (i.e. $\R^2$).

This process is often referred to as [parameterization]() because the
two-dimensional coordinate system of the flattened mesh can now be interpreted
as a parameterization of the 3D surface.

> Surface in 3D --> parameterization in 2D --> parameterized surface in 3d


### Whitelist

 - `igl::boundary_list.h`
