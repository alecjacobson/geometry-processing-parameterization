# Geometry Processing – Parameterization

> **To get started:** Clone this repository with
> 
>     git clone --recursive http://github.com/alecjacobson/geometry-processing-parameterization.git
>

## Installation, Layout, and Compilation

See
[introduction](http://github.com/alecjacobson/geometry-processing-introduction).

## Execution

Once built, you can execute the assignment from inside the `build/` by running
on a given mesh:

    ./parameterization [path to mesh with boundary.obj]

![When your app is implemented correctly you'll be able to cycle through
viewing the mesh, viewing the parameterization, toggling the checkerboard
pattern and switching between parameterization
methods.](images/beetle-cycle-screenshots.gif)

## Background

In this assignment we will explore how to _flatten_ a surface
[embedded](https://en.wikipedia.org/wiki/Embedding) (or even just
[immersed](https://en.wikipedia.org/wiki/Immersion)) in <img src="./tex/d03c1e146df015e061405cc425738d83.svg?invert_in_darkmode" align=middle width=18.424726649999986pt height=26.76175259999998pt/> to the flat
plane (i.e. <img src="./tex/433badc501d4f8a183b14684b47f305e.svg?invert_in_darkmode" align=middle width=18.424726649999986pt height=26.76175259999998pt/>).

This process is often referred to as
[parameterization](https://en.wikipedia.org/wiki/Parametrization#Parametrization_techniques)
because the two-dimensional coordinate system of the flattened mesh can now be
interpreted as a parameterization of the 3D surface.

![A triangle mesh of a [VW
Beetle](https://en.wikipedia.org/wiki/Volkswagen_Beetle) is _parameterized_ by
flattening the mesh to the <img src="./tex/ee5f11272c9cd93256bbf7ba019c3953.svg?invert_in_darkmode" align=middle width=17.96813369999999pt height=14.15524440000002pt/>-plane. There the <img src="./tex/6dbb78540bd76da3f1625782d42d6d16.svg?invert_in_darkmode" align=middle width=9.41027339999999pt height=14.15524440000002pt/>- and <img src="./tex/6c4adbc36120d62b98deef2a20d5d303.svg?invert_in_darkmode" align=middle width=8.55786029999999pt height=14.15524440000002pt/>- coordinates
(orange and white lines) can be directly interpreted as a parameterization of
the surface. ](images/beetle-uv-parameterization-low-res.png)

In this assignment, we are given a representation of the surface in 3D as a
[triangle mesh](https://en.wikipedia.org/wiki/Triangle_mesh) with a list of
<img src="./tex/55a049b8f161ae7cfeb0197d75aff967.svg?invert_in_darkmode" align=middle width=9.86687624999999pt height=14.15524440000002pt/> vertex positions <img src="./tex/375e8c9ba3d3b259a673c6015274f7fe.svg?invert_in_darkmode" align=middle width=71.47064264999999pt height=26.76175259999998pt/>, so then the goal is to assign <img src="./tex/ee5f11272c9cd93256bbf7ba019c3953.svg?invert_in_darkmode" align=middle width=17.96813369999999pt height=14.15524440000002pt/>
coordinates to each vertex <img src="./tex/1b1d3e72f12695d11aa8bd4d01b7e32e.svg?invert_in_darkmode" align=middle width=71.45921639999999pt height=26.76175259999998pt/>. 

In general, a 3D surface cannot be flattened onto the plane without
_**distortion**_. Some parts of the surface will have to be stretched and other
squished. Surfaces with [topological
handles](https://en.wikipedia.org/wiki/Handle_decomposition) or without [a
boundary](https://en.wikipedia.org/wiki/Surface_(topology)#Closed_surfaces) must be
cut. 

### Mass-spring methods

If we view our triangle mesh surface as a simple
[graph](https://en.wikipedia.org/wiki/Graph_(discrete_mathematics)) then the
surface flattening problem reduces to [graph
drawing](https://en.wikipedia.org/wiki/Graph_drawing). Distortion can be
measured in terms of the relative change in lengths between neighboring
vertices.

We can pose the graph drawing problem as an optimization over node locations so
that the lengths between neighboring vertices are minimized:

<p align="center"><img src="./tex/3e21fa6f1989f96bf48dfaf5659e9898.svg?invert_in_darkmode" align=middle width=161.4557736pt height=40.548151049999994pt/></p>

where <img src="./tex/cd58ffd8803a3e8e5065d6293ab76e5f.svg?invert_in_darkmode" align=middle width=127.65767849999997pt height=27.91243950000002pt/> holds a list of edge indices into <img src="./tex/26eb59da31fb48cb17abfe4c6dc80375.svg?invert_in_darkmode" align=middle width=14.554737449999989pt height=22.55708729999998pt/>. This
energy has a physical interpretation as the [potential
energy](https://en.wikipedia.org/wiki/Potential_energy) of
[mass-spring](https://en.wikipedia.org/wiki/Simple_harmonic_motion#Mass_on_a_spring)
system. Each edge represents a spring with zero rest length, all springs have
uniform [stiffness](https://en.wikipedia.org/wiki/Hooke's_law#Spring_energy),
and all vertices have equivalent (unit) mass.

Without additional constraints, this minimization has a trivial solution: map
all vertices to the same point, e.g., <img src="./tex/e6157f8faddf86037cbc50411258d3a0.svg?invert_in_darkmode" align=middle width=100.17695159999998pt height=24.65753399999998pt/>.

We can avoid this by fixing the mapping of certain vertices. If we choose these
fixed vertices arbitrarily we will in general get overlaps in the flattening.
For graph drawing this means that edges cross each other; for surface
parameterization this means that multiple triangles cover the same patch of the
<img src="./tex/ee5f11272c9cd93256bbf7ba019c3953.svg?invert_in_darkmode" align=middle width=17.96813369999999pt height=14.15524440000002pt/>-plane and some of those triangles are upside down. This problem is often
referred to as _fold overs_ or lack of
[injectivity](https://en.wikipedia.org/wiki/Injective_function).

![The face of [Max Planck](https://en.wikipedia.org/wiki/Max_Planck) is
parameterized using a mass-spring system. More and more vertices are fixed
explicitly along the boundary. With only a few fixed vertices there are severe
overlaps and degeneracies in the interior. When the entire boundary is fixed to
the circle, there are no overlaps.](images/max-tutte-boundary.gif)

In 1963, [Tutte showed](https://en.wikipedia.org/wiki/Tutte_embedding) that if
the boundary of a disk-topology mesh is fixed to a [convex
polygon](https://en.wikipedia.org/wiki/Convex_polygon) (and all spring
stiffness are positive) then minimizing the energy above will result in an
injective (i.e., fold-over-free) flattening.

While avoiding fold-overs is important, Tutte-style mappings suffer from a
couple problems.

If uniform spring stiffness are used, then the mapping in the <img src="./tex/ee5f11272c9cd93256bbf7ba019c3953.svg?invert_in_darkmode" align=middle width=17.96813369999999pt height=14.15524440000002pt/> domain will
try to make all edges the same length. Combined with the boundary constraints,
the flattened mesh will have smoothly varying edge-lengths and near-equilateral
triangles _regardless_ of the triangles shapes and sizes on the surface mesh.

![The Tutte embedding of the 3D Ogre mesh leads to severe distortion.
Overlaying a checkerboard pattern on the 2D domain and visualizing it on the 3D
surface shows the wobbliness of the non-smooth mapping and
stretching.](images/keenan-ogre-tutte.jpg)

We can try to remedy this by introducing a non-uniform weight or spring
stiffness  <img src="./tex/64e70e84545b2941bed8aa7fe2211cde.svg?invert_in_darkmode" align=middle width=22.523917349999987pt height=14.15524440000002pt/> for each edge <img src="./tex/45cec34bd24a46fe74066f1fef04d815.svg?invert_in_darkmode" align=middle width=37.11794624999999pt height=24.65753399999998pt/>:
<p align="center"><img src="./tex/a92d747bc2194259cbdb6367e631c1fb.svg?invert_in_darkmode" align=middle width=184.80156419999997pt height=40.548151049999994pt/></p>


For example, we could weigh the distortion of shorter edges (on the 3D mesh)
more than longer ones: <img src="./tex/98cfd14aaed04f2386a2169b68c44d6c.svg?invert_in_darkmode" align=middle width=130.58488244999998pt height=24.65753399999998pt/>. See ["Parametrization and
smooth approximation of surface triangulations" [Floater 1996]](papers/Floater97.pdf). This will at
best help tame _**length distortion**_. The "shapes" (i.e., aspect ratios) of
triangles will only be indirectly preserved. We need a way to discourage _area
distortion_ and _angle distortion_.

To do this, let's write the energy minimization problem above in matrix form:

<p align="center"><img src="./tex/4c50bfb852a530ece7a3d0f4c46303c5.svg?invert_in_darkmode" align=middle width=129.8754798pt height=32.990165999999995pt/></p>

where <img src="./tex/e9a4432c2ea64a0aa58762f307040216.svg?invert_in_darkmode" align=middle width=69.85918335pt height=26.17730939999998pt/> is a sparse matrix with:
<p align="center"><img src="./tex/d17d7ae3aaf4ec7a1ba7cf1262a3d47a.svg?invert_in_darkmode" align=middle width=320.4661317pt height=78.90491235pt/></p>


> #### What's up with the <img src="./tex/7f255cacf9e27c987ea2c249653a3e6c.svg?invert_in_darkmode" align=middle width=28.35616739999999pt height=24.65753399999998pt/> in the energy?
>
> The degrees of freedom in our optimization are a collected in the _matrix_
> <img src="./tex/b59b3daaef5011363b567f9ba7c31bf6.svg?invert_in_darkmode" align=middle width=71.45921639999999pt height=26.76175259999998pt/> with two columns. The energy is written as the
> [trace](https://en.wikipedia.org/wiki/Trace_(linear_algebra)) of the
> quadratic form (a.k.a. matrix) <img src="./tex/3e1f41ee3f8955763490ada25e7706e2.svg?invert_in_darkmode" align=middle width=72.69021704999999pt height=26.17730939999998pt/> applied to <img src="./tex/35531be55273dc37ee90083451d089ff.svg?invert_in_darkmode" align=middle width=14.54330789999999pt height=22.55708729999998pt/>. In effect,
> this is really applying <img src="./tex/61ccc6d099c3b104d8de703a10b20230.svg?invert_in_darkmode" align=middle width=14.20083224999999pt height=22.55708729999998pt/> to each column of <img src="./tex/35531be55273dc37ee90083451d089ff.svg?invert_in_darkmode" align=middle width=14.54330789999999pt height=22.55708729999998pt/> independently and summing
> the result:
>

<p align="center"><img src="./tex/05052801ffe8e68655d9f3edda2424fa.svg?invert_in_darkmode" align=middle width=298.24741155pt height=166.0315833pt/></p>

>
> The benefits of energies written as the trace of a quadratic form applied to a
> matrix include: 1) each column can be optimized _independently_ (assuming
> constraints are also separable by column), and this is often the case when
> columns correspond to coordinates (u, v, etc.); and 2) the quadratic form for
> each columns is the same (the same <img src="./tex/61ccc6d099c3b104d8de703a10b20230.svg?invert_in_darkmode" align=middle width=14.20083224999999pt height=22.55708729999998pt/>). For quadratic energy minimization,
> this means that we can precompute work (e.g., [Cholesky
> facotorization](https://en.wikipedia.org/wiki/Cholesky_decomposition)) on
> <img src="./tex/61ccc6d099c3b104d8de703a10b20230.svg?invert_in_darkmode" align=middle width=14.20083224999999pt height=22.55708729999998pt/> and take advantage of it for solving with <img src="./tex/ea9f4d50b9b585000b4579d3d9e39a01.svg?invert_in_darkmode" align=middle width=21.09585554999999pt height=22.55708729999998pt/> and <img src="./tex/4dc30bac56c3bd7d970c8369cbe6b99a.svg?invert_in_darkmode" align=middle width=21.09585554999999pt height=22.55708729999998pt/> and we might
> even solve [in parallel](https://en.wikipedia.org/wiki/SIMD).

#### Dirichlet energy

We should immediately recognize this sparsity structure from the discrete
Laplacians considered in the previous assignments. If <img src="./tex/e85343fc9fd54d37ce80267ae3fa3764.svg?invert_in_darkmode" align=middle width=53.482629749999994pt height=21.18721440000001pt/>, then <img src="./tex/80637df1ca7533740cc7b3fdd1ab540b.svg?invert_in_darkmode" align=middle width=11.36979854999999pt height=22.55708729999998pt/>
is the _uniform Laplacian_ (a.k.a., [graph
Laplacian](https://en.wikipedia.org/wiki/Laplacian_matrix)). If <img src="./tex/64e70e84545b2941bed8aa7fe2211cde.svg?invert_in_darkmode" align=middle width=22.523917349999987pt height=14.15524440000002pt/> is
based on edge-lengths, then <img src="./tex/80637df1ca7533740cc7b3fdd1ab540b.svg?invert_in_darkmode" align=middle width=11.36979854999999pt height=22.55708729999998pt/> corresponds to a physical static equilibrium
problem for a linear spring system. 

But we have more information than edges: we know that our graph is really a
discrete representation of a two-dimensional surface. Wobbliness distortions in
the parameterization correspond to high _variation_ in the <img src="./tex/6dbb78540bd76da3f1625782d42d6d16.svg?invert_in_darkmode" align=middle width=9.41027339999999pt height=14.15524440000002pt/> and <img src="./tex/6c4adbc36120d62b98deef2a20d5d303.svg?invert_in_darkmode" align=middle width=8.55786029999999pt height=14.15524440000002pt/>
functions over the surface.

We can model the problem of parametrization as an energy minimization of
the variation in the <img src="./tex/6dbb78540bd76da3f1625782d42d6d16.svg?invert_in_darkmode" align=middle width=9.41027339999999pt height=14.15524440000002pt/>- and <img src="./tex/6c4adbc36120d62b98deef2a20d5d303.svg?invert_in_darkmode" align=middle width=8.55786029999999pt height=14.15524440000002pt/>-coordinate functions over the surface <img src="./tex/4870d18d47ab6d0e32510c4b1ccf4927.svg?invert_in_darkmode" align=middle width=10.502226899999991pt height=22.55708729999998pt/>:
<p align="center"><img src="./tex/a6c2a7c023e320d7c3b5ff87f567c0a5.svg?invert_in_darkmode" align=middle width=195.10869344999998pt height=37.3519608pt/></p>
 
This familiar energy is called the
[Dirichlet energy](https://en.wikipedia.org/wiki/Dirichlet's_energy).

We may discretize this problem immediately using [piecewise linear
functions](https://en.wikipedia.org/wiki/Piecewise_linear_function) spanned by
<img src="./tex/6dbb78540bd76da3f1625782d42d6d16.svg?invert_in_darkmode" align=middle width=9.41027339999999pt height=14.15524440000002pt/> and <img src="./tex/6c4adbc36120d62b98deef2a20d5d303.svg?invert_in_darkmode" align=middle width=8.55786029999999pt height=14.15524440000002pt/> values at vertices. This corresponds to using the _cotangent
Laplacian_ as <img src="./tex/80637df1ca7533740cc7b3fdd1ab540b.svg?invert_in_darkmode" align=middle width=11.36979854999999pt height=22.55708729999998pt/> in the discrete minimization problem above.

> In the smooth setting, minimizing the variation of <img src="./tex/6dbb78540bd76da3f1625782d42d6d16.svg?invert_in_darkmode" align=middle width=9.41027339999999pt height=14.15524440000002pt/> and <img src="./tex/6c4adbc36120d62b98deef2a20d5d303.svg?invert_in_darkmode" align=middle width=8.55786029999999pt height=14.15524440000002pt/> will lead to
> an injective mapping if the boundary is constrained to a closed convex curve.
> In the discrete setting, poor triangle shapes in the original 3D mesh could
> lead to _negative_ cotangent weights <img src="./tex/64e70e84545b2941bed8aa7fe2211cde.svg?invert_in_darkmode" align=middle width=22.523917349999987pt height=14.15524440000002pt/> so the positive stiffness
> weight assumption of [Tutte's
> theorem](https://en.wikipedia.org/wiki/Tutte_embedding) is broken and
> fold-overs _might_ occur. Keep in mind that positive weights are a
> _sufficient_ condition for injectivity, but this [does not
> imply](https://en.wikipedia.org/wiki/Denying_the_antecedent) that having a
> few negative weights will necessarily cause a fold-over. Even so, Floater
> proposes an alternative discrete Laplacian in ["Mean value coordinates"](papers/Floater03.pdf) in
> 2003 that retains some nice shape-preserving properties without negative
> weights.

Modeling distortion as an integral of variation over the given 3D surface is
going in the right direction, but so far we are treating <img src="./tex/6dbb78540bd76da3f1625782d42d6d16.svg?invert_in_darkmode" align=middle width=9.41027339999999pt height=14.15524440000002pt/> and <img src="./tex/6c4adbc36120d62b98deef2a20d5d303.svg?invert_in_darkmode" align=middle width=8.55786029999999pt height=14.15524440000002pt/>
_separately_. Intuitively <img src="./tex/6dbb78540bd76da3f1625782d42d6d16.svg?invert_in_darkmode" align=middle width=9.41027339999999pt height=14.15524440000002pt/> and <img src="./tex/6c4adbc36120d62b98deef2a20d5d303.svg?invert_in_darkmode" align=middle width=8.55786029999999pt height=14.15524440000002pt/> cannot "talk" to one-another during
optimization. There's no reason to expect that they will be able to minimize
_area distortation_ and _angle distortion_ directly. For that we will need to
consider <img src="./tex/6dbb78540bd76da3f1625782d42d6d16.svg?invert_in_darkmode" align=middle width=9.41027339999999pt height=14.15524440000002pt/> and <img src="./tex/6c4adbc36120d62b98deef2a20d5d303.svg?invert_in_darkmode" align=middle width=8.55786029999999pt height=14.15524440000002pt/> simultaneously.

### Least Squares Conformal Mappings

We can reason about distortion in terms of differential quantities of the
mapping from <img src="./tex/4870d18d47ab6d0e32510c4b1ccf4927.svg?invert_in_darkmode" align=middle width=10.502226899999991pt height=22.55708729999998pt/> to <img src="./tex/3177e934cf575c08431076a1a5479ba5.svg?invert_in_darkmode" align=middle width=18.424726649999986pt height=26.76175259999998pt/>. Now, ultimately we are trying to parametrize <img src="./tex/4870d18d47ab6d0e32510c4b1ccf4927.svg?invert_in_darkmode" align=middle width=10.502226899999991pt height=22.55708729999998pt/>
using the <img src="./tex/6dbb78540bd76da3f1625782d42d6d16.svg?invert_in_darkmode" align=middle width=9.41027339999999pt height=14.15524440000002pt/> and <img src="./tex/6c4adbc36120d62b98deef2a20d5d303.svg?invert_in_darkmode" align=middle width=8.55786029999999pt height=14.15524440000002pt/> coordinate functions spanning <img src="./tex/433badc501d4f8a183b14684b47f305e.svg?invert_in_darkmode" align=middle width=18.424726649999986pt height=26.76175259999998pt/>, but in order to
describe energies over these unknown functions we will assume (without
explicitly using) that we have a parameterization of <img src="./tex/4870d18d47ab6d0e32510c4b1ccf4927.svg?invert_in_darkmode" align=middle width=10.502226899999991pt height=22.55708729999998pt/> (e.g., with
coordinates <img src="./tex/332cc365a4987aacce0ead01b8bdcc0b.svg?invert_in_darkmode" align=middle width=9.39498779999999pt height=14.15524440000002pt/> and <img src="./tex/deceeaf6940a8c7a5a02373728002b0f.svg?invert_in_darkmode" align=middle width=8.649225749999989pt height=14.15524440000002pt/>). This way we can write about small changes in the
mapping function <img src="./tex/6dbb78540bd76da3f1625782d42d6d16.svg?invert_in_darkmode" align=middle width=9.41027339999999pt height=14.15524440000002pt/> with respect to moving a small amount on the surface of
<img src="./tex/4870d18d47ab6d0e32510c4b1ccf4927.svg?invert_in_darkmode" align=middle width=10.502226899999991pt height=22.55708729999998pt/> (small changes in <img src="./tex/332cc365a4987aacce0ead01b8bdcc0b.svg?invert_in_darkmode" align=middle width=9.39498779999999pt height=14.15524440000002pt/> and <img src="./tex/deceeaf6940a8c7a5a02373728002b0f.svg?invert_in_darkmode" align=middle width=8.649225749999989pt height=14.15524440000002pt/>).

#### Area distortion

We would like that regions on <img src="./tex/4870d18d47ab6d0e32510c4b1ccf4927.svg?invert_in_darkmode" align=middle width=10.502226899999991pt height=22.55708729999998pt/> have a proportionally similarly sized region
under the <img src="./tex/6dbb78540bd76da3f1625782d42d6d16.svg?invert_in_darkmode" align=middle width=9.41027339999999pt height=14.15524440000002pt/>, <img src="./tex/6c4adbc36120d62b98deef2a20d5d303.svg?invert_in_darkmode" align=middle width=8.55786029999999pt height=14.15524440000002pt/> mapping to <img src="./tex/3177e934cf575c08431076a1a5479ba5.svg?invert_in_darkmode" align=middle width=18.424726649999986pt height=26.76175259999998pt/>. On an infinitesimal scale, a small change
in <img src="./tex/332cc365a4987aacce0ead01b8bdcc0b.svg?invert_in_darkmode" align=middle width=9.39498779999999pt height=14.15524440000002pt/> and <img src="./tex/deceeaf6940a8c7a5a02373728002b0f.svg?invert_in_darkmode" align=middle width=8.649225749999989pt height=14.15524440000002pt/> on <img src="./tex/4870d18d47ab6d0e32510c4b1ccf4927.svg?invert_in_darkmode" align=middle width=10.502226899999991pt height=22.55708729999998pt/> should incur an equally small change in <img src="./tex/6dbb78540bd76da3f1625782d42d6d16.svg?invert_in_darkmode" align=middle width=9.41027339999999pt height=14.15524440000002pt/> and <img src="./tex/6c4adbc36120d62b98deef2a20d5d303.svg?invert_in_darkmode" align=middle width=8.55786029999999pt height=14.15524440000002pt/>. In
other words, the [determinant of the
Jacobian](https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant) of the
mapping should be one:

<p align="center"><img src="./tex/4f2ed257d20f2b84a9ab41042013cb8a.svg?invert_in_darkmode" align=middle width=117.42881039999999pt height=49.315569599999996pt/></p>

where <img src="./tex/d61117ce545c578957bf7e3f0413b249.svg?invert_in_darkmode" align=middle width=85.20516344999999pt height=24.65753399999998pt/> for a square matrix <img src="./tex/d05b996d2c08252f77613c25205a0f04.svg?invert_in_darkmode" align=middle width=14.29216634999999pt height=22.55708729999998pt/>.

> The determinant of the Jacobian of a mapping corresponds to the scale factor
> by which local area expands or shrinks. This quantity also appears during
> [integration by
> substitution](https://en.wikipedia.org/wiki/Integration_by_substitution) when
> multivariate functions are involved.

It is tempting to try to throw this equality into a least squares energy an
minimize it. Unfortunately the determinant is already a quadratic function of
<img src="./tex/6dbb78540bd76da3f1625782d42d6d16.svg?invert_in_darkmode" align=middle width=9.41027339999999pt height=14.15524440000002pt/> and <img src="./tex/6c4adbc36120d62b98deef2a20d5d303.svg?invert_in_darkmode" align=middle width=8.55786029999999pt height=14.15524440000002pt/> so a least-squares energy would be quartic and minimizing it would
be non-trivial. In the deformation assignment, we consider area-distortion
directly when looking for "as-rigid-as-possible" mappings. This idea can also be
applied to parameterization (also leading to a non-linear optimization).
But for this assignment, let us put aside area distortion and focus instead on
angle or aspect-ratio distortion. 

#### Angle distortion

We would also like that local regions on <img src="./tex/4870d18d47ab6d0e32510c4b1ccf4927.svg?invert_in_darkmode" align=middle width=10.502226899999991pt height=22.55708729999998pt/> are parameterized without
[shearing](https://en.wikipedia.org/wiki/Shear_mapping). This ensures that two
orthogonal directions on the surface <img src="./tex/4870d18d47ab6d0e32510c4b1ccf4927.svg?invert_in_darkmode" align=middle width=10.502226899999991pt height=22.55708729999998pt/> correspond to two orthogonal
directions on the parameteric plane <img src="./tex/3177e934cf575c08431076a1a5479ba5.svg?invert_in_darkmode" align=middle width=18.424726649999986pt height=26.76175259999998pt/>. We can capture this by requiring
that a small change in the <img src="./tex/332cc365a4987aacce0ead01b8bdcc0b.svg?invert_in_darkmode" align=middle width=9.39498779999999pt height=14.15524440000002pt/> and <img src="./tex/deceeaf6940a8c7a5a02373728002b0f.svg?invert_in_darkmode" align=middle width=8.649225749999989pt height=14.15524440000002pt/> directions on <img src="./tex/4870d18d47ab6d0e32510c4b1ccf4927.svg?invert_in_darkmode" align=middle width=10.502226899999991pt height=22.55708729999998pt/> corresponds to equal
magnitude, small changes in <img src="./tex/6dbb78540bd76da3f1625782d42d6d16.svg?invert_in_darkmode" align=middle width=9.41027339999999pt height=14.15524440000002pt/> and <img src="./tex/6c4adbc36120d62b98deef2a20d5d303.svg?invert_in_darkmode" align=middle width=8.55786029999999pt height=14.15524440000002pt/> in perpendicular directions:

<p align="center"><img src="./tex/fd159d21efcb590ede4c20639e90f164.svg?invert_in_darkmode" align=middle width=107.9787621pt height=94.67909055pt/></p>

where <img src="./tex/b9faa7ee40e6de9421dd0f8d81df1683.svg?invert_in_darkmode" align=middle width=20.25113804999999pt height=27.91243950000002pt/> indicates the vector <img src="./tex/b0ea07dc5c00127344a1cad40467b8de.svg?invert_in_darkmode" align=middle width=9.97711604999999pt height=14.611878600000017pt/> rotated by <img src="./tex/cce24ebf2e55323325713fbf54e7da16.svg?invert_in_darkmode" align=middle width=23.17361309999999pt height=22.63850490000001pt/>.


> By enlisting [complex
> analysis](https://en.wikipedia.org/wiki/Complex_analysis), we can reinterpret
> the mapping to the real plane <img src="./tex/3177e934cf575c08431076a1a5479ba5.svg?invert_in_darkmode" align=middle width=18.424726649999986pt height=26.76175259999998pt/> as a mapping to the [complex
> plane](https://en.wikipedia.org/wiki/Complex_plane) <img src="./tex/81324f07e9ffb7920321df72cc0bee1b.svg?invert_in_darkmode" align=middle width=11.87217899999999pt height=22.648391699999998pt/>. The angle
> preservation equality above corresponds to the [Cauchy-Riemann
> equations](https://en.wikipedia.org/wiki/Cauchy–Riemann_equations). Complex
> functions that satisfy these equations are called [_**conformal
> maps**_](https://en.wikipedia.org/wiki/Conformal_map).

This equality is linear in <img src="./tex/6dbb78540bd76da3f1625782d42d6d16.svg?invert_in_darkmode" align=middle width=9.41027339999999pt height=14.15524440000002pt/> and <img src="./tex/6c4adbc36120d62b98deef2a20d5d303.svg?invert_in_darkmode" align=middle width=8.55786029999999pt height=14.15524440000002pt/>. We can immediately build a quadratic
energy that minimizes deviation from satisfying this equation over the surface
<img src="./tex/4870d18d47ab6d0e32510c4b1ccf4927.svg?invert_in_darkmode" align=middle width=10.502226899999991pt height=22.55708729999998pt/> in a [least squares sense](https://en.wikipedia.org/wiki/Least_squares):

<p align="center"><img src="./tex/4b54ec32bf62c8d2f78f0c6a2d2079ec.svg?invert_in_darkmode" align=middle width=197.29578825pt height=37.3519608pt/></p>


This energy was employed for surface parameterization of triangle meshes as
early as ["Intrinsic parameterizations of surface meshes" [Desbrun et al. 2002]](papers/desbrun02.pdf)
and ["Least squares conformal maps for automatic texture atlas generation" [Lévy
et al. 2002]](papers/Levy02.pdf). Written in this form, it's perhaps not obvious how we can discretize
this over a triangle mesh. Let us massage the equations a bit, starting by
expanding the squared term:


<p align="center"><img src="./tex/22bcf8bcdc8039a2265aeaa859c64a5b.svg?invert_in_darkmode" align=middle width=304.66563105pt height=39.452455349999994pt/></p>


We should recognize the first two terms as the [Dirichlet
energies](https://en.wikipedia.org/wiki/Dirichlet's_energy) of <img src="./tex/6dbb78540bd76da3f1625782d42d6d16.svg?invert_in_darkmode" align=middle width=9.41027339999999pt height=14.15524440000002pt/> and <img src="./tex/6c4adbc36120d62b98deef2a20d5d303.svg?invert_in_darkmode" align=middle width=8.55786029999999pt height=14.15524440000002pt/>.  The third term is
at first glance not familiar. Let's massage it a bit by expanding the gradient
and dot product:

<p align="center"><img src="./tex/a6aff6fb96e244c4fe8b81505517d08d.svg?invert_in_darkmode" align=middle width=340.4940165pt height=157.78318875pt/></p>

where we end up with the integrated [determinant of the
Jacobian](https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant) of the
<img src="./tex/6dbb78540bd76da3f1625782d42d6d16.svg?invert_in_darkmode" align=middle width=9.41027339999999pt height=14.15524440000002pt/> and <img src="./tex/6c4adbc36120d62b98deef2a20d5d303.svg?invert_in_darkmode" align=middle width=8.55786029999999pt height=14.15524440000002pt/> mapping over <img src="./tex/4870d18d47ab6d0e32510c4b1ccf4927.svg?invert_in_darkmode" align=middle width=10.502226899999991pt height=22.55708729999998pt/>. By the rules of [integration by
substitution](https://en.wikipedia.org/wiki/Integration_by_substitution), this
is equivalent to integrating the unit density function over the image of the
mapping, i.e., the **_signed_** area of the flattened surface. If we apply
[Stoke's theorem](https://en.wikipedia.org/wiki/Stokes%27_theorem) we can
convert this area integral into a boundary integral:

<p align="center"><img src="./tex/e90b242e84abc3ff6f4b7341edba76e3.svg?invert_in_darkmode" align=middle width=349.00928205pt height=58.137570149999995pt/></p>

where <img src="./tex/b56595d2a30a0af329086562ca12d521.svg?invert_in_darkmode" align=middle width=10.502226899999991pt height=14.611878600000017pt/> is the unit vector pointing in the outward direction along the
boundary of the image of the mapping.

If we discretize <img src="./tex/6dbb78540bd76da3f1625782d42d6d16.svg?invert_in_darkmode" align=middle width=9.41027339999999pt height=14.15524440000002pt/> and <img src="./tex/6c4adbc36120d62b98deef2a20d5d303.svg?invert_in_darkmode" align=middle width=8.55786029999999pt height=14.15524440000002pt/> using piecewise-linear functions then the boundary
of the mapping will also be piecewise linear and the boundary integral for the
[**vector area**](https://en.wikipedia.org/wiki/Vector_area) is given by the sum over all
boundary edges of the integral of the position vector <img src="./tex/129c5b884ff47d80be4d6261a476e9f1.svg?invert_in_darkmode" align=middle width=10.502226899999991pt height=14.611878600000017pt/> dotted with that
edge's unit normal vector:

<p align="center"><img src="./tex/b8821167b301f1c2d023ea6588620b6b.svg?invert_in_darkmode" align=middle width=592.03065735pt height=216.7015653pt/></p>

where finally we have a simply quadratic expression: sum over all boundary
edges the determinant of the matrix with vertex positions as columns. This
quadratic form can be written as <img src="./tex/2c9cdcd8cbfa6581f86ca865a03c775e.svg?invert_in_darkmode" align=middle width=52.02617639999999pt height=27.91243950000002pt/> with the _vectorized_
<img src="./tex/6dbb78540bd76da3f1625782d42d6d16.svg?invert_in_darkmode" align=middle width=9.41027339999999pt height=14.15524440000002pt/>- and <img src="./tex/6c4adbc36120d62b98deef2a20d5d303.svg?invert_in_darkmode" align=middle width=8.55786029999999pt height=14.15524440000002pt/>-coordinates of the mapping in <img src="./tex/d648c2d7231e16448b889f709af72e9c.svg?invert_in_darkmode" align=middle width=61.18519604999999pt height=26.76175259999998pt/> and <img src="./tex/d5bf70f4feb9ba841052f71ac4f2768b.svg?invert_in_darkmode" align=middle width=85.88664314999998pt height=26.76175259999998pt/> a sparse matrix involving only values for vertices on the boundary of
<img src="./tex/4870d18d47ab6d0e32510c4b1ccf4927.svg?invert_in_darkmode" align=middle width=10.502226899999991pt height=22.55708729999998pt/>. 


**_Achtung!_** A naive implementation of <img src="./tex/68dfb51ab9104a4b4efc5eaadca1c161.svg?invert_in_darkmode" align=middle width=115.00084365000001pt height=36.71249670000002pt/> into matrix form <img src="./tex/2c9cdcd8cbfa6581f86ca865a03c775e.svg?invert_in_darkmode" align=middle width=52.02617639999999pt height=27.91243950000002pt/> will likely produce an
_asymmetric_ matrix <img src="./tex/96458543dc5abd380904d95cae6aa2bc.svg?invert_in_darkmode" align=middle width=14.29216634999999pt height=22.55708729999998pt/>. From a theoretical point of view, this is fine.
<img src="./tex/96458543dc5abd380904d95cae6aa2bc.svg?invert_in_darkmode" align=middle width=14.29216634999999pt height=22.55708729999998pt/> just needs to compute the signed area of the flattened mesh. However, from
a numerical methods point of view we will almost always need our quadratic
coefficients matrix to be
[_symmetric_](https://en.wikipedia.org/wiki/Symmetric_matrix). Fortunately,
when a matrix is acting as a [quadratic
form](https://en.wikipedia.org/wiki/Quadratic_form) it is trivial to
_symmetrize_. Consider we have some asymmetric matrix <img src="./tex/8b05770d43df4bb11931bc6e9b4b6ce9.svg?invert_in_darkmode" align=middle width=14.29216634999999pt height=30.358891200000016pt/> defining a
quadratic form: <img src="./tex/85f90f61b04a79da97de1c834537976a.svg?invert_in_darkmode" align=middle width=43.41889694999999pt height=30.358891200000016pt/>. The output of a quadratic form is
just a scalar, so it's equal to its transpose: 
<p align="center"><img src="./tex/6b691235cfb408b50a85dffd0c79b8eb.svg?invert_in_darkmode" align=middle width=122.49415365pt height=15.179447249999999pt/></p>

These are also equal to their average:
<p align="center"><img src="./tex/3507582c0f936b4e36777d55f96fb1e3.svg?invert_in_darkmode" align=middle width=248.0769522pt height=55.460441849999995pt/></p>


Putting this together with the Dirichlet energy terms, we can write the
discrete _least squares conformal mappings_ minimization problem as:

<p align="center"><img src="./tex/1ddae5641c712b206471331e8da4c010.svg?invert_in_darkmode" align=middle width=231.60324659999998pt height=64.16014605pt/></p>

where <img src="./tex/e9a4432c2ea64a0aa58762f307040216.svg?invert_in_darkmode" align=middle width=69.85918335pt height=26.17730939999998pt/> is the Dirichlet energy quadratic form (a.k.a.
cotangent Laplacian) and <img src="./tex/37cc2b8abb912e7f757ec6669bb421e1.svg?invert_in_darkmode" align=middle width=85.79531069999999pt height=26.76175259999998pt/> is the resulting (sparse)
quadratic form.

#### Free boundary

Similar to the mass-spring methods above, without constraints the least squares
conformal mapping energy will also have a trivial solution: set <img src="./tex/35531be55273dc37ee90083451d089ff.svg?invert_in_darkmode" align=middle width=14.54330789999999pt height=22.55708729999998pt/> to a
single point.

To avoid this solution, we _could_ "fix two vertices" (as originally suggested by
both [Desbrun et al. 2002] and [Lévy et al. 2002]). However, this will
introduce _bias_. Depending on the two vertices we choose we will get a different
solution.  If we're really unlucky, then we might choose two vertices that the
energy would rather like to place near each other and so placing them at
arbitrary positions will introduce unnecessary distortion (i.e., high energy).

Instead we would like [natural boundary
conditions](https://en.wikipedia.org/wiki/Natural_boundary_condition) (not to
be confused with [Neumann boundary
conditions](https://en.wikipedia.org/wiki/Neumann_boundary_condition)). Natural
boundary conditions minimize the given energy in the absence of explicit (or
_essential_) boundary conditions. Natural boundary conditions are convenient if
we discretize the energy _before_ differentiating to find the minimum. If our
discretization is "good", then natural boundary conditions will fall out for
free (_natural_ indeed!).

To obtain natural boundary conditions without bias _and_ avoid the trivial
solution, we can require that the solution:

  1. minimizes the given energy,
  2. has non-zero norm, and
  3. is [orthogonal](https://en.wikipedia.org/wiki/Orthogonality) to trivial
  solutions.

Let's break these down. The first requirement simply ensures that we're still
minimizing the given energy without [monkeying around with
it](https://en.wikipedia.org/wiki/Regularization_(mathematics)) in any way.

The second requirement adds the constraint that the solution <img src="./tex/35531be55273dc37ee90083451d089ff.svg?invert_in_darkmode" align=middle width=14.54330789999999pt height=22.55708729999998pt/> has _unit
norm_:

<p align="center"><img src="./tex/3ba1f27a3c43560f3956f5e4a2b45cb1.svg?invert_in_darkmode" align=middle width=116.3412294pt height=37.3519608pt/></p>


In our discrete case this corresponds to:

<p align="center"><img src="./tex/8cce898c3286e5afa75569736629d653.svg?invert_in_darkmode" align=middle width=171.40946955pt height=61.9227015pt/></p>

where <img src="./tex/e81a1a11e934fd54074805470b4f9ad1.svg?invert_in_darkmode" align=middle width=76.43450429999999pt height=26.17730939999998pt/> is the _mass matrix_ for a piecewise-linear triangle
mesh and <img src="./tex/b1be205a79507444f2c720851d0a98f3.svg?invert_in_darkmode" align=middle width=85.04189759999998pt height=26.76175259999998pt/> is the sparse, square constraint matrix

This is a _quadratic constraint_. Normally that would be [bad
news](https://en.wikipedia.org/wiki/The_Bad_News_Bears), but this type of
constraint results in a well-studied [generalized Eigen value
problem](https://en.wikipedia.org/wiki/Eigendecomposition_of_a_matrix#Generalized_eigenvalue_problem).

> ##### Generalized Eigenvalue problem 
>
> Consider a discrete quadratic minimization problem in <img src="./tex/dc9365440bfbf500cb54700abd8ecd7b.svg?invert_in_darkmode" align=middle width=50.32902434999999pt height=22.648391699999998pt/>:
>
> <p align="center"><img src="./tex/2a7b083cf8285198ee19aa8e3bc5f4f0.svg?invert_in_darkmode" align=middle width=246.7876158pt height=32.990165999999995pt/></p>
>
> where <img src="./tex/67498efca0dc8b3900e026123be5186e.svg?invert_in_darkmode" align=middle width=93.53485349999998pt height=26.17730939999998pt/> are [positive
> semi-definite](https://en.wikipedia.org/wiki/Positive-definite_matrix#Positive-semidefinite)
> matrices.
>
> We can enforce this constraint via the [Lagrange multiplier
> method](https://en.wikipedia.org/wiki/Lagrange_multiplier) by introducing the
> scalar Lagrange multiplier <img src="./tex/1b109d8b4484cf614f27126d788c510e.svg?invert_in_darkmode" align=middle width=9.58908224999999pt height=22.831056599999986pt/> and looking for the saddle-point of the
> _Lagrangian_:
> <p align="center"><img src="./tex/297ea4b66f9bc0380b23d1ee240ebf25.svg?invert_in_darkmode" align=middle width=247.72572164999997pt height=32.990165999999995pt/></p>
>
> This occurs when <img src="./tex/e25ac18ec7c8f297cb770d63b921eca6.svg?invert_in_darkmode" align=middle width=79.21448864999998pt height=24.65753399999998pt/> and <img src="./tex/eb957d3f1d02134b3e03a0935d5b07b3.svg?invert_in_darkmode" align=middle width=78.56388704999999pt height=24.65753399999998pt/>:
>
> <p align="center"><img src="./tex/4fbd6d670e048123c9242bcb85ca2c31.svg?invert_in_darkmode" align=middle width=217.89852975pt height=40.0567299pt/></p>
>
> This is the canonical form of the 
> [generalized Eigen value
> problem](https://en.wikipedia.org/wiki/Eigendecomposition_of_a_matrix#Generalized_eigenvalue_problem)
> for which there are available numerical algorithms.

Finally, our third constraint is that the solution is orthogonal to the trivial
solutions. There are _two_ trivial solutions. They correspond to mapping all
<img src="./tex/6dbb78540bd76da3f1625782d42d6d16.svg?invert_in_darkmode" align=middle width=9.41027339999999pt height=14.15524440000002pt/> values to zero and all <img src="./tex/6c4adbc36120d62b98deef2a20d5d303.svg?invert_in_darkmode" align=middle width=8.55786029999999pt height=14.15524440000002pt/> values to a constant-but-non-zero value and
[_vice-versa_](https://en.wikipedia.org/wiki/List_of_Latin_phrases_(V)#vice_versa).
These solutions will have _zero_ energy and thus their corresponding
eigenvalues <img src="./tex/1b109d8b4484cf614f27126d788c510e.svg?invert_in_darkmode" align=middle width=9.58908224999999pt height=22.831056599999986pt/> will be zero. The _next_ eigenmode (with next smallest
eigenvalue) will satisfy all of our criteria. See ["Spectral conformal
parameterization" [Mullen et al. 2008]](papers/Mullen08.pdf).

> This eigenvector is sometimes called the [Fiedler
> vector](https://en.wikipedia.org/wiki/Algebraic_connectivity#Fiedler_vector).

#### Canonical rotation

The least squares conformal mapping energy is _invariant_ to translation and
rotation. The eigen decomposition process described above will naturally take
care of "picking" a canonical translation by pulling the solution <img src="./tex/35531be55273dc37ee90083451d089ff.svg?invert_in_darkmode" align=middle width=14.54330789999999pt height=22.55708729999998pt/> toward
the origin. The rotation it returns, however, will be arbitrary.

We can try to find a canonical rotation by using [principle component
analysis](https://en.wikipedia.org/wiki/Principal_component_analysis) on the
returned <img src="./tex/ee5f11272c9cd93256bbf7ba019c3953.svg?invert_in_darkmode" align=middle width=17.96813369999999pt height=14.15524440000002pt/> coordinates in <img src="./tex/1b1d3e72f12695d11aa8bd4d01b7e32e.svg?invert_in_darkmode" align=middle width=71.45921639999999pt height=26.76175259999998pt/> (where now <img src="./tex/35531be55273dc37ee90083451d089ff.svg?invert_in_darkmode" align=middle width=14.54330789999999pt height=22.55708729999998pt/> places all <img src="./tex/6dbb78540bd76da3f1625782d42d6d16.svg?invert_in_darkmode" align=middle width=9.41027339999999pt height=14.15524440000002pt/>
coordinates in the first column and <img src="./tex/6c4adbc36120d62b98deef2a20d5d303.svg?invert_in_darkmode" align=middle width=8.55786029999999pt height=14.15524440000002pt/> coordinates in the second column). 

For mappings with strong reflectional symmetry then singular value
decomposition on the [covariance
matrix](https://en.wikipedia.org/wiki/Covariance) <img src="./tex/6751d46c9e987165f1c42b89df40c53f.svg?invert_in_darkmode" align=middle width=93.60155474999998pt height=27.91243950000002pt/> will produce a rotation that aligns the principle direction of <img src="./tex/35531be55273dc37ee90083451d089ff.svg?invert_in_darkmode" align=middle width=14.54330789999999pt height=22.55708729999998pt/> with
the "<img src="./tex/332cc365a4987aacce0ead01b8bdcc0b.svg?invert_in_darkmode" align=middle width=9.39498779999999pt height=14.15524440000002pt/>"-axis of the parametric domain.

![The least squares conformal mapping of the 3D Ogre mesh with natural boundary
conditions produces a more smooth, less distorted and canonically aligned
parameterization than the Tutte embedding above.](images/keenan-ogre-lscm.jpg)

##### Why is everything squished up in the interior?

![The entire camel head is parameterized _inside_ the neck boundary. The area
distortion for the face is extreme: in the parametric domain the face is tiny;
the checkerboard on 3D face has enormous
stretching.](images/camel-head-lscm.jpg)

If the surface has only a small boundary then _all_ of the surface will have to
be packed inside the interior. We're not directly punishing _area_ distortion
so in order to satisfy the _angle_ distortion. _Freeing_ the boundary helps a
little, but ultimately the only way to mitigate this is to: 1) trade area
distortion for angle distortion or 2) _cut_ (a.k.a. "interrupt") the mapping
with discontinuities (see, e.g., [Goode homolosine projection used for maps of
Earth](https://en.wikipedia.org/wiki/Goode_homolosine_projection)).

Cutting new boundaries is _always_ necessary for parameterizing _closed_
surfaces. There has been much work in the last few years on choosing cuts
automatically, many of which are good candidates for a final implementation
project.

![The animal model is a seemingly closed surface, but it has been cut so that
its parameterization is possible and has relatively low area
distortion.](images/animal-lscm.jpg)


## Tasks

### Blacklist

 - `igl::harmonic`
 - `igl::lscm`
 - `igl::vector_area_matrix`

### Whitelist

 - `igl::boundary_loop`
 - `igl::cotmatrix` (or your previous implementation)
 - `igl::eigs` (Use the `igl::EIGS_TYPE_SM` type)
 - `igl::map_vertices_to_circle`
 - `igl::massmatrix` (or your previous implementation)
 - `igl::min_quad_with_fixed` (for minimizing a quadratic energy subject to
     fixed value constraints)
 - `igl::repdiag`

### `src/tutte.cpp`

Given a 3D mesh (`V`,`F`) with a disk topology (i.e., a manifold with single
boundary), compute a 2D parameterization according to Tutte's mapping inside
the unit disk. All boundary vertices should be mapped to the unit circle and
interior vertices mapped inside the disk _without_ flips.

### `src/vector_area_matrix.cpp`

Constructs the symmetric area matrix `A`, s.t.  `[V.col(0)' V.col(1)'] * A *
[V.col(0); V.col(1)]` is the **vector area** of the mesh (`V`,`F`).

### `src/lscm.cpp`

Given a 3D mesh (`V`,`F`) with boundary compute a 2D parameterization that
minimizes the "least squares conformal" energy:

<p align="center"><img src="./tex/7b9d91347f27ec2cb353403dbad762ac.svg?invert_in_darkmode" align=middle width=165.04016099999998pt height=37.3519608pt/></p>


where <img src="./tex/6dbb78540bd76da3f1625782d42d6d16.svg?invert_in_darkmode" align=middle width=9.41027339999999pt height=14.15524440000002pt/> and <img src="./tex/6c4adbc36120d62b98deef2a20d5d303.svg?invert_in_darkmode" align=middle width=8.55786029999999pt height=14.15524440000002pt/> are the unknown (output) coordinates in the parametric domain
`U`.

Use eigen-decomposition to find an un-biased, non-trivial minimizer. Then use
singular value decomposition to find a canonical rotation to line the principle
axis of <img src="./tex/35531be55273dc37ee90083451d089ff.svg?invert_in_darkmode" align=middle width=14.54330789999999pt height=22.55708729999998pt/> with the <img src="./tex/332cc365a4987aacce0ead01b8bdcc0b.svg?invert_in_darkmode" align=middle width=9.39498779999999pt height=14.15524440000002pt/>-axis of the parametric domain.
