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

    ./parameterization [path to mesh with boundary.obj]

![When you're app is implemented correctly you'll be able to cycle through
viewing the mesh, viewing the parameterization, toggling the checkerboard
pattern and switching between parameterization
methods.](images/beetle-cycle-screenshots.gif)

## Background

In this assignment we will explore how to _flatten_ a surface
[embedded](https://en.wikipedia.org/wiki/Embedding) (or even just
[immersed](https://en.wikipedia.org/wiki/Immersion)) in $ℝ^3$ to the flat
plane (i.e. $ℝ^2$).

This process is often referred to as
[parameterization](https://en.wikipedia.org/wiki/Parametrization#Parametrization_techniques)
because the two-dimensional coordinate system of the flattened mesh can now be
interpreted as a parameterization of the 3D surface.

![A triangle mesh of a [VW
Beetle](https://en.wikipedia.org/wiki/Volkswagen_Beetle) is _parameterized_ by
flattening the mesh to the $uv$-plane. There the $u$- and $v$- coordinates
(orange and white lines) can be directly interpreted as a parameterization of
the surface. ](images/beetle-uv-parameterization-low-res.png)

In this assignment, we are given a representation of the surface in 3D as a
[triangle mesh](https://en.wikipedia.org/wiki/Triangle_mesh) with a list of
$n$ vertex positions $\V ∈ ℝ^{n × 3}$, so then the goal is to assign $uv$
coordinates to each vertex $\U ∈ ℝ^{n × 2}$. 

In general, a 3D surface cannot be flattened onto the plane without
_**distortion**_. Some parts of the surface will have to be stretched and other
squished. Surfaces with [topological
handles](https://en.wikipedia.org/wiki/Handle_decomposition) or without [a
boundary](https://en.wikipedia.org/wiki/Surface_(topology)#Closed_surfaces) 

## Mass-spring methods

If we view our triangle mesh surface as a simple
[graph](https://en.wikipedia.org/wiki/Graph_(discrete_mathematics)) then the
surface flattening problem reduces to [graph
drawing](https://en.wikipedia.org/wiki/Graph_drawing). Distortion can be
measured in terms of the relative change in lengths between neighboring
vertices.

We can pose the graph drawing problem as an optimization over node locations so
that the lengths between neighboring vertices are minimized:

\\[
\min_\U  ∑\limits_{\{i,j\} ∈ \E} ‖\u_i - \u_j‖²,
\\]
where $\E∈[1,n]^{k × 2}$ holds a list of edge indices into $\V$. This energy
has a physical interpretation as the [potential
energy](https://en.wikipedia.org/wiki/Potential_energy) of
[mass-spring](https://en.wikipedia.org/wiki/Simple_harmonic_motion#Mass_on_a_spring)
system. Each edge represents a spring with zero rest length, all springs have
uniform [stiffness](https://en.wikipedia.org/wiki/Hooke's_law#Spring_energy),
and all vertices have equivalent (unit) mass.

Without additional constraints, this minimization has a trivial solution: map
all vertices to the same point, e.g., $\u_i = (0\ 0),\ ∀ i$.

We can avoid this by fixing the mapping of certain vertices. If we choose these
fixed vertices arbitrarily we will in general get overlaps in the flattening.
For graph drawing this means that edges cross each other; for surface
parameterization this means that multiple triangles cover the same patch of the
$uv$-plane and some of those triangles are upside down. This problem is often
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
injective (i.e., foldover-free) flattening.

While avoiding foldovers is important, Tutte-style mappings suffer from a
couple problems.

If uniform spring stiffness are used, then the mapping in the $uv$ domain will
try to make all edges the same length. Combined with the boundary constraints,
the flattened mesh will have smoothly varying edge-lengths and near-equilateral
triangles _regardless_ of the triangles shapes and sizes on the surface mesh.

![The Tutte embedding of the 3D Ogre mesh leads to severe distortion.
Overlaying a checkerboard pattern on the 2D domain and visualizing it on the 3D
surface shows the wobbliness of the non-smooth mapping and
stretching.](images/keenan-ogre-tutte.jpg)

We can try to remedy this by introducing a non-uniform weight or spring
stiffness for each edge $w_{ij}$:
\\[
\min_\U  ∑\limits_{\{i,j\} ∈ \E} w_{ij} ‖\u_i - \u_j‖².
\\]

For example, we could weigh the distortion of shorter edges (on the 3D mesh)
more than longer ones: $w_{ij} = 1/‖\v_i - \v_j‖$. See "Parametrization and
smooth approximation of surface triangulations" [Floater 1996]. This will at
best help tame _**length distortion**_. The "shapes" (i.e., aspect ratios) of
triangles will only be indirectly preserved. We need a way to discourage _area
distortion_ and _angle distortion_.

To do this, let's write the energy minimization problem above in matrix form:

\\[
\min_\U = ½ \tr{\U^transpose \L \U},
\\]
where $\L ∈ \R^{n × n}$ is a sparse matrix with:
\\[
L_{ij} = \begin{cases}
w_{ij} & \text{ if $i≠j$ and $∃ \{ij\} ∈ \E$, }\\
-∑\limits_{\ell≠i} L_{i\ell} & \text{ if $i = j$, or } \\
0 & \text{ otherwise}
\end{cases}.
\\]

### Dirichlet energy

We should immediately recognize this sparsity structure from the discrete
Laplacians considered in the previous assignments. If $w_{ij} = 1$, then $\L$
is the _uniform Laplacian_ (a.k.a., [graph
Laplacian](https://en.wikipedia.org/wiki/Laplacian_matrix)). If $w_{ij}$ is
based on edge-lengths, then $\L$ corresponds to a physical spring system. 

But we have more information then edges, since we know that our graph is really
a discrete representation of a surface. If we model distortion as _variation_
then we can model the parametrization as an energy minimization of the
variation in the $u$- and $v$-coordinate functions over the surface $\S$:
\\[
\min_{u,v} ∫_\S ‖∇u‖² + ‖∇v‖² \ dA.
\\]
This familiar energy is called the [Dirichelt
energy](https://en.wikipedia.org/wiki/Dirichlet's_energy).

We may discretize this problem immediately using [piecewise linear
functions](https://en.wikipedia.org/wiki/Piecewise_linear_function) spanned by
$u$ and $v$ values at vertices. This corresponds to using the _cotangent
Laplacian_ as $\L$ in the discrete minimization problem above.

> In the smooth setting, minimizing the variation of $u$ and $v$ will lead to
> an injective mapping if the boundary is constrained to a closed convex curve.
> In the discrete setting, poor triangle shapes in the original 3D mesh could
> lead to _negative_ cotangent weights $w_{ij}$ so the positive stiffness
> weight of [Tutte's theorem](https://en.wikipedia.org/wiki/Tutte_embedding) is
> broken and foldovers _might_ occur. Keep in mind that positive weights are a
> _sufficient_ condition for injectivity, but this [does not
> imply](https://en.wikipedia.org/wiki/Denying_the_antecedent) that having a
> few negative weights will necessarily cause a foldover. Even so, Floater
> proposes an alternative discrete Laplacian in "Mean value coordinates" 2003,
> that retains some nice shape-preserving properties without negative weights.

Modeling distortion as an integral of variation over the given 3D surface is
going in the right direction, but so far we are treating $u$ and $v$
_separately_. Intuitively $u$ and $v$ cannot "talk" to one-another during
optimization. There's no reason to expect that they will be able to minimize
_area distortation_ and _angle distortion_ directly. For that we will need to
consider $u$ and $v$ simultaneously.

## Least Squares Conformal Mappings

We can reason about distortion in terms of differential quantities of the
mapping from $\S$ to $\R²$. Now, ultimately we are trying to parametrize $\S$
using the $u$ and $v$ coordinate functions spanning $\R^2$, but in order to
describe energies over these unknown functions we will assume (without
explicitly using) that we have a parameterization of $\S$ (e.g., with
coordinates $x$ and $y$). This way we can write about small changes in the
mapping function $u$ with respect to moving a small amount on the surface of
$\S$ (small changes in $x$ and $y$).

### Area distortion

We would like that regions on $\S$ have a proportionally similarly sized region
under the $u$, $v$ mapping to $\R²$. On an infinitesimal scale, a small change
in $x$ and $y$ on $\S$ should incur an equally small change in $u$ and $v$. In
other words, the [determinant of the
Jacobian](https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant) of the
mapping should be one:

\\[
\left|
\begin{array}{cc}
\frac{∂u}{∂x} & \frac{∂u}{∂y} \\
\frac{∂v}{∂x} & \frac{∂v}{∂y} 
\end{array}
\right|
= 1.
\\]

> The determinant of the Jacobian of a mapping corresponds to the scale factor
> by which local area expands or shrinks. This quantity also appears during
> [integration by
> substitution](https://en.wikipedia.org/wiki/Integration_by_substitution) when
> multivariate functions are involved.

It is tempting to try to throw this equality into a least squares energy an
minimize it. Unfortunately the determinant is already a quadratic function of
$u$ and $v$ so a least-squares energy would be quartic and minimizing it would
be non-trivial. We will reinvestigate this _later_ when we look into surface
deformation energies. But for now, let us put aside area distortion and focus
instead on angle or aspect ratio distortion.

### Angle distortion

We would also like that local regions on $\S$ are parameterized without
[shearing](https://en.wikipedia.org/wiki/Shear_mapping). This ensures that two
orthogonal directions on the surface $\S$ correspond to two orthogonal
directions on the parameteric plane $\R²$. We can capture this by requiring
that a small change in the $x$ and $y$ directions on $\S$ corresponds to equal
magnitude, small changes in $u$ and $v$ in perpendicular directions:

\\[
∇u = 
\left( \begin{array}{r}
\frac{∂u}{∂x} \\
\frac{∂u}{∂y} \\
\end{array} \right)
=
\left( \begin{array}{r}
 \frac{∂v}{∂y} \\
-\frac{∂v}{∂x} \\
\end{array} \right)
=
\left( \begin{array}{r}
 \frac{∂v}{∂x} \\
 \frac{∂v}{∂y} \\
\end{array} \right)^⊥
=
∇v^⊥,
\\]
where $\x^⊥$ indicates the vector $\x$ rotated by by 90°.


> If we consider for a moment that we are not trying to flatten a surface in 3D
> to 2D, but rather we're simply trying to a region $\S$ of the 2D plane to
> itself: $M : \R² → \R²$. Here we can reason about _area distortation_ and
> _angle distortion_ in terms of the differential properties of this mapping
> $M$. 
> 
> By enlisting [complex
> analysis](https://en.wikipedia.org/wiki/Complex_analysis), we can reinterpret
> the mapping to the real plane $\R²$ as a mapping to the [complex
> plane](https://en.wikipedia.org/wiki/Complex_plane) $\mathbb{C}$. The angle
> preservation equality above corresponds to the [Cauchy-Riemann
> equations](https://en.wikipedia.org/wiki/Cauchy–Riemann_equations). Complex
> functions that satisfy these equations are called _**conformal functions**_.

This equality is linear in $u$ and $v$. We can immediately build a quadratic
energy that minimizes deviation from satisfying this equation over the surface
$\S$ in a [least squares sense](https://en.wikipedia.org/wiki/Least_squares):

\\[
\min_{u,v} ½ ∫_\S ‖∇u - ∇v^⊥‖² \ dA.
\\]

This energy was employed for surface parameterization of triangle meshes as
early as "Intrinsic parameterizations of surface meshes" [Desbrun et al. 2002]
and "Least squares conformal maps for automatic texture atlas generation" [Lévy
2002]. Written in this form, it's perhaps not obvious how we can discretize
this over a triangle mesh. Let us massage the equations a bit, starting by
expanding the squared term:


\\[
½ ∫_\S ‖∇u - ∇v^⊥‖² \ dA = \\
∫_\S \left(½ ‖∇u‖² + ½ ‖∇v‖² - ∇u ⋅ ∇v^⊥ \right)\ dA.
\\]

We should recognize the first two terms as the [Dirichlet
energy](https://en.wikipedia.org/wiki/Dirichlet's_energy).  The third term is
at first glance not familiar. Let's massage it a bit by expanding the gradient
and dot product:

\\[
\begin{align}
∫_\S  ∇u ⋅ ∇v^⊥ \ dA &= \\
∫_\S \left( \begin{array}{r}
\frac{∂u}{∂x} \\
\frac{∂u}{∂y} \\
\end{array}\right)
\left( \begin{array}{r}
 \frac{∂v}{∂y} \\
-\frac{∂v}{∂x} \\
\end{array} \right) \ dA &= \\
∫_\S \left|
\begin{array}{cc}
\frac{∂u}{∂x} & \frac{∂u}{∂y} \\
\frac{∂v}{∂x} & \frac{∂v}{∂y} 
\end{array}
\right| \ dA &= 
∫_{\left(\begin{array}{c}u(\S) \\ v(\S) \end{array}\right)} 1 \ dA
,
\end{align}
\\]
where we end up with the integrated [determinant of the
Jacobian](https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant) of the
$u$ and $v$ mapping over $\S$. By the rules of [integration by
substitution](https://en.wikipedia.org/wiki/Integration_by_substitution), this
is equivalent to integrating the unit density function over the image of the
mapping, i.e., the **_signed_** area of the flattened surface. If we apply
[Stoke's theorem](https://en.wikipedia.org/wiki/Stokes%27_theorem) we can
convert this area integral into a boundary integral:

\\[
∫_{\left(\begin{array}{c}u(\S) \\ v(\S) \end{array}\right)} 1 \ dA
= ∮_{∂\left(\begin{array}{c}u(\S) \\ v(\S) \end{array}\right)} \x⋅\n \ ds,
\\]
where $\n$ is the unit vector pointing in the outward direction along the
boundary of the image of the mapping.

If we discretize $u$ and $v$ using piecewise-linear functions then the boundary
of the mapping will also be piecewise linear and the boundary integral for the
[**vector area**](https://en.wikipedia.org/wiki/Vector_area) is given by the sum over all
boundary edges of the integral of the position vector $\u$ dotted with that
edge's unit normal vector:

\\[
∮_{∂(\u(\S))} \u(s)⋅\n(s) \ ds = \\
  ∑\limits_{\{i,j\} ∈ ∂\S} ∫_0^1 
    (\u_i + t(\u_j - \u_i))⋅\frac{(\u_j-\u_i)^⊥}{‖\u_j - \u_i‖} \ dt = \\
  ∑\limits_{\{i,j\} ∈ ∂\S} (\u_j-\u_i)⋅(\u_j-\u_i)^⊥ = \\
  ∑\limits_{\{i,j\} ∈ ∂\S} | \u_i\  \u_j |,
\\]
where finally we have a simply quadratic expression: sum over all boundary
edges the determinant of the matrix with vertex positions as columns. This
quadratic form can be written as $\U^\transpose \A \U$ with the _vectorized_
$u$- and $v$-coordinates of the mapping in $\U ∈ \R^{2n}$ and $\A ∈ \R^{2n ×
2n}$ a sparse matrix involving only values for vertices on the boundary of
$\S$. 

Putting this together with the Dirichlet energy terms, we can write the
discrete _least squares conformal mappings_ minimization problem as:

\\[
\min_{\U ∈ \R^{2n}} 
  \U^\transpose 
  \underbrace{
  \left(
  \left(
  \begin{array}{rr}
    \L & 0 \\
    0 & \L
  \end{array}
  \right)
  + \A
  \right)
  }_{\Q}
  \U,
\\]
where $\L ∈ \R^{n × n}$ is the Dirichlet energy quadratic form (a.k.a.
cotangent Laplacian) and $\Q ∈ \R^{2n × 2n}$ is the resulting (sparse)
quadratic form.

### Free boundary

Similar to the mass-spring methods above, without constraints the least squares
conformal mapping energy will also have a trivial solution: set $\U$ to a
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

To obtain natural boundary conditions without bias by picking two arbitrary
fixed vertices _and_ avoid the trivial solution, we can require that the
solution:

  1. minimizes the given energy,
  2. has non-zero norm,
  3. is [orthogonal](https://en.wikipedia.org/wiki/Orthogonality) to trivial
  solutions.

Let's break these down. The first requirement simply ensures that we're still
minimizing the given energy without [monkeying around with
it](https://en.wikipedia.org/wiki/Regularization_(mathematics)) in any way.

The second requirement adds the constraint that the solution $\U$ has _unit
norm_:

\\[
∫_\S ‖\u‖² \ dA = 1.
\\]

In our discrete case this corresponds to:

\\[
\U^\transpose 
  \underbrace{\left(\begin{array}{cc}\M & 0 \\ 0 & \M \end{array}\right)}_{\B} 
  \U = 1,
\\]
where $\M ∈ \R^{n × n}$ is the _mass matrix_ for a piecewise-linear triangle
mesh and $\B ∈ \R^{2n × 2n}$ is the sparse, square constraint matrix

This is a _quadratic constraint_. Normally that would be [bad
news](https://en.wikipedia.org/wiki/The_Bad_News_Bears), but this type of
constraint results in a well-studied [generalized Eigen value
problem](https://en.wikipedia.org/wiki/Eigendecomposition_of_a_matrix#Generalized_eigenvalue_problem).

> #### Generalized Eigenvalue problem 
>
> Consider a discrete quadratic minimization problem in $\v ∈ \R^n$:
>
> \\[
> \min_{\v} ½ \v^\transpose \A \v \text{ subject to } \v^\transpose \B \v = 1,
> \\]
>
> where $\A,\B ∈ \R^{n × n}$ are [positive
> semi-definite](https://en.wikipedia.org/wiki/Positive-definite_matrix#Positive-semidefinite)
> matrices.
>
> We can enforce this constraint via the [Lagrange multiplier
> method](https://en.wikipedia.org/wiki/Lagrange_multiplier) by introducing the
> scalar Lagrange multiplier $λ$ and looking for the saddle-point of the
> _Lagrangian_:
>
> \\[
> \mathcal{L}(\v,λ) = ½ \v^\transpose \A \v + λ (1 - \v^\transpose \B \v ).
> \\]
>
> This occurs when $∂\mathcal{L}/∂\v = 0$ and $∂\mathcal{L}/∂λ = 0$:
>
> \\[
> \begin{align}
> \A \v - λ \B \v  = 0 & → \A \v = λ \B \v,\\
> 1 - \v^\transpose \B \v = 0 &→  \v^\transpose \B \v = 1.
> \end{align}
> \\]
>
> This is the canonical form of the 
> [generalized Eigen value
> problem](https://en.wikipedia.org/wiki/Eigendecomposition_of_a_matrix#Generalized_eigenvalue_problem)
> for which there are available numerical algorithms.

Finally, our third constraint is that the solution is orthogonal to the trivial
solutions. There are _two_ trivial solutions. They correspond to mapping all
$u$ values to zero and all $v$ values to a constant-but-non-zero value and
[_vice-versa_](https://en.wikipedia.org/wiki/List_of_Latin_phrases_(V)#vice_versa).
These solutions will have _zero_ energy and thus their corresponding
eigenvalues $λ$ will be zero. The _next_ eigenmode (with next smallest
eigenvalue) will satisfy all of our criteria. See "Spectral conformal
parameterization" [Mullen et al. 2008].

> This eigenvector is sometimes called the [Fiedler
> vector](https://en.wikipedia.org/wiki/Algebraic_connectivity#Fiedler_vector).

### Canonical rotation

The least squares conformal mapping energy is _invariant_ to translation and
rotation. The eigen decomposition process described above will naturally take
care of "picking" a canonical translation by pulling the solution $\U$ toward
the origin. The rotation it returns, however, will be arbitrary.

We can try to find a canonical rotation by using [principle component
analysis](https://en.wikipedia.org/wiki/Principal_component_analysis) on the
returned $uv$ coordinates in $U ∈ \R^{n × 2}$ (where now $U$ places all $u$
coordinates in the first column and $v$ coordinates in the second column). 

For mappings with strong reflectional symmetry then singular value
decomposition on the [covariance
matrix](https://en.wikipedia.org/wiki/Covariance) $\U^transpose \U ∈ \R^{2 ×
2}$ will produce a rotation that aligns the principle direction of $\U$ with
the "$x$"-axis of the parametric domain.

![The least squares conformal mapping of the 3D Ogre mesh with natural boundary
conditions produces a more smooth, less distorted and canonically aligned
parameterization than the Tutte embedding above.](images/keenan-ogre-lscm.jpg)

#### Why is everything squished up in the interior?

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
surfaces. It is (still in 2017) difficult to choose the cuts in an automatic
way. _Good opportunity for a final project ;-)_.

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
 - `igl::boundary_loop`
 - `igl::cotmatrix` (or your previous implementation)
 - `igl::eigs`
 - `igl::map_vertices_to_circle`
 - `igl::massmatrix`(or your previous implementation)
 - `igl::repdiag`
