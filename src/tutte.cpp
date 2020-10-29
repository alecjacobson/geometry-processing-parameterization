#include "tutte.h"
#include <igl/boundary_loop.h>
#include <vector>
#include <igl/map_vertices_to_circle.h>
#include <igl/edges.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <igl/min_quad_with_fixed.h>


void tutte(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  // Replace with your code

  // first find the border vertices of the given mesh
  // assuming the mesh is a manifold, each edge can be 
  // shared by utmost 2 triangle faces. 
  // the edge represents "border" iff it belongs to 1 triangle.


/* 
First tried mapping *ALL* the boundary vertices to unit disk edges. 
But, that would result in a distorted mapping since a car has multiple boundaries (windows, etc).
The correct way is to only pick the largest boundary!

  // use igl::boundary_loop to get the boundaries
  std::vector<std::vector<int> > Loops;
  igl::boundary_loop(F, Loops);

  // extract all the boundary vertices from the boundary loops
  Eigen::VectorXi boundary_vertices_idx(10);
  
  
  int tot_vertices = 0;
  int vertices_arr_counter = 0;

  // over boundary loops
  for (int loop=0; loop<Loops.size(); loop++) {
    tot_vertices += Loops[loop].size();
    boundary_vertices_idx.conservativeResize(tot_vertices);

    // over the vertices per boundary loop
    for (int bnd_vertices=0; bnd_vertices<Loops[loop].size(); bnd_vertices++) {
      boundary_vertices_idx[vertices_arr_counter] = Loops[loop][bnd_vertices];
      vertices_arr_counter += 1;
    }
  }
*/


  // need to pick only the largest boundary for mapping it to the unit 
  // disk's boundary in the parametrized space
  Eigen::VectorXi boundary_vertices_idx;
  igl::boundary_loop(F, boundary_vertices_idx);



  // map the boundary points list boundary_vertices_idx to a 
  // unit circle in the uv-space 
  // to be used as a constraint while solving 

  Eigen::MatrixXd uv_constraints;
  igl::map_vertices_to_circle(V, boundary_vertices_idx, uv_constraints);

  // Prepare the L matrix for Tutte (positive stiffness case of mass-spring energy)
  // Get all the edges using igl::edges
  Eigen::MatrixXi edges_idx;

  igl::edges(F, edges_idx);

  // set triplets for creating the sparse L matrix with positive stiffness
  std::vector<Eigen::Triplet<double> > L_vector;

  double wij = 0; // tutte stifness
  int i, j;

  // iterate over the unq edges
  for (int edg=0; edg<edges_idx.rows(); edg++) {
    i = edges_idx.row(edg)(0);
    j = edges_idx.row(edg)(1);

    // positive weights
    wij = 1.0 / (V.row(i) - V.row(j)).norm();

    // take care of the i!=j case
    L_vector.push_back(Eigen::Triplet<double>(i, j, wij));
    L_vector.push_back(Eigen::Triplet<double>(j, i, wij));

    // take care of the i=j case
    // exploit the property that setFromTriplets adds across duplicates
    L_vector.push_back(Eigen::Triplet<double>(i, i, -1.0 * wij));
    L_vector.push_back(Eigen::Triplet<double>(j, j, -1.0 * wij));

  }

  Eigen::SparseMatrix<double> tutte_L(V.rows(), V.rows());
  tutte_L.setFromTriplets(L_vector.begin(), L_vector.end());
  


  // dummy
  Eigen::SparseMatrix<double> Aeq;

  // data
  igl::min_quad_with_fixed_data<double> data;

  // Solve the problem min_U [0.5 * tr(U'*tutte_L*U)] subject to boundary conditions
  // trace( 0.5*Z'*A*Z + Z'*B + constant )
  igl::min_quad_with_fixed_precompute(tutte_L, boundary_vertices_idx, Aeq, false, data);

  // dummy
  Eigen::MatrixXd B, Beq;
  B = Eigen::MatrixXd::Zero(tutte_L.rows(), 2);

  // Solve!
  igl::min_quad_with_fixed_solve(data, B, uv_constraints, Beq, U);

}

