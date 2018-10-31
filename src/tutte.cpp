#include "tutte.h"
#include <igl/map_vertices_to_circle.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/boundary_loop.h>
#include <igl/cotmatrix.h>

void tutte(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  // find which vertices make up the boundary loop
  Eigen::VectorXi boundary_indices;
  igl::boundary_loop(F, boundary_indices);
  // pin these to a unit circle
  Eigen::MatrixXd boundary_on_circle;
  igl::map_vertices_to_circle(V, boundary_indices, boundary_on_circle);
  // get the laplacian (L)
  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V, F, L);
  // solve for U! given L and the fixed points on the unit circle
  igl::min_quad_with_fixed_data<double> data;
  igl::min_quad_with_fixed_precompute(L, boundary_indices, Eigen::SparseMatrix<double>(), false, data);
  Eigen::MatrixXd B = Eigen::MatrixXd::Zero(data.n, 2);
  igl::min_quad_with_fixed_solve(data, B, boundary_on_circle, Eigen::MatrixXd(), U);
}

