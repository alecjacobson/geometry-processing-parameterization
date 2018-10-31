#include "tutte.h"
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix.h>

void tutte(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  // Replace with your code

  // Get boundary loop
  Eigen::VectorXi bloop;
  igl::boundary_loop(F, bloop);

  // Map boundary to unit circle
  Eigen::MatrixXd UV;
  igl::map_vertices_to_circle(V, bloop, UV);

  // Compute cotmatrix
  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V, F, L);

  // Solve
  igl::min_quad_with_fixed_data<double> data;
  bool success = 
    igl::min_quad_with_fixed_precompute(L, bloop, Eigen::SparseMatrix<double>(), false, data);

  if (!success)
    throw std::runtime_error("[tutte] Precompute returned false.");

  success = min_quad_with_fixed_solve(
    data, Eigen::VectorXd::Zero(data.n), UV, Eigen::MatrixXd(), U);

  if (!success)
    throw std::runtime_error("[tutte] Solve returned false.");
}

