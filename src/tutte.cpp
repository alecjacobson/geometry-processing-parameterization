#include "tutte.h"
#include <Eigen/Sparse>
#include <igl/boundary_loop.h>
#include <igl/cotmatrix.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/min_quad_with_fixed.h>

void tutte(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  Eigen::VectorXi bL;
  igl::boundary_loop(F, bL);

  Eigen::MatrixXd UV;
  igl::map_vertices_to_circle(V, bL, UV);

  Eigen::SparseMatrix<double> L(V.rows(), V.rows());
  igl::cotmatrix(V, F, L);

  igl::min_quad_with_fixed_data<double> data;
  igl::min_quad_with_fixed_precompute(L, bL, Eigen::SparseMatrix<double>(), false, data);

  Eigen::VectorXd B = Eigen::VectorXd::Zero(data.n, 1);
  igl::min_quad_with_fixed_solve(data, B, UV, Eigen::MatrixXd(), U);
  U.col(0) = -U.col(0);
}

