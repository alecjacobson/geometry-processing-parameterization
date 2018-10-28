#include "tutte.h"
#include <igl/cotmatrix.h>
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/min_quad_with_fixed.h>

using namespace Eigen;

void tutte(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & U)
{
  // compute the cotangent matrix L
  SparseMatrix<double> L;
  igl::cotmatrix(V, F, L);

  // construct boundary constraints: convex polygon UV
  VectorXi b;
  igl::boundary_loop(F, b);
  MatrixXd UV;
  igl::map_vertices_to_circle(V, b, UV);

  // compute U
  SparseMatrix<double> Aeq;
  igl::min_quad_with_fixed_data<double> data;
  igl::min_quad_with_fixed_precompute(L, b, Aeq, false, data);
  MatrixXd Beq, B;
  B = MatrixXd::Zero(V.rows(), 2);
  igl::min_quad_with_fixed_solve(data, B, UV, Beq, U);
}
